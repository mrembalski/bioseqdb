#include "mmseq2.h"

#include <utility>
#include <thread>
#include <cmath>

namespace {
    // from hsp to bit score, using K, lambda
    double evalBitScore(double rawScore) {
        static double K = 3.0, lambda = 0.5;
        return (lambda * rawScore - std::log(K)) / std::log(2);
    }

    // number of expected hits of similar score
    double evalEValue(double bitScore, uint32_t m, uint32_t n) {
        return (double)m * (double)n * (std::pow(2.0, -bitScore));
    }
}

void processQueries(const common::InputParams::InputParamsPtr &inputParams, std::mutex *mtx, uint32_t *nextQuery, std::mutex *resMtx, const common::VecResPtr &resultPtr);

void processSingleQuery(DB::DBconn& dbconn, uint64_t qId, common::InputParams::StrPtr queryStr, common::InputParams::InputParamsPtr inputParams, std::mutex *resMtx, const common::VecResPtr &resultPtr);

common::VecRes mmseq2::MMSeq2(common::InputParams inputParams) {
    common::InputParams::InputParamsPtr inputParamsPtr = std::make_shared<common::InputParams>(inputParams);
    std::vector<std::thread> workers{};
    std::mutex mtx, resMtx;
    uint32_t nextQuery = 0;
    common::VecResPtr resultPtr = std::make_shared<common::VecRes>();

    for (uint32_t i = 0; i < inputParamsPtr.get()->getThreadNumber(); ++i) {
        workers.emplace_back(std::thread(
            processQueries, inputParamsPtr, &mtx, &nextQuery, &resMtx, resultPtr));
    }

    for (std::thread &worker : workers) {
        worker.join();
    }

    return *resultPtr;
}

// Returns the smalles id of query that was not processed yet
uint32_t getNextQuery(uint32_t q_len, std::mutex *mtx, uint32_t *nextQuery) {
    uint32_t res = -1;

    mtx->lock();
    uint32_t tmpNextQuery = *nextQuery;

    if (tmpNextQuery < q_len) {
        res = tmpNextQuery;
        *nextQuery = tmpNextQuery + 1;
    }

    mtx->unlock();

    return res;
}

void processQueries(const common::InputParams::InputParamsPtr &inputParams, std::mutex *mtx, uint32_t *nextQuery, std::mutex *resMtx, const common::VecResPtr &resultPtr) {
    //TODO: names as variables, possibly in InputParams
    DB::DBconn dbconn("my_table", "my_column");

    uint32_t tmpNextQuery;

    while ((tmpNextQuery = getNextQuery(inputParams.get()->getQLen(), mtx, nextQuery)) != -1) {
        processSingleQuery(dbconn, inputParams.get()->getQIds()->at(tmpNextQuery), inputParams.get()->getQueries()->at(tmpNextQuery), inputParams, resMtx, resultPtr);
    }

    dbconn.CloseConnection();
}

void processSingleQuery(DB::DBconn& dbconn, uint64_t qId, common::InputParams::StrPtr queryStr, common::InputParams::InputParamsPtr inputParams, std::mutex *resMtx, const common::VecResPtr &resultPtr) {

    mmseq2::Query query{qId, std::move(queryStr), std::move(inputParams)};

    query.findPrefilterKmerStageResults(dbconn);

    query.executeAlignment(dbconn, resMtx, resultPtr);
}

void mmseq2::Query::findPrefilterKmerStageResults(DB::DBconn& dbconn) {
    int32_t SMaxSuf = 0;

    for (uint32_t i = 0; i < this->kMerLength - 1; ++i) {
        SMaxSuf += AminoAcid::getPenalty(this->substitutionMatrixId, this->sequence.get()->at(i), this->sequence.get()->at(i));
    }

    for (uint32_t kMerPos = 0; kMerPos + this->kMerLength <= this->sequence.get()->length(); ++kMerPos) {
        SMaxSuf += AminoAcid::getPenalty(this->substitutionMatrixId, this->sequence.get()->at(kMerPos + this->kMerLength - 1), this->sequence.get()->at(kMerPos + this->kMerLength - 1));

        std::string kMer = this->sequence.get()->substr(kMerPos, this->kMerLength);

        this->processSimilarKMers(dbconn, kMerPos, kMer, SMaxSuf);

        SMaxSuf -= AminoAcid::getPenalty(this->substitutionMatrixId, this->sequence.get()->at(kMerPos), this->sequence.get()->at(kMerPos));
    }
}

void mmseq2::Query::processSimilarKMers(DB::DBconn& dbconn, uint32_t kMerPos, std::string &kMer, int32_t SMaxSuf,
                                        int32_t Spref, uint32_t indx) {
    if (indx == 0 && evalBitScore(Spref + SMaxSuf) >= this->kMerGenThreshold) {
        processSingleKmer(dbconn, kMerPos, kMer);
    }

    if (indx >= kMer.size()) {
        return;
    }

    const char currentAA = kMer[indx];

    const uint32_t currentAAId = AminoAcid::charToId(currentAA);

    SMaxSuf -= AminoAcid::getPenalty(this->substitutionMatrixId, currentAAId, currentAAId);

    for (uint32_t aaId = 0; aaId < AminoAcid::getAlphabetSize(); ++aaId) {
        Spref += AminoAcid::getPenalty(this->substitutionMatrixId, currentAAId, aaId);

        if (evalBitScore(Spref + SMaxSuf) >= this->kMerGenThreshold) {
            kMer[indx] = AminoAcid::idToChar(aaId);

            if (currentAAId != aaId) {
                processSingleKmer(dbconn, kMerPos, kMer);
            }

            processSimilarKMers(dbconn, kMerPos, kMer, SMaxSuf, Spref, indx + 1);
        }

        Spref -= AminoAcid::getPenalty(this->substitutionMatrixId, currentAAId, aaId);
    }

    kMer[indx] = currentAA;
}

void mmseq2::Query::processSingleKmer(DB::DBconn& dbconn, uint32_t kMerPos, std::string &kMer) {
    // uint32_t n = mock::get_indexes(targetTableName.get()->c_str(), kMer.c_str(), getKMerLength());

    // for (uint32_t i = 0; i < n; ++i) {
    uint32_t i = 0; 
    while (true) {
        try {
            uint64_t target_id;
            uint32_t position;

            // added kmer for new interface
            dbconn.GetIthIndex(kMer.c_str(), (int32_t)i, &target_id, &position);
            
            // mock::get_ith_index((int32_t)i, &target_id, &position, kMer.c_str(), getKMerLength());
            int32_t diagonal = (int32_t)position - (int32_t)kMerPos;

            if (diagonalPreVVisited[target_id] && diagonalPrev[target_id] == diagonal) {
                addMatch(target_id, diagonal);
            }

            diagonalPrev[target_id] = diagonal;
            diagonalPreVVisited[target_id] = true;
        }
        catch (const std::exception& e) {
            // std::cout << "processSingleKmer: got exception at " << i << std::endl;
            break;
        }
        i++;
    }
}

double mmseq2::Query::ungappedAlignment(const StrPtr &querySequence, const StrPtr &targetSequence, int32_t diagonal) const {
    uint32_t querySequenceLength = querySequence.get()->size(), targetSequenceLength = targetSequence.get()->size();
    int32_t queryPosition = 0, queryLastPosition = (int32_t)querySequenceLength - 1;
    int32_t targetPosition = queryPosition + diagonal;

    if (targetPosition < 0) {
        queryPosition -= targetPosition;
        targetPosition = 0;
    }
    if (queryLastPosition + diagonal >= targetSequenceLength) {
        queryLastPosition = (int32_t)targetSequenceLength - diagonal - 1;
    }

    int32_t score = 0, maxScore = 0;

    while (queryPosition <= queryLastPosition) {
        score += AminoAcid::getPenalty(getSubstitutionMatrixId(), querySequence.get()->at(queryPosition), targetSequence.get()->at(targetPosition));
        score = std::max(score, 0);
        maxScore = std::max(maxScore, score);

        queryPosition++;
        targetPosition++;
    }

    return evalBitScore(maxScore);
}

void mmseq2::Query::gappedAlignment(const StrPtr &querySequence, const StrPtr &targetSequence, common::MmseqResult &result) const {
    uint32_t qSeqLen = querySequence.get()->size(), tSeqLen = targetSequence.get()->size();
    int32_t costOp = this->gapOpenCost, costEx = this->costGapExtended;
    // E - gap in row, F - gap in column, H - best score
    std::vector<std::vector<int32_t> > E(qSeqLen, std::vector<int>(tSeqLen, -costOp));
    std::vector<std::vector<int32_t> > F(qSeqLen, std::vector<int>(tSeqLen, -costOp));
    std::vector<std::vector<int32_t> > H(qSeqLen, std::vector<int>(tSeqLen, 0));

    int32_t bestScore = 0;
    uint32_t qPos = -1, tPos = -1;

    // compute matrix
    for (uint32_t qInd = 0; qInd < qSeqLen; qInd++) {
        for (uint32_t tInd = 0; tInd < tSeqLen; tInd++) {

            if (tInd > 0) {
                E[qInd][tInd] = std::max(E[qInd][tInd - 1] - costEx, H[qInd][tInd - 1] - costOp);
                H[qInd][tInd] = std::max(H[qInd][tInd], E[qInd][tInd]);
            }
            if (qInd > 0) {
                F[qInd][tInd] = std::max(F[qInd - 1][tInd] - costEx, H[qInd - 1][tInd] - costOp);
                H[qInd][tInd] = std::max(H[qInd][tInd], F[qInd][tInd]);
            }

            int32_t match = AminoAcid::getPenalty(this->substitutionMatrixId, querySequence.get()->at(qInd), targetSequence.get()->at(tInd));
            H[qInd][tInd] = std::max(H[qInd][tInd], match);
            if (qInd > 0 && tInd > 0) {
                H[qInd][tInd] = std::max(H[qInd][tInd], H[qInd - 1][tInd - 1] + match);
            }

            if (H[qInd][tInd] > bestScore) {
                bestScore = H[qInd][tInd];
                qPos = qInd;
                tPos = tInd;
            }
        }
    }

    result.setQEnd(qPos);
    result.setTEnd(tPos);
    std::string qAl, tAl, cigar;
    int lastAction = -1, identicalMatches = 0;
    auto qInd = (int32_t)qPos, tInd = (int32_t)tPos;

    // backtrace
    while (qInd >= 0 && tInd >= 0 && H[qInd][tInd] > 0) {
        int32_t match = AminoAcid::getPenalty(this->substitutionMatrixId,
                                              querySequence.get()->at(qInd), targetSequence.get()->at(tInd));
        if ((qInd > 0 && tInd > 0 && H[qInd][tInd] == H[qInd - 1][tInd - 1] + match) || (H[qInd][tInd] == match)) {
            if (querySequence.get()->at(qInd) != targetSequence.get()->at(tInd)) {
                result.incrMismatch();
            } else {
                identicalMatches += 1;
            }
            cigar += 'M';
            qAl += querySequence.get()->at(qInd);
            tAl += targetSequence.get()->at(tInd);
            qInd--;
            tInd--;
            lastAction = 0;
        } else if (H[qInd][tInd] == E[qInd][tInd]) {
            result.incrGapOpen();
            while (tInd > 0 && E[qInd][tInd] == E[qInd][tInd - 1] - costEx) {
                cigar += 'D';
                qAl += ' ';
                tAl += targetSequence.get()->at(tInd);
                tInd--;
            }
            if (tInd > 0 && E[qInd][tInd] == H[qInd][tInd - 1] - costOp) {
                cigar += 'D';
                qAl += ' ';
                tAl += targetSequence.get()->at(tInd);
                tInd--;
            }
            lastAction = 1;
        } else if (H[qInd][tInd] == F[qInd][tInd]) {
            result.incrGapOpen();
            while (qInd > 0 && F[qInd][tInd] == F[qInd - 1][tInd] - costEx) {
                cigar += 'I';
                tAl += ' ';
                qAl += querySequence.get()->at(qInd);
                qInd--;
            }
            if (qInd > 0 && F[qInd][tInd] == H[qInd - 1][tInd] - costOp) {
                cigar += 'I';
                tAl += ' ';
                qAl += querySequence.get()->at(qInd);
                qInd--;
            }
            lastAction = 2;
        }
    }

    std::reverse(qAl.begin(), qAl.end());
    std::reverse(tAl.begin(), tAl.end());
    std::reverse(cigar.begin(), cigar.end());

    result.setQStart(qInd + (lastAction % 2 == 0 ? 1 : 0));
    result.setTStart(tInd + (lastAction < 2 ? 1 : 0));
    result.setQAln(qAl);
    result.setTAln(tAl);
    result.setAlnLen(qAl.size());

    result.setRawScore(bestScore);
    result.setBitScore(evalBitScore(result.getRawScore()));
    result.setEValue(evalEValue(result.getBitScore(), result.getQLen(), result.getTLen()));

    result.setPident((double)identicalMatches / (double)qAl.size());
    result.setCigar(cigar);
}

void mmseq2::Query::executeAlignment(DB::DBconn& dbconn, std::mutex *resMtx, const common::VecResPtr &mmseqResult) {
    const PrefilterKmerStageResults &kmerStageResults = mmseq2::Query::getPrefilterKmerStageResults();
    for (uint32_t i = 0; i < kmerStageResults.getTargetsNumber(); i++) {
        const int32_t diagonal = kmerStageResults.getDiagonal((int)i);
        const uint32_t targetId = kmerStageResults.getTargetId((int)i);

        StrPtr querySequence = this->sequence;
        // const std::string &targetSequence = "DDDDDAAGGGGG";
        // StrPtr targetSequence = mock::get_sequence(this->targetColumnName.get()->c_str(), targetId);
        StrPtr targetSequence = dbconn.GetTargetById(targetId);
        
        if (ungappedAlignment(querySequence, targetSequence, diagonal) >= this->ungappedAlignmentScore && filteredTargetIds.find(targetId) == filteredTargetIds.end()) {
            filteredTargetIds.insert(targetId);

            common::MmseqResult result(this->queryId, targetId);
            result.setQLen(querySequence->size());
            result.setTLen(targetSequence->size());
            gappedAlignment(querySequence, targetSequence, result);

            if (result.getEValue() <= this->evalTreshold) {
                resMtx->lock();
                mmseqResult->push_back(result);
                resMtx->unlock();
            }
        }
    }
}
