#include "mmseq2.h"
#include "mock_structures.h"
#include <iostream>
#include <algorithm>
#include <mutex>

void processQueries(uint32_t q_len, uint32_t t_len,
                    uint64_t *q_ids, uint64_t *t_ids,
                    char **queries,
                    char *target_table_name, char *target_column_name,
                    std::mutex *mtx, uint32_t *nextQuery);

void processSingleQuery(uint64_t q_id, char *query,
                        uint32_t t_len, uint64_t *t_ids,
                        char *target_table_name, char *target_column_name);

void mmseq2::cpp_mmseq2(uint32_t q_len, uint32_t t_len,
                        uint64_t *q_ids, uint64_t *t_ids,
                        char **queries,
                        char *target_table_name, char *target_column_name)
{
    std::vector<std::thread> workers{};
    std::mutex mtx;
    uint32_t nextQuery = 0;

    for (uint32_t i = 0; i < mock::threadNumber; ++i)
    {
        workers.emplace_back(
            processQueries, q_len, t_len, q_ids, t_ids, queries, target_table_name, target_column_name, &mtx, &nextQuery);
    }

    for (std::thread &worker : workers)
    {
        worker.join();
    }
}

// Returns the smalles id of query that was not processed yet
uint32_t getNextQuery(uint32_t q_len, std::mutex *mtx, uint32_t *nextQuery)
{
    uint32_t res = -1;

    mtx->lock();
    uint32_t tmpNextQuery = *nextQuery;

    if (tmpNextQuery < q_len)
    {
        res = tmpNextQuery;
        *nextQuery = tmpNextQuery + 1;
    }

    mtx->unlock();

    return res;
}

void processQueries(uint32_t q_len, uint32_t t_len,
                    uint64_t *q_ids, uint64_t *t_ids,
                    char **queries,
                    char *target_table_name, char *target_column_name,
                    std::mutex *mtx, uint32_t *nextQuery)
{
    uint32_t tmpNextQuery;

    while ((tmpNextQuery = getNextQuery(q_len, mtx, nextQuery)) != -1)
    {
        processSingleQuery(q_ids[tmpNextQuery], queries[tmpNextQuery], t_len, t_ids, target_table_name, target_column_name);
    }
}

void processSingleQuery(uint64_t q_id, char *query_str,
                        uint32_t t_len, uint64_t *t_ids,
                        char *target_table_name, char *target_column_name)
{

    mmseq2::Query query{q_id, query_str, t_len, t_ids, target_table_name, target_column_name};

    query.findPrefilterKmerStageResults();

    query.executeAlignment();
}

void mmseq2::Query::findPrefilterKmerStageResults()
{
    int32_t SMaxSuf = 0;

    for (uint32_t i = 0; i < mock::kMerSize - 1; ++i)
    {
        SMaxSuf += mock::vtml80[mock::get_aa_id(this->sequence[i])][mock::get_aa_id(this->sequence[i])];
    }

    for (uint32_t kMerPos = 0; kMerPos + mock::kMerSize <= this->sequence.size(); ++kMerPos)
    {
        SMaxSuf += mock::vtml80[mock::get_aa_id(this->sequence[kMerPos + mock::kMerSize - 1])][mock::get_aa_id(this->sequence[kMerPos + mock::kMerSize - 1])];

        std::string kMer = this->sequence.substr(kMerPos, mock::kMerSize);

        this->processSimilarKMers(kMerPos, kMer, SMaxSuf);

        SMaxSuf -= mock::vtml80[mock::get_aa_id(this->sequence[kMerPos])][mock::get_aa_id(this->sequence[kMerPos])];
    }
}

void mmseq2::Query::processSimilarKMers(uint32_t kMerPos, std::string &kMer, int32_t SMaxSuf,
                                        int32_t Spref, uint32_t indx)
{
    if (indx == 0 && Spref + SMaxSuf >= mock::Smin)
    {
        processSingleKmer(kMerPos, kMer);
    }

    if (indx >= kMer.size())
    {
        return;
    }

    const char currentAA = kMer[indx];

    const uint32_t currentAAId = mock::get_aa_id(currentAA);

    SMaxSuf -= mock::vtml80[currentAAId][currentAAId];

    for (uint32_t aaId = 0; aaId < mock::aa_number; ++aaId)
    {
        Spref += mock::vtml80[currentAAId][aaId];

        if (Spref + SMaxSuf >= mock::Smin)
        {
            kMer[indx] = mock::get_aa_by_id(aaId);

            if (currentAAId != aaId)
            {
                processSingleKmer(kMerPos, kMer);
            }

            processSimilarKMers(kMerPos, kMer, SMaxSuf, Spref, indx + 1);
        }

        Spref -= mock::vtml80[currentAAId][aaId];
    }

    kMer[indx] = currentAA;
}

void mmseq2::Query::processSingleKmer(uint32_t kMerPos, std::string &kMer)
{
    std::vector<std::pair<uint64_t, uint32_t>> indexes_vec;
    if (kMer == "DDDDDAA")
        indexes_vec.push_back({1, 0});
    else if (kMer == "DDDDAAG")
        indexes_vec.push_back({1, 1});
    else if (kMer == "DDDAAGG")
        indexes_vec.push_back({1, 2});
    else if (kMer == "DDAAGGG")
        indexes_vec.push_back({1, 3});
    else if (kMer == "DAAGGGG")
        indexes_vec.push_back({1, 4});
    else if (kMer == "AAGGGGG")
        indexes_vec.push_back({1, 5});
    uint32_t n = indexes_vec.size();

    for (uint32_t i = 0; i < n; ++i)
    {
        uint64_t target_id = indexes_vec[i].first;
        uint32_t position = indexes_vec[i].second;

        // added kmer for new interface
        int32_t diagonal = (int32_t)position - (int32_t)kMerPos;

        if (diagonalPreVVisited[target_id] && diagonalPrev[target_id] == diagonal)
        {
            addMatch(target_id, diagonal);
        }

        diagonalPrev[target_id] = diagonal;
        diagonalPreVVisited[target_id] = true;
    }
}

int32_t mmseq2::Query::ungappedAlignment(const std::string &querySequence, const std::string &targetSequence, int32_t diagonal)
{
    uint32_t querySequenceLength = querySequence.size(), targetSequenceLength = targetSequence.size();
    int32_t queryPosition = 0, queryLastPosition = (int32_t)querySequenceLength - 1;
    int32_t targetPosition = queryPosition + diagonal;

    if (targetPosition < 0)
    {
        queryPosition -= targetPosition;
        targetPosition = 0;
    }
    if (queryLastPosition + diagonal >= targetSequenceLength)
    {
        queryLastPosition = (int32_t)targetSequenceLength - diagonal - 1;
    }

    int32_t score = 0, maxScore = 0;

    while (queryPosition <= queryLastPosition)
    {
        score += mock::vtml80[mock::get_aa_id(querySequence[queryPosition])][mock::get_aa_id(targetSequence[targetPosition])];
        score = std::max(score, 0);
        maxScore = std::max(maxScore, score);

        queryPosition++;
        targetPosition++;
    }

    return maxScore;
}

std::string mmseq2::Query::gappedAlignment(const std::string &querySequence, const std::string &targetSequence)
{
    uint32_t qSeqLen = querySequence.size(), tSeqLen = targetSequence.size();
    int32_t costOp = mock::costGapOpen, costEx = mock::costGapExtend;
    // E - gap in row, F - gap in column, H - best score
    std::vector<std::vector<int32_t>> E(qSeqLen, std::vector<int>(tSeqLen, -costOp));
    std::vector<std::vector<int32_t>> F(qSeqLen, std::vector<int>(tSeqLen, -costOp));
    std::vector<std::vector<int32_t>> H(qSeqLen, std::vector<int>(tSeqLen, 0));

    int32_t bestScore = 0;
    uint32_t qPos = -1, tPos = -1;

    // compute matrix
    for (uint32_t qInd = 0; qInd < qSeqLen; qInd++)
    {
        for (uint32_t tInd = 0; tInd < tSeqLen; tInd++)
        {

            if (tInd > 0)
            {
                E[qInd][tInd] = std::max(E[qInd][tInd - 1] - costEx, H[qInd][tInd - 1] - costOp);
                H[qInd][tInd] = std::max(H[qInd][tInd], E[qInd][tInd]);
            }
            if (qInd > 0)
            {
                F[qInd][tInd] = std::max(F[qInd - 1][tInd] - costEx, H[qInd - 1][tInd] - costOp);
                H[qInd][tInd] = std::max(H[qInd][tInd], F[qInd][tInd]);
            }

            int32_t match = mock::vtml80[mock::get_aa_id(querySequence[qInd])][mock::get_aa_id(targetSequence[tInd])];
            H[qInd][tInd] = std::max(H[qInd][tInd], match);
            if (qInd > 0 && tInd > 0)
            {
                H[qInd][tInd] = std::max(H[qInd][tInd], H[qInd - 1][tInd - 1] + match);
            }

            if (H[qInd][tInd] > bestScore)
            {
                bestScore = H[qInd][tInd];
                qPos = qInd;
                tPos = tInd;
            }

            // std::cout << E[qInd][tInd] << " " << F[qInd][tInd] << " " << H[qInd][tInd] << std::endl;
        }
    }

    std::string qAl, tAl;
    auto qInd = (int32_t)qPos, tInd = (int32_t)tPos;

    // backtrace
    while (qInd >= 0 && tInd >= 0 && H[qInd][tInd] > 0)
    {
        int32_t match = mock::vtml80[mock::get_aa_id(querySequence[qInd])][mock::get_aa_id(targetSequence[tInd])];
        if ((qInd > 0 && tInd > 0 && H[qInd][tInd] == H[qInd - 1][tInd - 1] + match) || (H[qInd][tInd] == match))
        {
            qAl += querySequence[qInd];
            tAl += targetSequence[tInd];
            qInd--;
            tInd--;
        }
        else if (H[qInd][tInd] == E[qInd][tInd])
        {
            while (tInd > 0 && E[qInd][tInd] == E[qInd][tInd - 1] - costEx)
            {
                qAl += ' ';
                tAl += targetSequence[tInd];
                tInd--;
            }
            if (tInd > 0 && E[qInd][tInd] == H[qInd][tInd - 1] - costOp)
            {
                qAl += ' ';
                tAl += targetSequence[tInd];
                tInd--;
            }
        }
        else if (H[qInd][tInd] == F[qInd][tInd])
        {
            while (qInd > 0 && F[qInd][tInd] == F[qInd - 1][tInd] - costEx)
            {
                tAl += ' ';
                qAl += querySequence[qInd];
                qInd--;
            }
            if (qInd > 0 && F[qInd][tInd] == H[qInd - 1][tInd] - costOp)
            {
                tAl += ' ';
                qAl += querySequence[qInd];
                qInd--;
            }
        }
    }

    std::reverse(qAl.begin(), qAl.end());
    std::reverse(tAl.begin(), tAl.end());
    return qAl.append("\n").append(tAl);
}

void mmseq2::Query::executeAlignment()
{
    const PrefilterKmerStageResults &kmerStageResults = mmseq2::Query::getPrefilterKmerStageResults();
    for (uint32_t i = 0; i < kmerStageResults.getTargetsNumber(); i++)
    {
        const int32_t diagonal = kmerStageResults.getDiagonal((int)i);
        const uint32_t targetId = kmerStageResults.getTargetId((int)i);

        const std::string &querySequence = this->sequence;
        const std::string &targetSequence = "DDDDDAAGGGGG";

        if (ungappedAlignment(querySequence, targetSequence, diagonal) >= mock::minUngappedScore && filteredTargetIds.find(targetId) == filteredTargetIds.end())
        {
            filteredTargetIds.insert(targetId);
            std::string swResult = gappedAlignment(querySequence, targetSequence);

            // postgres - stdout
            std::string swOut = "SW alignment for qId: " + std::to_string(this->queryId) +
                                ", tId: " + std::to_string(targetId) + "\n" + swResult + "\n";
            // log_from_cpp(swOut.c_str());
        }
    }
}
