#include "gtest/gtest.h"
#include "mmseq2.h"
#include "data/example.h"

namespace run {
    struct runInput {
        std::string matrixName;
        uint32_t kMerLength;
        int32_t kMerGenThreshold, ungappedAlignmentScore;
        double evalTreshold;
        int32_t gapOpenCost, gapPenaltyCost;
        uint32_t threadNumber;
    };

    mmseq2::InputParams::InputParamsPtr prepareInput(const runInput& input) {

        uint32_t qLen = mock::querySequences.size();
        uint32_t tLen = mock::targetSequences.size();

        mmseq2::InputParams::Vec64Ptr qIds = std::make_shared<std::vector<uint64_t>>(qLen, 0);
        for (uint32_t i = 0; i < qLen; i++) {
            (*qIds)[i] = (uint64_t)i;
        }
        mmseq2::InputParams::Vec64Ptr tIds = std::make_shared<std::vector<uint64_t>>(tLen, 0);
        for (uint32_t i = 0; i < tLen; i++) {
            (*tIds)[i] = (uint64_t)i;
        }

        mmseq2::InputParams::VecStrPtr queries = std::make_shared<std::vector<mmseq2::InputParams::StrPtr>>(qLen, nullptr);
        for (uint32_t i = 0; i < qLen; i++) {
            (*queries)[i] = std::make_shared<std::string>(mock::querySequences[i]);
        }

        mmseq2::InputParams::StrPtr targetTableName = std::make_shared<std::string>("TNAME");
        mmseq2::InputParams::StrPtr targetColumnName = std::make_shared<std::string>("CNAME");
        mmseq2::InputParams::StrPtr substitutionMatrixName = std::make_shared<std::string>(input.matrixName);

        return std::make_shared<mmseq2::InputParams>(
                qLen, tLen, qIds, tIds, queries, targetTableName, targetColumnName,
                substitutionMatrixName, input.kMerLength, input.kMerGenThreshold, input.ungappedAlignmentScore,
                input.evalTreshold, input.gapOpenCost, input.gapPenaltyCost, input.threadNumber);
    }

    mmseq2::VecRes MMSeq2res(std::vector<std::string> &&querySequences, std::vector<std::string> &&targetSequences, const runInput& input) {
        mock::init_mock_test(querySequences, targetSequences, input.kMerLength);
        mmseq2::InputParams::InputParamsPtr inputParams = prepareInput(input);
        return mmseq2::MMSeq2(inputParams);
    }
}

namespace egData {
    using VecPairIntType = std::vector<std::pair<uint64_t, uint64_t>>;

    // get matches from example with int ids
    VecPairIntType getMatchesWithIntIds() {
        VecPairIntType matches;
        auto query = getQueryIds();
        auto target = getTargetIds();

        for (const auto& el : getMatchesId()) {
            uint64_t qId = std::find(query.begin(), query.end(), el.first) - query.begin();
            uint64_t tId = std::find(target.begin(), target.end(), el.second) - target.begin();
            matches.push_back({qId, tId});
        }
        return matches;
    }
}

// look at data/example.h
TEST(testLegalAlignments, similarityWithExampleMmseqs2) {
    std::vector<run::runInput> simTestInputs {
        {"blosum62", 5, 15, 15, 1, 11, 1, 1},
        {"blosum45", 7, 35, 10, 0.5, 12, 2, 4}
    };
    auto matches = egData::getMatchesWithIntIds();

    for (const auto &input : simTestInputs) {
        auto ress = run::MMSeq2res(getQuerySequences(), getTargetSequences(), input);

        for (const auto &res : ress) {
            std::pair<uint64_t, uint64_t> match = {res.getQueryId(), res.getTargetId()};
            EXPECT_NE(std::find(matches.begin(), matches.end(), match), matches.end());
        }
    }
}

// checked sense of ids, cigar, starts, ends, lens, gapopen, mismatch, pident, rawScore
TEST(testLegalAlignments, correctnessOfResult) {
    VecStr query{"AACCTTGG", "ACTGACTGACTG", "TACTCAT"};
    VecStr target{"TACGGTAGCTTACTGA", "CTAGCTTACGATGCAAG", "CTTACAGCATACAGCATCGAT"};
    int32_t go = 4, ge = 1;
    auto ress = run::MMSeq2res(
            {"AACCTTGG", "ACTGACTGACTG", "TACTCAT"},
            {"TACGGTAGCTTACTGA", "CTAGCTTACGATGCAAG", "CTTACAGCATACAGCATCGAT"},
            {"blosum62", 5, 10, 10, 1, go, ge, 5});

    for (const auto& res : ress) {
        auto qId = res.getQueryId(), tId = res.getTargetId();
        auto cigar = res.getCigar();

        EXPECT_LE(qId, 2);
        EXPECT_LE(tId, 2);

        EXPECT_EQ(res.getQLen(), query[res.getQueryId()].size());
        EXPECT_EQ(res.getTLen(), target[res.getTargetId()].size());

        uint32_t alnLen = cigar.size();
        EXPECT_EQ(alnLen, res.getAlnLen());

        EXPECT_LE(res.getQStart(), res.getQEnd());
        EXPECT_LE(res.getTStart(), res.getTEnd());
        uint32_t qAlnLen = res.getQEnd() - res.getQStart() + 1;
        uint32_t tAlnLen = res.getTEnd() - res.getTStart() + 1;

        std::string qAln = res.getQAln(), tAln = res.getTAln();
        std::string qAlnWithoutSpace = qAln, tAlnWithoutSpace = tAln;
        qAlnWithoutSpace.erase(remove(qAlnWithoutSpace.begin(), qAlnWithoutSpace.end(), ' '), qAlnWithoutSpace.end());
        tAlnWithoutSpace.erase(remove(tAlnWithoutSpace.begin(), tAlnWithoutSpace.end(), ' '), tAlnWithoutSpace.end());

        EXPECT_EQ(qAlnWithoutSpace, query[res.getQueryId()].substr(res.getQStart(), qAlnLen));
        EXPECT_EQ(tAlnWithoutSpace, target[res.getTargetId()].substr(res.getTStart(), tAlnLen));

        uint32_t gapOpen = 0, allMatch = 0, raw = 0, pidentNum = 0;
        for (uint32_t i = 0; i < alnLen; i++) {
            if (cigar[i] == 'I') {
                if (i == 0 || cigar[i - 1] != 'I') {
                    gapOpen++;
                    raw -= go;
                } else {
                    raw -= ge;
                }
                EXPECT_EQ(' ', tAln[i]);
            }
            if (cigar[i] == 'D') {
                if (i == 0 || cigar[i - 1] != 'D') {
                    gapOpen++;
                    raw -= go;
                } else {
                    raw -= ge;
                }
                EXPECT_EQ(' ', qAln[i]);
            }
            if (cigar[i] == 'M') {
                allMatch++;
                raw += mmseq2::AminoAcid::getPenalty(2, qAln[i], tAln[i]); // blosum62
                if (qAln[i] == tAln[i]) {
                    pidentNum++;
                }
            }
        }
        EXPECT_EQ(gapOpen, res.getGapOpen());
        EXPECT_LE(res.getMismatch(), allMatch);
        EXPECT_EQ(res.getPident(), (double)pidentNum / (double)alnLen);
        EXPECT_EQ(raw, res.getRawScore());
    }
}

