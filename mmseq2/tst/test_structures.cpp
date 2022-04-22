#include "gtest/gtest.h"
#include "mmseq2.h"

namespace {
    mmseq2::InputParams simpleInputParams(uint32_t qLen = 0, uint32_t tLen = 0, uint32_t kMerLength = 5,
                                          int32_t kMerGenThreshold = 15, int32_t ungappedAlignmentScore = 15,
                                          double evalTreshold = 0.001, std::vector<std::string> querySequences = {}) {

        std::string matrixName = "blosum62";
        uint32_t threadNumber = 4;
        int32_t gapOpenCost = 11, gapPenaltyCost = 1;

        mmseq2::InputParams::StrPtr targetTableName = std::make_shared<std::string>("TEST TARGET NAME");
        mmseq2::InputParams::StrPtr targetColumnName = std::make_shared<std::string>("TEST COLUMN NAME");
        mmseq2::InputParams::StrPtr substitutionMatrixName = std::make_shared<std::string>(matrixName);
        mmseq2::InputParams::Vec64Ptr qIds = std::make_shared<std::vector<uint64_t>>(qLen, 0);
        mmseq2::InputParams::Vec64Ptr tIds = std::make_shared<std::vector<uint64_t>>(tLen, 0);
        mmseq2::InputParams::VecStrPtr queries = std::make_shared<std::vector<mmseq2::InputParams::StrPtr>>(qLen, nullptr);

        for (uint32_t i = 0; i < qLen; i++) {
            (*qIds)[i] = (uint64_t)i;
        }
        for (uint32_t i = 0; i < tLen; i++) {
            (*tIds)[i] = (uint64_t)i;
        }
        for (uint32_t i = 0; i < qLen; i++) {
            (*queries)[i] = std::make_shared<std::string>(querySequences[i]);
        }

        return {qLen, tLen, qIds, tIds, queries,
                targetTableName, targetColumnName, substitutionMatrixName, kMerLength,
                kMerGenThreshold, ungappedAlignmentScore, evalTreshold, gapOpenCost,
                gapPenaltyCost, threadNumber};
    }
}

TEST(testAminoAcidClass, charId) {
    std::vector<char> aminoAcids(256, 0);
    std::vector<uint32_t> aminoAcidPositions(256, 0);

    // few exampless, results for char should be 0 <= res <= 24
    for (uint32_t i = 0; i <= 255; i++) {
        uint32_t aaId = mmseq2::AminoAcid::charToId((char)i);
        EXPECT_TRUE(0 <= aaId && aaId <= 24);

        aminoAcids[i] = (aaId == 24 ? '*' : (char)i);
        aminoAcidPositions[i] = aaId;
    }

    // in reverse the result should be in table below, without U, O
    std::string charIds = "ARNDCQEGHILKMFPSTWYVBJZX*";
    for (uint32_t i = 0; i <= 255; i++) {
        EXPECT_TRUE(charIds.find(mmseq2::AminoAcid::idToChar(i)) != std::string::npos);
    }

    // moreover it should be the same as before -> aminoAcids
    for (uint32_t i = 0; i <= 255; i++) {
        EXPECT_TRUE(aminoAcids[i] == mmseq2::AminoAcid::idToChar(aminoAcidPositions[i]));
    }
}

TEST(testAminoAcidClass, getters) {
    // unique symbols of aplhabet size
    std::set<char> alph;
    for (uint32_t i = 0; i < mmseq2::AminoAcid::getAlphabetSize(); i++) {
        alph.insert(mmseq2::AminoAcid::idToChar(i));
    }
    EXPECT_EQ(alph.size(), mmseq2::AminoAcid::getAlphabetSize());

    // blosum 62, few examples
    EXPECT_EQ(4, mmseq2::AminoAcid::getPenalty(2, (uint32_t)0, (uint32_t)0));
    EXPECT_EQ(4, mmseq2::AminoAcid::getPenalty(2, 'A', 'A'));
    EXPECT_EQ(0, mmseq2::AminoAcid::getPenalty(2, 'E', 'R'));
    EXPECT_EQ(0, mmseq2::AminoAcid::getPenalty(2, 'V', 'A'));
    EXPECT_EQ(-2, mmseq2::AminoAcid::getPenalty(2, 'W', 'H'));
    EXPECT_EQ(1, mmseq2::AminoAcid::getPenalty(2, 'K', 'E'));
    EXPECT_EQ(-4, mmseq2::AminoAcid::getPenalty(2, 'I', 'G'));
}

TEST(testAminoAcidClass, blosumIdtoMatrixId) {
    EXPECT_EQ(mmseq2::AminoAcid::blosumIdToMatrixId(62), 2);
    EXPECT_EXIT(mmseq2::AminoAcid::blosumIdToMatrixId(100), testing::ExitedWithCode(1), "");
}

TEST(testPrefilterKmerStageResultsClass, getters) {
    mmseq2::PrefilterKmerStageResults kmerStageResults(0);
    kmerStageResults.addDiagonal(0, 0);
    kmerStageResults.addDiagonal(1, 2);

    EXPECT_EQ(kmerStageResults.getTargetsNumber(), 2);

    EXPECT_EQ(kmerStageResults.getTargetId(0), 0);
    EXPECT_EQ(kmerStageResults.getTargetId(1), 1);
    EXPECT_EQ(kmerStageResults.getDiagonal(0), 0);
    EXPECT_EQ(kmerStageResults.getDiagonal(1), 2);

    kmerStageResults.addDiagonal(100, -300);
    kmerStageResults.addDiagonal(3, 7);

    EXPECT_EQ(kmerStageResults.getTargetsNumber(), 4);

    EXPECT_EQ(kmerStageResults.getTargetId(0), 0);
    EXPECT_EQ(kmerStageResults.getTargetId(1), 1);
    EXPECT_EQ(kmerStageResults.getDiagonal(0), 0);
    EXPECT_EQ(kmerStageResults.getDiagonal(1), 2);

    EXPECT_EQ(kmerStageResults.getTargetId(2), 100);
    EXPECT_EQ(kmerStageResults.getTargetId(3), 3);
    EXPECT_EQ(kmerStageResults.getDiagonal(2), -300);
    EXPECT_EQ(kmerStageResults.getDiagonal(3), 7);
}

TEST(testInputParams, getters) {

    uint32_t qLen = 2, tLen = 3;
    uint32_t kMerLength = 5;
    int32_t kMerGenThreshold = 20, ungappedAlignmentScore = 15;
    double evalTreshold = 0.001;
    mmseq2::InputParams inputParams = simpleInputParams(qLen, tLen, kMerLength, kMerGenThreshold,
                                                        ungappedAlignmentScore, evalTreshold,
                                                        {"ABCABC", "TTTTT"});

    EXPECT_EQ(inputParams.getQLen(), 2);
    EXPECT_EQ(inputParams.getTLen(), 3);
    EXPECT_EQ(inputParams.getQIds().get()->size(), 2);
    EXPECT_EQ(inputParams.getQIds().get()->at(0), 0);
    EXPECT_EQ(inputParams.getQIds().get()->at(1), 1);
    EXPECT_EQ(inputParams.getTIds().get()->size(), 3);
    EXPECT_EQ(inputParams.getTIds().get()->at(0), 0);
    EXPECT_EQ(inputParams.getTIds().get()->at(1), 1);
    EXPECT_EQ(inputParams.getTIds().get()->at(2), 2);

    EXPECT_TRUE(*inputParams.getQueries().get()->at(0).get() == "ABCABC");
    EXPECT_TRUE(*inputParams.getQueries().get()->at(1).get() == "TTTTT");
    EXPECT_TRUE(*inputParams.getTargetTableName().get() == "TEST TARGET NAME");
    EXPECT_TRUE(*inputParams.getTargetColumnName().get() == "TEST COLUMN NAME");

    EXPECT_EQ(inputParams.getKMerLength(), kMerLength);
    EXPECT_EQ(inputParams.getKMerGenThreshold(), kMerGenThreshold);
    EXPECT_EQ(inputParams.getUngappedAlignmentScore(), ungappedAlignmentScore);
    EXPECT_EQ(inputParams.getEvalThreshold(), evalTreshold);
    EXPECT_EQ(inputParams.getGapOpenCost(), 11);
    EXPECT_EQ(inputParams.getGapPenaltyCost(), 1);
    EXPECT_EQ(inputParams.getThreadNumber(), 4);
    EXPECT_EQ(inputParams.getSubstitutionMatrixId(), 2);
}

TEST(testQuery, getters) {

    mmseq2::InputParams::StrPtr queryStrPtr = nullptr;
    mmseq2::InputParams inputParams = simpleInputParams();
    mmseq2::InputParams::InputParamsPtr inputParamsPtr = std::make_shared<mmseq2::InputParams>(inputParams);
    mmseq2::Query query(0, queryStrPtr, inputParamsPtr);

    EXPECT_EQ(query.getPrefilterKmerStageResults().getTargetsNumber(), 0);
    query.addMatch(0, 1);
    query.addMatch(2, 3);
    EXPECT_EQ(query.getPrefilterKmerStageResults().getTargetsNumber(), 2);
    EXPECT_EQ(query.getPrefilterKmerStageResults().getTargetId(0), 0);
    EXPECT_EQ(query.getPrefilterKmerStageResults().getDiagonal(1), 3);

    EXPECT_EQ(query.getSubstitutionMatrixId(), 2);
    EXPECT_EQ(query.getKMerLength(), 5);
}

TEST(testQuery, kMerStage) {
    mock::init_mock_test({"ALA"}, {"KAT", "MA", "ALE"}, 3);
    mmseq2::InputParams::StrPtr queryStrPtr = std::make_shared<std::string>("ALA");
    mmseq2::InputParams::InputParamsPtr inputParamsPtr = std::make_shared<mmseq2::InputParams>(simpleInputParams(
            1, 3, 3, /* kmer gen thold */ 0, 0, 0.001, {"ALA"}));
    mmseq2::Query queryTest(0, queryStrPtr, inputParamsPtr);
    queryTest.findPrefilterKmerStageResults();
    // no double hits
    EXPECT_EQ(queryTest.getPrefilterKmerStageResults().getTargetsNumber(), 0);

    mock::init_mock_test({"AB"}, {"**AB"}, 1);
    queryStrPtr = std::make_shared<std::string>("AB");
    inputParamsPtr = std::make_shared<mmseq2::InputParams>(simpleInputParams(1, 3, 1, 0, 0, 0.001, {"AAAAAAABBBBBBB"}));
    queryTest = mmseq2::Query(0, queryStrPtr, inputParamsPtr);
    queryTest.findPrefilterKmerStageResults();
    // raw: A-A (4), B-B (4), A-B (-2), A,B-* (-4), so result should > 0 for diag = 2 in other matches / results < 0
    EXPECT_EQ(queryTest.getPrefilterKmerStageResults().getTargetsNumber(), 1);
    EXPECT_EQ(queryTest.getPrefilterKmerStageResults().getTargetId(0), 0);
    EXPECT_EQ(queryTest.getPrefilterKmerStageResults().getDiagonal(0), 2);
}

TEST(testQuery, executeAlignment) {
    // one double hit with diag -2
    mock::init_mock_test({"**A**B"}, {"A**B**"}, 1);
    mmseq2::InputParams::StrPtr queryStrPtr = std::make_shared<std::string>("**A**B");
    mmseq2::InputParams::InputParamsPtr inputParamsPtr = std::make_shared<mmseq2::InputParams>(simpleInputParams(1, 1, 1, 0, 1000000, 0.001, {"**A**B"}));
    mmseq2::Query queryTest(0, queryStrPtr, inputParamsPtr);
    queryTest.findPrefilterKmerStageResults();

    // ungapped alignment filter out result
    std::mutex mtx;
    std::shared_ptr<std::vector<mmseq2::MmseqResult>> mmseqResult = std::make_shared<std::vector<mmseq2::MmseqResult>>();
    queryTest.executeAlignment(&mtx, mmseqResult);
    EXPECT_EQ(mmseqResult.get()->size(), 0);

    // evalue to big
    inputParamsPtr = std::make_shared<mmseq2::InputParams>(simpleInputParams(1, 1, 1, 0, 0, 0.001, {"**A**B"}));
    queryTest = mmseq2::Query(0, queryStrPtr, inputParamsPtr);
    queryTest.findPrefilterKmerStageResults();
    mmseqResult->clear();
    queryTest.executeAlignment(&mtx, mmseqResult);
    EXPECT_EQ(mmseqResult.get()->size(), 0);

    // pass
    inputParamsPtr = std::make_shared<mmseq2::InputParams>(simpleInputParams(1, 1, 1, 0, 0, 1000, {"**A**B"}));
    queryTest = mmseq2::Query(0, queryStrPtr, inputParamsPtr);
    queryTest.findPrefilterKmerStageResults();
    mmseqResult->clear();
    queryTest.executeAlignment(&mtx, mmseqResult);
    EXPECT_EQ(mmseqResult.get()->size(), 1);
}

