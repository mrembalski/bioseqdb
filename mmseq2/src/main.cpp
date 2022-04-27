#include "mmseq2.h"

clock_t startTime, endTime;

namespace {
    mmseq2::InputParams::InputParamsPtr prepareInput(
            std::string matrixName, uint32_t kMerLength,
            int32_t kMerGenThreshold, int32_t ungappedAlignmentScore, double evalTreshold,
            int32_t gapOpenCost, int32_t gapPenaltyCost, uint32_t threadNumber) {

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
        mmseq2::InputParams::StrPtr substitutionMatrixName = std::make_shared<std::string>(matrixName);

        return std::make_shared<mmseq2::InputParams>(
                qLen, tLen, qIds, tIds, queries, targetTableName, targetColumnName,
                substitutionMatrixName, kMerLength, kMerGenThreshold, ungappedAlignmentScore,
                evalTreshold, gapOpenCost, gapPenaltyCost, threadNumber);
    }

    void runMMSeq2(std::vector<std::string> &&querySequences, std::vector<std::string> &&targetSequences,
                 std::string&& matrixName, uint32_t kMerLength,
                 int32_t kMerGenThreshold, int32_t ungappedAlignmentScore, double evalTreshold,
                 int32_t gapOpenCost, int32_t gapPenaltyCost, uint32_t threadNumber) {

        mock::init_mock_test(querySequences, targetSequences, kMerLength);

        mmseq2::InputParams::InputParamsPtr inputParams = prepareInput(
                matrixName, kMerLength, kMerGenThreshold, ungappedAlignmentScore,
                evalTreshold, gapOpenCost, gapPenaltyCost, threadNumber);

        startTime = clock();
        auto res = mmseq2::MMSeq2(inputParams);
        endTime = clock();

//        for (auto& el : res) {
//            std::cout << "qId: " << el.getQueryId() << ", ";
//            std::cout << "tId: " << el.getTargetId() << "\n";
//            std::cout << "rawScore: " << el.getRawScore() << ", ";
//            std::cout << "bitScore: " << el.getBitScore() << ", ";
//            std::cout << "eValue: " << el.getEValue() << "\n";
//            std::cout << "qStart: " << el.getQStart() << ", ";
//            std::cout << "qEnd: " << el.getQEnd() << ", ";
//            std::cout << "qLen: " << el.getQLen() << "\n";
//            std::cout << "tStart: " << el.getTStart() << ", ";
//            std::cout << "tEnd: " << el.getTEnd() << ", ";
//            std::cout << "tLen: " << el.getTLen() << "\n";
//            std::cout << "qAln: " << el.getQAln() << "\n";
//            std::cout << "tAln: " << el.getTAln() << "\n";
//            std::cout << "cigar: " << el.getCigar() << "\n";
//            std::cout << "alnLen: " << el.getAlnLen() << "\n";
//            std::cout << "mismatch: " << el.getMismatch() << "\n";
//            std::cout << "gapOpen: " << el.getGapOpen() << "\n";
//            std::cout << "pident: " << el.getPident() << "\n";
//            std::cout << "\n";
//        }
    }
}

void runMmseq(int i) {
    switch (i) {
        case 0:
            runMMSeq2({"DDDDDDDAAGGGGGGG"},
                      {"AADDDDDDDCCGGGGGGGAA"},
                      "blosum62", 7, 20, 0, 1, 4, 1, 1);
            break;
        case 1:
            runMMSeq2({"DDDDDDDDDCCGGGGGGGAA", "AAADDDDDDDCCGGGGGGGDD"},
                      {"DDDDDDDAAGGGGGGG"},
                      "blosum62", 7, 22, 0, 1, 4, 1, 2);
            break;
        case 2:
            runMMSeq2({"AAADDDDDDDCCGGGGGGGDD"},
                      {"DDDDDDDAAGGGGGGG", "DDDDDDDDDCCGGGGGGGAA"},
                      "blosum62", 7, 22, 0, 1, 4, 1, 1);
            break;
        default:
            runMMSeq2({"AACCTTGG", "ACTGACTGACTG", "TACTCAT"},
                      {"TACGGTAGCTTACTGA", "CTAGCTTACGATGCAAG", "CTTACAGCATACAGCATCGAT"},
                      "blosum62", 7, 15, 10, 1, 4, 1, 5);
            break;
    }
}

int main() {

    double times[4] = {0.0, 0.0, 0.0, 0.0};

    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 4; j++) {
            runMmseq(j);
            double elapsed = (double) (endTime - startTime) / CLOCKS_PER_SEC;
            times[j] += elapsed;
        }
    }

    for (double & time : times) {
        time /= 10.0;
        printf("Time: %.5f\n", time);
    }

    return 0;
}
