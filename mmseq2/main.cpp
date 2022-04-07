#include "mmseq2.h"

namespace {
    mmseq2::InputParams::InputParamsPtr prepareInput(
            std::string matrixName, uint32_t kMerLength,
            int32_t kMerGenThreshold, int32_t ungappedAlignmentScore, int32_t evalTreshold,
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
                 int32_t kMerGenThreshold, int32_t ungappedAlignmentScore, int32_t evalTreshold,
                 int32_t gapOpenCost, int32_t gapPenaltyCost, uint32_t threadNumber) {

        mock::querySequences = querySequences;
        mock::targetSequences = targetSequences;

        mmseq2::InputParams::InputParamsPtr inputParams = prepareInput(
                matrixName, kMerLength, kMerGenThreshold, ungappedAlignmentScore,
                evalTreshold, gapOpenCost, gapPenaltyCost, threadNumber);

        auto res = mmseq2::MMSeq2(std::move(inputParams));
        std::cout << "[Next result]\n";
        for (auto& el : res) {
            std::cout << "qId: " << el.queryId << ", ";
            std::cout << "tId: " << el.targetId << "\n";
            std::cout << "rawScore: " << el.rawScore << ", ";
            std::cout << "bitScore: " << el.bitScore << ", ";
            std::cout << "evaLue: " << el.eValue << "\n";
            std::cout << "qStart: " << el.qStart << ", ";
            std::cout << "qEnd: " << el.qEnd << ", ";
            std::cout << "qLen: " << el.qLen << "\n";
            std::cout << "tStart: " << el.tStart << ", ";
            std::cout << "tEnd: " << el.tEnd << ", ";
            std::cout << "tLen: " << el.tLen << "\n";
            std::cout << "qAln: " << el.qAln << "\n";
            std::cout << "tAln: " << el.tAln << "\n";
            std::cout << "cigar: " << el.cigar << "\n";
            std::cout << "alnLen: " << el.alnLen << "\n";
            std::cout << "mismatch: " << el.mismatch << "\n";
            std::cout << "gapOpen: " << el.gapOpen << "\n";
            std::cout << "pident: " << el.pident << "\n";
            std::cout << "\n";
        }
    }
};

int main() {

    runMMSeq2({"DDDDDAAGGGGG"},
            {"AADDDDDCCGGGGGAA"},
            "blosum62", 5, 30, 0, 1, 4, 1, 1);

//    runMMSeq2({"DDDDDDDDDCCGGGGGGGAA", "AAADDDDDDDCCGGGGGGGDD"},
//            {"DDDDDDDAAGGGGGGG"},
//            "blosum62", 7, 42, 0, 1, 4, 1, 2);

//    runMMSeq2({"AAADDDDDDDCCGGGGGGGDD"},
//            {"DDDDDDDAAGGGGGGG", "DDDDDDDDDCCGGGGGGGAA"},
//            "blosum62", 7, 42, 0, 1, 4, 1, 1);

//    runMMSeq2({"AACCTTGG", "ACTGACTGACTG", "TACTCAT"},
//            {"TACGGTAGCTTACTGA", "CTAGCTTACGATGCAAG", "CTTACAGCATACAGCATCGAT"},
//            "blosum62", 5, 10, 10, 1, 4, 1, 5);

//    std::cout << "[Get all hits - example]\n";
//    std::string kmer = "TAC";
//    mock::targetSequences = {"TACGGTAGCTTACTGA", "CTAGCTTACGATGCAAG", "CTTACAGCATACAGCATCGAT"};
//    for (size_t i = 0; i < mock::get_indexes(nullptr, kmer.c_str(), 3); i++) {
//        uint64_t target_id;
//        uint32_t position;
//        mock::get_ith_index((int32_t)i, &target_id, &position, kmer.c_str(), 3);
//        std::cout << "id: " << target_id << ", pos: " << position << "\n";
//    }

//    std::cout << "[Ungapped alignment - example]\n"; // move ungapped to public
//    mock::querySequences = {"AACCTTGG"};
//    mock::targetSequences = {"AAAACCTTGG"};
//    auto querySequencePtr = std::make_shared<std::string>(mock::querySequences[0]);
//    auto targetSequencePtr = std::make_shared<std::string>(mock::targetSequences[0]);
//    mmseq2::InputParams::InputParamsPtr inputParams = prepareInput("blosum62", 0, 0, 0, 0, 0, 0, 0);
//    auto query = mmseq2::Query(0, querySequencePtr, inputParams);
//    for (int i = -5; i <= 5; i++) {
//        std::cout << "diagonal: " << i << ", ";
//        std::cout << "score: " << query.mmseq2::Query::ungappedAlignment(querySequencePtr, targetSequencePtr, i) << '\n';
//    }

//    std::cout << "[Gapped alignment - example]\n"; // move gapped to public
//    std::vector<std::pair<std::string, std::string>> seqs = {
//            {"AATTCCGG", "AAAAATTGGCC"},
//            {"TACTCAT", "CTAGCTTACGATGCAAG"},
//            {"AAAAAACCCCCCTTTTTTGGGGGG", "GGGGGAAACCCCAAGGGGTTGGGGGAAA"},
//            {"AAAAAACCCCCCTTTTTTGGGGGG", "CTTACAGCATACAGCATCGAT"}
//    };
//    for (const auto & seq : seqs) {
//        std::cout << "q: " << seq.first << ", t: " << seq.second << "\n";
//        mock::querySequences = {seq.first};
//        mock::targetSequences = {seq.second};
//        auto querySequencePtr = std::make_shared<std::string>(seq.first);
//        auto targetSequencePtr = std::make_shared<std::string>(seq.second);
//        mmseq2::InputParams::InputParamsPtr inputParams = prepareInput("blosum62", 0, 0, 0, 0, 4, 1, 0);
//        auto query = mmseq2::Query(0, querySequencePtr, inputParams);
//
//        mmseq2::MmseqResult result(0, 0);
//        query.mmseq2::Query::gappedAlignment(querySequencePtr, targetSequencePtr, result);
//        std::cout << result.qAln << "\n" << result.tAln << "\n";
//    }

    return 0;
}