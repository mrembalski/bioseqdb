#include "mmseq2.h"
#include "mock_structures.h"

namespace {
    void runTest(mock::TestsParameter&& test) {

        test.setGlobalParameteres();

        uint32_t q_len = mock::querySequences.size();
        uint32_t t_len = mock::targetSequences.size();
        auto *q_ids = new uint64_t[q_len];
        for (uint32_t i = 0; i < q_len; i++) {
            q_ids[i] = (uint64_t)i;
        }
        auto *t_ids = new uint64_t[t_len];
        for (uint32_t i = 0; i < t_len; i++) {
            t_ids[i] = (uint64_t)i;
        }
        char **queries = new char*[q_len];
        for (uint32_t i = 0; i < q_len; i++) {
            queries[i] = new char[test.querySequences[i].size()];
            strcpy(queries[i], test.querySequences[i].c_str());
        }
        static char target_table_name[] = {'T', 'N', 'A', 'M', 'E'};
        static char target_column_name[] = {'C', 'N', 'A', 'M', 'E'};

        mmseq2::cpp_mmseq2(q_len, t_len, q_ids, t_ids,
                queries, target_table_name, target_column_name);

        delete[] q_ids;
        delete[] t_ids;
        for (uint32_t i = 0; i < q_len; i++) {
            delete[] queries[i];
        }
        delete[] queries;
    }
};

int main() {

    runTest(mock::TestsParameter(7, 45, 15, 11, 1,
                                 {"AAAAAACCCCCCTTTTTTGGGGGG"},
                                 {"GGGGGAAACCCCAAGGGGTTGGGGGAAA"}));
    runTest(mock::TestsParameter(3, 5, 10, 5, 1,
                                 {"AAAAAACCCCCCTTTTTTGGGGGG"},
                                 {"GGGGGAAACCCCAAGGGGTTGGGGGAAA"}));
    runTest(mock::TestsParameter(5, 10, 10, 4, 1,
                                 {"AACCTTGG", "ACTGACTGACTG", "TACTCAT"},
                                 {"TACGGTAGCTTACTGA", "CTAGCTTACGATGCAAG", "CTTACAGCATACAGCATCGAT"}));

//    std::cout << "[Ungapped alignment - example]\n"; // move ungapped to public
//    const std::string qSeq = "AACCTTGG", tSeq = "AAAACCTTGG";
//    for (int i = -5; i <= 5; i++) {
//        std::cout << "diagonal: " << i << ", ";
//        std::cout << "score: " << mmseq2::Query::ungappedAlignment(qSeq, tSeq, i) << '\n';
//    }

//    std::cout << "[Gapped alignment - example]\n"; // move gapped to public
//    mock::TestsParameter test(0, 0, 0, 4, 1, {}, {});
//    test.setGlobalParameteres();
//    const std::vector<std::pair<const std::string, const std::string>> seqs = {
//            {"AATTCCGG", "AAAAATTGGCC"},
//            {"TACTCAT", "CTAGCTTACGATGCAAG"},
//            {"AAAAAACCCCCCTTTTTTGGGGGG", "GGGGGAAACCCCAAGGGGTTGGGGGAAA"},
//            {"AAAAAACCCCCCTTTTTTGGGGGG", "CTTACAGCATACAGCATCGAT"}
//    };
//    for (const auto & seq : seqs) {
//        std::cout << "q: " << seq.first << ", t: " << seq.second << "\n";
//        std::cout << mmseq2::Query::gappedAlignment(seq.first, seq.second) << "\n";
//    }

    return 0;
}