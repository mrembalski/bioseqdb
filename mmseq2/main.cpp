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
        char* target_table_name = nullptr;
        char* target_column_name = nullptr;

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

    runTest(mock::TestsParameter(7, 10, 15, 11, 1,
                                 {"AAAAAACCCCCCTTTTTTGGGGGG"},
                                 {"GGGGGAAACCCCAAGGGGTTGGGGGAAA"}));
    runTest(mock::TestsParameter(3, 5, 10, 5, 1,
                                 {"AAAAAACCCCCCTTTTTTGGGGGG"},
                                 {"GGGGGAAACCCCAAGGGGTTGGGGGAAA"}));
    runTest(mock::TestsParameter(5, 10, 10, 4, 1,
                                 {"AACCTTGG", "ACTGACTGACTG", "TACTCAT"},
                                 {"TACGGTAGCTTACTGA", "CTAGCTTACGATGCAAG", "CTTACAGCATACAGCATCGAT"}));
    return 0;
}
