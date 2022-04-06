#ifndef BIOSEQDB_MOCK_STRUCTURES_H
#define BIOSEQDB_MOCK_STRUCTURES_H

#include <exception>
#include <string>
#include <vector>

namespace mock {
//    static constexpr int threadNumber = 5;

//    static constexpr int aa_number = 21;

//    uint32_t get_aa_id(char aa);

//    char get_aa_by_id(uint32_t aa_id);

//    static int32_t vtml80[aa_number][aa_number] = {
//            {5,       0,       -2,      -2,      -4,      -1,      -3,      -3,      -2,      -3,      -2,      -2,      -1,      -2,      -3,       1,       0,       0,       -5,      -4,      -1},
//            {0,       10,      -7,      -7,      -6,      -3,      -3,      -2,      -6,      -5,      -1,      -4,      -4,      -6,      -4,       0,      -2,       0,       -8,      -1,      -1},
//            {-2,      -7,       7,       2,      -9,      -2,      -1,      -7,      -2,      -8,      -5,       1,      -2,      -1,      -5,      -1,      -2,      -5,       -7,      -7,      -1},
//            {-2,      -7,       2,       6,      -7,      -3,      -2,      -5,       0,      -5,      -4,      -1,      -2,       2,      -3,      -1,      -2,      -4,       -8,      -4,      -1},
//            {-4,      -6,      -9,      -7,       8,      -6,      -1,      -1,      -7,       0,       0,      -5,      -5,      -4,      -6,      -3,      -4,      -2,        1,       3,      -1},
//            {-1,      -3,      -2,      -3,      -6,       7,      -3,      -8,      -3,      -7,      -6,      -1,      -4,      -4,      -3,      -1,      -4,      -5,       -5,      -6,      -1},
//            {-3,      -3,      -1,      -2,      -1,      -3,       9,      -5,      -1,      -3,      -5,       0,      -3,       1,       0,      -1,      -2,      -4,       -2,       1,      -1},
//            {-3,      -2,      -7,      -5,      -1,      -8,      -5,       6,      -5,       1,       1,      -5,      -6,      -5,      -5,      -5,      -2,       3,       -3,      -3,      -1},
//            {-2,      -6,      -2,       0,      -7,      -3,      -1,      -5,       6,      -4,      -2,       0,      -2,       1,       3,      -2,      -1,      -4,       -5,      -4,      -1},
//            {-3,      -5,      -8,      -5,       0,      -7,      -3,       1,      -4,       5,       2,      -5,      -4,      -3,      -4,      -4,      -3,       0,       -2,      -2,      -1},
//            {-2,      -1,      -5,      -4,       0,      -6,      -5,       1,      -2,       2,       8,      -4,      -5,      -2,      -3,      -4,      -1,       0,       -6,      -4,      -1},
//            {-2,      -4,       1,      -1,      -5,      -1,       0,      -5,       0,      -5,      -4,       7,      -4,      -1,      -2,       1,      -1,      -5,       -6,      -2,      -1},
//            {-1,      -4,      -2,      -2,      -5,      -4,      -3,      -6,      -2,      -4,      -5,      -4,       8,      -2,      -3,      -1,      -2,      -4,       -5,      -7,      -1},
//            {-2,      -6,      -1,       2,      -4,      -4,       1,      -5,       1,      -3,      -2,      -1,      -2,       7,       1,      -1,      -2,      -3,       -8,      -5,      -1},
//            {-3,      -4,      -5,      -3,      -6,      -3,       0,      -5,       3,      -4,      -3,      -2,      -3,       1,       7,      -2,      -3,      -5,       -4,      -3,      -1},
//            {1,        0,      -1,      -1,      -3,      -1,      -1,      -5,      -2,      -4,      -4,       1,      -1,      -1,      -2,       5,       1,      -3,       -4,      -3,      -1},
//            {0,       -2,      -2,      -2,      -4,      -4,      -2,      -2,      -1,      -3,      -1,      -1,      -2,      -2,      -3,       1,       6,      -1,       -7,      -4,      -1},
//            {0,        0,      -5,      -4,      -2,      -5,      -4,       3,      -4,       0,       0,      -5,      -4,      -3,      -5,      -3,      -1,       5,       -6,      -4,      -1},
//            {-5,      -8,      -7,      -8,       1,      -5,      -2,      -3,      -5,      -2,      -6,      -6,      -5,      -8,      -4,      -4,      -7,      -6,       11,       1,      -1},
//            {-4,      -1,      -7,      -4,       3,      -6,       1,      -3,      -4,      -2,      -4,      -2,      -7,      -5,      -3,      -3,      -4,      -4,        1,       8,      -1},
//            {-1,      -1,      -1,      -1,      -1,      -1,      -1,      -1,      -1,      -1,      -1,      -1,      -1,      -1,      -1,      -1,      -1,      -1,       -1,      -1,      -1}};

    extern std::vector<std::string> querySequences;
    extern std::vector<std::string> targetSequences;

    // Mock POSTRES structures

    // For logging in Postgres.
    void log_from_cpp(const char *str);

    // For fetching targets (or queries if you'd like to).
    // Overwrites SPI_tuptable.
    // changed interface for const structures in mock
    std::shared_ptr<std::string> get_sequence(const char *table_name, uint64_t sequence_id);

    // Fetches indexes for a given kmer into SPI_tuptable.
    // To access them from C++ use get_ith_index() but you have to do so
    // before calling get_indexes() or get_sequence() again.
    uint32_t get_indexes(const char *table_name, const char *kmer, uint32_t kMerLength);

    // Fetches i-th index from SPI_tuptable (assuming SPI_tuptable contains indexes).
    // changed interface
    void get_ith_index(int32_t i, uint64_t *target_id, uint32_t *position, const char *kmer, uint32_t kMerLength);

    // for tests
    struct TestsParameter {
        const std::vector<std::string> querySequences, targetSequences;

        TestsParameter(std::vector<std::string>&&, std::vector<std::string>&&);
        void setGlobalParameteres();
    };
};
#endif //BIOSEQDB_MOCK_STRUCTURES_H
