#ifndef BIOSEQDB_MOCK_STRUCTURES_H
#define BIOSEQDB_MOCK_STRUCTURES_H

#include <exception>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <mutex>
#include <algorithm>
#include <memory>

namespace mock {
    extern int c;
    extern std::vector<std::string> querySequences;
    extern std::vector<std::string> targetSequences;
    extern std::unordered_map<std::string, std::vector<std::pair<uint32_t, uint32_t>>> hits; // id, diag

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

    void init_mock_test(std::vector<std::string> querySeqs,
                        std::vector<std::string> targetSeqs,
                        uint32_t kMerLength);
}
#endif //BIOSEQDB_MOCK_STRUCTURES_H
