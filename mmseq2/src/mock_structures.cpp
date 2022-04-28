#include "mock_structures.h"
#include <iostream>
#include <utility>

int mock::c;
std::vector<std::string> mock::querySequences;
std::vector<std::string> mock::targetSequences;
std::unordered_map<std::string, std::vector<std::pair<uint32_t, uint32_t>>> mock::hits; // id, diag

void mock::log_from_cpp(const char *str) {
    std::cout << str;
}

std::shared_ptr<std::string> mock::get_sequence(const char *table_name, uint64_t sequence_id) {
    return std::make_shared<std::string>(mock::targetSequences[sequence_id].c_str());
}

uint32_t mock::get_indexes(const char *table_name, const char *kmer, uint32_t kMerLength) {
    if (hits.find(kmer) == hits.end()) {
        return 0;
    }
    return hits[kmer].size();
}

// we assume that grt ith index will with safe kmer, i
void mock::get_ith_index(int32_t i, uint64_t *target_id, uint32_t *position, const char *kmer, uint32_t kMerLength) {
    auto hit = hits[kmer][i];
    *target_id = hit.first;
    *position = hit.second;
}

void mock::init_mock_test(std::vector<std::string> querySeqs, std::vector<std::string> targetSeqs, uint32_t kMerLength) {
    querySequences = std::move(querySeqs);
    targetSequences = std::move(targetSeqs);
    hits.clear();

    for (uint32_t i = 0; i < targetSequences.size(); i++) {
        std::string target = mock::targetSequences[i];
        if (target.size() < kMerLength) {
            continue;
        }
        for (uint j = 0; j <= target.size() - kMerLength; j++) {
            std::string kMer = target.substr(j, kMerLength);
            mock::hits[kMer].push_back({i, j});
        }
    }
}

