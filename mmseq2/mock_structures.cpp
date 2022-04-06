#include "mock_structures.h"
#include <string>
#include <vector>
#include <iostream>

std::vector<std::string> mock::querySequences;
std::vector<std::string> mock::targetSequences;

// static char aa_to_id[mock::aa_number] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X'};

//uint32_t mock::get_aa_id(char aa) {
//    for (uint32_t id = 0; id < aa_number; ++id) {
//        if (aa_to_id[id] == aa) {
//            return id;
//        }
//    }
//
//    throw invalid_aa_ex;
//}

//char mock::get_aa_by_id(uint32_t aa_id) {
//    if (aa_id >= 0 && aa_id < aa_number) {
//        return aa_to_id[aa_id];
//    }
//
//    throw invalid_aa_ex;
//}

void mock::log_from_cpp(const char *str) {
    std::cout << str;
}

std::shared_ptr<std::string> mock::get_sequence(const char *table_name, uint64_t sequence_id) {
    if (std::string(table_name) == "QUERY") {
        return std::make_shared<std::string>(mock::querySequences[sequence_id].c_str());
    } else {
        return std::make_shared<std::string>(mock::targetSequences[sequence_id].c_str());
    }
}

// get_ith_index need to know somehow about result of get_indexes
// hits can't be global bcs of threads, can't be hold by thread bcs
// this is mock not mmseq function
std::vector<std::pair<uint32_t, int32_t>> kmerHits(const char *kmer, uint32_t kMerLength) {
    std::vector<std::pair<uint32_t, int32_t>> hits;
    std::string kmerPattern(kmer);
    uint32_t targetId = 0;

    for (const auto &targetSequence : mock::targetSequences) {
        std::string targetSeq(targetSequence);
        for (uint32_t kmerPos = 0; kmerPos + kMerLength <= targetSeq.size(); kmerPos++) {
            std::string proposedKmer(targetSeq.begin() + kmerPos, targetSeq.begin() + kmerPos + kMerLength);
            if (kmerPattern == proposedKmer) {
                hits.emplace_back(targetId, kmerPos);
            }
        }
        targetId++;
    }
    return hits;
}

uint32_t mock::get_indexes(const char *table_name, const char *kmer, uint32_t kMerLength) {
    return kmerHits(kmer, kMerLength).size();
}

// added par kmer bcs we don't have any information about kmer in get_indexes
void mock::get_ith_index(int32_t i, uint64_t *target_id, uint32_t *position, const char *kmer, uint32_t kMerLength) {
    auto hits = kmerHits(kmer, kMerLength);
    *target_id = hits[i].first;
    *position = hits[i].second;
}

mock::TestsParameter::TestsParameter(std::vector<std::string> &&querySequences,
                                     std::vector<std::string> &&targetSequences)
                                     : querySequences{querySequences}, targetSequences{targetSequences} {}

void mock::TestsParameter::setGlobalParameteres() {
    mock::querySequences = this->querySequences;
    mock::targetSequences = this->targetSequences;
}
