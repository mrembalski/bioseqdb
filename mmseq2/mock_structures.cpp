//
// Created by user on 13/03/2022.
//
#include "mock_structures.h"


class mock::invalid_aa_exception: public std::exception
{
    virtual const char* what() const throw()
    {
        return "Given amino acid does not exist";
    }
} invalid_aa_ex;

static char aa_to_id[mock::aa_number] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X'};

uint32_t mock::get_aa_id(char aa) {
    for (uint32_t id = 0; id < aa_number; ++id) {
        if (aa_to_id[id] == aa) {
            return id;
        }
    }

    throw invalid_aa_ex;
}

char mock::get_aa_by_id(uint32_t aa_id) {
    if (aa_id > 0 && aa_id < aa_number) {
        return aa_to_id[aa_id];
    }

    throw invalid_aa_ex;
}

uint32_t mock::get_indexes(const char *table_name, const char *kmer) {
    return 5; // TODO: Marcin uzupelnij o dane testowe
}

void mock::get_ith_index(int i, uint64_t *target_id, uint32_t *position) {
    return; //TODO: Marcin uzupelnij o dane testowe
}