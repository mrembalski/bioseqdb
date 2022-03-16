//
// Created by user on 16/03/2022.
//
#include "mmseq2.h"
#include "mock_structures.h"


void processQueries(uint32_t q_len, uint32_t t_len,
                    uint64_t *q_ids, uint64_t *t_ids,
                    char **queries,
                    char* target_table_name, char* target_column_name,
                    std::mutex *mtx);

void processSingleQuery(uint64_t q_id, char *query,
                        uint32_t t_len, uint64_t *t_ids,
                        char* target_table_name, char* target_column_name);

void mmseq2::cpp_mmseq2(uint32_t q_len, uint32_t t_len,
                uint64_t *q_ids, uint64_t *t_ids,
                char **queries,
                char* target_table_name, char* target_column_name) {
    std::vector<std::thread> workers{};
    std::mutex mtx;

    for (uint32_t i = 0; i < mock::threadNumber; ++i) {
        workers.push_back(std::thread(
            processQueries, q_len, t_len, q_ids, t_ids, queries, target_table_name, target_column_name, &mtx
        ));
    }

    for (std::thread &worker : workers) {
        worker.join();
    }
}

// Returns the smalles id of query that was not processed yet
uint32_t getNextQuery(uint32_t q_len, std::mutex *mtx) {
    uint32_t res = -1;

    mtx->lock();
    static uint32_t nextQuery = 0;

    if (nextQuery < q_len) {
        res = nextQuery++;
    }

    mtx->unlock();

    return res;
}

void processQueries(uint32_t q_len, uint32_t t_len,
                        uint64_t *q_ids, uint64_t *t_ids,
                        char **queries,
                        char* target_table_name, char* target_column_name,
                        std::mutex *mtx) {
    uint32_t nextQuery;

    while ((nextQuery = getNextQuery(q_len, mtx)) != -1) {
        processSingleQuery(q_ids[nextQuery], queries[nextQuery], t_len, t_ids, target_table_name, target_column_name);
    }
}

void processSingleQuery(uint64_t q_id, char *query_str,
                            uint32_t t_len, uint64_t *t_ids,
                            char* target_table_name, char* target_column_name) {

    mmseq2::Query query{q_id, query_str, t_len, t_ids, target_table_name, target_column_name};

    query.findPrefilterKmerStageResults();

    // TODO Marcin: Tutaj masz ready matche.
}

void mmseq2::Query::findPrefilterKmerStageResults() {
    int32_t SMaxSuf = 0;

    for (uint32_t i = 0; i < mock::kMerSize - 1; ++i) {
        SMaxSuf += mock::vtml80[mock::get_aa_id(this->sequence[i])][mock::get_aa_id(this->sequence[i])];
    }

    for (uint32_t diagonalNumber = 0; diagonalNumber + mock::kMerSize <= this->sequence.size(); ++diagonalNumber) {
        SMaxSuf += this->sequence[diagonalNumber + mock::kMerSize];

        std::string kMer = this->sequence.substr(diagonalNumber, mock::kMerSize);

        this->processSimilarKMers(diagonalNumber, kMer, SMaxSuf);

        SMaxSuf -= mock::vtml80[mock::get_aa_id(this->sequence[diagonalNumber])][mock::get_aa_id(this->sequence[diagonalNumber])];
    }
}

void mmseq2::Query::processSimilarKMers(uint32_t diagonalNumber, std::string &kMer, int32_t SMaxSuf,
                                        int32_t Spref, uint32_t indx) {
    if (indx > kMer.size()) {
        return;
    }

    const char currentAA = kMer[indx];

    const uint32_t currentAAId = mock::get_aa_id(currentAA);

    SMaxSuf -= mock::vtml80[currentAAId][currentAAId];

    for (uint32_t aaId = 0; aaId < mock::aa_number; ++aaId) {
        Spref += mock::vtml80[currentAAId][aaId];

        if (Spref + SMaxSuf >= mock::Smin) {
            kMer[indx] = mock::get_aa_by_id(aaId);

            processSingleKmer(diagonalNumber, kMer);

            processSimilarKMers(diagonalNumber, kMer, SMaxSuf, Spref, indx + 1);
        }

        Spref -= mock::vtml80[this->sequence[indx]][this->sequence[aaId]];
    }

    kMer[indx] = currentAA;

    return;
}

void mmseq2::Query::processSingleKmer(uint32_t diagonal, std::string &kMer) {
    uint32_t n = mock::get_indexes(targetTableName.c_str(), kMer.c_str());

    for (uint32_t i = 0; i < n; ++i) {
        uint64_t target_id;
        uint32_t position;

        mock::get_ith_index(i, &target_id, &position);

        if (diagonalPreVVisited[target_id]
        && diagonalPreVVisited[target_id] == position - diagonal) {
            addMatch(target_id, position - diagonal);
        }

        diagonalPrev[target_id] = position - diagonal;
        diagonalPreVVisited[target_id] = true;
    }

}