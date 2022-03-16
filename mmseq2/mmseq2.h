//
// Created by user on 13/03/2022.
//

#ifndef BIOSEQDB_MMSEQ2_H
#define BIOSEQDB_MMSEQ2_H
#include <exception>
#include <vector>
#include <string>
#include <thread>

namespace mmseq2 {
    void cpp_mmseq2(uint32_t q_len, uint32_t t_len,
                    uint64_t *q_ids, uint64_t *t_ids,
                    char **queries,
                    char* target_table_name, char* target_column_name);

    class PrefilterKmerStageResults {
    private:
        uint32_t queryId;
        std::vector <uint32_t> targetIds;
        std::vector <uint32_t> diagonals;

    public:
        PrefilterKmerStageResults(uint32_t queryId) : queryId{ queryId }, targetIds{ std::vector < uint32_t > {}},
        diagonals{ std::vector < uint32_t > {}}  {}

        PrefilterKmerStageResults() = delete;

        void addDiagonal(uint32_t targetId, uint32_t diagonal) {
            targetIds.push_back(targetId);
            diagonals.push_back(diagonal);
        }

        uint32_t getDiagonal(int index) const {
            return diagonals[index];
        }

        uint32_t getTargetId(int index) const {
            return targetIds[index];
        }

        uint32_t getTargetsNumber() const {
            return targetIds.size();
        }
    };

    class Query {
    private:
        uint64_t queryId;
        std::string sequence;

        uint32_t targetLength;
        uint64_t *targetIds;
        std::string targetTableName;
        std::string targetColumnName;

        PrefilterKmerStageResults prefilterKmerStageResults;

        std::vector<uint32_t> diagonalPrev;

        std::vector<bool> diagonalPreVVisited;

        void processSimilarKMers(uint32_t diagonalNumber, std::string &kMer, int32_t SMaxSuf,
                                 int32_t Spref = 0, uint32_t indx = 0);

        void processSingleKmer(uint32_t diagonal, std::string &kMer);

    public:
        Query() = delete;

        Query(uint64_t queryId, char *query,
              uint32_t t_len, uint64_t *t_ids,
              char *target_table_name, char *target_column_name) : queryId{queryId}, sequence{std::string(query)},
                                                                   targetLength{t_len}, targetIds{t_ids},
                                                                   prefilterKmerStageResults{PrefilterKmerStageResults(queryId)},
                                                                   diagonalPreVVisited{std::vector<bool>(t_len, false)},
                                                                   diagonalPrev{std::vector{t_len}},
                                                                   targetTableName{std::string(target_table_name)},
                                                                   targetColumnName{std::string(target_column_name)} { }

        const PrefilterKmerStageResults& getPrefilterKmerStageResults() const {
            return prefilterKmerStageResults;
        }

        void addMatch(uint32_t targetId, uint32_t diagonal) {
            prefilterKmerStageResults.addDiagonal(targetId, diagonal);
        }

        void findPrefilterKmerStageResults();

        void executeVectorizedUngappedAlignment(); // TODO Marcin

        void executeVectorizedGappedAlignment(); // TODO Marcin

        void displayResults(); // TODO Marcin: to na razie mock funkcja wypisujaca na stdout
    };
}

#endif //BIOSEQDB_MMSEQ2_H
