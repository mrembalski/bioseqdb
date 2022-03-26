#ifndef BIOSEQDB_MMSEQ2_H
#define BIOSEQDB_MMSEQ2_H
#include <exception>
#include <vector>
#include <string>
#include <thread>
#include <set>

namespace mmseq2 {
    void cpp_mmseq2(uint32_t q_len, uint32_t t_len,
                    uint64_t *q_ids, uint64_t *t_ids,
                    char **queries,
                    char* target_table_name, char* target_column_name);

    class PrefilterKmerStageResults {
    private:
        uint32_t queryId;
        std::vector <uint32_t> targetIds;
        std::vector <int32_t> diagonals;

    public:
        PrefilterKmerStageResults(uint32_t queryId) : queryId{ queryId }, targetIds{ std::vector<uint32_t> {}},
        diagonals{ std::vector<int32_t> {}}  {}

        PrefilterKmerStageResults() = delete;

        void addDiagonal(uint32_t targetId, int32_t diagonal) {
            targetIds.push_back(targetId);
            diagonals.push_back(diagonal);
        }

        int32_t getDiagonal(int index) const {
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

        std::vector<int32_t> diagonalPrev;

        std::vector<bool> diagonalPreVVisited;

        std::set<uint32_t> filteredTargetIds;

        void processSimilarKMers(uint32_t diagonalNumber, std::string &kMer, int32_t SMaxSuf,
                                 int32_t Spref = 0, uint32_t indx = 0);

        void processSingleKmer(uint32_t diagonal, std::string &kMer);

        // returns best score on diagonal
        static int32_t ungappedAlignment(const std::string& querySequence, const std::string& targetSequence, int32_t diagonal);

        // returns the best alignment
        static std::string gappedAlignment(const std::string& querySequence, const std::string& targetSequence);

    public:
        Query() = delete;

        Query(uint64_t queryId, char *query,
              uint32_t t_len, uint64_t *t_ids,
              char *target_table_name, char *target_column_name) : queryId{queryId}, sequence{std::string(query)},
                                                                   targetLength{t_len}, targetIds{t_ids},
                                                                   prefilterKmerStageResults{PrefilterKmerStageResults(queryId)},
                                                                   diagonalPreVVisited{std::vector<bool>(t_len, false)},
                                                                   diagonalPrev{std::vector<int32_t>(t_len, 0)},
                                                                   targetTableName{std::string(target_table_name)},
                                                                   targetColumnName{std::string(target_column_name)} { }

        const PrefilterKmerStageResults& getPrefilterKmerStageResults() const {
            return prefilterKmerStageResults;
        }

        void addMatch(uint32_t targetId, int32_t diagonal) {
            prefilterKmerStageResults.addDiagonal(targetId, diagonal);
        }

        void findPrefilterKmerStageResults();

        void executeAlignment();
    };
}

#endif //BIOSEQDB_MMSEQ2_H
