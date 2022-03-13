//
// Created by user on 26/01/2022.
//

#ifndef MMSEQ2_TSEQUENCE_H
#define MMSEQ2_TSEQUENCE_H

#include <map>
#include <utility>
#include <vector>
#include <string>
#include <utility>
#include "mock_structures.h"
#include "blosum.h"

class TSequences {
public:
    TSequences() = delete;

    using Seq = std::string;

    using Seqs = std::vector<Seq>;

    using SeqsPtr = std::shared_ptr<Seqs>;

    using Num = int;

    using KMerPos = std::pair<Num, Num>;

    using KMerPosVec = std::vector<KMerPos>;

    using KMerPosVecPtr = std::shared_ptr<std::vector<std::pair<Num, Num>>>;

    using IndexTableType = std::map<Seq, KMerPosVecPtr>;

    TSequences(SeqsPtr sequences) {
        this->sequences = std::move(sequences);

        this->indexTable = IndexTableType();
    }

    inline KMerPosVecPtr getKMerPosVec(const Seq &seq) {
        if (indexTable.find(seq) != indexTable.end()) {
            return indexTable[seq];
        }
        else {
            return std::make_shared<std::vector<std::pair<Num, Num>>>();
        }
    }

    void precomputeIndexTable();

    inline Seq getSequence(int indx) const {
        return sequences->at(indx);
    }

private:
    SeqsPtr sequences;

    IndexTableType indexTable;

    void updateIndexTable(Seq &kMer, KMerPos kMerPos);
};
#endif //MMSEQ2_TSEQUENCE_H
