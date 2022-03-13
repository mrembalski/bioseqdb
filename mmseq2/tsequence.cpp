//
// Created by user on 26/01/2022.
//

#include <string>
#include "tsequence.h"

TSequences::Seq getFirstKVal(TSequences::Seq &seq) {
    TSequences::Seq kMer;

    if (seq.size() < kVal) {
        return std::move(kMer);
    }

    for (int i = 0;i < kVal; ++i) {
        kMer.push_back(seq[i]);
    }

    return std::move(kMer);
}

void TSequences::updateIndexTable(Seq &kMer, KMerPos kMerPos) {
    auto kMerPosVecPtr = indexTable.find(kMer);

    if (kMerPosVecPtr == indexTable.end()) {
        indexTable[kMer] = std::make_shared<KMerPosVec>();
    }

    indexTable[kMer]->push_back(kMerPos);
}

void TSequences::precomputeIndexTable() { // This is complete mock function. It will be implemented differently later
    auto seqId = -1;

    for (auto &sequence : *sequences) {
        ++seqId;

        Seq kMer = getFirstKVal(sequence);

        if (kMer.empty()) {
            return;
        }

        updateIndexTable(kMer, std::move(KMerPos(seqId, 0)));

        int nextLastPos = kVal;

        while (nextLastPos < sequence.size()) {
            kMer.erase(0, 1);

            kMer.push_back(sequence[nextLastPos]);

            nextLastPos++;

            updateIndexTable(kMer, std::move(KMerPos(seqId, nextLastPos - kVal)));
        }
    }
}