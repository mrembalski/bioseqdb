//
// Created by user on 26/01/2022.
//

#include "qsequence.h"
#include <string>
#include <utility>

QSequence::SeqPtr QSequence::getFirstKMer() {
    QSequence::SeqPtr kMer = std::make_shared<std::string>("");

    if (sequence->size() < kVal) {
        return std::move(kMer);
    }

    for (int i = 0;i < kVal; ++i) {
        kMer->push_back(sequence->at(i));
    }

    return std::move(kMer);
}

// Make diagonal size be at least k
void QSequence::resizeDiagonal(int k) {
    while (diagonals.size() <= k) {
        diagonals.push_back(EmptyDiag);
    }

    while (lastFound.size() <= k) {
        lastFound.push_back(EmptyDiag);
    }
}

void QSequence::addResult(int tSequenceId, int diagId, int posQ, int lastPosQ) {
    results.emplace_back(tSequenceId, diagId, posQ, lastPosQ);
}

void QSequence::fastKMerMatchStage(const TSequencesPtr& tSequencesPtr) {
    SeqPtr kMer = getFirstKMer();

    if (kMer->empty()) {
        return;
    }

    for (auto i = 0; i + kVal - 1 < sequence->size(); ++i) {
        // loop 2
        SeqsPtr simKMers = Blosum62SimilarKmers(kMer);

        // loop 3
        for (const auto& simKMer : *simKMers) {

            TSequences::KMerPosVecPtr tSequencesKMerPosVecPtr = tSequencesPtr->getKMerPosVec(*simKMer);

            if (tSequencesKMerPosVecPtr->empty()) {
                continue;
            }

            // loop 4
            for (auto tSequencesKMerPosPtr : *tSequencesKMerPosVecPtr) {
                int tSequenceId = tSequencesKMerPosPtr.first;
                int pos = tSequencesKMerPosPtr.second;

                resizeDiagonal(tSequenceId);

                // printf("kapibara %d %s\n", pos - i, simKMer->c_str());
                if (diagonals[tSequenceId] == pos - i) {
                    addResult(tSequenceId, pos - i,  i, lastFound[tSequenceId]);
                }

                lastFound[tSequenceId] = i;
                diagonals[tSequenceId] = pos - i;
            }
        }


        if (i + kVal < sequence->size()) {
            kMer->erase(0, 1);
            kMer->push_back(sequence->at(i + kVal));
        }
    }
}

void QSequence::vectorizedUngappedAlignment() {
    // TODO
}

void QSequence::vectorizedGappedAlignment() {
    // TODO
}

const Result QSequence::getResult(int i) {
    return results[i];
}

int QSequence::getResultSize() const {
    return results.size();
}