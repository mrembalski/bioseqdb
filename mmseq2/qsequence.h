//
// Created by user on 26/01/2022.
//

#ifndef MMSEQ2_QSEQUENCE_H
#define MMSEQ2_QSEQUENCE_H

#include <utility>
#include <vector>
#include <memory>
#include <string>
#include <utility>
#include "tsequence.h"

class QSequence {
public:
    QSequence() = delete;

    using SeqPtr = std::shared_ptr<std::string>;

    using SeqsPtr = std::shared_ptr<std::vector<SeqPtr>>;

    using TSequencesPtr = std::shared_ptr<TSequences>;

    using Num = int;

    explicit QSequence(SeqPtr sequence) { //TODO: change for a data from the database
        this->sequence = std::move(sequence);
    }

    void fastKMerMatchStage(const TSequencesPtr& tSequencesPtr);

    void addResult(int tSequenceId, int diagId, int posQ, int lastPosQ);

    void vectorizedUngappedAlignment();

    void vectorizedGappedAlignment();

    const Result getResult(int i);

    inline SeqPtr getSequence() const {
        return sequence;
    }

    int getResultSize() const;

private:
    SeqPtr sequence;

    std::vector<Num> diagonals;

    std::vector<Num> lastFound;

    const int EmptyDiag = -1e9;

    std::vector<Result> results;

    SeqPtr getFirstKMer();

    void resizeDiagonal(int k);
};
#endif //MMSEQ2_QSEQUENCE_H
