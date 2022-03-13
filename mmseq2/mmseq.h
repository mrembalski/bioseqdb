//
// Created by user on 26/01/2022.
//

#ifndef MMSEQ2_MMSEQ_H
#define MMSEQ2_MMSEQ_H

#include "qsequence.h"

#include <utility>
#include <memory>

class MMSeq2 {
public:

    using QSequences = std::shared_ptr<std::vector<QSequence>>;

    using TSequencesPtr = QSequence::TSequencesPtr;

    MMSeq2() = delete;

    MMSeq2(QSequences qSequences, TSequencesPtr tSequencesPtr) : qSequences{std::move(qSequences)}, tSequencesPtr{std::move(tSequencesPtr)} { }

    void fastKMerMatchStage(); //TODO: Get these two dickheads you already work with and do the matchin

    void vectorizedUngappedAlignment();

    void vectorizedGappedAlignment(); //SW

    void execute();

    void getResults(); //TODO: just write down debug info

    private:
        QSequences qSequences;
        TSequencesPtr tSequencesPtr;
};

#endif //MMSEQ2_MMSEQ_H
