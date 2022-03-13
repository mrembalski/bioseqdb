//
// Created by user on 26/01/2022.
//

#include "mmseq.h"

void MMSeq2::execute() {

    puts("tSequencesPtr->precomputeIndexTable");
    tSequencesPtr->precomputeIndexTable();

    puts("fastKMerMatchStage");
    fastKMerMatchStage();

    puts("vectorizedGappedAlignment");
    vectorizedGappedAlignment();

    puts("vectorizedUngappedAlignment");

    vectorizedUngappedAlignment();
}

// loop 1
void MMSeq2::fastKMerMatchStage() {
    for (auto &v : *qSequences) {
        puts("fastKMerMatchStage");
        v.fastKMerMatchStage(tSequencesPtr);
        puts("fastKMerMatchStage");
    }
}

void MMSeq2::vectorizedUngappedAlignment() {
    for (auto &v : *qSequences) {
        v.vectorizedUngappedAlignment();
    }
}

void MMSeq2::vectorizedGappedAlignment() {
    for (auto &v : *qSequences) {
        v.vectorizedGappedAlignment();
    }
}

void MMSeq2::getResults() {
    //TODO: NOW just generic bullshit
    int qSeqId = 0;
    for (auto &q : *qSequences) {
        printf("Results for qSequence %d:\n", qSeqId++);

        for (int i = 0; i < q.getResultSize(); ++i) {
            auto res = q.getResult(i);

            auto tSequenceId = res.getTSequenceId();
            auto posQ = res.getPosQ();
            auto lastPosQ = res.getLastPosQ();
            auto posT = posQ + res.getDiagId();
            auto lastPosT = lastPosQ + res.getDiagId();

            printf("tSequenceId: %d, Qpositions: (%d %d) Tpositions: (%d %d)",
                   tSequenceId,
                   lastPosQ, posQ,
                   lastPosT, posT
                   );

            for (auto c = 0; c < kVal; ++c) {
                printf("%c", q.getSequence()->at(lastPosQ + c));
            }

            puts("");

            for (auto c = 0; c < kVal; ++c) {
                printf("%c", q.getSequence()->at(posQ + c));
            }

            puts("");

            for (auto c = 0; c < kVal; ++c) {
                printf("%c", tSequencesPtr->getSequence(tSequenceId)[posT + c]);
            }

            puts("");

            for (auto c = 0; c < kVal; ++c) {
                printf("%c", tSequencesPtr->getSequence(tSequenceId)[lastPosT + c]);
            }

            puts("");
        }
    }
}