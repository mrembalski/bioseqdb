#include "mmseq.h"
#include <iostream>

int main() {
    MMSeq2::QSequences q = std::make_shared<std::vector<QSequence>>();
    //MMSeq2::TSequencesPtr t;
    TSequences::SeqsPtr tPtr = std::make_shared<TSequences::Seqs>();

    for (auto &i : Q) {
        q->push_back((QSequence(std::make_shared<std::string>(i))));
    }

    for (auto &i : T) {
        tPtr->push_back(i);
    }

    MMSeq2::TSequencesPtr t = std::make_shared<TSequences>(tPtr);

    MMSeq2 mmSeq2 = MMSeq2(q, t);

    puts("to execution");
    mmSeq2.execute();

    puts("xdd");

    mmSeq2.getResults();

    return 0;
}
