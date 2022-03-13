//
// Created by user on 27/01/2022.
//

#ifndef MMSEQ2_BLOSUM_H
#define MMSEQ2_BLOSUM_H

#include <string>
#include <vector>

using BlosumSeqPtr = std::shared_ptr<std::string>;
using BlosumSeqsPtr = std::shared_ptr<std::vector<BlosumSeqPtr>>;

inline BlosumSeqsPtr Blosum62SimilarKmers(const BlosumSeqPtr& seq) {
    BlosumSeqsPtr res = std::make_shared<std::vector<BlosumSeqPtr>>();

    res->push_back(seq);

    return res;
}

#endif //MMSEQ2_BLOSUM_H
