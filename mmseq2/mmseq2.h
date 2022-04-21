#ifndef BIOSEQDB_MMSEQ2_H
#define BIOSEQDB_MMSEQ2_H

#include <stdint.h>
#include <memory>
#include <exception>
#include <utility>
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <mutex>
#include "mock_structures.h"
#include "rpc/server.h"

namespace mmseq2
{
    class InputParams;
    class AminoAcid
    {
    public:
        static uint32_t charToId(char aa_id)
        {
            for (uint32_t id = 0; id < alphabetSize; id++)
            {
                if (charId[id] == aa_id)
                {
                    return id;
                }
            }
            return alphabetSize - 1;
        }

        static char idToChar(uint32_t aa_char)
        {
            if (aa_char >= alphabetSize)
            {
                return '*';
            }
            return charId[aa_char];
        }

        static int32_t getPenalty(uint32_t matrixId, uint32_t currentId, uint32_t replacementId)
        {
            return blosum[matrixId][currentId][replacementId];
        }

        static int32_t getPenalty(uint32_t matrixId, char currentId, char replacementId)
        {
            return blosum[matrixId][charToId(currentId)][charToId(replacementId)];
        }

        static uint32_t getAlphabetSize()
        {
            return alphabetSize;
        }

        static uint32_t blosumIdToMatrixId(uint32_t blosumId)
        {
            switch (blosumId)
            {
            case 45:
                return 0;
            case 50:
                return 1;
            case 62:
                return 2;
            case 80:
                return 3;
            case 90:
                return 4;
            default:
                std::cout << "Wrong blosum number" << std::endl;
                exit(1);
            }
        }

    private:
        static constexpr uint32_t alphabetSize = 25;

        static constexpr char charId[25] = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'J', 'Z', 'X', '*'};

        static constexpr int32_t blosum[5][25][25] = {
            {{5, -2, -1, -2, -1, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -2, -2, 0, -1, -1, -1, -1, -5},
             {-2, 7, 0, -1, -3, 1, 0, -2, 0, -3, -2, 3, -1, -2, -2, -1, -1, -2, -1, -2, -1, -3, 1, -1, -5},
             {-1, 0, 6, 2, -2, 0, 0, 0, 1, -2, -3, 0, -2, -2, -2, 1, 0, -4, -2, -3, 5, -3, 0, -1, -5},
             {-2, -1, 2, 7, -3, 0, 2, -1, 0, -4, -3, 0, -3, -4, -1, 0, -1, -4, -2, -3, 6, -3, 1, -1, -5},
             {-1, -3, -2, -3, 12, -3, -3, -3, -3, -3, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -2, -2, -3, -1, -5},
             {-1, 1, 0, 0, -3, 6, 2, -2, 1, -2, -2, 1, 0, -4, -1, 0, -1, -2, -1, -3, 0, -2, 4, -1, -5},
             {-1, 0, 0, 2, -3, 2, 6, -2, 0, -3, -2, 1, -2, -3, 0, 0, -1, -3, -2, -3, 1, -3, 5, -1, -5},
             {0, -2, 0, -1, -3, -2, -2, 7, -2, -4, -3, -2, -2, -3, -2, 0, -2, -2, -3, -3, -1, -4, -2, -1, -5},
             {-2, 0, 1, 0, -3, 1, 0, -2, 10, -3, -2, -1, 0, -2, -2, -1, -2, -3, 2, -3, 0, -2, 0, -1, -5},
             {-1, -3, -2, -4, -3, -2, -3, -4, -3, 5, 2, -3, 2, 0, -2, -2, -1, -2, 0, 3, -3, 4, -3, -1, -5},
             {-1, -2, -3, -3, -2, -2, -2, -3, -2, 2, 5, -3, 2, 1, -3, -3, -1, -2, 0, 1, -3, 4, -2, -1, -5},
             {-1, 3, 0, 0, -3, 1, 1, -2, -1, -3, -3, 5, -1, -3, -1, -1, -1, -2, -1, -2, 0, -3, 1, -1, -5},
             {-1, -1, -2, -3, -2, 0, -2, -2, 0, 2, 2, -1, 6, 0, -2, -2, -1, -2, 0, 1, -2, 2, -1, -1, -5},
             {-2, -2, -2, -4, -2, -4, -3, -3, -2, 0, 1, -3, 0, 8, -3, -2, -1, 1, 3, 0, -3, 1, -3, -1, -5},
             {-1, -2, -2, -1, -4, -1, 0, -2, -2, -2, -3, -1, -2, -3, 9, -1, -1, -3, -3, -3, -2, -3, -1, -1, -5},
             {1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -3, -1, -2, -2, -1, 4, 2, -4, -2, -1, 0, -2, 0, -1, -5},
             {0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -1, -1, 2, 5, -3, -1, 0, 0, -1, -1, -1, -5},
             {-2, -2, -4, -4, -5, -2, -3, -2, -3, -2, -2, -2, -2, 1, -3, -4, -3, 15, 3, -3, -4, -2, -2, -1, -5},
             {-2, -1, -2, -2, -3, -1, -2, -3, 2, 0, 0, -1, 0, 3, -3, -2, -1, 3, 8, -1, -2, 0, -2, -1, -5},
             {0, -2, -3, -3, -1, -3, -3, -3, -3, 3, 1, -2, 1, 0, -3, -1, 0, -3, -1, 5, -3, 2, -3, -1, -5},
             {-1, -1, 5, 6, -2, 0, 1, -1, 0, -3, -3, 0, -2, -3, -2, 0, 0, -4, -2, -3, 5, -3, 1, -1, -5},
             {-1, -3, -3, -3, -2, -2, -3, -4, -2, 4, 4, -3, 2, 1, -3, -2, -1, -2, 0, 2, -3, 4, -2, -1, -5},
             {-1, 1, 0, 1, -3, 4, 5, -2, 0, -3, -2, 1, -1, -3, -1, 0, -1, -2, -2, -3, 1, -2, 5, -1, -5},
             {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5},
             {-5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, 1}},
            {{5, -2, -1, -2, -1, -1, -1, 0, -2, -1, -2, -1, -1, -3, -1, 1, 0, -3, -2, 0, -2, -2, -1, -1, -5},
             {-2, 7, -1, -2, -4, 1, 0, -3, 0, -4, -3, 3, -2, -3, -3, -1, -1, -3, -1, -3, -1, -3, 0, -1, -5},
             {-1, -1, 7, 2, -2, 0, 0, 0, 1, -3, -4, 0, -2, -4, -2, 1, 0, -4, -2, -3, 5, -4, 0, -1, -5},
             {-2, -2, 2, 8, -4, 0, 2, -1, -1, -4, -4, -1, -4, -5, -1, 0, -1, -5, -3, -4, 6, -4, 1, -1, -5},
             {-1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -3, -2, -3, -1, -5},
             {-1, 1, 0, 0, -3, 7, 2, -2, 1, -3, -2, 2, 0, -4, -1, 0, -1, -1, -1, -3, 0, -3, 4, -1, -5},
             {-1, 0, 0, 2, -3, 2, 6, -3, 0, -4, -3, 1, -2, -3, -1, -1, -1, -3, -2, -3, 1, -3, 5, -1, -5},
             {0, -3, 0, -1, -3, -2, -3, 8, -2, -4, -4, -2, -3, -4, -2, 0, -2, -3, -3, -4, -1, -4, -2, -1, -5},
             {-2, 0, 1, -1, -3, 1, 0, -2, 10, -4, -3, 0, -1, -1, -2, -1, -2, -3, 2, -4, 0, -3, 0, -1, -5},
             {-1, -4, -3, -4, -2, -3, -4, -4, -4, 5, 2, -3, 2, 0, -3, -3, -1, -3, -1, 4, -4, 4, -3, -1, -5},
             {-2, -3, -4, -4, -2, -2, -3, -4, -3, 2, 5, -3, 3, 1, -4, -3, -1, -2, -1, 1, -4, 4, -3, -1, -5},
             {-1, 3, 0, -1, -3, 2, 1, -2, 0, -3, -3, 6, -2, -4, -1, 0, -1, -3, -2, -3, 0, -3, 1, -1, -5},
             {-1, -2, -2, -4, -2, 0, -2, -3, -1, 2, 3, -2, 7, 0, -3, -2, -1, -1, 0, 1, -3, 2, -1, -1, -5},
             {-3, -3, -4, -5, -2, -4, -3, -4, -1, 0, 1, -4, 0, 8, -4, -3, -2, 1, 4, -1, -4, 1, -4, -1, -5},
             {-1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3, -2, -3, -1, -1, -5},
             {1, -1, 1, 0, -1, 0, -1, 0, -1, -3, -3, 0, -2, -3, -1, 5, 2, -4, -2, -2, 0, -3, 0, -1, -5},
             {0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 2, 5, -3, -2, 0, 0, -1, -1, -1, -5},
             {-3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1, 1, -4, -4, -3, 15, 2, -3, -5, -2, -2, -1, -5},
             {-2, -1, -2, -3, -3, -1, -2, -3, 2, -1, -1, -2, 0, 4, -3, -2, -2, 2, 8, -1, -3, -1, -2, -1, -5},
             {0, -3, -3, -4, -1, -3, -3, -4, -4, 4, 1, -3, 1, -1, -3, -2, 0, -3, -1, 5, -3, 2, -3, -1, -5},
             {-2, -1, 5, 6, -3, 0, 1, -1, 0, -4, -4, 0, -3, -4, -2, 0, 0, -5, -3, -3, 6, -4, 1, -1, -5},
             {-2, -3, -4, -4, -2, -3, -3, -4, -3, 4, 4, -3, 2, 1, -3, -3, -1, -2, -1, 2, -4, 4, -3, -1, -5},
             {-1, 0, 0, 1, -3, 4, 5, -2, 0, -3, -3, 1, -1, -4, -1, 0, -1, -2, -2, -3, 1, -3, 5, -1, -5},
             {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5},
             {-5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, 1}},
            {{4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, -1, -1, -4},
             {-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, -2, 0, -1, -4},
             {-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 4, -3, 0, -1, -4},
             {-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, -3, 1, -1, -4},
             {0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -1, -3, -1, -4},
             {-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, -2, 4, -1, -4},
             {-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, -3, 4, -1, -4},
             {0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -4, -2, -1, -4},
             {-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, -3, 0, -1, -4},
             {-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, 3, -3, -1, -4},
             {-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, 3, -3, -1, -4},
             {-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, -3, 1, -1, -4},
             {-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, 2, -1, -1, -4},
             {-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, 0, -3, -1, -4},
             {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -3, -1, -1, -4},
             {1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, -2, 0, -1, -4},
             {0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, -1, -1, -4},
             {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -2, -2, -1, -4},
             {-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -1, -2, -1, -4},
             {0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, 2, -2, -1, -4},
             {-2, -1, 4, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, -3, 0, -1, -4},
             {-1, -2, -3, -3, -1, -2, -3, -4, -3, 3, 3, -3, 2, 0, -3, -2, -1, -2, -1, 2, -3, 3, -3, -1, -4},
             {-1, 0, 0, 1, -3, 4, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -2, -2, -2, 0, -3, 4, -1, -4},
             {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -4},
             {-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1}},
            {
                {5, -2, -2, -2, -1, -1, -1, 0, -2, -2, -2, -1, -1, -3, -1, 1, 0, -3, -2, 0, -2, -2, -1, -1, -6},
                {-2, 6, -1, -2, -4, 1, -1, -3, 0, -3, -3, 2, -2, -4, -2, -1, -1, -4, -3, -3, -1, -3, 0, -1, -6},
                {-2, -1, 6, 1, -3, 0, -1, -1, 0, -4, -4, 0, -3, -4, -3, 0, 0, -4, -3, -4, 5, -4, 0, -1, -6},
                {-2, -2, 1, 6, -4, -1, 1, -2, -2, -4, -5, -1, -4, -4, -2, -1, -1, -6, -4, -4, 5, -5, 1, -1, -6},
                {-1, -4, -3, -4, 9, -4, -5, -4, -4, -2, -2, -4, -2, -3, -4, -2, -1, -3, -3, -1, -4, -2, -4, -1,
                 -6},
                {-1, 1, 0, -1, -4, 6, 2, -2, 1, -3, -3, 1, 0, -4, -2, 0, -1, -3, -2, -3, 0, -3, 4, -1, -6},
                {-1, -1, -1, 1, -5, 2, 6, -3, 0, -4, -4, 1, -2, -4, -2, 0, -1, -4, -3, -3, 1, -4, 5, -1, -6},
                {0, -3, -1, -2, -4, -2, -3, 6, -3, -5, -4, -2, -4, -4, -3, -1, -2, -4, -4, -4, -1, -5, -3, -1,
                 -6},
                {-2, 0, 0, -2, -4, 1, 0, -3, 8, -4, -3, -1, -2, -2, -3, -1, -2, -3, 2, -4, -1, -4, 0, -1, -6},
                {-2, -3, -4, -4, -2, -3, -4, -5, -4, 5, 1, -3, 1, -1, -4, -3, -1, -3, -2, 3, -4, 3, -4, -1, -6},
                {-2, -3, -4, -5, -2, -3, -4, -4, -3, 1, 4, -3, 2, 0, -3, -3, -2, -2, -2, 1, -4, 3, -3, -1, -6},
                {-1, 2, 0, -1, -4, 1, 1, -2, -1, -3, -3, 5, -2, -4, -1, -1, -1, -4, -3, -3, -1, -3, 1, -1, -6},
                {-1, -2, -3, -4, -2, 0, -2, -4, -2, 1, 2, -2, 6, 0, -3, -2, -1, -2, -2, 1, -3, 2, -1, -1, -6},
                {-3, -4, -4, -4, -3, -4, -4, -4, -2, -1, 0, -4, 0, 6, -4, -3, -2, 0, 3, -1, -4, 0, -4, -1, -6},
                {-1, -2, -3, -2, -4, -2, -2, -3, -3, -4, -3, -1, -3, -4, 8, -1, -2, -5, -4, -3, -2, -4, -2, -1,
                 -6},
                {1, -1, 0, -1, -2, 0, 0, -1, -1, -3, -3, -1, -2, -3, -1, 5, 1, -4, -2, -2, 0, -3, 0, -1, -6},
                {0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -2, -1, -1, -2, -2, 1, 5, -4, -2, 0, -1, -1, -1, -1, -6},
                {-3, -4, -4, -6, -3, -3, -4, -4, -3, -3, -2, -4, -2, 0, -5, -4, -4, 11, 2, -3, -5, -3, -3, -1,
                 -6},
                {-2, -3, -3, -4, -3, -2, -3, -4, 2, -2, -2, -3, -2, 3, -4, -2, -2, 2, 7, -2, -3, -2, -3, -1,
                 -6},
                {0, -3, -4, -4, -1, -3, -3, -4, -4, 3, 1, -3, 1, -1, -3, -2, 0, -3, -2, 4, -4, 2, -3, -1, -6},
                {-2, -1, 5, 5, -4, 0, 1, -1, -1, -4, -4, -1, -3, -4, -2, 0, -1, -5, -3, -4, 5, -4, 0, -1, -6},
                {-2, -3, -4, -5, -2, -3, -4, -5, -4, 3, 3, -3, 2, 0, -4, -3, -1, -3, -2, 2, -4, 3, -3, -1, -6},
                {-1, 0, 0, 1, -4, 4, 5, -3, 0, -4, -3, 1, -1, -4, -2, 0, -1, -3, -3, -3, 0, -3, 5, -1, -6},
                {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                 -6},
                {-6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6,
                 1},
            },
            {{5, -2, -2, -3, -1, -1, -1, 0, -2, -2, -2, -1, -2, -3, -1, 1, 0, -4, -3, -1, -2, -2, -1, -1, -6},
             {-2, 6, -1, -3, -5, 1, -1, -3, 0, -4, -3, 2, -2, -4, -3, -1, -2, -4, -3, -3, -2, -3, 0, -1, -6},
             {-2, -1, 7, 1, -4, 0, -1, -1, 0, -4, -4, 0, -3, -4, -3, 0, 0, -5, -3, -4, 5, -4, -1, -1, -6},
             {-3, -3, 1, 7, -5, -1, 1, -2, -2, -5, -5, -1, -4, -5, -3, -1, -2, -6, -4, -5, 5, -5, 1, -1, -6},
             {-1, -5, -4, -5, 9, -4, -6, -4, -5, -2, -2, -4, -2, -3, -4, -2, -2, -4, -4, -2, -4, -2, -5, -1, -6},
             {-1, 1, 0, -1, -4, 7, 2, -3, 1, -4, -3, 1, 0, -4, -2, -1, -1, -3, -3, -3, -1, -3, 5, -1, -6},
             {-1, -1, -1, 1, -6, 2, 6, -3, -1, -4, -4, 0, -3, -5, -2, -1, -1, -5, -4, -3, 1, -4, 5, -1, -6},
             {0, -3, -1, -2, -4, -3, -3, 6, -3, -5, -5, -2, -4, -5, -3, -1, -3, -4, -5, -5, -2, -5, -3, -1, -6},
             {-2, 0, 0, -2, -5, 1, -1, -3, 8, -4, -4, -1, -3, -2, -3, -2, -2, -3, 1, -4, -1, -4, 0, -1, -6},
             {-2, -4, -4, -5, -2, -4, -4, -5, -4, 5, 1, -4, 1, -1, -4, -3, -1, -4, -2, 3, -5, 3, -4, -1, -6},
             {-2, -3, -4, -5, -2, -3, -4, -5, -4, 1, 5, -3, 2, 0, -4, -3, -2, -3, -2, 0, -5, 4, -4, -1, -6},
             {-1, 2, 0, -1, -4, 1, 0, -2, -1, -4, -3, 6, -2, -4, -2, -1, -1, -5, -3, -3, -1, -3, 1, -1, -6},
             {-2, -2, -3, -4, -2, 0, -3, -4, -3, 1, 2, -2, 7, -1, -3, -2, -1, -2, -2, 0, -4, 2, -2, -1, -6},
             {-3, -4, -4, -5, -3, -4, -5, -5, -2, -1, 0, -4, -1, 7, -4, -3, -3, 0, 3, -2, -4, 0, -4, -1, -6},
             {-1, -3, -3, -3, -4, -2, -2, -3, -3, -4, -4, -2, -3, -4, 8, -2, -2, -5, -4, -3, -3, -4, -2, -1, -6},
             {1, -1, 0, -1, -2, -1, -1, -1, -2, -3, -3, -1, -2, -3, -2, 5, 1, -4, -3, -2, 0, -3, -1, -1, -6},
             {0, -2, 0, -2, -2, -1, -1, -3, -2, -1, -2, -1, -1, -3, -2, 1, 6, -4, -2, -1, -1, -2, -1, -1, -6},
             {-4, -4, -5, -6, -4, -3, -5, -4, -3, -4, -3, -5, -2, 0, -5, -4, -4, 11, 2, -3, -6, -3, -4, -1, -6},
             {-3, -3, -3, -4, -4, -3, -4, -5, 1, -2, -2, -3, -2, 3, -4, -3, -2, 2, 8, -3, -4, -2, -3, -1, -6},
             {-1, -3, -4, -5, -2, -3, -3, -5, -4, 3, 0, -3, 0, -2, -3, -2, -1, -3, -3, 5, -4, 1, -3, -1, -6},
             {-2, -2, 5, 5, -4, -1, 1, -2, -1, -5, -5, -1, -4, -4, -3, 0, -1, -6, -4, -4, 5, -5, 0, -1, -6},
             {-2, -3, -4, -5, -2, -3, -4, -5, -4, 3, 4, -3, 2, 0, -4, -3, -2, -3, -2, 1, -5, 4, -4, -1, -6},
             {-1, 0, -1, 1, -5, 5, 5, -3, 0, -4, -4, 1, -2, -4, -2, -1, -1, -4, -3, -3, 0, -4, 5, -1, -6},
             {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -6},
             {-6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, 1}}};
    };

    class PrefilterKmerStageResults
    {
    private:
        uint64_t queryId;
        std::vector<uint32_t> targetIds;
        std::vector<int32_t> diagonals;

    public:
        explicit PrefilterKmerStageResults(uint64_t queryId) : queryId{queryId}, targetIds{std::vector<uint32_t>{}},
                                                               diagonals{std::vector<int32_t>{}} {}

        PrefilterKmerStageResults() = delete;

        void addDiagonal(uint32_t targetId, int32_t diagonal)
        {
            targetIds.push_back(targetId);
            diagonals.push_back(diagonal);
        }

        [[nodiscard]] int32_t getDiagonal(int index) const
        {
            return diagonals[index];
        }

        [[nodiscard]] uint32_t getTargetId(int index) const
        {
            return targetIds[index];
        }

        [[nodiscard]] uint32_t getTargetsNumber() const
        {
            return targetIds.size();
        }
    };

    class InputParams
    { // Jeremi - feel free to add your own constructors
    public:
        using Vec64Ptr = std::shared_ptr<std::vector<uint64_t>>;
        using StrPtr = std::shared_ptr<std::string>;
        using VecStrPtr = std::shared_ptr<std::vector<StrPtr>>;
        using InputParamsPtr = std::shared_ptr<InputParams>;

        InputParams() {}

        InputParams(uint32_t qLen, uint32_t tLen, Vec64Ptr qIds, Vec64Ptr tIds, VecStrPtr queries,
                    StrPtr targetTableName, StrPtr targetColumnName, const StrPtr &substitutionMatrixName, uint32_t kMerLength,
                    int32_t kMerGenThreshold, int32_t ungappedAlignmentScore, double evalTreshold, int32_t gapOpenCost,
                    int32_t gapPenaltyCost, uint32_t threadNumber);

        [[nodiscard]] uint32_t getQLen() const
        {
            return qLen;
        }

        [[nodiscard]] uint32_t getTLen() const
        {
            return tLen;
        }

        Vec64Ptr getQIds()
        {
            return qIds;
        }

        [[nodiscard]] Vec64Ptr getTIds() const
        {
            return tIds;
        }

        [[nodiscard]] VecStrPtr getQueries() const
        {
            return queries;
        }

        [[nodiscard]] StrPtr getTargetTableName() const
        {
            return targetTableName;
        }

        [[nodiscard]] StrPtr getTargetColumnName() const
        {
            return targetColumnName;
        }

        [[nodiscard]] uint32_t getKMerLength() const
        {
            return kMerLength;
        }

        [[nodiscard]] int32_t getKMerGenThreshold() const
        {
            return kMerGenThreshold;
        }

        [[nodiscard]] int32_t getUngappedAlignmentScore() const
        {
            return ungappedAlignmentScore;
        }

        [[nodiscard]] double getEvalThreshold() const
        {
            return evalTreshold;
        }

        [[nodiscard]] int32_t getGapOpenCost() const
        {
            return gapOpenCost;
        }

        [[nodiscard]] int32_t getGapPenaltyCost() const
        {
            return gapPenaltyCost;
        }

        [[nodiscard]] uint32_t getThreadNumber() const
        {
            return threadNumber;
        }

        [[nodiscard]] uint32_t getSubstitutionMatrixId() const
        {
            return substitutionMatrixId;
        }

        MSGPACK_DEFINE_MAP(qLen, tLen, qIds, tIds,
                           queries, targetTableName, targetColumnName, substitutionMatrixName,
                           kMerLength, kMerGenThreshold, ungappedAlignmentScore, evalTreshold, gapOpenCost,
                           gapPenaltyCost, threadNumber, substitutionMatrixId);

    private:
        uint32_t qLen;
        uint32_t tLen;
        Vec64Ptr qIds;
        Vec64Ptr tIds;

        VecStrPtr queries;
        StrPtr targetTableName;
        StrPtr targetColumnName;
        StrPtr substitutionMatrixName;

        uint32_t kMerLength;
        int32_t kMerGenThreshold;
        int32_t ungappedAlignmentScore;
        double evalTreshold;
        int32_t gapOpenCost;
        int32_t gapPenaltyCost;
        uint32_t threadNumber;
        uint32_t substitutionMatrixId;
    };

    class MmseqResult
    {
    private:
        uint64_t queryId, targetId;
        double rawScore = 0.0, bitScore = 0.0, eValue = 0.0;
        uint32_t qStart = 0, qEnd = 0, qLen = 0, tStart = 0, tEnd = 0, tLen = 0;
        std::string qAln, tAln, cigar;
        uint32_t alnLen = 0, mismatch = 0, gapOpen = 0;
        double pident = 0.0;

    public:
        MmseqResult(uint64_t queryId, uint64_t targetId) : queryId{queryId}, targetId{targetId} {}

        uint64_t getQueryId() const
        {
            return queryId;
        }

        uint64_t getTargetId() const
        {
            return targetId;
        }

        double getRawScore() const
        {
            return rawScore;
        }

        double getBitScore() const
        {
            return bitScore;
        }

        double getEValue() const
        {
            return eValue;
        }

        uint32_t getQStart() const
        {
            return qStart;
        }

        uint32_t getQEnd() const
        {
            return qEnd;
        }

        uint32_t getQLen() const
        {
            return qLen;
        }

        uint32_t getTStart() const
        {
            return tStart;
        }

        uint32_t getTEnd() const
        {
            return tEnd;
        }

        uint32_t getTLen() const
        {
            return tLen;
        }

        std::string getQAln() const
        {
            return qAln;
        }

        std::string getTAln() const
        {
            return tAln;
        }

        std::string getCigar() const
        {
            return cigar;
        }

        uint32_t getAlnLen() const
        {
            return alnLen;
        }

        uint32_t getMismatch() const
        {
            return mismatch;
        }

        uint32_t getGapOpen() const
        {
            return gapOpen;
        }

        double getPident() const
        {
            return pident;
        }

        void setQLen(uint32_t qLen)
        {
            MmseqResult::qLen = qLen;
        }

        void setTLen(uint32_t tLen)
        {
            MmseqResult::tLen = tLen;
        }

        void setQEnd(uint32_t qEnd)
        {
            this->qEnd = qEnd;
        }

        void setTEnd(uint32_t tEnd)
        {
            this->tEnd = tEnd;
        }

        void incrMismatch()
        {
            this->mismatch++;
        }

        void incrGapOpen()
        {
            this->gapOpen++;
        }

        void setQStart(uint32_t qStart)
        {
            this->qStart = qStart;
        }

        void setRawScore(double rawScore)
        {
            MmseqResult::rawScore = rawScore;
        }

        void setBitScore(double bitScore)
        {
            MmseqResult::bitScore = bitScore;
        }

        void setEValue(double eValue)
        {
            MmseqResult::eValue = eValue;
        }

        void setTStart(uint32_t tStart)
        {
            MmseqResult::tStart = tStart;
        }

        void setQAln(const std::string &qAln)
        {
            MmseqResult::qAln = qAln;
        }

        void setTAln(const std::string &tAln)
        {
            MmseqResult::tAln = tAln;
        }

        void setCigar(const std::string &cigar)
        {
            MmseqResult::cigar = cigar;
        }

        void setAlnLen(uint32_t alnLen)
        {
            MmseqResult::alnLen = alnLen;
        }

        void setPident(double pident)
        {
            MmseqResult::pident = pident;
        }

        MSGPACK_DEFINE_MAP(queryId, targetId,
                           rawScore, bitScore, eValue, qStart, qEnd, qLen, tStart, tEnd, tLen,
                           qAln, tAln, cigar,
                           alnLen, mismatch, gapOpen,
                           pident);
    };

    using VecRes = std::vector<mmseq2::MmseqResult>;
    using VecResPtr = std::shared_ptr<VecRes>;
    VecRes MMSeq2(mmseq2::InputParams::InputParamsPtr &inputParams);

    class Query
    {
    public:
        using StrPtr = mmseq2::InputParams::StrPtr;
        using VecStrPtr = mmseq2::InputParams::VecStrPtr;

        Query() = delete;

        Query(uint64_t queryId, StrPtr query, mmseq2::InputParams::InputParamsPtr inputParams) : queryId{queryId}, sequence{std::move(query)},
                                                                                                 kMerLength{inputParams.get()->getKMerLength()},
                                                                                                 substitutionMatrixId{inputParams.get()->getSubstitutionMatrixId()},
                                                                                                 kMerGenThreshold{inputParams.get()->getKMerGenThreshold()},
                                                                                                 ungappedAlignmentScore{inputParams.get()->getUngappedAlignmentScore()},
                                                                                                 targetColumnName{inputParams.get()->getTargetColumnName()},
                                                                                                 evalTreshold{inputParams.get()->getEvalThreshold()},
                                                                                                 gapOpenCost{inputParams.get()->getGapOpenCost()},
                                                                                                 costGapExtended{inputParams.get()->getGapPenaltyCost()},
                                                                                                 targetTableName{inputParams.get()->getTargetTableName()},
                                                                                                 prefilterKmerStageResults{queryId},
                                                                                                 diagonalPreVVisited{std::vector<bool>(inputParams.get()->getTLen(), false)},
                                                                                                 diagonalPrev{std::vector<int32_t>(inputParams.get()->getTLen(), 0)} {}

        [[nodiscard]] const PrefilterKmerStageResults &getPrefilterKmerStageResults() const
        {
            return prefilterKmerStageResults;
        }

        void addMatch(uint32_t targetId, int32_t diagonal)
        {
            prefilterKmerStageResults.addDiagonal(targetId, diagonal);
        }

        void findPrefilterKmerStageResults();

        void executeAlignment(std::mutex *resMtx, const VecResPtr &mmseqResult);

        [[nodiscard]] uint32_t getSubstitutionMatrixId() const
        {
            return substitutionMatrixId;
        }

        [[nodiscard]] uint32_t getKMerLength() const
        {
            return kMerLength;
        }

    private:
        uint64_t queryId;
        StrPtr sequence;
        StrPtr targetTableName;

        uint32_t kMerLength;
        uint32_t substitutionMatrixId;
        int32_t kMerGenThreshold;
        int32_t ungappedAlignmentScore;
        StrPtr targetColumnName;
        double evalTreshold;
        int32_t gapOpenCost;
        int32_t costGapExtended;

        PrefilterKmerStageResults prefilterKmerStageResults;

        std::vector<int32_t> diagonalPrev;

        std::vector<bool> diagonalPreVVisited;

        std::set<uint32_t> filteredTargetIds;

        void processSimilarKMers(uint32_t diagonalNumber, std::string &kMer, int32_t SMaxSuf,
                                 int32_t Spref = 0, uint32_t indx = 0);

        void processSingleKmer(uint32_t diagonal, std::string &kMer);

        [[nodiscard]] double ungappedAlignment(const StrPtr &querySequence, const StrPtr &targetSequence, int32_t diagonal) const;

        void gappedAlignment(const StrPtr &querySequence, const StrPtr &targetSequence, MmseqResult &mmseqResult) const;
    };
}

#endif // BIOSEQDB_MMSEQ2_H
