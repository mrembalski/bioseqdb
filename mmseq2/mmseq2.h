#ifndef BIOSEQDB_MMSEQ2_H
#define BIOSEQDB_MMSEQ2_H
#include <exception>
#include <utility>
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include "mock_structures.h"

namespace mmseq2 {
    class InputParams;
    class AminoAcid {
    public:
        static uint32_t charToId(char aa_id) {
            switch(aa_id) {
                case 'A':
                    return 0;
                case 'R':
                    return 1;
                case 'N':
                    return 2;
                case 'D':
                    return 3;
                case 'C':
                    return 4;
                case 'Q':
                    return 5;
                case 'E':
                    return 6;
                case 'G':
                    return 7;
                case 'H':
                    return 8;
                case 'I':
                    return 9;
                case 'L':
                    return 10;
                case 'K':
                    return 11;
                case 'M':
                    return 12;
                case 'F':
                    return 13;
                case 'P':
                    return 14;
                case 'S':
                    return 15;
                case 'T':
                    return 16;
                case 'W':
                    return 17;
                case 'Y':
                    return 18;
                case 'V':
                    return 19;
                case 'B':
                    return 20;
                case 'J':
                    return 21;
                case 'Z':
                    return 22;
                case 'X':
                    return 23;
                default: // case '*'
                    return 24;
            }
        }

        static char idToChar(uint32_t aa_char) {
            switch(aa_char) {
                case 0:
                    return 'A';
                case 1:
                    return 'R';
                case 2:
                    return 'N';
                case 3:
                    return 'D';
                case 4:
                    return 'C';
                case 5:
                    return 'Q';
                case 6:
                    return 'E';
                case 7:
                    return 'G';
                case 8:
                    return 'H';
                case 9:
                    return 'I';
                case 10:
                    return 'L';
                case 11:
                    return 'K';
                case 12:
                    return 'M';
                case 13:
                    return 'F';
                case 14:
                    return 'P';
                case 15:
                    return 'S';
                case 16:
                    return 'T';
                case 17:
                    return 'W';
                case 18:
                    return 'Y';
                case 19:
                    return 'V';
                case 20:
                    return 'B';
                case 21:
                    return 'J';
                case 22:
                    return 'Z';
                case 23:
                    return 'X';
                default: // case 24
                    return '*';
            }
        }

        static int32_t getPenalty(uint32_t matrixId, uint32_t currentId, uint32_t replacementId) {
            return blosum[matrixId][currentId][replacementId];
        }

        static int32_t getPenalty(uint32_t matrixId, char currentId, char replacementId) {
            return blosum[matrixId][charToId(currentId)][charToId(replacementId)];
        }

        static uint32_t getAlphabetSize() {
            return alphabetSize;
        }

        static uint32_t blosumIdToMatrixId(uint32_t blosumId) {
            switch(blosumId) {
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

        static constexpr int32_t blosum[5][25][25] = {
                {
                        {5, -2, -1, -2, -1, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -2, -2, 0, -1, -1, -1, -1, -5},
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
                        {-5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, 1}
                },
                {
                        {5, -2, -1, -2, -1, -1, -1, 0, -2, -1, -2, -1, -1, -3, -1, 1, 0, -3, -2, 0, -2, -2, -1, -1, -5},
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
                        {-5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, 1}
                },
                {
                        {4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, -1, -1, -4},
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
                        {-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1}
                },
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
                {
                        {5, -2, -2, -3, -1, -1, -1, 0, -2, -2, -2, -1, -2, -3, -1, 1, 0, -4, -3, -1, -2, -2, -1, -1, -6},
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
                        {-6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, 1}
                }
        };
    };

    class PrefilterKmerStageResults {
    private:
        uint64_t queryId;
        std::vector <uint32_t> targetIds;
        std::vector <int32_t> diagonals;

    public:
        explicit PrefilterKmerStageResults(uint64_t queryId) : queryId{ queryId }, targetIds{ std::vector<uint32_t> {}},
        diagonals{ std::vector<int32_t> {}}  {}

        PrefilterKmerStageResults() = delete;

        void addDiagonal(uint32_t targetId, int32_t diagonal) {
            targetIds.push_back(targetId);
            diagonals.push_back(diagonal);
        }

        [[nodiscard]] int32_t getDiagonal(int index) const {
            return diagonals[index];
        }

        [[nodiscard]] uint32_t getTargetId(int index) const {
            return targetIds[index];
        }

        [[nodiscard]] uint32_t getTargetsNumber() const {
            return targetIds.size();
        }
    };

    class InputParams { // Jeremi - feel free to add your own constructors
    public:
        using Vec64Ptr = std::shared_ptr<std::vector<uint64_t>>;
        using StrPtr = std::shared_ptr<std::string>;
        using VecStrPtr = std::shared_ptr<std::vector<StrPtr>>;
        using InputParamsPtr = std::shared_ptr<InputParams>;

        InputParams() = delete;

        InputParams(uint32_t qLen, uint32_t tLen, Vec64Ptr qIds, Vec64Ptr tIds, VecStrPtr queries,
                    StrPtr targetTableName, StrPtr targetColumnName, const StrPtr& substitutionMatrixName, uint32_t kMerLength,
                    int32_t kMerGenThreshold, int32_t ungappedAlignmentScore, int32_t evalTreshold, int32_t gapOpenCost,
                    int32_t gapPenaltyCost, uint32_t threadNumber);

        [[nodiscard]] uint32_t getQLen() const {
            return qLen;
        }

        [[nodiscard]] uint32_t getTLen() const {
            return tLen;
        }

        Vec64Ptr getQIds() {
            return qIds;
        }

        [[nodiscard]] Vec64Ptr getTIds() const {
            return tIds;
        }

        [[nodiscard]] VecStrPtr getQueries() const {
            return queries;
        }

        [[nodiscard]] StrPtr getTargetTableName() const {
            return targetTableName;
        }

        [[nodiscard]] StrPtr getTargetColumnName() const {
            return targetColumnName;
        }

        [[nodiscard]] uint32_t getKMerLength() const {
            return kMerLength;
        }

        [[nodiscard]] int32_t getKMerGenThreshold() const {
            return kMerGenThreshold;
        }

        [[nodiscard]] int32_t getUngappedAlignmentScore() const {
            return ungappedAlignmentScore;
        }

        [[nodiscard]] int32_t getEvalThreshold() const {
            return evalTreshold;
        }

        [[nodiscard]] int32_t getGapOpenCost() const {
            return gapOpenCost;
        }

        [[nodiscard]] int32_t getGapPenaltyCost() const {
            return gapPenaltyCost;
        }

        [[nodiscard]] uint32_t getThreadNumber() const {
            return threadNumber;
        }

        [[nodiscard]] uint32_t getSubstitutionMatrixId() const {
            return substitutionMatrixId;
        }

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
        int32_t evalTreshold;
        int32_t gapOpenCost;
        int32_t gapPenaltyCost;
        uint32_t threadNumber;
        uint32_t substitutionMatrixId;
    };

    struct MmseqResult {
        uint64_t queryId, targetId;
        double rawScore = 0.0, bitScore = 0.0, eValue = 0.0;
        uint32_t qStart = 0, qEnd = 0, qLen = 0, tStart = 0, tEnd = 0, tLen = 0;
        std::string qAln, tAln, cigar;
        uint32_t alnLen = 0, mismatch = 0, gapOpen = 0;
        double pident = 0.0;

        MmseqResult(uint64_t queryId, uint64_t targetId) : queryId{queryId}, targetId{targetId} {}
    };

    using VecRes = std::vector<mmseq2::MmseqResult>;
    using VecResPtr = std::shared_ptr<VecRes>;
    VecRes MMSeq2(mmseq2::InputParams::InputParamsPtr inputParams);

    class Query {
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
                                                                                 diagonalPrev{std::vector<int32_t>(inputParams.get()->getTLen(), 0)} { }

        [[nodiscard]] const PrefilterKmerStageResults& getPrefilterKmerStageResults() const {
            return prefilterKmerStageResults;
        }

        void addMatch(uint32_t targetId, int32_t diagonal) {
            prefilterKmerStageResults.addDiagonal(targetId, diagonal);
        }

        void findPrefilterKmerStageResults();

        void executeAlignment(std::mutex *resMtx, const VecResPtr& mmseqResult);

        [[nodiscard]] uint32_t getSubstitutionMatrixId() const {
            return substitutionMatrixId;
        }

        [[nodiscard]] uint32_t getKMerLength() const {
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
        int32_t evalTreshold;
        int32_t gapOpenCost;
        int32_t costGapExtended;

        PrefilterKmerStageResults prefilterKmerStageResults;

        std::vector<int32_t> diagonalPrev;

        std::vector<bool> diagonalPreVVisited;

        std::set<uint32_t> filteredTargetIds;

        void processSimilarKMers(uint32_t diagonalNumber, std::string &kMer, int32_t SMaxSuf,
                                 int32_t Spref = 0, uint32_t indx = 0);

        void processSingleKmer(uint32_t diagonal, std::string &kMer);

        [[nodiscard]] int32_t ungappedAlignment(const StrPtr& querySequence, const StrPtr& targetSequence, int32_t diagonal) const;

        void gappedAlignment(const StrPtr& querySequence, const StrPtr& targetSequence, MmseqResult& mmseqResult) const;
    };
}

#endif //BIOSEQDB_MMSEQ2_H
