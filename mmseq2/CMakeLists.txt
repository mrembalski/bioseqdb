cmake_minimum_required(VERSION 3.17)
project(MMSeq2)

set(CMAKE_CXX_STANDARD 20)


add_executable(MMSeq2 main.cpp mmseq.cpp mmseq.h qsequence.cpp qsequence.h tsequence.cpp tsequence.h blosum.cpp blosum.h mock_structures.cpp mock_structures.h)
install(TARGETS DESTINATION .)

