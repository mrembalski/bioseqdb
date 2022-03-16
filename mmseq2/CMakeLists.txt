cmake_minimum_required(VERSION 3.17)
project(MMSeq2)

set(CMAKE_CXX_STANDARD 20)


add_executable(MMSeq2 main.cpp mmseq2.h mmseq2.cpp mock_structures.h mock_structures.cpp)
install(TARGETS DESTINATION .)
