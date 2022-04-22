set(BINARY ${CMAKE_PROJECT_NAME})

set(SOURCES
        main.cpp
        mmseq2.cpp
        mmseq2.h
        mock_structures.cpp
        mock_structures.h)

add_executable(${BINARY}_run ${SOURCES})

add_library(${BINARY}_lib STATIC ${SOURCES})