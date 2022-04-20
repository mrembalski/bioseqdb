
# KOMENTARZ OD MATEUSZA:
# W PONIZSZYCH KOMENTARZACH STARY MMSEQ2 DLA ARCHITEKTURY LACZONEGO POSTGRES Z MMSEQ2
# NOWY CMAKE ODPOWIADA NOWEJ ARCHITEKTURZE TEGO SERWISU

cmake_minimum_required(VERSION 3.17)
project(mmseq2)

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(rpclib)
include_directories(rpclib/include)

add_executable(mmseq2 main.cpp mmseq2.h mmseq2.cpp mock_structures.h mock_structures.cpp)

target_link_libraries(mmseq2 PRIVATE rpc)

#cmake_minimum_required(VERSION 3.17)
#project(mmseq2 LANGUAGES CXX)

#set(CMAKE_CXX_STANDARD 20)
#set(THREADS_PREFER_PTHREAD_FLAG ON)

##find_package(PostgreSQL REQUIRED)
#find_package(Threads REQUIRED)

# add_library(mmseq2 SHARED main.cpp mmseq2.h mmseq2.cpp mock_structures.h mock_structures.cpp)

#target_link_libraries(mmseq2 PRIVATE Threads::Threads)
# target_link_libraries(mmseq2 PRIVATE PostgreSQL::PostgreSQL)

#find_program(PG_CONFIG pg_config)
#execute_process(COMMAND ${PG_CONFIG} --pkglibdir OUTPUT_VARIABLE PG_CONFIG_PKGLIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
#execute_process(COMMAND ${PG_CONFIG} --sharedir OUTPUT_VARIABLE PG_CONFIG_SHAREDIR OUTPUT_STRIP_TRAILING_WHITESPACE)

#install(TARGETS mmseq2 DESTINATION ${PG_CONFIG_PKGLIBDIR})
#install(FILES mmseq2.control DESTINATION ${PG_CONFIG_SHAREDIR}/extension)
#install(FILES mmseq2--0.0.0.sql DESTINATION ${PG_CONFIG_SHAREDIR}/extension)
