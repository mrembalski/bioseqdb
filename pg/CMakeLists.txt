
cmake_minimum_required(VERSION 3.17)
project(mmseq2 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)

find_package(PostgreSQL REQUIRED)

add_subdirectory(rpclib)
include_directories(rpclib/include)

add_library(mmseq2 SHARED extension.cpp)

target_link_libraries(mmseq2 PRIVATE PostgreSQL::PostgreSQL)

find_program(PG_CONFIG pg_config)

execute_process(COMMAND ${PG_CONFIG} --pkglibdir OUTPUT_VARIABLE PG_CONFIG_PKGLIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${PG_CONFIG} --sharedir OUTPUT_VARIABLE PG_CONFIG_SHAREDIR OUTPUT_STRIP_TRAILING_WHITESPACE)

install(TARGETS mmseq2 DESTINATION ${PG_CONFIG_PKGLIBDIR})
install(FILES mmseq2.control DESTINATION ${PG_CONFIG_SHAREDIR}/extension)
install(FILES mmseq2--0.0.0.sql DESTINATION ${PG_CONFIG_SHAREDIR}/extension)
