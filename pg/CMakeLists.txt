
cmake_minimum_required(VERSION 3.17)
project(mmseq2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

add_library(mmseq2 SHARED 
    extension.cpp 
    ../mmseq2/mmseq2.h 
    ../mmseq2/mmseq2.cpp 
    ../mmseq2/mock_structures.h 
    ../mmseq2/mock_structures.cpp)

# Add Postgres support
find_package(PostgreSQL REQUIRED)
target_link_libraries(mmseq2 PRIVATE PostgreSQL::PostgreSQL)

find_program(PG_CONFIG pg_config)
execute_process(COMMAND ${PG_CONFIG} --pkglibdir OUTPUT_VARIABLE PG_CONFIG_PKGLIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${PG_CONFIG} --sharedir OUTPUT_VARIABLE PG_CONFIG_SHAREDIR OUTPUT_STRIP_TRAILING_WHITESPACE)

# Add RPC support
add_subdirectory(rpclib)
include_directories(rpclib/include)
target_link_libraries(mmseq2 PRIVATE rpc)

# Actually install files
install(TARGETS mmseq2 DESTINATION ${PG_CONFIG_PKGLIBDIR})
install(FILES mmseq2.control DESTINATION ${PG_CONFIG_SHAREDIR}/extension)
install(FILES mmseq2--0.0.0.sql DESTINATION ${PG_CONFIG_SHAREDIR}/extension)
