#include <string.h>

#include "postgres.h"
#include "executor/spi.h"

#include "psql_interface.h"

void log_from_cpp(char *str)
{
    elog(WARNING, "%s", str);
}

uint64_t get_indexes(char *column_name, char *kmer)
{
    char query[200];
    memset(query, 0, 200);
    sprintf(query, "%s%s%s%s%s",
            "SELECT * FROM ",
            column_name,
            "__index WHERE kmer=\'",
            kmer,
            "\';");
    SPI_exec(query, 0);
    return SPI_processed;
}

void get_ith_index(int i, uint64_t *target_id, uint32_t *position)
{
    TupleDesc spi_tupdesc = SPI_tuptable->tupdesc;
    int ncols = spi_tupdesc->natts;
    HeapTuple spi_tuple = SPI_tuptable->vals[i];
    char *col_one = SPI_getvalue(spi_tuple, spi_tupdesc, 2);
    char *col_two = SPI_getvalue(spi_tuple, spi_tupdesc, 3);

    *position = atoi(col_one);
    *target_id = atoi(col_two);
}
