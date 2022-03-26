#include "postgres.h"

#include "access/htup_details.h"
#include "catalog/pg_type.h"
#include "executor/spi.h"
#include "funcapi.h"
#include "lib/stringinfo.h"
#include "miscadmin.h"
#include "executor/tablefunc.h"
#include "utils/builtins.h"

#include "mmseq2.h"

PG_MODULE_MAGIC;

PG_FUNCTION_INFO_V1(seq_search_mmseqs);

// 3 first arguments - columns from our input table containing genomic data
// 4th arguemnt - number of rows of the input table
// 5th argument - name of the target table for accessing indexes via SPI
// For now returns an empty table
Datum seq_search_mmseqs(PG_FUNCTION_ARGS)
{
    // Get all arguments
    int64 sequence_id = PG_GETARG_INT64(0);
    bool is_query = PG_GETARG_BOOL(1);
    char *sequence = text_to_cstring(PG_GETARG_TEXT_PP(2));
    int32 no_of_rows = PG_GETARG_INT32(3);
    char *target_name = text_to_cstring(PG_GETARG_TEXT_PP(4));

    // Get rsinfo and per_query_ctx
    ReturnSetInfo *rsinfo = (ReturnSetInfo *)fcinfo->resultinfo;
    MemoryContext per_query_ctx = rsinfo->econtext->ecxt_per_query_memory;
    MemoryContext oldcontext;

    // Tuple descriptor
    TupleDesc tupdesc;
    get_call_result_type(fcinfo, NULL, &tupdesc);
    AttInMetadata *attinmeta = TupleDescGetAttInMetadata(tupdesc);

    // Switch to per_query memory context
    oldcontext = MemoryContextSwitchTo(per_query_ctx);

    // Initialize the output table as empty
    tupdesc = CreateTupleDescCopy(tupdesc);
    Tuplestorestate *tupstore = tuplestore_begin_heap(
        rsinfo->allowedModes & SFRM_Materialize_Random,
        false,
        work_mem);

    // Allocate the memory for the input data
    const size_t row_size = sizeof(int64) + sizeof(bool) + sizeof(void *);
    if (fcinfo->flinfo->fn_extra == NULL)
        fcinfo->flinfo->fn_extra = palloc0(sizeof(int32) + sizeof(char *) * no_of_rows);
    char *memory = fcinfo->flinfo->fn_extra;

    // n - current row's number
    int32 n = (*((int32 *)memory))++;

    // Put the input data in memory
    char *row_ptr = memory + sizeof(int32) + n * sizeof(char *);
    *((int64 *)row_ptr) = sequence_id;
    *((bool *)(row_ptr + sizeof(int64))) = is_query;
    *((char **)row_ptr) = palloc(sizeof(char) * strlen(sequence) + 1);
    strcpy(*((char **)row_ptr), sequence);

    // Switch to the old memory context
    MemoryContextSwitchTo(oldcontext);

    // Only for the last row
    if (n + 1 == no_of_rows)
    {
        SPI_connect();

        uint32_t q_len = 2;
        uint32_t t_len = 1;
        uint64_t *q_ids = (uint64_t *)palloc(sizeof(uint64_t) * q_len);
        uint64_t *t_ids = (uint64_t *)palloc(sizeof(uint64_t) * t_len);
        char **queries = (char **)(memory + sizeof(int32));
        char *target_table_name = "taco";
        char *target_column_name = "sequence";
        t_ids[0] = 1;
        q_ids[0] = 1;
        q_ids[1] = 2;
        cpp_mmseq2(q_len, t_len, q_ids, t_ids, queries, target_table_name, target_column_name);

        SPI_finish();

        for (int i = 0; i < no_of_rows; i++)
        {
            char **values = (char **)palloc0(sizeof(char *));
            values[0] = *((char **)(memory + sizeof(int32) + i * sizeof(char *)));
            pfree(*((char **)(memory + sizeof(int32) + i * sizeof(char *))));

            HeapTuple tuple;
            tuple = BuildTupleFromCStrings(attinmeta, values);
            tuplestore_puttuple(tupstore, tuple);
            heap_freetuple(tuple);
            pfree(values);
        }

        pfree(memory); // For now only free the memory
    }

    // Set the return mode to materialize and return tupstore
    rsinfo->returnMode = SFRM_Materialize;
    rsinfo->setResult = tupstore;
    rsinfo->setDesc = tupdesc;
    return (Datum)0;
}
