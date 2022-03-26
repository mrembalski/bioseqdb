#include "postgres.h"

#include "access/htup_details.h"
#include "catalog/pg_type.h"
#include "executor/spi.h"
#include "funcapi.h"
#include "lib/stringinfo.h"
#include "miscadmin.h"
#include "executor/tablefunc.h"
#include "utils/builtins.h"

#include "mock_function.h"

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
    void *sequence = PG_GETARG_TEXT_PP(2);
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
        fcinfo->flinfo->fn_extra = palloc0(sizeof(int32) + (sizeof(int64) + no_of_rows * row_size));
    char *memory = fcinfo->flinfo->fn_extra;

    // n - current row's number
    int32 n = (*((int32 *)memory))++;

    // Put the input data in memory
    char *row_ptr = memory + sizeof(int32) + n * row_size;
    *((int64 *)row_ptr) = sequence_id;
    *((bool *)(row_ptr + sizeof(int64))) = is_query;
    *((void **)(row_ptr + sizeof(int64) + sizeof(bool))) = sequence;

    // Switch to the old memory context
    MemoryContextSwitchTo(oldcontext);

    // Only for the last row
    if (n + 1 == no_of_rows)
    {
        SPI_connect();
        mock_function(memory, target_name);
        SPI_finish();
        pfree(memory); // For now only free the memory
    }

    // Set the return mode to materialize and return tupstore
    rsinfo->returnMode = SFRM_Materialize;
    rsinfo->setResult = tupstore;
    rsinfo->setDesc = tupdesc;
    return (Datum)0;
}
