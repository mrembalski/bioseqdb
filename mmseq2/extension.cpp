extern "C" {
#include "postgres.h"
#include "access/htup_details.h"
#include "catalog/pg_type.h"
#include "executor/spi.h"
#include "funcapi.h"
#include "lib/stringinfo.h"
#include "miscadmin.h"
#include "executor/tablefunc.h"
#include "utils/builtins.h"
}

#include "mmseq2.h"

extern "C" {
    PG_MODULE_MAGIC;

    PG_FUNCTION_INFO_V1(seq_search_mmseqs);

    Datum seq_search_mmseqs(PG_FUNCTION_ARGS)
    {
        // Get all arguments
        char  *query_tblname = text_to_cstring(PG_GETARG_TEXT_PP(0));
        char  *query_colname = text_to_cstring(PG_GETARG_TEXT_PP(1));
        char *target_tblname = text_to_cstring(PG_GETARG_TEXT_PP(2));
        char *target_colname = text_to_cstring(PG_GETARG_TEXT_PP(3));

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

        // Switch to the old memory context
        MemoryContextSwitchTo(oldcontext);

        // Initialize arguments for the C++ function
        /* uint32_t q_len = 2;
        uint32_t t_len = 1;
        uint64_t *q_ids = (uint64_t *)palloc(sizeof(uint64_t) * q_len);
        uint64_t *t_ids = (uint64_t *)palloc(sizeof(uint64_t) * t_len);
        char **queries = (char **)(memory + sizeof(int32));
        char *target_table_name = "taco";
        char *target_column_name = "sequence";
        t_ids[0] = 1;
        q_ids[0] = 1;
        q_ids[1] = 2;
        SPI_connect();
        cpp_mmseq2(q_len, t_len, q_ids, t_ids, queries, target_table_name, target_column_name);
        SPI_finish(); */

        // Build the output table
        for (int i = 0; i < 12; i++)
        {
            char **values = (char **)palloc0(sizeof(char *) * 4);
            values[0] = "asdf";
            values[1] = "asdf";
            values[2] = "asdf";
            values[3] = "asdf";

            HeapTuple tuple;
            tuple = BuildTupleFromCStrings(attinmeta, values);
            tuplestore_puttuple(tupstore, tuple);
            heap_freetuple(tuple);
            pfree(values);
        }

        // Set the return mode to materialize and return tupstore
        rsinfo->returnMode = SFRM_Materialize;
        rsinfo->setResult = tupstore;
        rsinfo->setDesc = tupdesc;
        return (Datum)0;
    }
}
