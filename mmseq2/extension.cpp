#include <vector>
#include <string>

extern "C"
{
#include <postgres.h>
#include <fmgr.h>
#include <access/htup_details.h>
#include <catalog/pg_type.h>
#include <executor/spi.h>
#include <funcapi.h>
#include <lib/stringinfo.h>
#include <miscadmin.h>
#include <executor/tablefunc.h>
#include <utils/array.h>
#include <utils/builtins.h>
#include <utils/lsyscache.h>
}

#include "mmseq2.h"

constexpr uint32_t OUT_TUPLE_ARITY = 18;

extern "C"
{
    PG_MODULE_MAGIC;

    PG_FUNCTION_INFO_V1(seq_search_mmseqs);

    Datum seq_search_mmseqs(PG_FUNCTION_ARGS)
    {
        // Table and column names
        char *query_tblname = text_to_cstring(PG_GETARG_TEXT_PP(0));
        char *query_colname = text_to_cstring(PG_GETARG_TEXT_PP(1));
        char *target_tblname = text_to_cstring(PG_GETARG_TEXT_PP(2));
        char *target_colname = text_to_cstring(PG_GETARG_TEXT_PP(3));

        // Determine the array element types.
        ArrayType *queries_array_type = PG_GETARG_ARRAYTYPE_P(4);
        ArrayType *targets_array_type = PG_GETARG_ARRAYTYPE_P(5);
        Oid element_type = ARR_ELEMTYPE(queries_array_type);
        int16 elem_type_width;
        bool elem_type_by_val;
        char elem_type_align;
        get_typlenbyvalalign(element_type, &elem_type_width, &elem_type_by_val, &elem_type_align);

        uint64_t *queries_array, *targets_array;
        int queries_num, targets_num;
        deconstruct_array(queries_array_type, element_type, elem_type_width, elem_type_by_val, elem_type_align,
            &queries_array, NULL, &queries_num);
        deconstruct_array(targets_array_type, element_type, elem_type_width, elem_type_by_val, elem_type_align,
            &targets_array, NULL, &targets_num);

        for (int i = 0; i < queries_num; i++)
            elog(WARNING, "%d", queries_array[i]);
        for (int i = 0; i < targets_num; i++)
            elog(WARNING, "%d", targets_array[i]);

        // Optional parameters
        char *substitution_matrix_name = text_to_cstring(PG_GETARG_TEXT_PP(6));
        uint32_t kmer_length = PG_GETARG_INT32(7);
        uint32_t kmer_gen_threshold = PG_GETARG_INT32(8);
		uint32_t ungapped_alignment_score = PG_GETARG_INT32(9);
		double eval_threshold = PG_GETARG_FLOAT8(10);
		uint32_t gap_open_cost = PG_GETARG_INT32(11);
		uint32_t gap_penalty_cost = PG_GETARG_INT32(12);

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

        // Create the vectors to pass to cpp_mmseq2
        mmseq2::Vec64Ptr qIds(new std::vector<uint64_t>{});
        mmseq2::Vec64Ptr tIds(new std::vector<uint64_t>{});
        mmseq2::VecStrPtr queries(new std::vector<mmseq2::StrPtr>{});

        std::string getQueriesQuery =
            std::string("SELECT id, ") +
            query_colname +
            std::string(" FROM ") +
            query_tblname +
            std::string(";");
        std::string getTargetsQuery =
            std::string("SELECT id FROM ") +
            target_tblname +
            std::string(";");

        SPI_connect();

        // Queries
        SPI_exec(getQueriesQuery.data(), 0);
        TupleDesc spi_tupdesc = SPI_tuptable->tupdesc;
        uint64_t processed = SPI_processed;
        for (uint32_t i = 0; i < processed; i++)
        {
            HeapTuple spi_tuple = SPI_tuptable->vals[i];
            std::string id_str = SPI_getvalue(spi_tuple, spi_tupdesc, 1);
            std::string sequence = SPI_getvalue(spi_tuple, spi_tupdesc, 2);
            qIds->push_back(std::stol(id_str));
            queries->push_back(std::make_shared<std::string>(sequence));
        }

        // Targets
        SPI_exec(getTargetsQuery.data(), 0);
        spi_tupdesc = SPI_tuptable->tupdesc;
        processed = SPI_processed;
        for (uint32_t i = 0; i < processed; i++)
        {
            HeapTuple spi_tuple = SPI_tuptable->vals[i];
            std::string id_str = SPI_getvalue(spi_tuple, spi_tupdesc, 1);
            tIds->push_back(std::stol(id_str));
        }

        SPI_finish();

        mmseq2::MMseqOutTuple tuple;
        tuple.queryId = 12;
        tuple.targetId = 34;
        tuple.rawScore = 5.34;
        tuple.bitScore = 11.2323;

        std::vector<mmseq2::MMseqOutTuple> mmseq_result;
        mmseq_result.push_back(tuple);
        uint32_t n = mmseq_result.size();

        // Build the output table
        for (uint32_t i = 0; i < n; i++)
        {
            mmseq2::MMseqOutTuple t = mmseq_result[i];
            char **values = (char **)palloc0(sizeof(char *) * OUT_TUPLE_ARITY);
            values[0] = std::to_string(t.queryId).data();
            values[1] = std::to_string(t.targetId).data();
            values[2] = std::to_string(t.rawScore).data();
            values[3] = std::to_string(t.bitScore).data();
            values[4] = std::to_string(t.eValue).data();
            values[5] = std::to_string(t.qStart).data();
            values[6] = std::to_string(t.qEnd).data();
            values[7] = std::to_string(t.qLen).data();
            values[8] = std::to_string(t.tStart).data();
            values[9] = std::to_string(t.tEnd).data();
            values[10] = std::to_string(t.tLen).data();
            values[11] = t.qAln.data();
            values[12] = t.tAln.data();
            values[13] = t.cigar.data();
            values[14] = std::to_string(t.alnLen).data();
            values[15] = std::to_string(t.mismatch).data();
            values[16] = std::to_string(t.gapOpen).data();
            values[17] = std::to_string(t.pident).data();

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
