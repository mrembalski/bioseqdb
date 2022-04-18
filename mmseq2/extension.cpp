// Postgres C libraries
extern "C"
{
#include <postgres.h>
#include <fmgr.h>
#include <executor/spi.h>    // SPI_*
#include <funcapi.h>         // get_call_result_type, AttInMetadata, TupleDescGetAttInMetadata, BuildTupleFromCStrings
#include <miscadmin.h>       // work_mem
#include <utils/array.h>     // deconstruct_array, ARR_ELEMTYPE, PG_GETARG_ARRAYTYPE_P
#include <utils/builtins.h>  // text_to_cstring
#include <utils/lsyscache.h> // get_typlenbyvalalign
}

// C++ STL
#include <vector>
#include <string>

// C++ MMSeq microservice libraries
#include "mmseq2.h"

namespace
{
    constexpr uint32_t OUT_TUPLE_ARITY = 18;

    void add_one_sequence(const std::string &sequence,
                          const mmseq2::Vec64Ptr &ids, const mmseq2::VecStrPtr &sequences)
    {
        ids->push_back(1);
        sequences->push_back(std::make_shared<std::string>(sequence));
    }

    void add_array_of_sequences(ArrayType *const sequences_array_type,
                                const mmseq2::Vec64Ptr &ids, const mmseq2::VecStrPtr &sequences)
    {
        // Get the element type of the array of queries
        const Oid element_type = ARR_ELEMTYPE(sequences_array_type);
        int16 elem_type_width;
        bool elem_type_by_val;
        char elem_type_align;

        // Deconstruct the array
        Datum *sequences_array;
        int sequences_num;
        get_typlenbyvalalign(element_type, &elem_type_width, &elem_type_by_val, &elem_type_align);
        deconstruct_array(sequences_array_type, element_type, elem_type_width, elem_type_by_val, elem_type_align,
                          &sequences_array, NULL, &sequences_num);

        for (int i = 0; i < sequences_num; i++)
        {
            ids->push_back(i + 1);
            const std::string sequence = TextDatumGetCString(sequences_array[i]);
            sequences->push_back(std::make_shared<std::string>(sequence));
        }
    }

    void add_queries_all(const std::string &query_tblname, const std::string &query_colname,
                         const mmseq2::Vec64Ptr &qIds, const mmseq2::VecStrPtr &queries)
    {
        // Perform a query via SPI
        const std::string getQueriesQuery =
            std::string("SELECT id, ") +
            query_colname +
            std::string(" FROM ") +
            query_tblname +
            std::string(";");
        SPI_exec(getQueriesQuery.data(), 0);

        // Move data from SPI_tuptable to std::vectors qIds and queries
        const TupleDesc spi_tupdesc = SPI_tuptable->tupdesc;
        const uint64_t processed = SPI_processed;
        for (uint32_t i = 0; i < processed; i++)
        {
            const HeapTuple spi_tuple = SPI_tuptable->vals[i];
            const std::string id_str = SPI_getvalue(spi_tuple, spi_tupdesc, 1);
            const std::string sequence = SPI_getvalue(spi_tuple, spi_tupdesc, 2);
            qIds->push_back(std::stol(id_str));
            queries->push_back(std::make_shared<std::string>(sequence));
        }
    }

    void add_queries_with_ids(const std::string &query_tblname, const std::string &query_colname,
                              const mmseq2::Vec64Ptr &qIds, const mmseq2::VecStrPtr &queries,
                              ArrayType *const queries_array_type)
    {
        // Get the element type of the array of queries
        const Oid element_type = ARR_ELEMTYPE(queries_array_type);
        int16 elem_type_width;
        bool elem_type_by_val;
        char elem_type_align;

        // Deconstruct the array
        uint64_t *queries_array;
        int queries_num;
        get_typlenbyvalalign(element_type, &elem_type_width, &elem_type_by_val, &elem_type_align);
        deconstruct_array(queries_array_type, element_type, elem_type_width, elem_type_by_val, elem_type_align,
                          &queries_array, NULL, &queries_num);

        for (int i = 0; i < queries_num; i++)
        {
            qIds->push_back(queries_array[i]);
            const std::string getOneQueryQuery =
                std::string("SELECT ") +
                query_colname +
                std::string(" FROM ") +
                query_tblname +
                std::string(" WHERE id = ") +
                std::to_string(queries_array[i]) +
                std::string(";");
            SPI_exec(getOneQueryQuery.data(), 0);
            const TupleDesc spi_tupdesc = SPI_tuptable->tupdesc;
            const HeapTuple spi_tuple = SPI_tuptable->vals[0];
            const std::string sequence = SPI_getvalue(spi_tuple, spi_tupdesc, 1);
            queries->push_back(std::make_shared<std::string>(sequence));
        }

        // Free the allocated memory
        pfree(queries_array);
    }

    void add_targets_all(const std::string &target_tblname, const std::string &target_colname,
                         const mmseq2::Vec64Ptr &tIds)
    {
        const std::string getTargetsQuery =
            std::string("SELECT id FROM ") +
            target_tblname +
            std::string(";");
        SPI_exec(getTargetsQuery.data(), 0);

        const TupleDesc spi_tupdesc = SPI_tuptable->tupdesc;
        const uint64_t processed = SPI_processed;
        for (uint32_t i = 0; i < processed; i++)
        {
            const HeapTuple spi_tuple = SPI_tuptable->vals[i];
            const std::string id_str = SPI_getvalue(spi_tuple, spi_tupdesc, 1);
            tIds->push_back(std::stol(id_str));
        }
    }

    void add_targets_with_ids(const std::string &target_tblname, const std::string &target_colname,
                              const mmseq2::Vec64Ptr &tIds,
                              ArrayType *const targets_array_type)
    {
        const Oid element_type = ARR_ELEMTYPE(targets_array_type);
        int16 elem_type_width;
        bool elem_type_by_val;
        char elem_type_align;

        uint64_t *targets_array;
        int targets_num;
        get_typlenbyvalalign(element_type, &elem_type_width, &elem_type_by_val, &elem_type_align);
        deconstruct_array(targets_array_type, element_type, elem_type_width, elem_type_by_val, elem_type_align,
                          &targets_array, NULL, &targets_num);

        for (int i = 0; i < targets_num; i++)
        {
            tIds->push_back(targets_array[i]);
        }
    }
}

extern "C"
{
    PG_MODULE_MAGIC;

    PG_FUNCTION_INFO_V1(seq_search_mmseqs_one_to_one);
    PG_FUNCTION_INFO_V1(seq_search_mmseqs_arr_to_arr);
    PG_FUNCTION_INFO_V1(seq_search_mmseqs_db_to_db);

    Datum seq_search_mmseqs_one_to_one(PG_FUNCTION_ARGS)
    {
        // Optional parameters
        char *substitution_matrix_name = text_to_cstring(PG_GETARG_TEXT_PP(2));
        uint32_t kmer_length = PG_GETARG_INT32(3);
        uint32_t kmer_gen_threshold = PG_GETARG_INT32(4);
        uint32_t ungapped_alignment_score = PG_GETARG_INT32(5);
        double eval_threshold = PG_GETARG_FLOAT8(6);
        uint32_t gap_open_cost = PG_GETARG_INT32(7);
        uint32_t gap_penalty_cost = PG_GETARG_INT32(8);

        // Construct vectors of ids and queries
        mmseq2::Vec64Ptr qIds(new std::vector<uint64_t>{});
        mmseq2::Vec64Ptr tIds(new std::vector<uint64_t>{});
        mmseq2::VecStrPtr queries(new std::vector<mmseq2::StrPtr>{});
        mmseq2::VecStrPtr targets(new std::vector<mmseq2::StrPtr>{});

        add_one_sequence(text_to_cstring(PG_GETARG_TEXT_PP(0)), qIds, queries);
        add_one_sequence(text_to_cstring(PG_GETARG_TEXT_PP(1)), tIds, targets);

        // Write out the contents of the vectors for testing purposes
        for (int i = 0; i < qIds->size(); i++)
            elog(WARNING, "%d", (*qIds)[i]);
        for (int i = 0; i < queries->size(); i++)
            elog(WARNING, "%s", ((*queries)[i])->data());
        for (int i = 0; i < tIds->size(); i++)
            elog(WARNING, "%d", (*tIds)[i]);
        for (int i = 0; i < targets->size(); i++)
            elog(WARNING, "%s", ((*targets)[i])->data());

        // Just return null for now
        PG_RETURN_NULL();
    }

    Datum seq_search_mmseqs_arr_to_arr(PG_FUNCTION_ARGS)
    {
        // Optional parameters
        char *substitution_matrix_name = text_to_cstring(PG_GETARG_TEXT_PP(2));
        uint32_t kmer_length = PG_GETARG_INT32(3);
        uint32_t kmer_gen_threshold = PG_GETARG_INT32(4);
        uint32_t ungapped_alignment_score = PG_GETARG_INT32(5);
        double eval_threshold = PG_GETARG_FLOAT8(6);
        uint32_t gap_open_cost = PG_GETARG_INT32(7);
        uint32_t gap_penalty_cost = PG_GETARG_INT32(8);

        // Construct vectors of ids and queries
        mmseq2::Vec64Ptr qIds(new std::vector<uint64_t>{});
        mmseq2::Vec64Ptr tIds(new std::vector<uint64_t>{});
        mmseq2::VecStrPtr queries(new std::vector<mmseq2::StrPtr>{});
        mmseq2::VecStrPtr targets(new std::vector<mmseq2::StrPtr>{});

        add_array_of_sequences(PG_GETARG_ARRAYTYPE_P(0), qIds, queries);
        add_array_of_sequences(PG_GETARG_ARRAYTYPE_P(1), tIds, targets);

        // Write out the contents of the vectors for testing purposes
        for (int i = 0; i < qIds->size(); i++)
            elog(WARNING, "%d", (*qIds)[i]);
        for (int i = 0; i < queries->size(); i++)
            elog(WARNING, "%s", ((*queries)[i])->data());
        for (int i = 0; i < tIds->size(); i++)
            elog(WARNING, "%d", (*tIds)[i]);
        for (int i = 0; i < targets->size(); i++)
            elog(WARNING, "%s", ((*targets)[i])->data());

        // Just return null for now
        PG_RETURN_NULL();
    }

    Datum seq_search_mmseqs_db_to_db(PG_FUNCTION_ARGS)
    {
        // Table and column names
        char *query_tblname = text_to_cstring(PG_GETARG_TEXT_PP(0));
        char *query_colname = text_to_cstring(PG_GETARG_TEXT_PP(1));
        char *target_tblname = text_to_cstring(PG_GETARG_TEXT_PP(2));
        char *target_colname = text_to_cstring(PG_GETARG_TEXT_PP(3));

        // Optional parameters
        char *substitution_matrix_name = text_to_cstring(PG_GETARG_TEXT_PP(6));
        uint32_t kmer_length = PG_GETARG_INT32(7);
        uint32_t kmer_gen_threshold = PG_GETARG_INT32(8);
        uint32_t ungapped_alignment_score = PG_GETARG_INT32(9);
        double eval_threshold = PG_GETARG_FLOAT8(10);
        uint32_t gap_open_cost = PG_GETARG_INT32(11);
        uint32_t gap_penalty_cost = PG_GETARG_INT32(12);

        // Construct vectors of ids and queries
        mmseq2::Vec64Ptr qIds(new std::vector<uint64_t>{});
        mmseq2::Vec64Ptr tIds(new std::vector<uint64_t>{});
        mmseq2::VecStrPtr queries(new std::vector<mmseq2::StrPtr>{});

        SPI_connect();

        // Prepare queries
        if (PG_ARGISNULL(4))
            add_queries_all(query_tblname, query_colname, qIds, queries);
        else
            add_queries_with_ids(query_tblname, query_colname, qIds, queries, PG_GETARG_ARRAYTYPE_P(4));

        // Prepare targets
        if (PG_ARGISNULL(5))
            add_targets_all(target_tblname, target_colname, tIds);
        else
            add_targets_with_ids(target_tblname, target_colname, tIds, PG_GETARG_ARRAYTYPE_P(5));

        SPI_finish();

        // Write out the contents of the vectors for testing purposes
        for (int i = 0; i < qIds->size(); i++)
            elog(WARNING, "%d", (*qIds)[i]);
        for (int i = 0; i < queries->size(); i++)
            elog(WARNING, "%s", ((*queries)[i])->data());
        for (int i = 0; i < tIds->size(); i++)
            elog(WARNING, "%d", (*tIds)[i]);

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

        std::vector<mmseq2::MMseqOutTuple> mmseq_result;
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

            const HeapTuple tuple = BuildTupleFromCStrings(attinmeta, values);
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
