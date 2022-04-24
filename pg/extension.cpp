#include <vector>
#include <string>
#include "../common/mmseq2lib.h"
#include "rpc/client.h"

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

namespace
{
    constexpr uint32_t OUT_TUPLE_ARITY = 18;

    void add_one_sequence(const std::string &sequence,
                          const common::InputParams::Vec64Ptr &ids, const common::InputParams::VecStrPtr &sequences)
    {
        ids->push_back(1);
        sequences->push_back(std::make_shared<std::string>(sequence));
    }

    void add_array_of_sequences(ArrayType *const sequences_array_type,
                                const common::InputParams::Vec64Ptr &ids, const common::InputParams::VecStrPtr &sequences)
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
                         const common::InputParams::Vec64Ptr &qIds, const common::InputParams::VecStrPtr &queries)
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
                              const common::InputParams::Vec64Ptr &qIds, const common::InputParams::VecStrPtr &queries,
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
                         const common::InputParams::Vec64Ptr &tIds)
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
                              const common::InputParams::Vec64Ptr &tIds,
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
    PG_FUNCTION_INFO_V1(seq_search_mmseqs_arr_to_one);
    PG_FUNCTION_INFO_V1(seq_search_mmseqs_db_to_one);
    PG_FUNCTION_INFO_V1(seq_search_mmseqs_one_to_arr);
    PG_FUNCTION_INFO_V1(seq_search_mmseqs_arr_to_arr);
    PG_FUNCTION_INFO_V1(seq_search_mmseqs_db_to_arr);
    PG_FUNCTION_INFO_V1(seq_search_mmseqs_one_to_db);
    PG_FUNCTION_INFO_V1(seq_search_mmseqs_arr_to_db);
    PG_FUNCTION_INFO_V1(seq_search_mmseqs_db_to_db);

    Datum seq_search_mmseqs_main(const std::optional<std::string> target_tblname,
                                 const std::optional<std::string> target_colname,
                                 const common::InputParams::Vec64Ptr &qIds,
                                 const common::InputParams::Vec64Ptr &tIds,
                                 const common::InputParams::VecStrPtr &queries,
                                 const common::InputParams::VecStrPtr &targets,
                                 const FunctionCallInfo &fcinfo,
                                 uint32_t fst_opt_param)
    {
        // Optional parameters
        uint32_t kmer_length = PG_GETARG_INT32(fst_opt_param);
        std::string substitution_matrix_name = text_to_cstring(PG_GETARG_TEXT_PP(fst_opt_param + 1));
        uint32_t kmer_gen_threshold = PG_GETARG_INT32(fst_opt_param + 2);
        uint32_t ungapped_alignment_score = PG_GETARG_INT32(fst_opt_param + 3);
        double eval_threshold = PG_GETARG_FLOAT8(fst_opt_param + 4);
        uint32_t gap_open_cost = PG_GETARG_INT32(fst_opt_param + 5);
        uint32_t gap_penalty_cost = PG_GETARG_INT32(fst_opt_param + 6);
        uint32_t thread_number = PG_GETARG_INT32(fst_opt_param + 7);

        // Write out the contents of the vectors
        elog(WARNING, "%s", "----- Passed sequences -----");
        elog(WARNING, "%s%ld%s", "Query ids (", qIds->size(), "):");
        for (int i = 0; i < qIds->size(); i++)
            elog(WARNING, "%ld", (*qIds)[i]);
        elog(WARNING, "%s%ld%s", "Queries (", queries->size(), "):");
        for (int i = 0; i < queries->size(); i++)
            elog(WARNING, "%s", ((*queries)[i])->data());
        elog(WARNING, "%s%ld%s", "Target ids (", tIds->size(), "):");
        for (int i = 0; i < tIds->size(); i++)
            elog(WARNING, "%ld", (*tIds)[i]);
        elog(WARNING, "%s%ld%s", "Targets (", targets->size(), "):");
        for (int i = 0; i < targets->size(); i++)
            elog(WARNING, "%s", ((*targets)[i])->data());

        // Write out the values of the optional parameters
        elog(WARNING, "%s", "----- Optional parameters -----");
        elog(WARNING, "%s%s", "Substitution matrix: ", substitution_matrix_name.data());
        elog(WARNING, "%s%d", "Kmer length: ", kmer_length);
        elog(WARNING, "%s%d", "Kmer gen threshold: ", kmer_gen_threshold);
        elog(WARNING, "%s%d", "Ungapped alignment score: ", ungapped_alignment_score);
        elog(WARNING, "%s%f", "Eval threshold: ", eval_threshold);
        elog(WARNING, "%s%d", "Gap open cost: ", gap_open_cost);
        elog(WARNING, "%s%d", "Gap penalty cost: ", gap_penalty_cost);

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

        // Input params
        common::InputParams input_params(qIds->size(), tIds->size(), qIds, tIds, queries,
                                         std::make_shared<std::string>(target_tblname.value()),
                                         std::make_shared<std::string>(target_colname.value()),
                                         std::make_shared<std::string>(substitution_matrix_name),
                                         kmer_length,
                                         kmer_gen_threshold,
                                         ungapped_alignment_score,
                                         eval_threshold,
                                         gap_open_cost,
                                         gap_penalty_cost,
                                         thread_number);

        // Client
        rpc::client c("localhost", 8080);
        common::VecRes mmseq_result(c.call("mmseq2", std::make_shared<common::InputParams>(input_params)).as<common::VecRes>());
        uint32_t n = mmseq_result.size();
        elog(WARNING, "%s%d", "Returned rows: ", n);

        // Build the output table
        for (uint32_t i = 0; i < n; i++)
        {
            common::MmseqResult t = mmseq_result[i];
            char **values = (char **)palloc0(sizeof(char *) * OUT_TUPLE_ARITY);
            values[0] = std::to_string(t.getQueryId()).data();
            values[1] = std::to_string(t.getTargetId()).data();
            values[2] = std::to_string(t.getRawScore()).data();
            values[3] = std::to_string(t.getBitScore()).data();
            values[4] = std::to_string(t.getEValue()).data();
            values[5] = std::to_string(t.getQStart()).data();
            values[6] = std::to_string(t.getQEnd()).data();
            values[7] = std::to_string(t.getQLen()).data();
            values[8] = std::to_string(t.getTStart()).data();
            values[9] = std::to_string(t.getTEnd()).data();
            values[10] = std::to_string(t.getTLen()).data();
            values[11] = t.getQAln().data();
            values[12] = t.getTAln().data();
            values[13] = t.getCigar().data();
            values[14] = std::to_string(t.getAlnLen()).data();
            values[15] = std::to_string(t.getMismatch()).data();
            values[16] = std::to_string(t.getGapOpen()).data();
            values[17] = std::to_string(t.getPident()).data();

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

    Datum seq_search_mmseqs_one_to_one(PG_FUNCTION_ARGS)
    {
        common::InputParams::Vec64Ptr qIds(new std::vector<uint64_t>{});
        common::InputParams::Vec64Ptr tIds(new std::vector<uint64_t>{});
        common::InputParams::VecStrPtr queries(new std::vector<common::InputParams::StrPtr>{});
        common::InputParams::VecStrPtr targets(new std::vector<common::InputParams::StrPtr>{});

        add_one_sequence(text_to_cstring(PG_GETARG_TEXT_PP(0)), qIds, queries);
        add_one_sequence(text_to_cstring(PG_GETARG_TEXT_PP(1)), tIds, targets);

        return seq_search_mmseqs_main(std::nullopt, std::nullopt,
                                      qIds, tIds, queries, targets,
                                      fcinfo, 2);
    }

    Datum seq_search_mmseqs_arr_to_one(PG_FUNCTION_ARGS)
    {
        common::InputParams::Vec64Ptr qIds(new std::vector<uint64_t>{});
        common::InputParams::Vec64Ptr tIds(new std::vector<uint64_t>{});
        common::InputParams::VecStrPtr queries(new std::vector<common::InputParams::StrPtr>{});
        common::InputParams::VecStrPtr targets(new std::vector<common::InputParams::StrPtr>{});

        add_array_of_sequences(PG_GETARG_ARRAYTYPE_P(0), qIds, queries);
        add_one_sequence(text_to_cstring(PG_GETARG_TEXT_PP(1)), tIds, targets);

        return seq_search_mmseqs_main(std::nullopt, std::nullopt,
                                      qIds, tIds, queries, targets,
                                      fcinfo, 2);
    }

    Datum seq_search_mmseqs_db_to_one(PG_FUNCTION_ARGS)
    {
        std::string query_tblname = text_to_cstring(PG_GETARG_TEXT_PP(0));
        std::string query_colname = text_to_cstring(PG_GETARG_TEXT_PP(1));

        common::InputParams::Vec64Ptr qIds(new std::vector<uint64_t>{});
        common::InputParams::Vec64Ptr tIds(new std::vector<uint64_t>{});
        common::InputParams::VecStrPtr queries(new std::vector<common::InputParams::StrPtr>{});
        common::InputParams::VecStrPtr targets(new std::vector<common::InputParams::StrPtr>{});

        SPI_connect();
        if (PG_ARGISNULL(3))
            add_queries_all(query_tblname, query_colname, qIds, queries);
        else
            add_queries_with_ids(query_tblname, query_colname, qIds, queries, PG_GETARG_ARRAYTYPE_P(3));
        SPI_finish();

        add_one_sequence(text_to_cstring(PG_GETARG_TEXT_PP(2)), tIds, targets);

        return seq_search_mmseqs_main(std::nullopt, std::nullopt,
                                      qIds, tIds, queries, targets,
                                      fcinfo, 4);
    }

    Datum seq_search_mmseqs_one_to_arr(PG_FUNCTION_ARGS)
    {
        common::InputParams::Vec64Ptr qIds(new std::vector<uint64_t>{});
        common::InputParams::Vec64Ptr tIds(new std::vector<uint64_t>{});
        common::InputParams::VecStrPtr queries(new std::vector<common::InputParams::StrPtr>{});
        common::InputParams::VecStrPtr targets(new std::vector<common::InputParams::StrPtr>{});

        add_one_sequence(text_to_cstring(PG_GETARG_TEXT_PP(0)), qIds, queries);
        add_array_of_sequences(PG_GETARG_ARRAYTYPE_P(1), tIds, targets);

        return seq_search_mmseqs_main(std::nullopt, std::nullopt,
                                      qIds, tIds, queries, targets,
                                      fcinfo, 2);
    }

    Datum seq_search_mmseqs_arr_to_arr(PG_FUNCTION_ARGS)
    {
        common::InputParams::Vec64Ptr qIds(new std::vector<uint64_t>{});
        common::InputParams::Vec64Ptr tIds(new std::vector<uint64_t>{});
        common::InputParams::VecStrPtr queries(new std::vector<common::InputParams::StrPtr>{});
        common::InputParams::VecStrPtr targets(new std::vector<common::InputParams::StrPtr>{});

        add_array_of_sequences(PG_GETARG_ARRAYTYPE_P(0), qIds, queries);
        add_array_of_sequences(PG_GETARG_ARRAYTYPE_P(1), tIds, targets);

        return seq_search_mmseqs_main(std::nullopt, std::nullopt,
                                      qIds, tIds, queries, targets,
                                      fcinfo, 2);
    }

    Datum seq_search_mmseqs_db_to_arr(PG_FUNCTION_ARGS)
    {
        std::string query_tblname = text_to_cstring(PG_GETARG_TEXT_PP(0));
        std::string query_colname = text_to_cstring(PG_GETARG_TEXT_PP(1));

        common::InputParams::Vec64Ptr qIds(new std::vector<uint64_t>{});
        common::InputParams::Vec64Ptr tIds(new std::vector<uint64_t>{});
        common::InputParams::VecStrPtr queries(new std::vector<common::InputParams::StrPtr>{});
        common::InputParams::VecStrPtr targets(new std::vector<common::InputParams::StrPtr>{});

        SPI_connect();
        if (PG_ARGISNULL(3))
            add_queries_all(query_tblname, query_colname, qIds, queries);
        else
            add_queries_with_ids(query_tblname, query_colname, qIds, queries, PG_GETARG_ARRAYTYPE_P(3));
        SPI_finish();

        add_array_of_sequences(PG_GETARG_ARRAYTYPE_P(2), tIds, targets);

        return seq_search_mmseqs_main(std::nullopt, std::nullopt,
                                      qIds, tIds, queries, targets,
                                      fcinfo, 4);
    }

    Datum seq_search_mmseqs_one_to_db(PG_FUNCTION_ARGS)
    {
        std::string target_tblname = text_to_cstring(PG_GETARG_TEXT_PP(1));
        std::string target_colname = text_to_cstring(PG_GETARG_TEXT_PP(2));

        common::InputParams::Vec64Ptr qIds(new std::vector<uint64_t>{});
        common::InputParams::Vec64Ptr tIds(new std::vector<uint64_t>{});
        common::InputParams::VecStrPtr queries(new std::vector<common::InputParams::StrPtr>{});
        common::InputParams::VecStrPtr targets(new std::vector<common::InputParams::StrPtr>{});

        add_one_sequence(text_to_cstring(PG_GETARG_TEXT_PP(0)), qIds, queries);

        SPI_connect();
        if (PG_ARGISNULL(3))
            add_targets_all(target_tblname, target_colname, tIds);
        else
            add_targets_with_ids(target_tblname, target_colname, tIds, PG_GETARG_ARRAYTYPE_P(3));
        SPI_finish();

        return seq_search_mmseqs_main(std::optional<std::string>{target_tblname},
                                      std::optional<std::string>{target_colname},
                                      qIds, tIds, queries, targets,
                                      fcinfo, 4);
    }

    Datum seq_search_mmseqs_arr_to_db(PG_FUNCTION_ARGS)
    {
        std::string target_tblname = text_to_cstring(PG_GETARG_TEXT_PP(1));
        std::string target_colname = text_to_cstring(PG_GETARG_TEXT_PP(2));

        common::InputParams::Vec64Ptr qIds(new std::vector<uint64_t>{});
        common::InputParams::Vec64Ptr tIds(new std::vector<uint64_t>{});
        common::InputParams::VecStrPtr queries(new std::vector<common::InputParams::StrPtr>{});
        common::InputParams::VecStrPtr targets(new std::vector<common::InputParams::StrPtr>{});

        add_array_of_sequences(PG_GETARG_ARRAYTYPE_P(0), qIds, queries);

        SPI_connect();
        if (PG_ARGISNULL(3))
            add_targets_all(target_tblname, target_colname, tIds);
        else
            add_targets_with_ids(target_tblname, target_colname, tIds, PG_GETARG_ARRAYTYPE_P(3));
        SPI_finish();

        return seq_search_mmseqs_main(std::optional<std::string>{target_tblname},
                                      std::optional<std::string>{target_colname},
                                      qIds, tIds, queries, targets,
                                      fcinfo, 4);
    }

    Datum seq_search_mmseqs_db_to_db(PG_FUNCTION_ARGS)
    {
        std::string query_tblname = text_to_cstring(PG_GETARG_TEXT_PP(0));
        std::string query_colname = text_to_cstring(PG_GETARG_TEXT_PP(1));
        std::string target_tblname = text_to_cstring(PG_GETARG_TEXT_PP(2));
        std::string target_colname = text_to_cstring(PG_GETARG_TEXT_PP(3));

        common::InputParams::Vec64Ptr qIds(new std::vector<uint64_t>{});
        common::InputParams::Vec64Ptr tIds(new std::vector<uint64_t>{});
        common::InputParams::VecStrPtr queries(new std::vector<common::InputParams::StrPtr>{});
        common::InputParams::VecStrPtr targets(new std::vector<common::InputParams::StrPtr>{});

        SPI_connect();

        if (PG_ARGISNULL(4))
            add_queries_all(query_tblname, query_colname, qIds, queries);
        else
            add_queries_with_ids(query_tblname, query_colname, qIds, queries, PG_GETARG_ARRAYTYPE_P(4));

        if (PG_ARGISNULL(5))
            add_targets_all(target_tblname, target_colname, tIds);
        else
            add_targets_with_ids(target_tblname, target_colname, tIds, PG_GETARG_ARRAYTYPE_P(5));

        SPI_finish();

        return seq_search_mmseqs_main(std::optional<std::string>{target_tblname},
                                      std::optional<std::string>{target_colname},
                                      qIds, tIds, queries, targets,
                                      fcinfo, 6);
    }
}
