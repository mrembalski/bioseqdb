#include "postgres.h"
#include "fmgr.h"

#include "executor/spi.h"
#include "utils/builtins.h"
#include "utils/array.h"
#include "utils/lsyscache.h"
#include "funcapi.h"

PG_MODULE_MAGIC;

PG_FUNCTION_INFO_V1(test_function);
PG_FUNCTION_INFO_V1(addTen);

Datum test_function(PG_FUNCTION_ARGS)
{
    SRF_RETURN_DONE(funcctx);
}