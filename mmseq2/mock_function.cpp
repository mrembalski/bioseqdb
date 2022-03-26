#include <cstdio>
#include <cstring>

#include "psql_interface.h"

#include "mock_function.h"

void mock_function(void *data, char *column_name)
{
    uint64_t num = get_indexes(column_name, "AAAAAAA");

    char answer[200];
    memset(answer, 0, 200);
    sprintf(answer, "%s%lu", "Number of AAAAAAA indexes: ", num);
    log_from_cpp(answer);

    for (int i = 0; i < num; i++)
    {
        uint64_t target_id;
        uint32_t position;
        get_ith_index(i, &target_id, &position);

        memset(answer, 0, 200);
        sprintf(answer, "%s%lu%s%d",
                "target_id = ",
                target_id,
                ", position = ",
                position);
        log_from_cpp(answer);
    }
}
