#ifndef BIOSEQDB_DBCONN_H
#define BIOSEQDB_DBCONN_H

#include <libpq-fe.h>
#include <string>

namespace DB
{
    class DBconn
    {
    public:
        DBconn(std::string, std::string);
        void GetIthIndex(std::string, uint32_t, uint64_t *, uint32_t *);

    private:
        PGconn *connection;
        std::string tableName;
        std::string columnName;
    };
}

#endif // BIOSEQDB_DBCONN_H
