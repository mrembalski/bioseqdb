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
        uint64_t GetIndexPosition(std::string, uint64_t);

    private:
        PGconn *connection;
        std::string tableName;
        std::string columnName;
    };
}

#endif // BIOSEQDB_DBCONN_H
