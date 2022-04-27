#ifndef BIOSEQDB_DBCONN_H
#define BIOSEQDB_DBCONN_H

extern "C" {
#include <libpq-fe.h>
}
#include <string>
#include <memory>

namespace DB
{
    class DBconn
    {
    public:
        DBconn(const std::string&, const std::string&);
        void GetIthIndex(std::string, uint32_t, uint64_t *, uint32_t *);
        void CloseConnection();
        std::shared_ptr<std::string> GetTargetById(uint64_t);

    private:
        PGconn *connection;
        std::string tableName;
        std::string columnName;
    };
}

#endif // BIOSEQDB_DBCONN_H
