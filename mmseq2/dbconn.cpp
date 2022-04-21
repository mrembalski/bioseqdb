#include "dbconn.h"
#include <libpq-fe.h>
#include <unistd.h>
#include <string>

// class DBConn
// {
// private:
//     PGconn *connection;
//     PGnotify *notify;

//     std::string tableName;
//     std::string columnName;

// public:
//     DBConn()
//     {
//         connection = PQconnectdb("host=localhost port=5433 dbname=bioseqdb user=postgres password=postgres");

//         if (PQstatus(connection) != CONNECTION_OK)
//         {
//             fprintf(stderr, "%s", PQerrorMessage(connection));

//             exit_nicely(connection);
//         }
//     }

//     uint64_t GetIndexPosition(std::string kmer, uint64_t id)
//     {
//         std::string getIndexQuery =
//             "SELECT starting_position FROM " +
//             tableName + "_" + columnName + "__index" +
//             "WHERE kmer=\'" + kmer + "\'" + ";";

//         PGresult *res = PQexec(connection, getIndexQuery.c_str());

//         if (PQresultStatus(res) != PGRES_COMMAND_OK)
//         {
//             fprintf(stderr, "%s", PQerrorMessage(connection));

//             PQclear(res);
//             exit_nicely(connection);
//         }

//         PQclear(res);
//     }
// };
