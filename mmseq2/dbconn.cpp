#include <postgres.h>
#include <libpq-fe.h>
#include <unistd.h>

static void exit_nicely(PGconn *conn)
{
    PQfinish(conn);
    exit(1);
}

class DBConn
{
private:
    PGconn *connection;
    PGnotify *notify;

    std::string tableName;
    std::string columnName;

public:
    DBConn()
    {
        connection = PQconnectdb("host=localhost port=5433 dbname=bioseqdb user=postgres password=postgres");

        if (PQstatus(connection) != CONNECTION_OK)
        {
            fprintf(stderr, "%s", PQerrorMessage(conn));

            exit_nicely(conn);
        }
    }

    uint64_t GetIndexPosition(std::string kmer, uint64_t id)
    {
        std::string getIndexQuery =
            "SELECT id FROM " +
            tableName + "_" + columnName + "__index" +
            "WHERE kmer=\'" + kmer + "\'" + ";";

        PGresult *res = PQexec(conn, getIndexQuery.cStr());

        if (PQresultStatus(res) != PGRES_COMMAND_OK)
        {
            fprintf(stderr, "%s", PQerrorMessage(conn));

            PQclear(res);
            exit_nicely(conn);
        }

        PQclear(res);
    }
};

int main()
{
}