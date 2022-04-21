#ifndef BIOSEQDB_DBCONN_H
#define BIOSEQDB_DBCONN_H

#include <libpq-fe.h>

namespace DB
{
    class DBconn
    {
    public:
        DBconn(std::string tableName, std::string columnName)
        {
            this->columnName = columnName;
            this->tableName = tableName;

            this->connection = PQconnectdb("host=localhost port=5433 dbname=bioseqdb user=postgres password=postgres");

            if (PQstatus(this->connection) != CONNECTION_OK)
            {
                fprintf(stderr, "%s", PQerrorMessage(this->connection));

                PQfinish(this->connection);
                exit(1);
            }
        }

        uint64_t GetIndexPosition(std::string kmer, uint64_t id)
        {
            std::string getIndexQuery =
                "SELECT starting_position FROM " +
                this->tableName + "_" + this->columnName + "__index" + " "
                                                                       "WHERE kmer=\'" +
                kmer + "\'" + " AND dna_sequence_id = " + std::to_string(id) +
                ";";

            PGresult *res = PQexec(connection, getIndexQuery.c_str());

            if (PQresultStatus(res) != PGRES_TUPLES_OK)
            {
                fprintf(stderr, "%s", PQerrorMessage(connection));

                PQclear(res);

                PQfinish(this->connection);
                exit(1);
            }

            int starting_position_fnum = PQfnumber(res, "starting_position");

            for (auto i = 0; i < PQntuples(res); i++)
            {
                auto iptr = PQgetvalue(res, i, starting_position_fnum);
                std::cout << iptr << std::endl;
            }

            PQclear(res);

            return 0;
        }

    private:
        PGconn *connection;
        std::string tableName;
        std::string columnName;
    };
}

#endif // BIOSEQDB_DBCONN_H
