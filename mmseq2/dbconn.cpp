#include "dbconn.h"
#include <libpq-fe.h>
#include <unistd.h>
#include <string>
#include <stdio.h>
#include <iostream>
#include <exception>

DB::DBconn::DBconn(std::string tableName, std::string columnName)
{
    this->columnName = columnName;
    this->tableName = tableName;

    this->connection = PQconnectdb("host=localhost port=5432 dbname=bioseqdb user=postgres password=postgres");

    if (PQstatus(this->connection) != CONNECTION_OK)
    {
        fprintf(stderr, "%s", PQerrorMessage(this->connection));

        PQfinish(this->connection);
        exit(1);
    }
}

void DB::DBconn::GetIthIndex(std::string kmer, uint32_t i, uint64_t *target_id, uint32_t *position)
{
    std::string getIndexQuery =
        "SELECT starting_position, dna_sequence_id FROM " +
        this->tableName + "_" + this->columnName + "__index" +
        " " + "WHERE kmer=\'" + kmer + "\'" +
        "OFFSET " + std::to_string(i) + "LIMIT 1;";

    PGresult *res = PQexec(connection, getIndexQuery.c_str());

    if (PQresultStatus(res) != PGRES_TUPLES_OK) {
        throw std::invalid_argument("not PGRES_TUPLES_OK");
    }

    int starting_position_fnum = PQfnumber(res, "starting_position");
    int dna_sequence_id_fnum = PQfnumber(res, "dna_sequence_id");

    /* TODO: do something when no value is returned */
    if (PQntuples(res) != 1) {
        throw std::invalid_argument("No ith index exists");
    }

    char *starting_position = PQgetvalue(res, 0, starting_position_fnum);
    char *dna_sequence_id = PQgetvalue(res, 0, dna_sequence_id_fnum);

    *position = strtoul(starting_position, NULL, 0);
    *target_id = strtoull(dna_sequence_id, NULL, 0);

    PQclear(res);
}


void DB::DBconn::CloseConnection()
{
    PQfinish(this->connection);
}

std::shared_ptr<std::string> DB::DBconn::GetTargetById(uint64_t id)
{
    std::string getTargetQuery =
        "SELECT " + this->columnName + " FROM " +
        this->tableName + " WHERE id=" + std::to_string(id) + ";";

    PGresult *res = PQexec(connection, getTargetQuery.c_str());

    if (PQresultStatus(res) != PGRES_TUPLES_OK) {
        throw std::invalid_argument("not PGRES_TUPLES_OK");
    }

    int target_seq_fnum = PQfnumber(res, this->columnName.c_str());

    /* TODO: do something when no value is returned */
    if (PQntuples(res) != 1) {        
        throw std::invalid_argument("No ith target exists");
    }

    char *target_seq = PQgetvalue(res, 0, target_seq_fnum);
    
    PQclear(res);

    return std::make_shared<std::string>(target_seq);
}
