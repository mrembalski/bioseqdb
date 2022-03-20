#!/bin/bash

cat ./indexes/create.sql | docker exec -i bioseqdb-db-1 psql -U postgres -d bioseqdb
