FROM postgres:14.2

ARG PG_DEV_VERSION=14

ARG POSTGRES_DB 
ARG POSTGRES_PASSWORD

COPY . .
 
RUN echo "$POSTGRES_DB"
RUN echo "$POSTGRES_PASSWORD"
RUN set -x \
    && apt-get update && apt-get install -y --no-install-recommends ca-certificates make postgresql-server-dev-$PG_DEV_VERSION gcc g++ libc6-dev libssl-dev libkrb5-dev \
    # todo: move to makefile
    && g++ -std=c++17 -fPIC -c -I /usr/include/postgresql/$PG_DEV_VERSION/server/ /mmseq2/extension.cpp -o /mmseq2/extension.o \
    && g++ -std=c++17 -fPIC -c /mmseq2/mock_structures.cpp -o /mmseq2/mock_structures.o \
    && g++ -std=c++17 -pthread -fPIC -c /mmseq2/mmseq2.cpp -o /mmseq2/mmseq2.o \
    && gcc -shared -o /mmseq2/extension.so /mmseq2/extension.o /mmseq2/mock_structures.o /mmseq2/mmseq2.o


# .sh or .sql files from docker-entrypoint-initdb.d
# are run right after db is created
COPY ./entrypoint/index.sql /docker-entrypoint-initdb.d/
