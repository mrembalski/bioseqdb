FROM postgres:14.2

ARG PG_DEV_VERSION=14

ARG POSTGRES_DB 
ARG POSTGRES_PASSWORD

COPY . .
 
RUN echo "$POSTGRES_DB"
RUN echo "$POSTGRES_PASSWORD"
RUN set -x \
    && apt-get update && apt-get install -y --no-install-recommends ca-certificates make postgresql-server-dev-$PG_DEV_VERSION gcc libc6-dev libssl-dev libkrb5-dev \
    # todo: move to makefile
    && gcc -fPIC -c -I /usr/include/postgresql/$PG_DEV_VERSION/server/ /mmseq2/test.c -o /mmseq2/test.o \
    && gcc -shared -o /mmseq2/test.so /mmseq2/test.o 


# .sh or .sql files from docker-entrypoint-initdb.d
# are run right after db is created
COPY ./entrypoint/index.sql /docker-entrypoint-initdb.d/
COPY ./entrypoint/mmseq.sql /docker-entrypoint-initdb.d/
