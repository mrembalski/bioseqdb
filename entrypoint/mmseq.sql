CREATE OR REPLACE FUNCTION test_function(integer)
    RETURNS SETOF integer
    AS '/mmseq2/test.so', 'test_function'
    LANGUAGE C STRICT;
    