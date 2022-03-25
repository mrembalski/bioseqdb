CREATE OR REPLACE FUNCTION seq_search_mmseqs(bigint, boolean, dna_sequence, int, text)
    RETURNS SETOF record
    AS '/mmseq2/test.so', 'seq_search_mmseqs'
    LANGUAGE C STRICT;

