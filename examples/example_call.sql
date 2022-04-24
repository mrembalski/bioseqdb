CREATE EXTENSION mmseq2;

CREATE TABLE query(
    id bigint,
    seq dna_sequence
);

CREATE TABLE target(
    id bigint,
    seq dna_sequence
);

INSERT INTO query VALUES (1, 'ACGT');
INSERT INTO query VALUES (2, 'TTTT');

INSERT INTO target VALUES (11, 'ATTA');
INSERT INTO target VALUES (12, 'GGGGG');

SELECT * FROM seq_search_mmseqs('query', 'seq', 'target', 'seq', ARRAY(SELECT id FROM query), ARRAY(SELECT id FROM target));
