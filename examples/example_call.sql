CREATE EXTENSION mmseq2;

CREATE TABLE query(
    id bigint,
    seq nucl_seq
);

CREATE TABLE target(
    id bigint,
    seq nucl_seq
);

INSERT INTO query VALUES (1, 'ACGTCACACACGAGGGGCGGTTTG');
INSERT INTO query VALUES (2, 'TTTTAGAGGAGACCCACCAGAGGAC');

INSERT INTO target VALUES (11, 'ATTAGCGAGAGCGCGTGTATTTT');
INSERT INTO target VALUES (12, 'GGGGGAGATATATTAGGGGACCCCCAGTTTAC');

SELECT * FROM nucl_search_db_to_db('query', 'seq', 'target', 'seq');
