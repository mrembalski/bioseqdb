CREATE OR REPLACE FUNCTION seq_search_mmseqs(text, text, text, text)
    RETURNS TABLE(
		queryId bigint,
		targetId bigint,
		rawScore double precision,
		bitScore double precision,
		eValue double precision,
		qStart integer,
		qEnd integer,
		qLen integer,
		tStart integer,
		tEnd integer,
		tLen integer,
		qAln text,
		tAln text,
		cigar text,
		alnLen integer,
		mismatch integer,
		gapOpen integer,
		pident double precision
	)
    AS 'MODULE_PATHNAME' LANGUAGE C IMMUTABLE STRICT;
