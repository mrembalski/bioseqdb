CREATE OR REPLACE FUNCTION seq_search_mmseqs(
		text,
		text,
		text,
		text,
		query_ids bigint[] = NULL,
		target_ids bigint[] = NULL,
		substitution_matrix_name text = 'BLOSUM62',
		kmer_length integer = 7,
		kmer_gen_threshold integer = 2147483647,
		ungapped_alignment_score integer = 15,
		eval_threshold double precision = 0.001,
		gap_open_cost integer = 11,
		gap_penalty_cost integer = 1,
		thread_number integer = 1 -- TODO: increase when we have multithreading
	)
    RETURNS TABLE(
		query_id bigint,
		target_id bigint,
		raw_score double precision,
		bit_score double precision,
		e_value double precision,
		q_start integer,
		q_end integer,
		q_len integer,
		t_start integer,
		t_end integer,
		t_len integer,
		q_aln text,
		t_aln text,
		cigar text,
		aln_len integer,
		mismatch integer,
		gap_open integer,
		pident double precision
	)
    AS 'MODULE_PATHNAME' LANGUAGE C;
