CREATE DOMAIN DNA_SEQUENCE AS TEXT
CHECK(
	-- Some PostgreSQL magic using \m and \M
   VALUE ~ '(\m[ACGT]*\M)'
);

-- TABLES STARTING WITH A PREFIX
CREATE OR REPLACE FUNCTION index_tables_for(
	preffix_t_name information_schema.sql_identifier
) RETURNS TABLE (
	"schema" information_schema.sql_identifier, 
	"tbl_name" information_schema.sql_identifier
) AS
$BODY$
BEGIN
	RETURN QUERY (
        SELECT
            "table_schema", "table_name"
        FROM
            information_schema.tables
        WHERE
            table_type = 'BASE TABLE'
        -- AND
        --    table_schema = changed_object.schema_name
        AND
            table_name LIKE (preffix_t_name || '%' || '__index')
	);
END;
$BODY$
LANGUAGE plpgsql;


CREATE OR REPLACE FUNCTION show_domain_usage(
	d_name text
) RETURNS SETOF information_schema.column_domain_usage AS
$BODY$
BEGIN
	RETURN QUERY (
		SELECT 
			*
		FROM 
			information_schema.column_domain_usage
		WHERE
			"domain_name" = d_name
	);
END;
$BODY$
LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION create_index_table(
	t_name information_schema.sql_identifier, 
	c_name information_schema.sql_identifier
) RETURNS void AS
$BODY$
BEGIN
	
	-- There could be a foreign key, but if someone wanted to delete 
	-- their table they would have to use 'CASCADE'
	-- FOREIGN KEY (dna_sequence_id) REFERENCES %I(id) 
	EXECUTE format(
	    'CREATE TABLE %I (
	        kmer	DNA_SEQUENCE NOT NULL,
	        starting_position INT NOT NULL, 
	        dna_sequence_id BIGINT NOT NULL
		)',
	   	(t_name || '_' || c_name || '__index'), 
	   	(t_name)
	);
END;
$BODY$
LANGUAGE plpgsql;

-- INSERT OCCURRANCE
CREATE OR REPLACE FUNCTION insert_index(
	t_name information_schema.sql_identifier, 
	c_name information_schema.sql_identifier, 

	kmer DNA_SEQUENCE,
	starting_position INT, 
	dna_sequence_id BIGINT
) RETURNS void AS
$BODY$
BEGIN 
	EXECUTE FORMAT('
		INSERT INTO %I("kmer", "starting_position", "dna_sequence_id") VALUES(%L, %L, %L)', 
	   (t_name || '_' || c_name || '__index'),
	   kmer, starting_position, dna_sequence_id
	);
END;
$BODY$
LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION insert_indexes(
	t_name information_schema.sql_identifier, 
	c_name information_schema.sql_identifier, 
	
	dna_sequence_id BIGINT,
	dna_sequence DNA_SEQUENCE 
) RETURNS void AS
$BODY$
DECLARE
	ch text;
	kmer text := '';
	ending_position SMALLINT := -1;

BEGIN
    FOR ch IN SELECT regexp_split_to_table(dna_sequence, '') LOOP
		ending_position := (ending_position + 1);
		kmer := (kmer || ch);
		
	 	IF ending_position >= 7 THEN 
			kmer := substr(kmer, 2); 
	  	END IF;
	
	  	IF ending_position >= 6 THEN 	  		
			EXECUTE insert_index(t_name, c_name, kmer, ending_position - 6, dna_sequence_id);
	  	END IF;

    END LOOP;
    RETURN;
END;
$BODY$
LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION update_index_trigger() RETURNS TRIGGER AS
$BODY$
DECLARE 
	c_name information_schema.sql_identifier;

	dna_sequence text;
BEGIN
	FOR c_name IN SELECT column_name FROM show_domain_usage('dna_sequence') WHERE table_name = TG_TABLE_NAME 
	LOOP 
		EXECUTE FORMAT('DELETE FROM %I WHERE "dna_sequence_id" = %L', 
			(TG_TABLE_NAME || '_' || c_name || '__index'), OLD.id);

		EXECUTE FORMAT ('SELECT %I FROM %I WHERE id = %L', 
			c_name, TG_TABLE_NAME, NEW.id) INTO dna_sequence;

		EXECUTE FORMAT ('SELECT insert_indexes(%L, %L, %s, %L)', 
			TG_TABLE_NAME, c_name, NEW.id, dna_sequence);
		
	END LOOP;

	RETURN NEW;
END;
$BODY$
LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION insert_index_trigger() RETURNS TRIGGER AS
$BODY$
DECLARE 
	c_name information_schema.sql_identifier;

	dna_sequence text;
BEGIN
	FOR c_name IN SELECT column_name FROM show_domain_usage('dna_sequence') WHERE table_name = TG_TABLE_NAME 
	LOOP 
		EXECUTE FORMAT ('SELECT %I FROM %I WHERE id = %L', 
			c_name, TG_TABLE_NAME, NEW.id) INTO dna_sequence;

		EXECUTE FORMAT ('SELECT insert_indexes(%L, %L, %s, %L)', 
			TG_TABLE_NAME, c_name, NEW.id, dna_sequence);
		
	END LOOP;

	RETURN NEW;
END;
$BODY$
LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION delete_index_trigger() RETURNS TRIGGER AS
$BODY$
DECLARE 
	c_name information_schema.sql_identifier;

	dna_sequence text;
BEGIN
	FOR c_name IN SELECT column_name FROM show_domain_usage('dna_sequence') WHERE table_name = TG_TABLE_NAME 
	LOOP 
		
		EXECUTE FORMAT('
			DELETE FROM %I WHERE "dna_sequence_id" = %L
		', 
		(TG_TABLE_NAME || '_' || c_name || '__index'), 
		OLD.id);

	END LOOP;

	RETURN NEW;
END;
$BODY$
LANGUAGE plpgsql;

-- CREATE create_mmseqs_index_tables_after_change
CREATE OR REPLACE FUNCTION create_mmseqs_index_tables_after_change()
        RETURNS event_trigger LANGUAGE plpgsql AS $$
DECLARE
    row record;
    changed_object record;

	t_name information_schema.sql_identifier; 
	c_name information_schema.sql_identifier;

	sequence_id BIGINT;
	dna_sequence DNA_SEQUENCE;
BEGIN
    FOR changed_object IN SELECT * FROM pg_event_trigger_ddl_commands()
    LOOP		
		-- Creating new tables with indexes for given table
    	IF changed_object.object_type = 'table' 
		AND changed_object.object_identity NOT LIKE '%__index' 
		THEN
			SELECT SPLIT_PART(changed_object.object_identity, '.', 2) INTO t_name;
			
			-- Deleting ALL previous tables
		    FOR row IN SELECT * FROM index_tables_for(SPLIT_PART(changed_object.object_identity, '.', 2))
		    LOOP EXECUTE format('DROP TABLE %I CASCADE', row.tbl_name); END LOOP;
			
			-- Creating new tables
			FOR c_name IN SELECT column_name FROM show_domain_usage('dna_sequence') 
			WHERE "table_schema" || '.' || "table_name" = changed_object.object_identity
			LOOP 
				
			   	EXECUTE create_index_table(SPLIT_PART(changed_object.object_identity, '.', 2), c_name);
												
				FOR sequence_id, dna_sequence IN EXECUTE format('SELECT id, %I as clmn from %I', c_name, t_name) LOOP
					EXECUTE insert_indexes(t_name, c_name, sequence_id, dna_sequence);
				END LOOP;
				
			END LOOP;

			EXECUTE FORMAT('
				DROP TRIGGER IF EXISTS %I
				ON %I;

				CREATE TRIGGER %I
			    AFTER UPDATE ON %I
			    FOR EACH ROW
			    WHEN (OLD IS DISTINCT FROM NEW)
			    EXECUTE PROCEDURE update_index_trigger()
			', 
			t_name || '__update_trigger', 
			t_name,
			t_name || '__update_trigger', 
			t_name
			);

			EXECUTE FORMAT('
				DROP TRIGGER IF EXISTS %I
				ON %I;

				CREATE TRIGGER %I
			    AFTER INSERT ON %I
			    FOR EACH ROW
			    EXECUTE PROCEDURE insert_index_trigger()
			', 
			t_name || '__insert_trigger', 
			t_name,
			t_name || '__insert_trigger', 
			t_name
			);

			EXECUTE FORMAT('
				DROP TRIGGER IF EXISTS %I
				ON %I;

				CREATE TRIGGER %I
			    AFTER DELETE ON %I
			    FOR EACH ROW
			    EXECUTE PROCEDURE delete_index_trigger()
			', 
			t_name || '__delete_trigger', 
			t_name,
			t_name || '__delete_trigger', 
			t_name
			);

		END IF;
    END LOOP;
END;
$$;

					
CREATE EVENT TRIGGER 
	indexise_after_change_trigger 
ON 
	ddl_command_end
EXECUTE PROCEDURE 
	create_mmseqs_index_tables_after_change();

CREATE OR REPLACE FUNCTION test_event_trigger_for_drops()
        RETURNS event_trigger LANGUAGE plpgsql AS $$
DECLARE
	row record;
    changed_object record;

	t_name information_schema.sql_identifier; 
BEGIN
    FOR changed_object IN SELECT * FROM pg_event_trigger_dropped_objects()
    LOOP
    	IF changed_object.object_type = 'table' AND changed_object.object_name NOT LIKE '%__index' THEN	                   
			SELECT SPLIT_PART(changed_object.object_identity, '.', 2) INTO t_name;
  
		    FOR row IN SELECT * FROM index_tables_for(SPLIT_PART(changed_object.object_identity, '.', 2))
		    LOOP EXECUTE format('DROP TABLE %I CASCADE', row.tbl_name); END LOOP;
	           
			EXECUTE FORMAT('DROP TRIGGER IF EXISTS %I ON %I', t_name || '__update_trigger', t_name);
			EXECUTE FORMAT('DROP TRIGGER IF EXISTS %I ON %I', t_name || '__insert_trigger', t_name);
			EXECUTE FORMAT('DROP TRIGGER IF EXISTS %I ON %I', t_name || '__delete_trigger', t_name);
		END IF;
    END LOOP;
END;
$$;

CREATE EVENT TRIGGER 
	test_event_trigger_for_drops_tr
ON 
	sql_drop
EXECUTE FUNCTION 
	test_event_trigger_for_drops();


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
    AS '/mmseq2/extension.so', 'seq_search_mmseqs'
    LANGUAGE C STRICT;

-- CREATE TABLE taco(id bigint, sequence dna_sequence);
-- INSERT INTO taco VALUES (1, 'DDDDDAAGGGGG');

-- CREATE TABLE que(id bigint, sequence dna_sequence);
-- INSERT INTO que VALUES (1, 'DDDDDDDCCGGGGGAA');
-- INSERT INTO que VALUES (2, 'AAADDDDDCCGGGGGDD');
