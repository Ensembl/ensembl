-- 
-- Created by SQL::Translator::Producer::SQLite
-- Created on Fri Jan 12 13:38:11 2018
-- 

BEGIN TRANSACTION;

--
-- Table: CAFE_gene_family
--
CREATE TABLE CAFE_gene_family (
  cafe_gene_family_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  root_id integer NOT NULL,
  lca_id integer NOT NULL,
  gene_tree_root_id integer NOT NULL,
  pvalue_avg double precision(5,4),
  lambdas varchar(100)
);

--
-- Table: CAFE_species_gene
--
CREATE TABLE CAFE_species_gene (
  cafe_gene_family_id integer NOT NULL,
  node_id integer NOT NULL,
  n_members integer NOT NULL,
  pvalue double precision(5,4),
  PRIMARY KEY (cafe_gene_family_id, node_id)
);

--
-- Table: conservation_score
--
CREATE TABLE conservation_score (
  genomic_align_block_id bigint NOT NULL,
  window_size smallint NOT NULL,
  position integer NOT NULL,
  expected_score blob,
  diff_score blob
);

--
-- Table: constrained_element
--
CREATE TABLE constrained_element (
  constrained_element_id bigint NOT NULL,
  dnafrag_id bigint NOT NULL,
  dnafrag_start integer NOT NULL,
  dnafrag_end integer NOT NULL,
  dnafrag_strand integer NOT NULL,
  method_link_species_set_id integer NOT NULL,
  p_value double precision NOT NULL DEFAULT 0,
  score double precision NOT NULL DEFAULT 0
);

--
-- Table: dnafrag
--
CREATE TABLE dnafrag (
  dnafrag_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  length integer NOT NULL DEFAULT 0,
  name varchar(255) NOT NULL DEFAULT '',
  genome_db_id integer NOT NULL,
  coord_system_name varchar(40) NOT NULL DEFAULT '',
  cellular_component enum NOT NULL DEFAULT 'NUC',
  is_reference tinyint NOT NULL DEFAULT 1,
  codon_table_id tinyint NOT NULL DEFAULT 1
);

CREATE UNIQUE INDEX name ON dnafrag (genome_db_id, name);

--
-- Table: dnafrag_region
--
CREATE TABLE dnafrag_region (
  synteny_region_id integer NOT NULL DEFAULT 0,
  dnafrag_id bigint NOT NULL DEFAULT 0,
  dnafrag_start integer NOT NULL DEFAULT 0,
  dnafrag_end integer NOT NULL DEFAULT 0,
  dnafrag_strand tinyint NOT NULL DEFAULT 0
);

--
-- Table: exon_boundaries
--
CREATE TABLE exon_boundaries (
  gene_member_id integer NOT NULL,
  seq_member_id integer NOT NULL,
  dnafrag_start integer NOT NULL,
  dnafrag_end integer NOT NULL,
  sequence_length integer NOT NULL,
  left_over tinyint NOT NULL DEFAULT 0
);

--
-- Table: external_db
--
CREATE TABLE external_db (
  external_db_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  db_name varchar(100) NOT NULL,
  db_release varchar(255),
  status enum NOT NULL,
  priority integer NOT NULL,
  db_display_name varchar(255),
  type enum,
  secondary_db_name varchar(255),
  secondary_db_table varchar(255),
  description text
);

CREATE UNIQUE INDEX db_name_db_release_idx ON external_db (db_name, db_release);

--
-- Table: family
--
CREATE TABLE family (
  family_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  stable_id varchar(40) NOT NULL,
  version integer NOT NULL,
  method_link_species_set_id integer NOT NULL,
  description text,
  description_score double precision
);

CREATE UNIQUE INDEX stable_id ON family (stable_id);

--
-- Table: family_member
--
CREATE TABLE family_member (
  family_id integer NOT NULL,
  seq_member_id integer NOT NULL,
  cigar_line mediumtext,
  PRIMARY KEY (family_id, seq_member_id)
);

--
-- Table: gene_align
--
CREATE TABLE gene_align (
  gene_align_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  seq_type varchar(40),
  aln_method varchar(40) NOT NULL DEFAULT '',
  aln_length integer NOT NULL DEFAULT 0
);

--
-- Table: gene_align_member
--
CREATE TABLE gene_align_member (
  gene_align_id integer NOT NULL,
  seq_member_id integer NOT NULL,
  cigar_line mediumtext,
  PRIMARY KEY (gene_align_id, seq_member_id)
);

--
-- Table: gene_member
--
CREATE TABLE gene_member (
  gene_member_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  stable_id varchar(128) NOT NULL,
  version integer DEFAULT 0,
  source_name enum NOT NULL,
  taxon_id integer NOT NULL,
  genome_db_id integer,
  biotype_group enum NOT NULL DEFAULT 'coding',
  canonical_member_id integer,
  description text,
  dnafrag_id bigint,
  dnafrag_start integer,
  dnafrag_end integer,
  dnafrag_strand tinyint,
  display_label varchar(128)
);

CREATE UNIQUE INDEX stable_id02 ON gene_member (stable_id);

--
-- Table: gene_member_hom_stats
--
CREATE TABLE gene_member_hom_stats (
  gene_member_id integer NOT NULL,
  collection varchar(40) NOT NULL,
  families integer NOT NULL DEFAULT 0,
  gene_trees tinyint NOT NULL DEFAULT 0,
  gene_gain_loss_trees tinyint NOT NULL DEFAULT 0,
  orthologues integer NOT NULL DEFAULT 0,
  paralogues integer NOT NULL DEFAULT 0,
  homoeologues integer NOT NULL DEFAULT 0,
  PRIMARY KEY (gene_member_id, collection)
);

--
-- Table: gene_member_qc
--
CREATE TABLE gene_member_qc (
  gene_member_stable_id varchar(128) NOT NULL,
  genome_db_id integer NOT NULL,
  seq_member_id integer,
  n_species integer,
  n_orth integer,
  avg_cov float,
  status varchar(50) NOT NULL
);

--
-- Table: gene_tree_node
--
CREATE TABLE gene_tree_node (
  node_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  parent_id integer,
  root_id integer,
  left_index integer NOT NULL DEFAULT 0,
  right_index integer NOT NULL DEFAULT 0,
  distance_to_parent double precision NOT NULL DEFAULT 1,
  seq_member_id integer
);

--
-- Table: gene_tree_node_attr
--
CREATE TABLE gene_tree_node_attr (
  node_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  node_type enum,
  species_tree_node_id integer,
  bootstrap tinyint,
  duplication_confidence_score double precision(5,4)
);

--
-- Table: gene_tree_node_tag
--
CREATE TABLE gene_tree_node_tag (
  node_id integer NOT NULL,
  tag varchar(50) NOT NULL,
  value mediumtext NOT NULL
);

--
-- Table: gene_tree_object_store
--
CREATE TABLE gene_tree_object_store (
  root_id integer NOT NULL,
  data_label varchar(255) NOT NULL,
  compressed_data mediumblob NOT NULL,
  PRIMARY KEY (root_id, data_label)
);

--
-- Table: gene_tree_root
--
CREATE TABLE gene_tree_root (
  root_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  member_type enum NOT NULL,
  tree_type enum NOT NULL,
  clusterset_id varchar(20) NOT NULL DEFAULT 'default',
  method_link_species_set_id integer NOT NULL,
  species_tree_root_id integer,
  gene_align_id integer,
  ref_root_id integer,
  stable_id varchar(40),
  version integer
);

CREATE UNIQUE INDEX stable_id03 ON gene_tree_root (stable_id);

--
-- Table: gene_tree_root_attr
--
CREATE TABLE gene_tree_root_attr (
  root_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  aln_after_filter_length integer,
  aln_length integer,
  aln_num_residues integer,
  aln_percent_identity float,
  best_fit_model_family varchar(10),
  best_fit_model_parameter varchar(5),
  gene_count integer,
  k_score float,
  k_score_rank integer,
  mcoffee_scores_gene_align_id integer,
  aln_n_removed_columns integer,
  aln_num_of_patterns integer,
  aln_shrinking_factor float,
  spec_count integer,
  tree_max_branch float,
  tree_max_length float,
  tree_num_dup_nodes integer,
  tree_num_leaves integer,
  tree_num_spec_nodes integer,
  lca_node_id integer,
  taxonomic_coverage float,
  ratio_species_genes float,
  model_name varchar(40),
  division varchar(10)
);

--
-- Table: gene_tree_root_tag
--
CREATE TABLE gene_tree_root_tag (
  root_id integer NOT NULL,
  tag varchar(50) NOT NULL,
  value mediumtext NOT NULL
);

--
-- Table: genome_db
--
CREATE TABLE genome_db (
  genome_db_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  taxon_id integer,
  name varchar(128) NOT NULL DEFAULT '',
  assembly varchar(100) NOT NULL DEFAULT '',
  genebuild varchar(100) NOT NULL DEFAULT '',
  has_karyotype tinyint NOT NULL DEFAULT 0,
  is_high_coverage tinyint NOT NULL DEFAULT 0,
  genome_component varchar(5),
  strain_name varchar(40),
  display_name varchar(255),
  locator varchar(400),
  first_release smallint,
  last_release smallint
);

CREATE UNIQUE INDEX name02 ON genome_db (name, assembly, genome_component);

--
-- Table: genomic_align
--
CREATE TABLE genomic_align (
  genomic_align_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  genomic_align_block_id bigint NOT NULL,
  method_link_species_set_id integer NOT NULL DEFAULT 0,
  dnafrag_id bigint NOT NULL DEFAULT 0,
  dnafrag_start integer NOT NULL DEFAULT 0,
  dnafrag_end integer NOT NULL DEFAULT 0,
  dnafrag_strand tinyint NOT NULL DEFAULT 0,
  cigar_line mediumtext NOT NULL,
  visible tinyint NOT NULL DEFAULT 1,
  node_id bigint
);

--
-- Table: genomic_align_block
--
CREATE TABLE genomic_align_block (
  genomic_align_block_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  method_link_species_set_id integer NOT NULL DEFAULT 0,
  score double precision,
  perc_id tinyint,
  length integer NOT NULL,
  group_id bigint,
  level_id tinyint NOT NULL DEFAULT 0
);

--
-- Table: genomic_align_tree
--
CREATE TABLE genomic_align_tree (
  node_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  parent_id bigint,
  root_id bigint NOT NULL DEFAULT 0,
  left_index integer NOT NULL DEFAULT 0,
  right_index integer NOT NULL DEFAULT 0,
  left_node_id bigint,
  right_node_id bigint,
  distance_to_parent double precision NOT NULL DEFAULT 1
);

--
-- Table: hmm_annot
--
CREATE TABLE hmm_annot (
  seq_member_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  model_id varchar(40),
  evalue float
);

--
-- Table: hmm_curated_annot
--
CREATE TABLE hmm_curated_annot (
  seq_member_stable_id varchar(40) NOT NULL,
  model_id varchar(40),
  library_version varchar(40) NOT NULL,
  annot_date timestamp NOT NULL DEFAULT current_timestamp,
  reason mediumtext,
  PRIMARY KEY (seq_member_stable_id)
);

--
-- Table: hmm_profile
--
CREATE TABLE hmm_profile (
  model_id varchar(40) NOT NULL,
  name varchar(40),
  type varchar(40) NOT NULL,
  compressed_profile mediumblob,
  consensus mediumtext,
  PRIMARY KEY (model_id, type)
);

--
-- Table: homology
--
CREATE TABLE homology (
  homology_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  method_link_species_set_id integer NOT NULL,
  description enum,
  is_tree_compliant tinyint NOT NULL DEFAULT 0,
  dn float(10,5),
  ds float(10,5),
  n float(10,1),
  s float(10,1),
  lnl float(10,3),
  species_tree_node_id integer,
  gene_tree_node_id integer,
  gene_tree_root_id integer,
  goc_score tinyint,
  wga_coverage decimal(5,2),
  is_high_confidence tinyint
);

--
-- Table: homology_member
--
CREATE TABLE homology_member (
  homology_id integer NOT NULL,
  gene_member_id integer NOT NULL,
  seq_member_id integer,
  cigar_line mediumtext,
  perc_cov float DEFAULT 0,
  perc_id float DEFAULT 0,
  perc_pos float DEFAULT 0,
  PRIMARY KEY (homology_id, gene_member_id)
);

--
-- Table: mapping_session
--
CREATE TABLE mapping_session (
  mapping_session_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  type enum,
  when_mapped timestamp NOT NULL DEFAULT current_timestamp,
  rel_from integer,
  rel_to integer,
  prefix char(4) NOT NULL
);

CREATE UNIQUE INDEX type ON mapping_session (type, rel_from, rel_to, prefix);

--
-- Table: member_xref
--
CREATE TABLE member_xref (
  gene_member_id integer NOT NULL,
  dbprimary_acc varchar(10) NOT NULL,
  external_db_id integer NOT NULL,
  PRIMARY KEY (gene_member_id, dbprimary_acc, external_db_id)
);

--
-- Table: meta
--
CREATE TABLE meta (
  meta_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  species_id integer DEFAULT 1,
  meta_key varchar(40) NOT NULL,
  meta_value text NOT NULL
);

CREATE UNIQUE INDEX species_key_value_idx ON meta (species_id, meta_key, meta_value);

--
-- Table: method_link
--
CREATE TABLE method_link (
  method_link_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  type varchar(50) NOT NULL DEFAULT '',
  class varchar(50) NOT NULL DEFAULT ''
);

CREATE UNIQUE INDEX type02 ON method_link (type);

--
-- Table: method_link_species_set
--
CREATE TABLE method_link_species_set (
  method_link_species_set_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  method_link_id integer NOT NULL,
  species_set_id integer NOT NULL,
  name varchar(255) NOT NULL DEFAULT '',
  source varchar(255) NOT NULL DEFAULT 'ensembl',
  url varchar(255) NOT NULL DEFAULT '',
  first_release smallint,
  last_release smallint
);

CREATE UNIQUE INDEX method_link_id ON method_link_species_set (method_link_id, species_set_id);

--
-- Table: method_link_species_set_attr
--
CREATE TABLE method_link_species_set_attr (
  method_link_species_set_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  n_goc_null integer,
  n_goc_0 integer,
  n_goc_25 integer,
  n_goc_50 integer,
  n_goc_75 integer,
  n_goc_100 integer,
  perc_orth_above_goc_thresh float,
  goc_quality_threshold integer,
  wga_quality_threshold integer,
  perc_orth_above_wga_thresh float,
  threshold_on_ds integer
);

--
-- Table: method_link_species_set_tag
--
CREATE TABLE method_link_species_set_tag (
  method_link_species_set_id integer NOT NULL,
  tag varchar(50) NOT NULL,
  value mediumtext NOT NULL,
  PRIMARY KEY (method_link_species_set_id, tag)
);

--
-- Table: ncbi_taxa_name
--
CREATE TABLE ncbi_taxa_name (
  taxon_id integer NOT NULL,
  name varchar(255) NOT NULL,
  name_class varchar(50) NOT NULL
);

--
-- Table: ncbi_taxa_node
--
CREATE TABLE ncbi_taxa_node (
  taxon_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  parent_id integer NOT NULL,
  rank char(32) NOT NULL DEFAULT '',
  genbank_hidden_flag tinyint NOT NULL DEFAULT 0,
  left_index integer NOT NULL DEFAULT 0,
  right_index integer NOT NULL DEFAULT 0,
  root_id integer NOT NULL DEFAULT 1
);

--
-- Table: other_member_sequence
--
CREATE TABLE other_member_sequence (
  seq_member_id integer NOT NULL,
  seq_type varchar(40) NOT NULL,
  length integer NOT NULL,
  sequence mediumtext NOT NULL,
  PRIMARY KEY (seq_member_id, seq_type)
);

--
-- Table: peptide_align_feature
--
CREATE TABLE peptide_align_feature (
  peptide_align_feature_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  qmember_id integer NOT NULL,
  hmember_id integer NOT NULL,
  qgenome_db_id integer,
  hgenome_db_id integer,
  qstart integer NOT NULL DEFAULT 0,
  qend integer NOT NULL DEFAULT 0,
  hstart integer NOT NULL DEFAULT 0,
  hend integer NOT NULL DEFAULT 0,
  score double precision(16,4) NOT NULL DEFAULT 0.0000,
  evalue double precision NOT NULL,
  align_length integer NOT NULL,
  identical_matches integer NOT NULL,
  perc_ident integer NOT NULL,
  positive_matches integer NOT NULL,
  perc_pos integer NOT NULL,
  hit_rank integer NOT NULL,
  cigar_line mediumtext
);

--
-- Table: seq_member
--
CREATE TABLE seq_member (
  seq_member_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  stable_id varchar(128) NOT NULL,
  version integer DEFAULT 0,
  source_name enum NOT NULL,
  taxon_id integer NOT NULL,
  genome_db_id integer,
  sequence_id integer,
  gene_member_id integer,
  has_transcript_edits tinyint NOT NULL DEFAULT 0,
  has_translation_edits tinyint NOT NULL DEFAULT 0,
  description text,
  dnafrag_id bigint,
  dnafrag_start integer,
  dnafrag_end integer,
  dnafrag_strand tinyint,
  display_label varchar(128)
);

CREATE UNIQUE INDEX stable_id04 ON seq_member (stable_id);

--
-- Table: seq_member_projection
--
CREATE TABLE seq_member_projection (
  target_seq_member_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  source_seq_member_id integer NOT NULL,
  identity float(5,2) NOT NULL
);

--
-- Table: seq_member_projection_stable_id
--
CREATE TABLE seq_member_projection_stable_id (
  target_seq_member_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  source_stable_id varchar(128) NOT NULL
);

--
-- Table: sequence
--
CREATE TABLE sequence (
  sequence_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  length integer NOT NULL,
  md5sum char(32) NOT NULL,
  sequence longtext NOT NULL
);

--
-- Table: species_set
--
CREATE TABLE species_set (
  species_set_id integer NOT NULL,
  genome_db_id integer NOT NULL,
  PRIMARY KEY (species_set_id, genome_db_id)
);

--
-- Table: species_set_header
--
CREATE TABLE species_set_header (
  species_set_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  name varchar(255) NOT NULL DEFAULT '',
  size integer NOT NULL,
  first_release smallint,
  last_release smallint
);

--
-- Table: species_set_tag
--
CREATE TABLE species_set_tag (
  species_set_id integer NOT NULL,
  tag varchar(50) NOT NULL,
  value mediumtext NOT NULL,
  PRIMARY KEY (species_set_id, tag)
);

--
-- Table: species_tree_node
--
CREATE TABLE species_tree_node (
  node_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  parent_id integer,
  root_id integer,
  left_index integer NOT NULL DEFAULT 0,
  right_index integer NOT NULL DEFAULT 0,
  distance_to_parent double precision DEFAULT 1,
  taxon_id integer,
  genome_db_id integer,
  node_name varchar(255)
);

--
-- Table: species_tree_node_attr
--
CREATE TABLE species_tree_node_attr (
  node_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  nb_long_genes integer,
  nb_short_genes integer,
  avg_dupscore float,
  avg_dupscore_nondub float,
  nb_dubious_nodes integer,
  nb_dup_nodes integer,
  nb_genes integer,
  nb_genes_in_tree integer,
  nb_genes_in_tree_multi_species integer,
  nb_genes_in_tree_single_species integer,
  nb_nodes integer,
  nb_orphan_genes integer,
  nb_seq integer,
  nb_spec_nodes integer,
  nb_gene_splits integer,
  nb_split_genes integer,
  root_avg_gene float,
  root_avg_gene_per_spec float,
  root_avg_spec float,
  root_max_gene integer,
  root_max_spec integer,
  root_min_gene integer,
  root_min_spec integer,
  root_nb_genes integer,
  root_nb_trees integer
);

--
-- Table: species_tree_node_tag
--
CREATE TABLE species_tree_node_tag (
  node_id integer NOT NULL,
  tag varchar(50) NOT NULL,
  value mediumtext NOT NULL
);

--
-- Table: species_tree_root
--
CREATE TABLE species_tree_root (
  root_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  method_link_species_set_id integer NOT NULL,
  label varchar(256) NOT NULL DEFAULT 'default'
);

CREATE UNIQUE INDEX method_link_species_set_id ON species_tree_root (method_link_species_set_id, label);

--
-- Table: stable_id_history
--
CREATE TABLE stable_id_history (
  mapping_session_id integer NOT NULL,
  stable_id_from varchar(40) NOT NULL DEFAULT '',
  version_from integer,
  stable_id_to varchar(40) NOT NULL DEFAULT '',
  version_to integer,
  contribution float,
  PRIMARY KEY (mapping_session_id, stable_id_from, stable_id_to)
);

--
-- Table: synteny_region
--
CREATE TABLE synteny_region (
  synteny_region_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  method_link_species_set_id integer NOT NULL
);

COMMIT;
