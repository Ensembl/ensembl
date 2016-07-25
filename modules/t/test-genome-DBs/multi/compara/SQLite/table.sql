-- 
-- Created by SQL::Translator::Producer::SQLite
-- Created on Thu Jul 21 10:40:50 2016
-- 

BEGIN TRANSACTION;

--
-- Table: CAFE_gene_family
--
CREATE TABLE CAFE_gene_family (
  cafe_gene_family_id INTEGER PRIMARY KEY NOT NULL,
  root_id int(10) NOT NULL,
  lca_id int(10) NOT NULL,
  gene_tree_root_id int(10) NOT NULL,
  pvalue_avg double(5,4) DEFAULT NULL,
  lambdas varchar(100) DEFAULT NULL
);

CREATE INDEX lca_id ON CAFE_gene_family (lca_id);

CREATE INDEX root_id ON CAFE_gene_family (root_id);

CREATE INDEX gene_tree_root_id ON CAFE_gene_family (gene_tree_root_id);

--
-- Table: CAFE_species_gene
--
CREATE TABLE CAFE_species_gene (
  cafe_gene_family_id int(10) NOT NULL,
  node_id int(10) NOT NULL,
  n_members int(4) NOT NULL,
  pvalue double(5,4) DEFAULT NULL,
  PRIMARY KEY (cafe_gene_family_id, node_id)
);

CREATE INDEX node_id ON CAFE_species_gene (node_id);

--
-- Table: conservation_score
--
CREATE TABLE conservation_score (
  genomic_align_block_id bigint(20) NOT NULL,
  window_size smallint(5) NOT NULL,
  position int(10) NOT NULL,
  expected_score blob,
  diff_score blob
);

CREATE INDEX genomic_align_block_id ON conservation_score (genomic_align_block_id, window_size);

--
-- Table: constrained_element
--
CREATE TABLE constrained_element (
  constrained_element_id bigint(20) NOT NULL,
  dnafrag_id bigint(20) NOT NULL,
  dnafrag_start int(12) NOT NULL,
  dnafrag_end int(12) NOT NULL,
  dnafrag_strand int(2) NOT NULL,
  method_link_species_set_id int(10) NOT NULL,
  p_value double(8,2) DEFAULT NULL,
  score double(8,2) NOT NULL DEFAULT 0
);

CREATE INDEX dnafrag_id ON constrained_element (dnafrag_id);

CREATE INDEX constrained_element_id_idx ON constrained_element (constrained_element_id);

CREATE INDEX mlssid_idx ON constrained_element (method_link_species_set_id);

CREATE INDEX mlssid_dfId_dfStart_dfEnd_idx ON constrained_element (method_link_species_set_id, dnafrag_id, dnafrag_start, dnafrag_end);

CREATE INDEX mlssid_dfId_idx ON constrained_element (method_link_species_set_id, dnafrag_id);

--
-- Table: dnafrag
--
CREATE TABLE dnafrag (
  dnafrag_id INTEGER PRIMARY KEY NOT NULL,
  length int(11) NOT NULL DEFAULT 0,
  name varchar(40) NOT NULL DEFAULT '',
  genome_db_id int(10) NOT NULL,
  coord_system_name varchar(40) NOT NULL DEFAULT '',
  is_reference tinyint(1) NOT NULL DEFAULT 1
);

CREATE UNIQUE INDEX name ON dnafrag (genome_db_id, name);

--
-- Table: dnafrag_region
--
CREATE TABLE dnafrag_region (
  synteny_region_id int(10) NOT NULL DEFAULT 0,
  dnafrag_id bigint(20) NOT NULL DEFAULT 0,
  dnafrag_start int(10) NOT NULL DEFAULT 0,
  dnafrag_end int(10) NOT NULL DEFAULT 0,
  dnafrag_strand tinyint(4) NOT NULL DEFAULT 0
);

CREATE INDEX synteny ON dnafrag_region (synteny_region_id, dnafrag_id);

CREATE INDEX synteny_reversed ON dnafrag_region (dnafrag_id, synteny_region_id);

--
-- Table: external_db
--
CREATE TABLE external_db (
  external_db_id INTEGER PRIMARY KEY NOT NULL,
  db_name varchar(100) NOT NULL,
  db_release varchar(255) DEFAULT NULL,
  status enum(9) NOT NULL,
  priority int(11) NOT NULL,
  db_display_name varchar(255) DEFAULT NULL,
  type enum(18) DEFAULT NULL,
  secondary_db_name varchar(255) DEFAULT NULL,
  secondary_db_table varchar(255) DEFAULT NULL,
  description text
);

CREATE UNIQUE INDEX db_name_db_release_idx ON external_db (db_name, db_release);

--
-- Table: family
--
CREATE TABLE family (
  family_id INTEGER PRIMARY KEY NOT NULL,
  stable_id varchar(40) NOT NULL,
  version int(10) NOT NULL,
  method_link_species_set_id int(10) NOT NULL,
  description text,
  description_score double(8,2) DEFAULT NULL
);

CREATE INDEX method_link_species_set_id ON family (method_link_species_set_id);

CREATE INDEX description ON family (description);

CREATE UNIQUE INDEX stable_id ON family (stable_id);

--
-- Table: family_member
--
CREATE TABLE family_member (
  family_id int(10) NOT NULL,
  seq_member_id int(10) NOT NULL,
  cigar_line mediumtext(16777215),
  PRIMARY KEY (family_id, seq_member_id)
);

CREATE INDEX family_id ON family_member (family_id);

CREATE INDEX seq_member_id ON family_member (seq_member_id);

--
-- Table: gene_align
--
CREATE TABLE gene_align (
  gene_align_id INTEGER PRIMARY KEY NOT NULL,
  seq_type varchar(40) DEFAULT NULL,
  aln_method varchar(40) NOT NULL DEFAULT '',
  aln_length int(10) NOT NULL DEFAULT 0
);

--
-- Table: gene_align_member
--
CREATE TABLE gene_align_member (
  gene_align_id int(10) NOT NULL,
  seq_member_id int(10) NOT NULL,
  cigar_line mediumtext(16777215),
  PRIMARY KEY (gene_align_id, seq_member_id)
);

CREATE INDEX seq_member_id02 ON gene_align_member (seq_member_id);

--
-- Table: gene_member
--
CREATE TABLE gene_member (
  gene_member_id INTEGER PRIMARY KEY NOT NULL,
  stable_id varchar(128) NOT NULL,
  version int(10) DEFAULT 0,
  source_name enum(12) NOT NULL,
  taxon_id int(10) NOT NULL,
  genome_db_id int(10) DEFAULT NULL,
  canonical_member_id int(10) DEFAULT NULL,
  description text,
  dnafrag_id bigint(20) DEFAULT NULL,
  dnafrag_start int(10) DEFAULT NULL,
  dnafrag_end int(10) DEFAULT NULL,
  dnafrag_strand tinyint(4) DEFAULT NULL,
  display_label varchar(128) DEFAULT NULL,
  families tinyint(1) DEFAULT 0,
  gene_trees tinyint(1) DEFAULT 0,
  gene_gain_loss_trees tinyint(1) DEFAULT 0,
  orthologues int(10) DEFAULT 0,
  paralogues int(10) DEFAULT 0,
  homoeologues int(10) DEFAULT 0
);

CREATE INDEX taxon_id ON gene_member (taxon_id);

CREATE INDEX genome_db_id ON gene_member (genome_db_id);

CREATE INDEX source_name ON gene_member (source_name);

CREATE INDEX canonical_member_id ON gene_member (canonical_member_id);

CREATE INDEX dnafrag_id_start ON gene_member (dnafrag_id, dnafrag_start);

CREATE INDEX dnafrag_id_end ON gene_member (dnafrag_id, dnafrag_end);

CREATE UNIQUE INDEX stable_id02 ON gene_member (stable_id);

--
-- Table: gene_tree_node
--
CREATE TABLE gene_tree_node (
  node_id INTEGER PRIMARY KEY NOT NULL,
  parent_id int(10) DEFAULT NULL,
  root_id int(10) DEFAULT NULL,
  left_index int(10) NOT NULL DEFAULT 0,
  right_index int(10) NOT NULL DEFAULT 0,
  distance_to_parent double(8,2) NOT NULL DEFAULT 1,
  seq_member_id int(10) DEFAULT NULL
);

CREATE INDEX parent_id ON gene_tree_node (parent_id);

CREATE INDEX seq_member_id03 ON gene_tree_node (seq_member_id);

CREATE INDEX root_id02 ON gene_tree_node (root_id);

CREATE INDEX root_id_left_index ON gene_tree_node (root_id, left_index);

--
-- Table: gene_tree_node_attr
--
CREATE TABLE gene_tree_node_attr (
  node_id INTEGER PRIMARY KEY NOT NULL,
  node_type enum(11) DEFAULT NULL,
  species_tree_node_id int(10) DEFAULT NULL,
  bootstrap tinyint(3) DEFAULT NULL,
  duplication_confidence_score double(5,4) DEFAULT NULL
);

CREATE INDEX species_tree_node_id ON gene_tree_node_attr (species_tree_node_id);

--
-- Table: gene_tree_node_tag
--
CREATE TABLE gene_tree_node_tag (
  node_id int(10) NOT NULL,
  tag varchar(50) NOT NULL,
  value mediumtext(16777215) NOT NULL
);

CREATE INDEX node_id_tag ON gene_tree_node_tag (node_id, tag);

CREATE INDEX tag ON gene_tree_node_tag (tag);

--
-- Table: gene_tree_object_store
--
CREATE TABLE gene_tree_object_store (
  root_id int(10) NOT NULL,
  data_label varchar(255) NOT NULL,
  compressed_data mediumblob(16777215) NOT NULL,
  PRIMARY KEY (root_id, data_label)
);

--
-- Table: gene_tree_root
--
CREATE TABLE gene_tree_root (
  root_id INTEGER PRIMARY KEY NOT NULL,
  member_type enum(7) NOT NULL,
  tree_type enum(10) NOT NULL,
  clusterset_id varchar(20) NOT NULL DEFAULT 'default',
  method_link_species_set_id int(10) NOT NULL,
  species_tree_root_id int(10) DEFAULT NULL,
  gene_align_id int(10) DEFAULT NULL,
  ref_root_id int(10) DEFAULT NULL,
  stable_id varchar(40) DEFAULT NULL,
  version int(10) DEFAULT NULL
);

CREATE INDEX method_link_species_set_id02 ON gene_tree_root (method_link_species_set_id);

CREATE INDEX gene_align_id ON gene_tree_root (gene_align_id);

CREATE INDEX ref_root_id ON gene_tree_root (ref_root_id);

CREATE INDEX tree_type ON gene_tree_root (tree_type);

CREATE UNIQUE INDEX stable_id03 ON gene_tree_root (stable_id);

--
-- Table: gene_tree_root_attr
--
CREATE TABLE gene_tree_root_attr (
  root_id INTEGER PRIMARY KEY NOT NULL,
  aln_after_filter_length int(10) DEFAULT NULL,
  aln_length int(10) DEFAULT NULL,
  aln_num_residues int(10) DEFAULT NULL,
  aln_percent_identity float(8,2) DEFAULT NULL,
  best_fit_model_family varchar(10) DEFAULT NULL,
  best_fit_model_parameter varchar(5) DEFAULT NULL,
  gene_count int(10) DEFAULT NULL,
  k_score float(8,2) DEFAULT NULL,
  k_score_rank int(10) DEFAULT NULL,
  mcoffee_scores_gene_align_id int(10) DEFAULT NULL,
  aln_n_removed_columns int(10) DEFAULT NULL,
  aln_num_of_patterns int(10) DEFAULT NULL,
  aln_shrinking_factor float(8,2) DEFAULT NULL,
  spec_count int(10) DEFAULT NULL,
  tree_max_branch float(8,2) DEFAULT NULL,
  tree_max_length float(8,2) DEFAULT NULL,
  tree_num_dup_nodes int(10) DEFAULT NULL,
  tree_num_leaves int(10) DEFAULT NULL,
  tree_num_spec_nodes int(10) DEFAULT NULL,
  lca_node_id int(10) DEFAULT NULL,
  taxonomic_coverage float(8,2) DEFAULT NULL,
  ratio_species_genes float(8,2) DEFAULT NULL,
  model_name varchar(40) DEFAULT NULL,
  division varchar(10) DEFAULT NULL
);

CREATE INDEX lca_node_id ON gene_tree_root_attr (lca_node_id);

--
-- Table: gene_tree_root_tag
--
CREATE TABLE gene_tree_root_tag (
  root_id int(10) NOT NULL,
  tag varchar(50) NOT NULL,
  value mediumtext(16777215) NOT NULL
);

CREATE INDEX root_id_tag ON gene_tree_root_tag (root_id, tag);

CREATE INDEX root_id03 ON gene_tree_root_tag (root_id);

CREATE INDEX tag02 ON gene_tree_root_tag (tag);

--
-- Table: genome_db
--
CREATE TABLE genome_db (
  genome_db_id INTEGER PRIMARY KEY NOT NULL,
  taxon_id int(10) DEFAULT NULL,
  name varchar(128) NOT NULL DEFAULT '',
  assembly varchar(100) NOT NULL DEFAULT '',
  genebuild varchar(100) NOT NULL DEFAULT '',
  has_karyotype tinyint(1) NOT NULL DEFAULT 0,
  is_high_coverage tinyint(1) NOT NULL DEFAULT 0,
  genome_component varchar(5) DEFAULT NULL,
  locator varchar(400) DEFAULT NULL,
  first_release smallint(6) DEFAULT NULL,
  last_release smallint(6) DEFAULT NULL
);

CREATE INDEX taxon_id02 ON genome_db (taxon_id);

CREATE UNIQUE INDEX name02 ON genome_db (name, assembly, genome_component);

--
-- Table: genomic_align
--
CREATE TABLE genomic_align (
  genomic_align_id INTEGER PRIMARY KEY NOT NULL,
  genomic_align_block_id bigint(20) NOT NULL,
  method_link_species_set_id int(10) NOT NULL DEFAULT 0,
  dnafrag_id bigint(20) NOT NULL DEFAULT 0,
  dnafrag_start int(10) NOT NULL DEFAULT 0,
  dnafrag_end int(10) NOT NULL DEFAULT 0,
  dnafrag_strand tinyint(4) NOT NULL DEFAULT 0,
  cigar_line mediumtext(16777215) NOT NULL,
  visible tinyint(2) NOT NULL DEFAULT 1,
  node_id bigint(20) DEFAULT NULL
);

CREATE INDEX genomic_align_block_id02 ON genomic_align (genomic_align_block_id);

CREATE INDEX method_link_species_set_id03 ON genomic_align (method_link_species_set_id);

CREATE INDEX dnafrag02 ON genomic_align (dnafrag_id, method_link_species_set_id, dnafrag_start, dnafrag_end);

CREATE INDEX node_id02 ON genomic_align (node_id);

--
-- Table: genomic_align_block
--
CREATE TABLE genomic_align_block (
  genomic_align_block_id INTEGER PRIMARY KEY NOT NULL,
  method_link_species_set_id int(10) NOT NULL DEFAULT 0,
  score double(8,2) DEFAULT NULL,
  perc_id tinyint(3) DEFAULT NULL,
  length int(10) NOT NULL,
  group_id bigint(20) DEFAULT NULL,
  level_id tinyint(2) NOT NULL DEFAULT 0
);

CREATE INDEX method_link_species_set_id04 ON genomic_align_block (method_link_species_set_id);

--
-- Table: genomic_align_tree
--
CREATE TABLE genomic_align_tree (
  node_id INTEGER PRIMARY KEY NOT NULL,
  parent_id bigint(20) NOT NULL DEFAULT 0,
  root_id bigint(20) NOT NULL DEFAULT 0,
  left_index int(10) NOT NULL DEFAULT 0,
  right_index int(10) NOT NULL DEFAULT 0,
  left_node_id bigint(10) NOT NULL DEFAULT 0,
  right_node_id bigint(10) NOT NULL DEFAULT 0,
  distance_to_parent double(8,2) NOT NULL DEFAULT 1
);

CREATE INDEX parent_id02 ON genomic_align_tree (parent_id);

CREATE INDEX root_id04 ON genomic_align_tree (root_id);

CREATE INDEX left_index ON genomic_align_tree (root_id, left_index);

--
-- Table: hmm_annot
--
CREATE TABLE hmm_annot (
  seq_member_id INTEGER PRIMARY KEY NOT NULL,
  model_id varchar(40) DEFAULT NULL,
  evalue float(8,2) DEFAULT NULL
);

CREATE INDEX model_id ON hmm_annot (model_id);

--
-- Table: hmm_curated_annot
--
CREATE TABLE hmm_curated_annot (
  seq_member_stable_id varchar(40) NOT NULL,
  model_id varchar(40) DEFAULT NULL,
  library_version varchar(40) NOT NULL,
  annot_date timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  reason mediumtext(16777215),
  PRIMARY KEY (seq_member_stable_id)
);

CREATE INDEX model_id02 ON hmm_curated_annot (model_id);

--
-- Table: hmm_profile
--
CREATE TABLE hmm_profile (
  model_id varchar(40) NOT NULL,
  name varchar(40) DEFAULT NULL,
  type varchar(40) NOT NULL,
  compressed_profile mediumblob(16777215),
  consensus mediumtext(16777215),
  PRIMARY KEY (model_id, type)
);

--
-- Table: homology
--
CREATE TABLE homology (
  homology_id INTEGER PRIMARY KEY NOT NULL,
  method_link_species_set_id int(10) NOT NULL,
  description enum(23) DEFAULT NULL,
  is_tree_compliant tinyint(1) NOT NULL DEFAULT 0,
  dn float(10,5) DEFAULT NULL,
  ds float(10,5) DEFAULT NULL,
  n float(10,1) DEFAULT NULL,
  s float(10,1) DEFAULT NULL,
  lnl float(10,3) DEFAULT NULL,
  species_tree_node_id int(10) DEFAULT NULL,
  gene_tree_node_id int(10) DEFAULT NULL,
  gene_tree_root_id int(10) DEFAULT NULL,
  goc_score tinyint(3) DEFAULT NULL,
  wga_coverage decimal(5,2) DEFAULT NULL
);

CREATE INDEX method_link_species_set_id05 ON homology (method_link_species_set_id);

CREATE INDEX species_tree_node_id02 ON homology (species_tree_node_id);

CREATE INDEX gene_tree_node_id ON homology (gene_tree_node_id);

CREATE INDEX gene_tree_root_id02 ON homology (gene_tree_root_id);

--
-- Table: homology_member
--
CREATE TABLE homology_member (
  homology_id int(10) NOT NULL,
  gene_member_id int(10) NOT NULL,
  seq_member_id int(10) DEFAULT NULL,
  cigar_line mediumtext(16777215),
  perc_cov float(8,2) DEFAULT 0,
  perc_id float(8,2) DEFAULT 0,
  perc_pos float(8,2) DEFAULT 0,
  PRIMARY KEY (homology_id, gene_member_id)
);

CREATE INDEX homology_id ON homology_member (homology_id);

CREATE INDEX gene_member_id ON homology_member (gene_member_id);

CREATE INDEX seq_member_id04 ON homology_member (seq_member_id);

--
-- Table: mapping_session
--
CREATE TABLE mapping_session (
  mapping_session_id INTEGER PRIMARY KEY NOT NULL,
  type enum(6) DEFAULT NULL,
  when_mapped timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  rel_from int(10) DEFAULT NULL,
  rel_to int(10) DEFAULT NULL,
  prefix char(4) NOT NULL
);

CREATE UNIQUE INDEX type ON mapping_session (type, rel_from, rel_to, prefix);

--
-- Table: member_xref
--
CREATE TABLE member_xref (
  gene_member_id int(10) NOT NULL,
  dbprimary_acc varchar(10) NOT NULL,
  external_db_id int(10) NOT NULL,
  PRIMARY KEY (gene_member_id, dbprimary_acc, external_db_id)
);

CREATE INDEX external_db_id ON member_xref (external_db_id);

--
-- Table: meta
--
CREATE TABLE meta (
  meta_id INTEGER PRIMARY KEY NOT NULL,
  species_id int(10) DEFAULT 1,
  meta_key varchar(40) NOT NULL,
  meta_value text NOT NULL
);

CREATE INDEX species_value_idx ON meta (species_id, meta_value);

CREATE UNIQUE INDEX species_key_value_idx ON meta (species_id, meta_key, meta_value);

--
-- Table: method_link
--
CREATE TABLE method_link (
  method_link_id INTEGER PRIMARY KEY NOT NULL,
  type varchar(50) NOT NULL DEFAULT '',
  class varchar(50) NOT NULL DEFAULT ''
);

CREATE UNIQUE INDEX type02 ON method_link (type);

--
-- Table: method_link_species_set
--
CREATE TABLE method_link_species_set (
  method_link_species_set_id INTEGER PRIMARY KEY NOT NULL,
  method_link_id int(10) NOT NULL,
  species_set_id int(10) NOT NULL,
  name varchar(255) NOT NULL DEFAULT '',
  source varchar(255) NOT NULL DEFAULT 'ensembl',
  url varchar(255) NOT NULL DEFAULT '',
  first_release smallint(6) DEFAULT NULL,
  last_release smallint(6) DEFAULT NULL
);

CREATE INDEX species_set_id ON method_link_species_set (species_set_id);

CREATE UNIQUE INDEX method_link_id ON method_link_species_set (method_link_id, species_set_id);

--
-- Table: method_link_species_set_attr
--
CREATE TABLE method_link_species_set_attr (
  method_link_species_set_id INTEGER PRIMARY KEY NOT NULL,
  n_goc_null int(11) DEFAULT NULL,
  n_goc_0 int(11) DEFAULT NULL,
  n_goc_25 int(11) DEFAULT NULL,
  n_goc_50 int(11) DEFAULT NULL,
  n_goc_75 int(11) DEFAULT NULL,
  n_goc_100 int(11) DEFAULT NULL,
  perc_orth_above_goc_thresh float(8,2) DEFAULT NULL,
  goc_quality_threshold int(11) DEFAULT NULL,
  wga_quality_threshold int(11) DEFAULT NULL,
  perc_orth_above_wga_thresh float(8,2) DEFAULT NULL,
  threshold_on_ds int(11) DEFAULT NULL
);

--
-- Table: method_link_species_set_tag
--
CREATE TABLE method_link_species_set_tag (
  method_link_species_set_id int(10) NOT NULL,
  tag varchar(50) NOT NULL,
  value mediumtext(16777215) NOT NULL,
  PRIMARY KEY (method_link_species_set_id, tag)
);

CREATE INDEX tag03 ON method_link_species_set_tag (tag);

--
-- Table: ncbi_taxa_name
--
CREATE TABLE ncbi_taxa_name (
  taxon_id int(10) NOT NULL,
  name varchar(255) NOT NULL,
  name_class varchar(50) NOT NULL
);

CREATE INDEX taxon_id03 ON ncbi_taxa_name (taxon_id);

CREATE INDEX name03 ON ncbi_taxa_name (name);

CREATE INDEX name_class ON ncbi_taxa_name (name_class);

--
-- Table: ncbi_taxa_node
--
CREATE TABLE ncbi_taxa_node (
  taxon_id INTEGER PRIMARY KEY NOT NULL,
  parent_id int(10) NOT NULL,
  rank char(32) NOT NULL DEFAULT '',
  genbank_hidden_flag tinyint(1) NOT NULL DEFAULT 0,
  left_index int(10) NOT NULL DEFAULT 0,
  right_index int(10) NOT NULL DEFAULT 0,
  root_id int(10) NOT NULL DEFAULT 1
);

CREATE INDEX parent_id03 ON ncbi_taxa_node (parent_id);

CREATE INDEX rank ON ncbi_taxa_node (rank);

CREATE INDEX left_index02 ON ncbi_taxa_node (left_index);

CREATE INDEX right_index ON ncbi_taxa_node (right_index);

--
-- Table: other_member_sequence
--
CREATE TABLE other_member_sequence (
  seq_member_id int(10) NOT NULL,
  seq_type varchar(40) NOT NULL,
  length int(10) NOT NULL,
  sequence mediumtext(16777215) NOT NULL,
  PRIMARY KEY (seq_member_id, seq_type)
);

--
-- Table: peptide_align_feature
--
CREATE TABLE peptide_align_feature (
  peptide_align_feature_id INTEGER PRIMARY KEY NOT NULL,
  qmember_id int(10) NOT NULL,
  hmember_id int(10) NOT NULL,
  qgenome_db_id int(10) DEFAULT NULL,
  hgenome_db_id int(10) DEFAULT NULL,
  qstart int(10) NOT NULL DEFAULT 0,
  qend int(10) NOT NULL DEFAULT 0,
  hstart int(11) NOT NULL DEFAULT 0,
  hend int(11) NOT NULL DEFAULT 0,
  score double(16,4) NOT NULL DEFAULT 0.0000,
  evalue double(8,2) NOT NULL,
  align_length int(10) NOT NULL,
  identical_matches int(10) NOT NULL,
  perc_ident int(10) NOT NULL,
  positive_matches int(10) NOT NULL,
  perc_pos int(10) NOT NULL,
  hit_rank int(10) NOT NULL,
  cigar_line mediumtext(16777215)
);

--
-- Table: seq_member
--
CREATE TABLE seq_member (
  seq_member_id INTEGER PRIMARY KEY NOT NULL,
  stable_id varchar(128) NOT NULL,
  version int(10) DEFAULT 0,
  source_name enum(17) NOT NULL,
  taxon_id int(10) NOT NULL,
  genome_db_id int(10) DEFAULT NULL,
  sequence_id int(10) DEFAULT NULL,
  gene_member_id int(10) DEFAULT NULL,
  description text,
  dnafrag_id bigint(20) DEFAULT NULL,
  dnafrag_start int(10) DEFAULT NULL,
  dnafrag_end int(10) DEFAULT NULL,
  dnafrag_strand tinyint(4) DEFAULT NULL,
  display_label varchar(128) DEFAULT NULL
);

CREATE INDEX taxon_id04 ON seq_member (taxon_id);

CREATE INDEX genome_db_id02 ON seq_member (genome_db_id);

CREATE INDEX source_name02 ON seq_member (source_name);

CREATE INDEX sequence_id ON seq_member (sequence_id);

CREATE INDEX gene_member_id02 ON seq_member (gene_member_id);

CREATE INDEX dnafrag_id_start02 ON seq_member (dnafrag_id, dnafrag_start);

CREATE INDEX dnafrag_id_end02 ON seq_member (dnafrag_id, dnafrag_end);

CREATE INDEX seq_member_gene_member_id_end ON seq_member (seq_member_id, gene_member_id);

CREATE UNIQUE INDEX stable_id04 ON seq_member (stable_id);

--
-- Table: sequence
--
CREATE TABLE sequence (
  sequence_id INTEGER PRIMARY KEY NOT NULL,
  length int(10) NOT NULL,
  md5sum char(32) NOT NULL,
  sequence longtext(4294967295) NOT NULL
);

CREATE INDEX md5sum ON sequence (md5sum);

--
-- Table: species_set
--
CREATE TABLE species_set (
  species_set_id int(10) NOT NULL,
  genome_db_id int(10) NOT NULL,
  PRIMARY KEY (species_set_id, genome_db_id)
);

CREATE INDEX genome_db_id03 ON species_set (genome_db_id);

--
-- Table: species_set_header
--
CREATE TABLE species_set_header (
  species_set_id INTEGER PRIMARY KEY NOT NULL,
  name varchar(255) NOT NULL DEFAULT '',
  size int(10) NOT NULL,
  first_release smallint(6) DEFAULT NULL,
  last_release smallint(6) DEFAULT NULL
);

--
-- Table: species_set_tag
--
CREATE TABLE species_set_tag (
  species_set_id int(10) NOT NULL,
  tag varchar(50) NOT NULL,
  value mediumtext(16777215) NOT NULL,
  PRIMARY KEY (species_set_id, tag)
);

CREATE INDEX tag04 ON species_set_tag (tag);

--
-- Table: species_tree_node
--
CREATE TABLE species_tree_node (
  node_id INTEGER PRIMARY KEY NOT NULL,
  parent_id int(10) DEFAULT NULL,
  root_id int(10) DEFAULT NULL,
  left_index int(10) NOT NULL DEFAULT 0,
  right_index int(10) NOT NULL DEFAULT 0,
  distance_to_parent double(8,2) DEFAULT 1,
  taxon_id int(10) DEFAULT NULL,
  genome_db_id int(10) DEFAULT NULL,
  node_name varchar(255) DEFAULT NULL
);

CREATE INDEX taxon_id05 ON species_tree_node (taxon_id);

CREATE INDEX genome_db_id04 ON species_tree_node (genome_db_id);

CREATE INDEX parent_id04 ON species_tree_node (parent_id);

CREATE INDEX root_id05 ON species_tree_node (root_id, left_index);

--
-- Table: species_tree_node_attr
--
CREATE TABLE species_tree_node_attr (
  node_id INTEGER PRIMARY KEY NOT NULL,
  nb_long_genes int(11) DEFAULT NULL,
  nb_short_genes int(11) DEFAULT NULL,
  avg_dupscore float(8,2) DEFAULT NULL,
  avg_dupscore_nondub float(8,2) DEFAULT NULL,
  nb_dubious_nodes int(11) DEFAULT NULL,
  nb_dup_nodes int(11) DEFAULT NULL,
  nb_genes int(11) DEFAULT NULL,
  nb_genes_in_tree int(11) DEFAULT NULL,
  nb_genes_in_tree_multi_species int(11) DEFAULT NULL,
  nb_genes_in_tree_single_species int(11) DEFAULT NULL,
  nb_nodes int(11) DEFAULT NULL,
  nb_orphan_genes int(11) DEFAULT NULL,
  nb_seq int(11) DEFAULT NULL,
  nb_spec_nodes int(11) DEFAULT NULL,
  nb_gene_splits int(11) DEFAULT NULL,
  nb_split_genes int(11) DEFAULT NULL,
  root_avg_gene float(8,2) DEFAULT NULL,
  root_avg_gene_per_spec float(8,2) DEFAULT NULL,
  root_avg_spec float(8,2) DEFAULT NULL,
  root_max_gene int(11) DEFAULT NULL,
  root_max_spec int(11) DEFAULT NULL,
  root_min_gene int(11) DEFAULT NULL,
  root_min_spec int(11) DEFAULT NULL,
  root_nb_genes int(11) DEFAULT NULL,
  root_nb_trees int(11) DEFAULT NULL
);

--
-- Table: species_tree_node_tag
--
CREATE TABLE species_tree_node_tag (
  node_id int(10) NOT NULL,
  tag varchar(50) NOT NULL,
  value mediumtext(16777215) NOT NULL
);

CREATE INDEX node_id_tag02 ON species_tree_node_tag (node_id, tag);

CREATE INDEX tag05 ON species_tree_node_tag (tag);

--
-- Table: species_tree_root
--
CREATE TABLE species_tree_root (
  root_id INTEGER PRIMARY KEY NOT NULL,
  method_link_species_set_id int(10) NOT NULL,
  label varchar(256) NOT NULL DEFAULT 'default'
);

CREATE UNIQUE INDEX method_link_species_set_id06 ON species_tree_root (method_link_species_set_id, label);

--
-- Table: stable_id_history
--
CREATE TABLE stable_id_history (
  mapping_session_id int(10) NOT NULL,
  stable_id_from varchar(40) NOT NULL DEFAULT '',
  version_from int(10) DEFAULT NULL,
  stable_id_to varchar(40) NOT NULL DEFAULT '',
  version_to int(10) DEFAULT NULL,
  contribution float(8,2) DEFAULT NULL,
  PRIMARY KEY (mapping_session_id, stable_id_from, stable_id_to)
);

--
-- Table: synteny_region
--
CREATE TABLE synteny_region (
  synteny_region_id INTEGER PRIMARY KEY NOT NULL,
  method_link_species_set_id int(10) NOT NULL
);

CREATE INDEX method_link_species_set_id07 ON synteny_region (method_link_species_set_id);

COMMIT;
