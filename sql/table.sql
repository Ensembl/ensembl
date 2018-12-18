-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2018] EMBL-European Bioinformatics Institute
--
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
--
--      http://www.apache.org/licenses/LICENSE-2.0
--
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

# Ensembl core table definitions
#

# Conventions:
#  - use lower case and underscores
#  - internal ids are integers named tablename_id
#  - same name is given in foreign key relations


/**
@header Assembly Tables
@colour #C70C09

*/


/**
@table assembly
@colour #C70C09
@desc The assembly table states, which parts of seq_regions are exactly equal. It enables to transform coordinates between seq_regions.
Typically this contains how chromosomes are made of contigs, clones out of contigs, and chromosomes out of supercontigs.
It allows you to artificially chunk chromosome sequence into smaller parts.

The data in this table defines the "static golden path", i.e. the best effort draft full genome sequence as determined by the UCSC or NCBI (depending which assembly you are using).
Each row represents a component, e.g. a contig,  (comp_seq_region_id, FK from seq_region table) at least part of which is present in the golden path.
The part of the component that is in the path is delimited by fields cmp_start and cmp_end (start < end), and the absolute position within the golden path chromosome (or other appropriate assembled structure) (asm_seq_region_id) is given by asm_start and asm_end.

@column asm_seq_region_id            Assembly sequence region id. Primary key, internal identifier. Foreign key references to the @link seq_region table.
@column cmp_seq_region_id            Component sequence region id. Foreign key references to the @link seq_region table.
@column asm_start                    Start absolute position within the golden path chromosome.
@column asm_end                      End absolute position within the golden path chromosome.
@column cmp_start                    Component start position within the golden path chromosome.
@column cmp_end                      Component start position within the golden path chromosome.
@column ori                          Orientation: 1 - sense; -1 - antisense.


@see seq_region
@see supercontigs

*/


CREATE TABLE assembly (

  asm_seq_region_id           INT(10) UNSIGNED NOT NULL,
  cmp_seq_region_id           INT(10) UNSIGNED NOT NULL,
  asm_start                   INT(10) NOT NULL,
  asm_end                     INT(10) NOT NULL,
  cmp_start                   INT(10) NOT NULL,
  cmp_end                     INT(10) NOT NULL,
  ori                         TINYINT  NOT NULL,

  KEY cmp_seq_region_idx (cmp_seq_region_id),
  KEY asm_seq_region_idx (asm_seq_region_id, asm_start),
  UNIQUE KEY all_idx (asm_seq_region_id, cmp_seq_region_id, asm_start, asm_end, cmp_start, cmp_end, ori)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table assembly_exception
@colour #C70C09
@desc Allows multiple sequence regions to point to the same sequence, analogous to a symbolic link in a filesystem pointing to the actual file.
This mechanism has been implemented specifically to support haplotypes and PARs, but may be useful for other similar structures in the future.

@column assembly_exception_id       Assembly exception sequence region id. Primary key, internal identifier.
@column seq_region_id               Sequence region id. Foreign key references to the @link seq_region table.
@column seq_region_start            Sequence start position.
@column seq_region_end              Sequence end position.
@column exc_type                    Exception type, e.g. PAR, HAP - haplotype.
@column exc_seq_region_id           Exception sequence region id. Foreign key references to the @link seq_region table.
@column exc_seq_region_start        Exception sequence start position.
@column exc_seq_region_end          Exception sequence end position.
@column ori                         Orientation: 1 - sense; -1 - antisense.


@see assembly

*/


CREATE TABLE assembly_exception (

  assembly_exception_id       INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  exc_type                    ENUM('HAP', 'PAR',
                                'PATCH_FIX', 'PATCH_NOVEL') NOT NULL,
  exc_seq_region_id           INT(10) UNSIGNED NOT NULL,
  exc_seq_region_start        INT(10) UNSIGNED NOT NULL,
  exc_seq_region_end          INT(10) UNSIGNED NOT NULL,
  ori                         INT NOT NULL,

  PRIMARY KEY (assembly_exception_id),
  KEY sr_idx (seq_region_id, seq_region_start),
  KEY ex_idx (exc_seq_region_id, exc_seq_region_start)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table coord_system
@colour #C70C09
@desc Stores information about the available co-ordinate systems for the species identified through the species_id field.
Note that for each species, there must be one co-ordinate system that has the attribute "top_level" and one that has the attribute "sequence_level".

@column coord_system_id      Primary key, internal identifier.
@column species_id           Indentifies the species for multi-species databases.
@column name                 Co-oridinate system name, e.g. 'chromosome', 'contig', 'scaffold' etc.
@column version              Assembly.
@column rank                 Co-oridinate system rank.
@column attrib               Co-oridinate system attrib (e.g. "top_level", "sequence_level").

@see seq_region
@see meta_coord
@see meta

*/


CREATE TABLE coord_system (

  coord_system_id             INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  species_id                  INT(10) UNSIGNED NOT NULL DEFAULT 1,
  name                        VARCHAR(40) NOT NULL,
  version                     VARCHAR(255) DEFAULT NULL,
  rank                        INT NOT NULL,
  attrib                      SET('default_version', 'sequence_level'),

  PRIMARY   KEY (coord_system_id),
  UNIQUE    KEY rank_idx (rank, species_id),
  UNIQUE    KEY name_idx (name, version, species_id),
            KEY species_idx (species_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table data_file
@colour #C70C09
@desc Allows the storage of flat file locations used to store large quanitities of data currently unsuitable in a traditional database table.

@column data_file_id      Auto-increment surrogate primary key
@column coord_system_id   Coordinate system this file is linked to. Used to decipher the assembly version it was mapped to
@column analysis_id       Analysis this file is linked to
@column name              Name of the file
@column version_lock      Indicates that this file is only compatible with the current Ensembl release version
@column absolute          Flags that the URL given is fully resolved and should be used without question
@column url               Optional path to the file (can be absolute or relative)
@column file_type         Type of file e.g. BAM, BIGBED, BIGWIG and VCF
*/

CREATE TABLE data_file (
  data_file_id      INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  coord_system_id   INT(10) UNSIGNED NOT NULL,
  analysis_id       SMALLINT UNSIGNED NOT NULL,
  name              VARCHAR(100) NOT NULL,
  version_lock      TINYINT(1) DEFAULT 0 NOT NULL,
  absolute          TINYINT(1) DEFAULT 0 NOT NULL,
  url               TEXT,
  file_type         ENUM('BAM','BAMCOV','BIGBED','BIGWIG','VCF'),

  PRIMARY KEY (data_file_id),
  UNIQUE KEY df_unq_idx (coord_system_id, analysis_id, name, file_type),
  KEY df_name_idx (name),
  KEY df_analysis_idx (analysis_id)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table dna
@colour #C70C09
@desc Contains DNA sequence. This table has a 1:1 relationship with the seq_region table.

@column seq_region_id           Primary key, internal identifier. Foreign key references to the @link seq_region table.
@column sequence                DNA sequence.

@see seq_region
*/


CREATE TABLE dna (

  seq_region_id       INT(10) UNSIGNED NOT NULL,
  sequence            LONGTEXT NOT NULL,

  PRIMARY KEY (seq_region_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM MAX_ROWS=750000 AVG_ROW_LENGTH=19000;



/**
@table genome_statistics
@colour #C70C09
@desc Contains genome and assembly related statistics
      These include but are not limited to: feature counts, sequence lengths

@column genome_statistics_id  Primary key, internal identifier.
@column statistic            Name of the statistics
@column value                 Corresponding value of the statistics (count/length)
@column species_id            Indentifies the species for multi-species databases.
@column attrib_type_id        To distinguish similar statistics for different cases
@column timestamp             Date the statistics was generated

*/


CREATE TABLE genome_statistics (

  genome_statistics_id     INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  statistic                VARCHAR(128) NOT NULL,
  value                    BIGINT(11) UNSIGNED DEFAULT '0' NOT NULL,
  species_id               INT UNSIGNED DEFAULT 1,
  attrib_type_id           INT(10) UNSIGNED DEFAULT NULL,
  timestamp                DATETIME DEFAULT NULL,

  PRIMARY KEY (genome_statistics_id),
  UNIQUE KEY stats_uniq (statistic, attrib_type_id, species_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;




/**
@table karyotype
@colour #C70C09
@desc Describes bands that can be stained on the chromosome.

@column karyotype_id            Primary key, internal identifier.
@column seq_region_id           Foreign key references to the @link seq_region table.
@column seq_region_start        Sequence start position.
@column seq_region_end          Sequence end position.
@column band                    Band.
@column stain                   Stain.

*/


CREATE TABLE karyotype (
  karyotype_id                INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  band                        VARCHAR(40) DEFAULT NULL,
  stain                       VARCHAR(40) DEFAULT NULL,

  PRIMARY KEY (karyotype_id),
  KEY region_band_idx (seq_region_id,band)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table meta
@colour #C70C09
@desc Stores data about the data in the current schema. Taxonomy information, version information and the default value for the type column in the assembly table are stored here.
Unlike other tables, data in the meta table is stored as key-value pairs. Also stores (via assembly.mapping keys) the relationships between co-ordinate systems in the assembly table.

The species_id field of the meta table is used in multi-species databases and makes it possible to have species-specific meta key-value pairs.
The species-specific meta key-value pairs needs to be repeated for each species_id.
Entries in the meta table that are not specific to any one species, such as the schema_version key and any other schema-related information must have their species_id field set to NULL.
The default species_id, and the only species_id value allowed in single-species databases, is 1.


@column meta_id                    Primary key, internal identifier.
@column species_id                 Indentifies the species for multi-species databases.
@column meta_key                   Name of the meta entry, e.g. "schema_version".
@column meta_value                 Corresponding value of the key, e.g. "61".

@see assembly
@see coord_system

*/


CREATE TABLE IF NOT EXISTS meta (

  meta_id                     INT NOT NULL AUTO_INCREMENT,
  species_id                  INT UNSIGNED DEFAULT 1,
  meta_key                    VARCHAR(40) NOT NULL,
  meta_value                  VARCHAR(255) NOT NULL,

  PRIMARY   KEY (meta_id),
  UNIQUE    KEY species_key_value_idx (species_id, meta_key, meta_value),
            KEY species_value_idx (species_id, meta_value)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


# Add schema type and schema version to the meta table.
INSERT INTO meta (species_id, meta_key, meta_value) VALUES
  (NULL, 'schema_type', 'core'),
  (NULL, 'schema_version', '96');

# Patches included in this schema file:
# NOTE: At start of release cycle, remove patch entries from last release.
# NOTE: Avoid line-breaks in values.
INSERT INTO meta (species_id, meta_key, meta_value)
  VALUES (NULL, 'patch', 'patch_95_96_a.sql|schema_version');


/**
@table meta_coord
@colour #C70C09
@desc Describes which co-ordinate systems the different feature tables use.

@column table_name              Ensembl database table name.
@column coord_system_id         Foreign key references to the @link coord_system table.
@column max_length              Longest sequence length.

@see coord_system

*/


CREATE TABLE meta_coord (

  table_name                  VARCHAR(40) NOT NULL,
  coord_system_id             INT(10) UNSIGNED NOT NULL,
  max_length                  INT,

  UNIQUE KEY cs_table_name_idx (coord_system_id, table_name)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table seq_region
@colour #C70C09
@desc Stores information about sequence regions. The primary key is used as a pointer into the dna table so that actual sequence can be obtained, and the coord_system_id allows sequence regions of multiple types to be stored.
Clones, contigs and chromosomes are all now stored in the seq_region table. Contigs are stored with the co-ordinate system 'contig'.
The relationship between contigs and clones is stored in the assembly table. The relationships between contigs and chromosomes, and between contigs and supercontigs, are stored in the assembly table.

@column seq_region_id             Primary key, internal identifier.
@column name                      Sequence region name.
@column coord_system_id           Foreign key references to the @link coord_system table.
@column length                    Sequence length.


@see dna
@see coord_system

*/


CREATE TABLE seq_region (

  seq_region_id               INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  name                        VARCHAR(255) NOT NULL,
  coord_system_id             INT(10) UNSIGNED NOT NULL,
  length                      INT(10) UNSIGNED NOT NULL,

  PRIMARY KEY (seq_region_id),
  UNIQUE KEY name_cs_idx (name, coord_system_id),
  KEY cs_idx (coord_system_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table seq_region_synonym
@colour #C70C09
@desc Allows for storing multiple names for sequence regions.

@column seq_region_synonym_id           Primary key, internal identifier.
@column seq_region_id                   Foreign key references to the @link seq_region table.
@column synonym                         Alternative name for sequence region.
@column external_db_id                  Foreign key references to the @link external_db table.

*/


CREATE TABLE seq_region_synonym (

  seq_region_synonym_id       INT UNSIGNED NOT NULL  AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  synonym                     VARCHAR(250) NOT NULL,
  external_db_id              INT UNSIGNED,

  PRIMARY KEY (seq_region_synonym_id),
  UNIQUE KEY syn_idx (synonym, seq_region_id),
  KEY seq_region_idx (seq_region_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table seq_region_attrib
@colour #C70C09
@desc Allows "attributes" to be defined for certain seq_regions. Provides a way of storing extra information about particular seq_regions without adding extra columns to the seq_region table. e.g.

@column seq_region_id       Foreign key references to the @link seq_region table.
@column attrib_type_id      Foreign key references to the @link attrib_type table.
@column value               Attribute value.

@see seq_region
@see attrib_type

*/


CREATE TABLE seq_region_attrib (

  seq_region_id               INT(10) UNSIGNED NOT NULL DEFAULT '0',
  attrib_type_id              SMALLINT(5) UNSIGNED NOT NULL DEFAULT '0',
  value                       TEXT NOT NULL,

  KEY type_val_idx (attrib_type_id, value(40)),
  KEY val_only_idx (value(40)),
  KEY seq_region_idx (seq_region_id),
  UNIQUE KEY region_attribx (seq_region_id, attrib_type_id, value(500))

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@header Fundamental Tables
@colour   #808000
*/


/**
@table alt_allele
@colour   #808000
@desc Stores information about genes on haplotypes that may be orthologous.

@column alt_allele_id          Primary key, internal identifier.
@column gene_id                Foreign key references to the @link gene table.
@column alt_allele_group_id    A group ID to show which alleles are related

*/

CREATE TABLE alt_allele (
        alt_allele_id INT UNSIGNED AUTO_INCREMENT,
        alt_allele_group_id INT UNSIGNED NOT NULL,
        gene_id INT UNSIGNED NOT NULL,

        PRIMARY KEY (alt_allele_id),
        UNIQUE KEY gene_idx (gene_id),
        KEY (gene_id,alt_allele_group_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

/**
@table alt_allele_attrib
@colour   #808000
@desc Holds all the different attributes assigned to individual alleles.

@column alt_allele_id           Primary key, internal identifier.
@column attrib                  Enum of attributes assigned to alternative alleles
*/


CREATE TABLE alt_allele_attrib (
         alt_allele_id INT UNSIGNED,
         attrib ENUM('IS_REPRESENTATIVE',
                     'IS_MOST_COMMON_ALLELE',
                     'IN_CORRECTED_ASSEMBLY',
                     'HAS_CODING_POTENTIAL',
                     'IN_ARTIFICIALLY_DUPLICATED_ASSEMBLY',
                     'IN_SYNTENIC_REGION',
                     'HAS_SAME_UNDERLYING_DNA_SEQUENCE',
                     'IN_BROKEN_ASSEMBLY_REGION',
                     'IS_VALID_ALTERNATE',
                     'SAME_AS_REPRESENTATIVE',
                     'SAME_AS_ANOTHER_ALLELE',
                     'MANUALLY_ASSIGNED',
                     'AUTOMATICALLY_ASSIGNED'),

        KEY aa_idx (alt_allele_id,attrib)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table alt_allele_group
@colour   #808000
@desc A minimal table used for tracking unique alt_allele_group_id's. MySQL does not allow multiple autoincrement fields. Further information about a group could be added here at a later date.

@column alt_allele_group_id     Primary key and only column.
*/

CREATE TABLE alt_allele_group (

         alt_allele_group_id INT UNSIGNED AUTO_INCREMENT,

         PRIMARY KEY (alt_allele_group_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table analysis
@colour   #808000
@desc Usually describes a program and some database that together are used to create a feature on a piece of sequence.
Each feature is marked with an analysis_id. The most important column is logic_name, which is used by the webteam to render a feature correctly on contigview (or even retrieve the right feature).
Logic_name is also used in the pipeline to identify the analysis which has to run in a given status of the pipeline.
The module column tells the pipeline which Perl module does the whole analysis, typically a RunnableDB module.

@column analysis_id                 Primary key, internal identifier.
@column created                     Date to distinguish newer and older versions off the same analysis.
@column logic_name                  String to identify the analysis. Used mainly inside pipeline.
@column db                          Database name.
@column db_version                  Database version.
@column db_file                     File system location of the database.
@column program                     The binary used to create a feature.
@column program_version             The binary version.
@column program_file                File system location of the binary.
@column parameters                  A parameter string which is processed by the perl module.
@column module                      Perl module names (RunnableDBS usually) executing this analysis.
@column module_version              Perl module version.
@column gff_source                  How to make a gff dump from features with this analysis.
@column gff_feature                 How to make a gff dump from features with this analysis.

@see analysis_description

*/

CREATE TABLE IF NOT EXISTS analysis (

  analysis_id                 SMALLINT UNSIGNED NOT NULL AUTO_INCREMENT,
  created                     datetime DEFAULT NULL,
  logic_name                  VARCHAR(128) NOT NULL,
  db                          VARCHAR(120),
  db_version                  VARCHAR(40),
  db_file                     VARCHAR(120),
  program                     VARCHAR(80),
  program_version             VARCHAR(40),
  program_file                VARCHAR(80),
  parameters                  TEXT,
  module                      VARCHAR(80),
  module_version              VARCHAR(40),
  gff_source                  VARCHAR(40),
  gff_feature                 VARCHAR(40),

  PRIMARY KEY (analysis_id),
  UNIQUE KEY logic_name_idx (logic_name)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table analysis_description
@colour   #808000
@desc Allows the storage of a textual description of the analysis, as well as a "display label", primarily for the EnsEMBL web site.

@column analysis_id            Primary key, internal identifier. Foreign key references to the @link analysis table.
@column description            Textual description of the analysis.
@column display_label          Display label for the EnsEMBL web site.
@column displayable            Flag indicating if the analysis description is to be displayed on the EnsEMBL web site.
@column web_data               Other data used by the EnsEMBL web site.

@see analysis

*/


CREATE TABLE IF NOT EXISTS analysis_description (

  analysis_id                  SMALLINT UNSIGNED NOT NULL,
  description                  TEXT,
  display_label                VARCHAR(255) NOT NULL,
  displayable                  TINYINT(1) NOT NULL DEFAULT 1,
  web_data                     TEXT,

  UNIQUE KEY analysis_idx (analysis_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table attrib_type
@colour   #808000
@desc Provides codes, names and desctriptions of attribute types.

@column attrib_type_id       Primary key, internal identifier.
@column code                 Attribute code, e.g. 'GapExons'.
@column name                 Attribute name, e.g. 'gap exons'.
@column description          Attribute description, e.g. 'number of gap exons'.

@see seq_region_attrib

*/


CREATE TABLE attrib_type (

  attrib_type_id              SMALLINT(5) UNSIGNED NOT NULL AUTO_INCREMENT,
  code                        VARCHAR(20) NOT NULL DEFAULT '',
  name                        VARCHAR(255) NOT NULL DEFAULT '',
  description                 TEXT,

  PRIMARY KEY (attrib_type_id),
  UNIQUE KEY code_idx (code)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table dna_align_feature
@colour   #808000
@desc Stores DNA sequence alignments generated from Blast (or Blast-like) comparisons.

@column dna_align_feature_id        Primary key, internal identifier.
@column seq_region_id               Foreign key references to the @link seq_region table.
@column seq_region_start            Sequence start position.
@column seq_region_end              Sequence end position.
@column seq_region_strand           Sequence region strand: 1 - forward; -1 - reverse.
@column hit_start                   Alignment hit start position.
@column hit_end                     Alignment hit end position.
@column hit_strand                  Alignment hit strand: 1 - forward; -1 - reverse.
@column hit_name                    Alignment hit name.
@column analysis_id                 Foreign key references to the @link analysis table.
@column score                       Alignment score.
@column evalue                      Alignment e-value.
@column perc_ident                  Alignment percentage identity.
@column cigar_line                  Used to encode gapped alignments.
@column external_db_id              Foreign key references to the @link external_db table.
@column hcoverage                   Hit coverage.
@column align_type                  Alignment string type used

@see cigar_line

*/

CREATE TABLE dna_align_feature (

  dna_align_feature_id        INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT(1) NOT NULL,
  hit_start                   INT NOT NULL,
  hit_end                     INT NOT NULL,
  hit_strand                  TINYINT(1) NOT NULL,
  hit_name                    VARCHAR(40) NOT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  score                       DOUBLE,
  evalue                      DOUBLE,
  perc_ident                  FLOAT,
  cigar_line                  TEXT,
  external_db_id              INT UNSIGNED,
  hcoverage                   DOUBLE,
  align_type                  ENUM('ensembl', 'cigar', 'vulgar', 'mdtag') DEFAULT 'ensembl',

  PRIMARY KEY (dna_align_feature_id),
  KEY seq_region_idx (seq_region_id, analysis_id, seq_region_start, score),
  KEY seq_region_idx_2 (seq_region_id, seq_region_start),
  KEY hit_idx (hit_name),
  KEY analysis_idx (analysis_id),
  KEY external_db_idx (external_db_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM MAX_ROWS=100000000 AVG_ROW_LENGTH=80;


/**
@table dna_align_feature_attrib
@colour   #808000
@desc Enables storage of attributes that relate to DNA sequence alignments.

@column dna_align_feature_id        Foreign key references to the @link dna_align_feature table.
@column attrib_type_id              Foreign key references to the @link attrib_type table.
@column value                       Attribute value.

@see dna_align_feature
*/

CREATE TABLE dna_align_feature_attrib (

  dna_align_feature_id        INT(10) UNSIGNED NOT NULL,
  attrib_type_id              SMALLINT(5) UNSIGNED NOT NULL,
  value                       TEXT NOT NULL,

  UNIQUE KEY dna_align_feature_attribx (dna_align_feature_id, attrib_type_id, value(500)),
  KEY dna_align_feature_idx (dna_align_feature_id),
  KEY type_val_idx (attrib_type_id, value(40)),
  KEY val_only_idx (value(40))

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table exon
@colour   #808000
@desc Stores data about exons. Associated with transcripts via exon_transcript. Allows access to contigs seq_regions.
Note seq_region_start is always less that seq_region_end, i.e. when the exon is on the other strand the seq_region_start is specifying the 3prime end of the exon.

@column exon_id                     Primary key, internal identifier.
@column seq_region_id               Foreign key references to the @link seq_region table.
@column seq_region_start            Sequence start position.
@column seq_region_end              Sequence end position.
@column seq_region_strand           Sequence region strand: 1 - forward; -1 - reverse.
@column phase                       The place where the intron lands inside the codon - 0 between codons, 1 between the 1st and second base, 2 between the second and 3rd base. Exons therefore have a start phase anda end phase, but introns have just one phase.
@column end_phase                   Usually, end_phase = (phase + exon_length)%3 but end_phase could be -1 if the exon is half-coding and its 3 prime end is UTR.
@column is_current                  1 - exon is current. Always set to 1 in ensembl dbs, but needed for otterlace dbs
@column is_constitutive             1 - exon is constitutive.
@column stable_id                   Release-independent stable identifier.
@column version                     Stable identifier version number.
@column created_date                Date created.
@column modified_date               Date modified.

@see exon_transcript

*/


CREATE TABLE exon (

  exon_id                     INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT(2) NOT NULL,

  phase                       TINYINT(2) NOT NULL,
  end_phase                   TINYINT(2) NOT NULL,

  is_current                  TINYINT(1) NOT NULL DEFAULT 1,
  is_constitutive             TINYINT(1) NOT NULL DEFAULT 0,

  stable_id                   VARCHAR(128) DEFAULT NULL,
  version                     SMALLINT UNSIGNED DEFAULT NULL,
  created_date                DATETIME DEFAULT NULL,
  modified_date               DATETIME DEFAULT NULL,

  PRIMARY KEY (exon_id),
  KEY seq_region_idx (seq_region_id, seq_region_start),
  KEY stable_id_idx (stable_id, version)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table exon_transcript
@colour   #808000
@desc Relationship table linking exons with transcripts. The rank column indicates the 5' to 3' position of the exon within the transcript, i.e. a rank of 1 means the exon is the 5' most within this transcript.

@column exon_id                Composite key. Foreign key references to the @link exon table.
@column transcript_id          Composite key. Foreign key references to the @link transcript table.
@column rank                   Composite key.

@see exon
@see transcript

*/

CREATE TABLE exon_transcript (

  exon_id                     INT(10) UNSIGNED NOT NULL,
  transcript_id               INT(10) UNSIGNED NOT NULL,
  rank                        INT(10) NOT NULL,

  PRIMARY KEY (exon_id,transcript_id,rank),
  KEY transcript (transcript_id),
  KEY exon (exon_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;



/**
@table gene
@colour   #808000
@desc Allows transcripts to be related to genes.

@column gene_id                     Primary key, internal identifier.
@column biotype                     Biotype, e.g. protein_coding.
@column analysis_id                 Foreign key references to the @link analysis table.
@column seq_region_id               Foreign key references to the @link seq_region table.
@column seq_region_start            Sequence start position.
@column seq_region_end              Sequence end position.
@column seq_region_strand           Sequence region strand: 1 - forward; -1 - reverse.
@column display_xref_id             External reference for EnsEMBL web site. Foreign key references to the @link xref table.
@column source                      e.g ensembl, havana etc.
@column description                 Gene description
@column is_current                  1 - gene is current. Always set to 1 in ensembl dbs, but needed for otterlace dbs
@column canonical_transcript_id     Foreign key references to the @link transcript table.
@column stable_id                   Release-independent stable identifier.
@column version                     Stable identifier version number.
@column created_date                Date created.
@column modified_date               Date modified.

@see transcript

*/


CREATE TABLE gene (

  gene_id                     INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  biotype                     VARCHAR(40) NOT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT(2) NOT NULL,
  display_xref_id             INT(10) UNSIGNED,
  source                      VARCHAR(40) NOT NULL,
  description                 TEXT,
  is_current                  TINYINT(1) NOT NULL DEFAULT 1,
  canonical_transcript_id     INT(10) UNSIGNED NOT NULL,
  stable_id                   VARCHAR(128) DEFAULT NULL,
  version                     SMALLINT UNSIGNED DEFAULT NULL,
  created_date                DATETIME DEFAULT NULL,
  modified_date               DATETIME DEFAULT NULL,

  PRIMARY KEY (gene_id),
  KEY seq_region_idx (seq_region_id, seq_region_start),
  KEY xref_id_index (display_xref_id),
  KEY analysis_idx (analysis_id),
  KEY stable_id_idx (stable_id, version),
  KEY canonical_transcript_id_idx (canonical_transcript_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table gene_attrib
@colour   #808000
@desc Enables storage of attributes that relate to genes.

@column gene_id             Foreign key references to the @link gene table.
@column attrib_type_id      Foreign key references to the @link attrib_type table.
@column value               Attribute value.


@see gene
*/


CREATE TABLE gene_attrib (

  gene_id                     INT(10) UNSIGNED NOT NULL DEFAULT '0',
  attrib_type_id              SMALLINT(5) UNSIGNED NOT NULL DEFAULT '0',
  value                       TEXT NOT NULL,

  KEY type_val_idx (attrib_type_id, value(40)),
  KEY val_only_idx (value(40)),
  KEY gene_idx (gene_id),
  UNIQUE KEY gene_attribx (gene_id, attrib_type_id, value(500))

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table protein_align_feature
@colour   #808000
@desc Stores translation alignments generated from Blast (or Blast-like) comparisons.

@column protein_align_feature_id    Primary key, internal identifier.
@column seq_region_id               Foreign key references to the @link seq_region table.
@column seq_region_start            Sequence start position.
@column seq_region_end              Sequence end position.
@column seq_region_strand           Sequence region strand: 1 - forward; -1 - reverse.
@column hit_start                   Alignment hit start position.
@column hit_end                     Alignment hit end position.
@column hit_name                    Alignment hit name.
@column analysis_id                 Foreign key references to the @link analysis table.
@column score                       Alignment score.
@column evalue                      Alignment e-value.
@column perc_ident                  Alignment percentage identity.
@column cigar_line                  Used to encode gapped alignments.
@column external_db_id              Foreign key references to the @link external_db table.
@column hcoverage                   Alignment hit coverage.
@column align_type                  Alignment string type used


@see cigar_line


*/


CREATE TABLE protein_align_feature (

  protein_align_feature_id    INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT(1) DEFAULT '1' NOT NULL,
  hit_start                   INT(10) NOT NULL,
  hit_end                     INT(10) NOT NULL,
  hit_name                    VARCHAR(40) NOT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  score                       DOUBLE,
  evalue                      DOUBLE,
  perc_ident                  FLOAT,
  cigar_line                  TEXT,
  external_db_id              INT UNSIGNED,
  hcoverage                   DOUBLE,
  align_type                  ENUM('ensembl', 'cigar', 'vulgar', 'mdtag') DEFAULT 'ensembl',

  PRIMARY KEY (protein_align_feature_id),
  KEY seq_region_idx (seq_region_id, analysis_id, seq_region_start, score),
  KEY seq_region_idx_2 (seq_region_id, seq_region_start),
  KEY hit_idx (hit_name),
  KEY analysis_idx (analysis_id),
  KEY external_db_idx (external_db_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM MAX_ROWS=100000000 AVG_ROW_LENGTH=80;


/**
@table protein_feature
@colour   #808000
@desc Describes features on the translations (as opposed to the DNA sequence itself), i.e. parts of the peptide. In peptide co-ordinates rather than contig co-ordinates.

@column protein_feature_id          Primary key, internal identifier.
@column translation_id              Foreign key references to the @link translation table.
@column seq_start                   Sequence start position.
@column seq_end                     Sequence end position.
@column hit_start                   Alignment hit start position.
@column hit_end                     Alignment hit end position.
@column hit_name                    Alignment hit name.
@column analysis_id                 Foreign key references to the @link analysis table.
@column score                       Alignment score.
@column evalue                      Alignment E-value.
@column perc_ident                  Alignment percentage identity.
@column external_data               External data for protein feature.
@column hit_description             Optional description of the hit. This can be a human readable name
@column cigar_line                  Used to encode gapped alignments.
@column align_type                  Alignment string type used

@see analysis

*/


CREATE TABLE protein_feature (

  protein_feature_id          INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  translation_id              INT(10) UNSIGNED NOT NULL,
  seq_start                   INT(10) NOT NULL,
  seq_end                     INT(10) NOT NULL,
  hit_start                   INT(10) NOT NULL,
  hit_end                     INT(10) NOT NULL,
  hit_name                    VARCHAR(40) CHARACTER SET latin1 COLLATE latin1_bin NOT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  score                       DOUBLE,
  evalue                      DOUBLE,
  perc_ident                  FLOAT,
  external_data               TEXT,
  hit_description             TEXT,
  cigar_line                  TEXT,
  align_type                  ENUM('ensembl', 'cigar', 'cigarplus', 'vulgar', 'mdtag') DEFAULT NULL,

  UNIQUE KEY aln_idx (translation_id,hit_name,seq_start,seq_end,hit_start,hit_end,analysis_id),
  PRIMARY KEY (protein_feature_id),
  KEY translation_idx (translation_id),
  KEY hitname_idx (hit_name),
  KEY analysis_idx (analysis_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table supporting_feature
@colour   #808000
@desc Describes the exon prediction process by linking exons to DNA or protein alignment features.
As in several other tables, the feature_id column is a foreign key; the feature_type column specifies which table feature_id refers to.

@column exon_id                    Foreign key references to the @link exon table.
@column feature_type               Feature type: 'dna_align_feature' or 'protein_align_feature'
@column feature_id                 Foreign key references to the @link dna_align_feature or @link protein_align_feature table depending on the feature type.


*/


CREATE TABLE supporting_feature (

  exon_id                     INT(10) UNSIGNED DEFAULT '0' NOT NULL,
  feature_type                ENUM('dna_align_feature','protein_align_feature'),
  feature_id                  INT(10) UNSIGNED DEFAULT '0' NOT NULL,

  UNIQUE KEY all_idx (exon_id,feature_type,feature_id),
  KEY feature_idx (feature_type,feature_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM MAX_ROWS=100000000 AVG_ROW_LENGTH=80;


/**
@table transcript
@colour   #808000
@desc Stores information about transcripts. Has seq_region_start, seq_region_end and seq_region_strand for faster retrieval and to allow storage independently of genes and exons.
Note that a transcript is usually associated with a translation, but may not be, e.g. in the case of pseudogenes and RNA genes (those that code for RNA molecules).

@column transcript_id               Primary key, internal identifier.
@column gene_id                     Foreign key references to the @link gene table.
@column analysis_id                 Foreign key references to the @link analysis table.
@column seq_region_id               Foreign key references to the @link seq_region table.
@column seq_region_start            Sequence start position.
@column seq_region_end              Sequence end position.
@column seq_region_strand           Sequence region strand: 1 - forward; -1 - reverse.
@column display_xref_id             External reference for EnsEMBL web site. Foreign key references to the @link xref table.
@column source                      e.g ensembl, havana etc.
@column biotype                     Biotype, e.g. protein_coding.
@column description                 Transcript description.
@column is_current                  Indicates a current transcript. Always set to 1 in ensembl dbs, but needed for otterlace dbs
@column canonical_translation_id    Foreign key references to the @link translation table.
@column stable_id                   Release-independent stable identifier.
@column version                     Stable identifier version number.
@column created_date                Date created.
@column modified_date               Date modified.

*/


CREATE TABLE transcript (

  transcript_id               INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  gene_id                     INT(10) UNSIGNED,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT(2) NOT NULL,
  display_xref_id             INT(10) UNSIGNED,
  source                      VARCHAR(40) NOT NULL default 'ensembl',
  biotype                     VARCHAR(40) NOT NULL,
  description                 TEXT,
  is_current                  TINYINT(1) NOT NULL DEFAULT 1,
  canonical_translation_id    INT(10) UNSIGNED,
  stable_id                   VARCHAR(128) DEFAULT NULL,
  version                     SMALLINT UNSIGNED DEFAULT NULL,
  created_date                DATETIME DEFAULT NULL,
  modified_date               DATETIME DEFAULT NULL,

  PRIMARY KEY (transcript_id),
  KEY seq_region_idx (seq_region_id, seq_region_start),
  KEY gene_index (gene_id),
  KEY xref_id_index (display_xref_id),
  KEY analysis_idx (analysis_id),
  UNIQUE KEY canonical_translation_idx (canonical_translation_id),
  KEY stable_id_idx (stable_id, version)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table transcript_attrib
@colour   #808000
@desc Enables storage of attributes that relate to transcripts.

@column transcript_id       Foreign key references to the @link transcript table.
@column attrib_type_id      Foreign key references to the @link attrib_type table.
@column value               Attribute value.

@see transcript

*/

CREATE TABLE transcript_attrib (

  transcript_id               INT(10) UNSIGNED NOT NULL DEFAULT '0',
  attrib_type_id              SMALLINT(5) UNSIGNED NOT NULL DEFAULT '0',
  value                       TEXT NOT NULL,

  KEY type_val_idx (attrib_type_id, value(40)),
  KEY val_only_idx (value(40)),
  KEY transcript_idx (transcript_id),
  UNIQUE KEY transcript_attribx (transcript_id, attrib_type_id, value(500))

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table transcript_supporting_feature
@colour   #808000
@desc Describes the exon prediction process by linking transcripts to DNA or protein alignment features.
As in several other tables, the feature_id column is a foreign key; the feature_type column specifies which table feature_id refers to.

@column transcript_id              Foreign key references to the @link transcript table.
@column feature_type               Feature type: 'dna_align_feature' or 'protein_align_feature'
@column feature_id                 Foreign key references to the @link dna_align_feature or @link protein_align_feature table depending on the feature type.

*/


CREATE TABLE transcript_supporting_feature (

  transcript_id               INT(10) UNSIGNED DEFAULT '0' NOT NULL,
  feature_type                ENUM('dna_align_feature','protein_align_feature'),
  feature_id                  INT(10) UNSIGNED DEFAULT '0' NOT NULL,

  UNIQUE KEY all_idx (transcript_id,feature_type,feature_id),
  KEY feature_idx (feature_type,feature_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM MAX_ROWS=100000000 AVG_ROW_LENGTH=80;


/**
@table translation
@colour   #808000
@desc Describes which parts of which exons are used in translation. The seq_start and seq_end columns are 1-based offsets into the relative coordinate system of start_exon_id and end_exon_id. i.e, if the translation starts at the first base of the exon, seq_start would be 1. Transcripts are related to translations by the transcript_id key in this table.

@column translation_id              Primary key, internal identifier.
@column transcript_id               Foreign key references to the @link transcript table.
@column seq_start                   1-based offset into the relative coordinate system of start_exon_id.
@column start_exon_id               Foreign key references to the @link exon table.
@column seq_end                     1-based offset into the relative coordinate system of end_exon_id.
@column end_exon_id                 Foreign key references to the @link exon table.
@column stable_id                   Release-independent stable identifier.
@column version                     Stable identifier version number.
@column created_date                Date created.
@column modified_date               Date modified.
*/


CREATE TABLE translation (

  translation_id              INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  transcript_id               INT(10) UNSIGNED NOT NULL,
  seq_start                   INT(10) NOT NULL,       # relative to exon start
  start_exon_id               INT(10) UNSIGNED NOT NULL,
  seq_end                     INT(10) NOT NULL,       # relative to exon start
  end_exon_id                 INT(10) UNSIGNED NOT NULL,
  stable_id                   VARCHAR(128) DEFAULT NULL,
  version                     SMALLINT UNSIGNED DEFAULT NULL,
  created_date                DATETIME DEFAULT NULL,
  modified_date               DATETIME DEFAULT NULL,

  PRIMARY KEY (translation_id),
  KEY transcript_idx (transcript_id),
  KEY stable_id_idx (stable_id, version)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table translation_attrib
@colour   #808000
@desc Enables storage of attributes that relate to translations.

@column translation_id      Foreign key references to the @link transcript table.
@column attrib_type_id      Foreign key references to the @link attrib_type table.
@column value               Attribute value.

@see translation

*/


CREATE TABLE translation_attrib (

  translation_id              INT(10) UNSIGNED NOT NULL DEFAULT '0',
  attrib_type_id              SMALLINT(5) UNSIGNED NOT NULL DEFAULT '0',
  value                       TEXT NOT NULL,

  KEY type_val_idx (attrib_type_id, value(40)),
  KEY val_only_idx (value(40)),
  KEY translation_idx (translation_id),
  UNIQUE KEY translation_attribx (translation_id, attrib_type_id, value(500))

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;



/**
@header Features
@colour   #3CB371
*/


/**
@table density_feature
@colour   #3CB371
@desc Describes features representing a density, or precentage coverage etc. in a given region.

@column density_feature_id        Primary key, internal identifier.
@column density_type_id           Density type. Foreign key references to the @link density_type table.
@column seq_region_id             Sequence region. Foreign key references to the @link seq_region table.
@column seq_region_start          Sequence start position.
@column seq_region_end            Sequence end position.
@column density_value             Density value.

@see density_type

*/


CREATE TABLE density_feature (

  density_feature_id    INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  density_type_id       INT(10) UNSIGNED NOT NULL,
  seq_region_id         INT(10) UNSIGNED NOT NULL,
  seq_region_start      INT(10) UNSIGNED NOT NULL,
  seq_region_end        INT(10) UNSIGNED NOT NULL,
  density_value         FLOAT NOT NULL,

  PRIMARY KEY (density_feature_id),
  KEY seq_region_idx (density_type_id, seq_region_id, seq_region_start),
  KEY seq_region_id_idx (seq_region_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table density_type
@colour   #3CB371
@desc Describes type representing a density, or percentage
coverage etc. in a given region.

@column density_type_id         Primary key, internal identifier.
@column analysis_id             Foreign key references to the @link analysis table.
@column block_size              Block size.
@column region_features         The number of features per sequence region inside this density type.
@column value_type              Value type, e.g. 'sum', 'ratio'.


@see density_feature

*/


CREATE TABLE density_type (

  density_type_id       INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  analysis_id           SMALLINT UNSIGNED NOT NULL,
  block_size            INT NOT NULL,
  region_features       INT NOT NULL,
  value_type            ENUM('sum','ratio') NOT NULL,

  PRIMARY KEY (density_type_id),
  UNIQUE KEY analysis_idx (analysis_id, block_size, region_features)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table ditag
@colour   #3CB371
@desc Represents a ditag object in the EnsEMBL database.
Corresponds to original tag containing the full sequence. This can be a single piece of sequence like CAGE tags or a ditag with concatenated sequence from 5' and 3' end like GIS or GSC tags.
This data is available as a DAS track in ContigView on the EnsEMBL web site.

@column ditag_id          Primary key, internal identifier.
@column name              Ditag name.
@column type              Ditag type.
@column tag_count         Tag count.
@column sequence          Sequence.

@see ditag_feature

*/

CREATE TABLE ditag (

       ditag_id          INT(10) UNSIGNED NOT NULL auto_increment,
       name              VARCHAR(30) NOT NULL,
       type              VARCHAR(30) NOT NULL,
       tag_count         smallint(6) UNSIGNED NOT NULL default 1,
       sequence          TINYTEXT NOT NULL,

       PRIMARY KEY (ditag_id)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table ditag_feature
@colour   #3CB371
@desc Describes where ditags hit on the genome. Represents a mapped ditag object in the EnsEMBL database. These are the original tags separated into start ("L") and end ("R") parts if applicable, successfully aligned to the genome.
Two DitagFeatures usually relate to one parent Ditag. Alternatively there are CAGE tags e.g. which only have a 5\'tag ("F").

@column ditag_feature_id          Primary key, internal identifier.
@column ditag_id                  Foreign key references to the @link ditag table.
@column ditag_pair_id             Ditag pair id.
@column seq_region_id             Foreign key references to the @link seq_region table.
@column seq_region_start          Sequence start position.
@column seq_region_end            Sequence end position.
@column seq_region_strand         Sequence region strand: 1 - forward; -1 - reverse.
@column analysis_id               Foreign key references to the @link analysis table.
@column hit_start                 Alignment hit start position.
@column hit_end                   Alignment hit end position.
@column hit_strand                Alignment hit strand: 1 - forward; -1 - reverse.
@column cigar_line                Used to encode gapped alignments.
@column ditag_side                Ditag side: L - start, R - end, F - 5\'tag only

@see ditag

*/


CREATE TABLE ditag_feature (

       ditag_feature_id   INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
       ditag_id           INT(10) UNSIGNED NOT NULL default '0',
       ditag_pair_id      INT(10) UNSIGNED NOT NULL default '0',
       seq_region_id      INT(10) UNSIGNED NOT NULL default '0',
       seq_region_start   INT(10) UNSIGNED NOT NULL default '0',
       seq_region_end     INT(10) UNSIGNED NOT NULL default '0',
       seq_region_strand  TINYINT(1) NOT NULL default '0',
       analysis_id        SMALLINT UNSIGNED NOT NULL default '0',
       hit_start          INT(10) UNSIGNED NOT NULL default '0',
       hit_end            INT(10) UNSIGNED NOT NULL default '0',
       hit_strand         TINYINT(1) NOT NULL default '0',
       cigar_line         TINYTEXT NOT NULL,
       ditag_side         ENUM('F', 'L', 'R') NOT NULL,

       PRIMARY KEY  (ditag_feature_id),
       KEY ditag_idx (ditag_id),
       KEY ditag_pair_idx (ditag_pair_id),
       KEY seq_region_idx (seq_region_id, seq_region_start, seq_region_end)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table intron_supporting_evidence
@colour   #3CB371
@desc Provides the evidence which we have used to declare an intronic region

@column intron_supporting_evidence_id Surrogate primary key
@column analysis_id                   Foreign key references to the @link analysis table.
@column seq_region_id                 Foreign key references to the @link seq_region table.
@column seq_region_start              Sequence start position.
@column seq_region_end                Sequence end position.
@column seq_region_strand             Sequence region strand: 1 - forward; -1 - reverse.
@column hit_name                      External entity name/identifier.
@column score                         Score supporting the intron
@column score_type                    The type of score e.g. NONE
@column is_splice_canonical           Indicates if the splice junction can be considered canonical i.e. behaves according to accepted rules

@see transcript_intron_supporting_evidence

*/


CREATE TABLE intron_supporting_evidence (
        intron_supporting_evidence_id INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
        analysis_id                   SMALLINT UNSIGNED NOT NULL,
        seq_region_id                 INT(10) UNSIGNED NOT NULL,
        seq_region_start              INT(10) UNSIGNED NOT NULL,
        seq_region_end                INT(10) UNSIGNED NOT NULL,
        seq_region_strand             TINYINT(2) NOT NULL,
        hit_name                      VARCHAR(100) NOT NULL,
        score                         DECIMAL(10,3),
        score_type                    ENUM('NONE', 'DEPTH') DEFAULT 'NONE',
        is_splice_canonical           TINYINT(1) NOT NULL DEFAULT 0,

        PRIMARY KEY (intron_supporting_evidence_id),

        UNIQUE KEY (analysis_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, hit_name),
        KEY seq_region_idx (seq_region_id, seq_region_start)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table map
@colour   #3CB371
@desc Stores the names of different genetic or radiation hybrid maps, for which there is marker map information.

@column map_id               Primary key, internal identifier.
@column map_name             Map name.


@see marker

*/

CREATE TABLE map (

  map_id                      INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  map_name                    VARCHAR(30) NOT NULL,

  PRIMARY KEY (map_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table marker
@colour   #3CB371
@desc Stores data about the marker itself. A marker in Ensembl consists of a pair of primer sequences, an expected product size and a set of associated identifiers known as synonyms.

@column marker_id                       Primary key, internal identifier.
@column display_marker_synonym_id       Marker synonym.
@column left_primer                     Left primer sequence.
@column right_primer                    Right primer sequence.
@column min_primer_dist                 Minimum primer distance.
@column max_primer_dist                 Maximum primer distance.
@column priority                        Priority.
@column type                            Type, e.g. 'est', 'microsatellite'.


@see marker_synonym
@see marker_map_location

*/


CREATE TABLE marker (

  marker_id                   INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  display_marker_synonym_id   INT(10) UNSIGNED,
  left_primer                 VARCHAR(100) NOT NULL,
  right_primer                VARCHAR(100) NOT NULL,
  min_primer_dist             INT(10) UNSIGNED NOT NULL,
  max_primer_dist             INT(10) UNSIGNED NOT NULL,
  priority                    INT,
  type                        ENUM('est', 'microsatellite'),

  PRIMARY KEY (marker_id),
  KEY marker_idx (marker_id, priority),
  KEY display_idx (display_marker_synonym_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table marker_feature
@colour   #3CB371
@desc Used to describe positions of markers on the assembly. Markers are placed on the genome electronically using an analysis program.

@column marker_feature_id       Primary key, internal identifier.
@column marker_id               Foreign key references to the @link marker table.
@column seq_region_id           Foreign key references to the @link seq_region table.
@column seq_region_start        Sequence start position.
@column seq_region_end          Sequence end position.
@column analysis_id             Foreign key references to the @link analysis table.
@column map_weight              The number of times that this marker has been mapped to the genome, e.g. a marker with map weight 3 has been mapped to 3 locations in the genome.


@see marker
@see marker_map_location
@see marker_synonym

*/


CREATE TABLE marker_feature (

  marker_feature_id           INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  marker_id                   INT(10) UNSIGNED NOT NULL,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  map_weight                  INT(10) UNSIGNED,

  PRIMARY KEY (marker_feature_id),
  KEY seq_region_idx (seq_region_id, seq_region_start),
  KEY analysis_idx (analysis_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table marker_map_location
@colour   #3CB371
@desc Stores map locations (genetic, radiation hybrid and in situ hybridization) for markers obtained from experimental evidence.

@column marker_id                Primary key, internal identifier.
@column map_id                   Foreign key references to the @link map table.
@column chromosome_name          Chromosome name
@column marker_synonym_id        Foreign key references to the @link marker_synonym table.
@column position                 Position of the map location.
@column lod_score                LOD score for map location.

@see marker
@see marker_feature


*/


CREATE TABLE marker_map_location (

  marker_id                   INT(10) UNSIGNED NOT NULL,
  map_id                      INT(10) UNSIGNED NOT NULL,
  chromosome_name             VARCHAR(15)  NOT NULL,
  marker_synonym_id           INT(10) UNSIGNED NOT NULL,
  position                    VARCHAR(15) NOT NULL,
  lod_score                   DOUBLE,

  PRIMARY KEY (marker_id, map_id),
  KEY map_idx (map_id, chromosome_name, position)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table marker_synonym
@colour   #3CB371
@desc Stores alternative names for markers, as well as their sources.

@column marker_synonym_id          Primary key, internal identifier.
@column marker_id                  Foreign key references to the @link marker table.
@column source                     Marker source.
@column name                       Alternative name for marker.


@see marker

*/


CREATE TABLE marker_synonym (

  marker_synonym_id           INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  marker_id                   INT(10) UNSIGNED NOT NULL,
  source                      VARCHAR(20),
  name                        VARCHAR(50),

  PRIMARY KEY (marker_synonym_id),
  KEY marker_synonym_idx (marker_synonym_id, name),
  KEY marker_idx (marker_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table misc_attrib
@colour   #3CB371
@desc Stores arbitrary attributes about the features in the misc_feature table.

@column misc_feature_id     Foreign key references to the @link misc_feature table.
@column attrib_type_id      Foreign key references to the @link attrib_type table.
@column value               Attribute value.

@see misc_feature

*/


CREATE TABLE misc_attrib (

  misc_feature_id             INT(10) UNSIGNED NOT NULL DEFAULT '0',
  attrib_type_id              SMALLINT(5) UNSIGNED NOT NULL DEFAULT '0',
  value                       TEXT NOT NULL,

  KEY type_val_idx (attrib_type_id, value(40)),
  KEY val_only_idx (value(40)),
  KEY misc_feature_idx (misc_feature_id),
  UNIQUE KEY misc_attribx (misc_feature_id, attrib_type_id, value(500))

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table misc_feature
@colour   #3CB371
@desc Allows for storage of arbitrary features.

@column misc_feature_id             Primary key, internal identifier.
@column seq_region_id               Foreign key references to the @link seq_region table.
@column seq_region_start            Sequence start position.
@column seq_region_end              Sequence end position.
@column seq_region_strand           Sequence region strand: 1 - forward; -1 - reverse.

@see misc_attrib

*/


CREATE TABLE misc_feature (

  misc_feature_id             INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL DEFAULT '0',
  seq_region_start            INT(10) UNSIGNED NOT NULL DEFAULT '0',
  seq_region_end              INT(10) UNSIGNED NOT NULL DEFAULT '0',
  seq_region_strand           TINYINT(4) NOT NULL DEFAULT '0',

  PRIMARY KEY (misc_feature_id),
  KEY seq_region_idx (seq_region_id, seq_region_start)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table misc_feature_misc_set
@colour   #3CB371
@desc This table classifies features into distinct sets.

@column misc_feature_id        Primary key, internal identifier. Foreign key references to the @link misc_feature table.
@column misc_set_id            Primary key, internal identifier. Foreign key references to the @link misc_feature table.

@see misc_feature
@see misc_set
*/


CREATE TABLE misc_feature_misc_set (

  misc_feature_id             INT(10) UNSIGNED NOT NULL DEFAULT '0',
  misc_set_id                 SMALLINT(5) UNSIGNED NOT NULL DEFAULT '0',

  PRIMARY KEY (misc_feature_id, misc_set_id),
  KEY reverse_idx (misc_set_id, misc_feature_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table misc_set
@colour   #3CB371
@desc Defines "sets" that the features held in the misc_feature table can be grouped into.

@column misc_set_id           Primary key, internal identifier.
@column code                  Set code, e.g. bac_map
@column name                  Code name, e.g. BAC map
@column description           Code description, e.g. Full list of mapped BAC clones
@column max_length            Longest feature, e.g. 500000

@see misc_feature_misc_set

*/


CREATE TABLE misc_set (

  misc_set_id                 SMALLINT(5) UNSIGNED NOT NULL AUTO_INCREMENT,
  code                        VARCHAR(25) NOT NULL DEFAULT '',
  name                        VARCHAR(255) NOT NULL DEFAULT '',
  description                 TEXT NOT NULL,
  max_length                  INT UNSIGNED NOT NULL,

  PRIMARY KEY (misc_set_id),
  UNIQUE KEY code_idx (code)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table prediction_exon
@colour   #3CB371
@desc Stores exons that are predicted by ab initio gene finder programs. Unlike EnsEMBL exons they are not supported by any evidence.

@column prediction_exon_id              Primary key, internal identifier.
@column prediction_transcript_id        Foreign key references to the @link prediction_transcript table.
@column exon_rank                       Exon rank
@column seq_region_id                   Foreign key references to the @link seq_region table.
@column seq_region_start                Sequence start position.
@column seq_region_end                  Sequence end position.
@column seq_region_strand               Sequence region strand: 1 - forward; -1 - reverse.
@column start_phase                     Exon start phase.
@column score                           Prediction score.
@column p_value                         Prediction p-value.

*/


CREATE TABLE prediction_exon (

  prediction_exon_id          INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  prediction_transcript_id    INT(10) UNSIGNED NOT NULL,
  exon_rank                   SMALLINT UNSIGNED NOT NULL,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT NOT NULL,
  start_phase                 TINYINT NOT NULL,
  score                       DOUBLE,
  p_value                     DOUBLE,

  PRIMARY KEY (prediction_exon_id),
  KEY transcript_idx (prediction_transcript_id),
  KEY seq_region_idx (seq_region_id, seq_region_start)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table prediction_transcript
@colour   #3CB371
@desc Stores transcripts that are predicted by ab initio gene finder programs (e.g. genscan, SNAP). Unlike EnsEMBL transcripts they are not supported by any evidence.

@column prediction_transcript_id        Primary key, internal identifier.
@column seq_region_id                   Foreign key references to the @link seq_region table.
@column seq_region_start                Sequence start position.
@column seq_region_end                  Sequence end position.
@column seq_region_strand               Sequence region strand: 1 - forward; -1 - reverse.
@column analysis_id                     Foreign key references to the @link analysis table.
@column display_label                   Display label for the EnsEMBL web site.

*/

CREATE TABLE prediction_transcript (

  prediction_transcript_id    INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT NOT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  display_label               VARCHAR(255),

  PRIMARY KEY (prediction_transcript_id),
  KEY seq_region_idx (seq_region_id, seq_region_start),
  KEY analysis_idx (analysis_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table repeat_consensus
@colour   #3CB371
@desc Stores consensus sequences obtained from analysing repeat features.

@column repeat_consensus_id              Primary key, internal identifier.
@column repeat_name                      Repeat name.
@column repeat_class                     E.g. 'Satellite', 'tRNA', 'LTR'.
@column repeat_type                      E.g. 'Satellite repeats', 'Tandem repeats', 'Low complexity regions'.
@column repeat_consensus                 Repeat consensus sequence.

*/


CREATE TABLE repeat_consensus (

  repeat_consensus_id         INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  repeat_name                 VARCHAR(255) NOT NULL,
  repeat_class                VARCHAR(100) NOT NULL,
  repeat_type                 VARCHAR(40) NOT NULL,
  repeat_consensus            TEXT,

  PRIMARY KEY (repeat_consensus_id),
  KEY name (repeat_name),
  KEY class (repeat_class),
  KEY consensus (repeat_consensus(10)),
  KEY type (repeat_type)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table repeat_feature
@colour   #3CB371
@desc Describes sequence repeat regions.

@column repeat_feature_id           Primary key, internal identifier.
@column seq_region_id               Foreign key references to the @link seq_region table.
@column seq_region_start            Sequence start position.
@column seq_region_end              Sequence end position.
@column seq_region_strand           Sequence region strand: 1 - forward; -1 - reverse.
@column repeat_start                Repeat sequence start.
@column repeat_end                  Repeat sequence end
@column repeat_consensus_id         Foreign key references to the @link repeat_consensus table.
@column analysis_id                 Foreign key references to the @link analysis table.
@column score                       Analysis score.

*/

CREATE TABLE repeat_feature (

  repeat_feature_id           INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT(1) DEFAULT '1' NOT NULL,
  repeat_start                INT(10) NOT NULL,
  repeat_end                  INT(10) NOT NULL,
  repeat_consensus_id         INT(10) UNSIGNED NOT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  score                       DOUBLE,

  PRIMARY KEY (repeat_feature_id),
  KEY seq_region_idx (seq_region_id, seq_region_start),
  KEY repeat_idx (repeat_consensus_id),
  KEY analysis_idx (analysis_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM MAX_ROWS=100000000 AVG_ROW_LENGTH=80;


/**
@table simple_feature
@colour   #3CB371
@desc Describes general genomic features that don't fit into any of the more specific feature tables.

@column simple_feature_id       Primary key, internal identifier.
@column seq_region_id           Foreign key references to the @link seq_region table.
@column seq_region_start        Sequence start position.
@column seq_region_end          Sequence end position.
@column seq_region_strand       Sequence region strand: 1 - forward; -1 - reverse.
@column display_label           Display label for the EnsEMBL web site.
@column analysis_id             Foreign key references to the @link analysis table.
@column score                   Analysis score.

*/

CREATE TABLE simple_feature (

  simple_feature_id           INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id               INT(10) UNSIGNED NOT NULL,
  seq_region_start            INT(10) UNSIGNED NOT NULL,
  seq_region_end              INT(10) UNSIGNED NOT NULL,
  seq_region_strand           TINYINT(1) NOT NULL,
  display_label               VARCHAR(255) NOT NULL,
  analysis_id                 SMALLINT UNSIGNED NOT NULL,
  score                       DOUBLE,

  PRIMARY KEY (simple_feature_id),
  KEY seq_region_idx (seq_region_id, seq_region_start),
  KEY analysis_idx (analysis_id),
  KEY hit_idx (display_label)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table transcript_intron_supporting_evidence
@colour   #3CB371
@desc Links intronic evidence to a pair of exons used within a transcript and to resolve the m:m relationship between introns and transcripts

@column intron_supporting_evidence_id Foreign key references to the @link intron_supporting_evidence table
@column transcript_id                 Foreign key references to the @link transcript table.
@column previous_exon_id              Foreign key to @link exon indicating the left hand flanking exon of the intron (assume forward strand)
@column next_exon_id                  Foreign key to @link exon indicating the right hand flanking exon of the intron (assume forward strand)

@see intron_supporting_evidence
@see transcript
@see exon

*/

CREATE TABLE transcript_intron_supporting_evidence (
transcript_id                 INT(10) UNSIGNED NOT NULL,
intron_supporting_evidence_id INT(10) UNSIGNED NOT NULL,
previous_exon_id              INT(10) UNSIGNED NOT NULL,
next_exon_id                  INT(10) UNSIGNED NOT NULL,
PRIMARY KEY (intron_supporting_evidence_id, transcript_id),
KEY transcript_idx (transcript_id)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@header ID Mapping
@colour #FF6666
*/

/**
@table gene_archive
@colour #FF6666
@desc Contains a snapshot of the stable IDs associated with genes deleted or changed between releases. Includes gene, transcript and translation stable IDs.

@column gene_stable_id              Stable ID of retired gene.
@column gene_version                Last live gene stable ID version.
@column transcript_stable_id        Stable ID of associated transcript.
@column transcript_version          Last live transcript stable ID version.
@column translation_stable_id       Stable ID of associated translation.
@column translation_version         Last live translation stable ID.
@column peptide_archive_id          Foreign key references to the @link peptide archive table.
@column mapping_session_id          Foreign key references to the @link mapping_session table.

@see gene

*/


CREATE TABLE gene_archive (

  gene_stable_id              VARCHAR(128) NOT NULL,
  gene_version                SMALLINT NOT NULL DEFAULT 1,
  transcript_stable_id        VARCHAR(128) NOT NULL,
  transcript_version          SMALLINT NOT NULL DEFAULT 1,
  translation_stable_id       VARCHAR(128),
  translation_version         SMALLINT NOT NULL DEFAULT 1,
  peptide_archive_id          INT(10) UNSIGNED,
  mapping_session_id          INT(10) UNSIGNED NOT NULL,

  KEY gene_idx (gene_stable_id, gene_version),
  KEY transcript_idx (transcript_stable_id, transcript_version),
  KEY translation_idx (translation_stable_id, translation_version),
  KEY peptide_archive_id_idx (peptide_archive_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table mapping_session
@colour #FF6666
@desc Stores details of ID mapping sessions - a mapping session represents the session when stable IDs where mapped from one database to another. Details of the "old" and "new" databases are stored.

@column mapping_session_id          Primary key, internal identifier.
@column old_db_name                 Old Ensembl database name.
@column new_db_name                 New Ensembl database name.
@column old_release                 Old Ensembl database release.
@column new_release                 New Ensembl database release.
@column old_assembly                Old assembly.
@column new_assembly                New assembly.
@column created                     Date created.

@see stable_id_event
@see stable_id

*/


CREATE TABLE mapping_session (

  mapping_session_id          INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  old_db_name                 VARCHAR(80) NOT NULL DEFAULT '',
  new_db_name                 VARCHAR(80) NOT NULL DEFAULT '',
  old_release                 VARCHAR(5) NOT NULL DEFAULT '',
  new_release                 VARCHAR(5) NOT NULL DEFAULT '',
  old_assembly                VARCHAR(20) NOT NULL DEFAULT '',
  new_assembly                VARCHAR(20) NOT NULL DEFAULT '',
  created                     DATETIME NOT NULL,

  PRIMARY KEY (mapping_session_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table peptide_archive
@colour #FF6666
@desc Contains the peptides for deleted or changed translations.

@column peptide_archive_id         Primary key, internal identifier.
@column md5_checksum               MD5 checksum hexadecimal digest of the peptide sequence.
@column peptide_seq                Peptide sequence of retired translation.


*/


CREATE TABLE peptide_archive (

  peptide_archive_id         INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  md5_checksum               VARCHAR(32),
  peptide_seq                MEDIUMTEXT NOT NULL,

  PRIMARY KEY (peptide_archive_id),
  KEY checksum (md5_checksum)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table mapping_set
@colour #FF6666
@desc Table structure for seq_region mapping between releases.

@column mapping_set_id            Primary key, internal identifier.
@column internal_schema_build     Schema version of the current database (eg 72_37)
@column external_schema_build     Schema version of the database the comparison was run against (eg 71_37)

*/


CREATE TABLE mapping_set (

        mapping_set_id  INT(10) UNSIGNED NOT NULL,
        internal_schema_build    VARCHAR(20) NOT NULL,
        external_schema_build    VARCHAR(20) NOT NULL,

        PRIMARY KEY (mapping_set_id),
        UNIQUE KEY mapping_idx (internal_schema_build, external_schema_build)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table stable_id_event
@colour #FF6666
@desc Represents what happened to all gene, transcript and translation stable IDs during a mapping session.
This includes which IDs where deleted, created and related to each other. Each event is represented by one or more rows in the table.

@column old_stable_id             Gene/transcript/translation stable id for the previous release.
@column old_version               Stable id version.
@column new_stable_id             Gene/transcript/translation stable id for the current release.
@column new_version               Stable id version.
@column mapping_session_id        Foreign key references to the @link mapping_session table.
@column type                      ENUM('gene', 'transcript', 'translation') NOT NULL,
@column score                     Combined mapping score.

@see mapping_session

*/


CREATE TABLE stable_id_event (

  old_stable_id             VARCHAR(128),
  old_version               SMALLINT,
  new_stable_id             VARCHAR(128),
  new_version               SMALLINT,
  mapping_session_id        INT(10) UNSIGNED NOT NULL DEFAULT '0',
  type                      ENUM('gene', 'transcript', 'translation') NOT NULL,
  score                     FLOAT NOT NULL DEFAULT 0,

  UNIQUE KEY uni_idx (mapping_session_id, old_stable_id, new_stable_id, type),

  KEY new_idx (new_stable_id),
  KEY old_idx (old_stable_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table seq_region_mapping
@colour #FF6666
@desc Describes how the core seq_region_id have changed from release to release.

@column external_seq_region_id            Foreign key references to the @link seq_region table.
@column internal_seq_region_id            Foreign key references to the @link seq_region table.
@column mapping_set_id                    Foreign key references to the @link mapping_set table.


*/


CREATE TABLE seq_region_mapping (

        external_seq_region_id  INT(10) UNSIGNED NOT NULL,
        internal_seq_region_id  INT(10) UNSIGNED NOT NULL,
        mapping_set_id          INT(10) UNSIGNED NOT NULL,

        KEY mapping_set_idx (mapping_set_id)

) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@header External References
@colour   #1E90FF
*/


/**
@table  associated_group
@colour   #1E90FF
@desc   Groups together xref associations under a single description. Used when more than one associated xref term must be used to describe a condition

@column associated_group_id Associated group id. Primary key, internal identifier
@column description Optional description for this group

@see associated_xref
*/


CREATE TABLE associated_group (

  associated_group_id            INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  description                    VARCHAR(128) DEFAULT NULL,

  PRIMARY KEY (associated_group_id)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table  associated_xref
@colour   #1E90FF
@desc   This table associates extra associated annotations with a given ontology xref evidence and source under a specific condition.   For GO this allows qualifiers (with/from) or annotation extensions to be added to a given ontology annotation.

@column associated_xref_id Associated xref id. Primary key, internal identifier
@column object_xref_id Object xref id this associated xref is linked to. Foreign key linked to the @link object_xref table
@column xref_id Xref which is the associated term. Foreign key linked to the @link xref table
@column source_xref_id  Xref which is source of this association. Foreign key linked to the @link xref table
@column condition_type The type of condition this link occurs in e.g. evidence, from, residue or assigned_by
@column associated_group_id Foreign key to allow for @link associated_group
@column rank The rank in which the association occurs within an @link associated_group

@see object_xref
@see associated_group
@see xref
*/


CREATE TABLE associated_xref (

  associated_xref_id             INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  object_xref_id                 INT(10) UNSIGNED DEFAULT '0' NOT NULL,
  xref_id                        INT(10) UNSIGNED DEFAULT '0' NOT NULL,
  source_xref_id                 INT(10) UNSIGNED DEFAULT NULL,
  condition_type                 VARCHAR(128) DEFAULT NULL,
  associated_group_id            INT(10) UNSIGNED DEFAULT NULL,
  rank                           INT(10) UNSIGNED DEFAULT '0',

  PRIMARY KEY (associated_xref_id),
  KEY associated_source_idx (source_xref_id),
  KEY associated_object_idx (object_xref_id),
  KEY associated_idx (xref_id),
  KEY associated_group_idx (associated_group_id),
  UNIQUE KEY object_associated_source_type_idx (object_xref_id, xref_id, source_xref_id, condition_type, associated_group_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table dependent_xref
@colour   #1E90FF
@desc Describes dependent external references which can't be directly mapped to Ensembl entities.
They are linked to primary external references instead.

@column object_xref_id          Primary key, internal identifier. Foreign key references to the @link object_xref table.
@column master_xref_id          Foreign key references to the @link xref table.
@column dependent_xref_id       Foreign key references to the @link xref table.

@see xref
@see object_xref

*/


CREATE TABLE dependent_xref (

  object_xref_id         INT(10) UNSIGNED NOT NULL,
  master_xref_id         INT(10) UNSIGNED NOT NULL,
  dependent_xref_id      INT(10) UNSIGNED NOT NULL,

  PRIMARY KEY (object_xref_id),
  KEY dependent (dependent_xref_id),
  KEY master_idx (master_xref_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table external_db
@colour   #1E90FF
@desc Stores data about the external databases in which the objects described in the xref table are stored.

@column external_db_id              Primary key, internal identifier.
@column db_name                     Database name.
@column db_release                  Database release.
@column status                      Status, e.g. 'KNOWNXREF','KNOWN','XREF','PRED','ORTH','PSEUDO'.
@column priority                    Determines which one of the xrefs will be used as the gene name.
@column db_display_name             Database display name.
@column type                        Type, e.g. 'ARRAY', 'ALT_TRANS', 'ALT_GENE', 'MISC', 'LIT', 'PRIMARY_DB_SYNONYM', 'ENSEMBL'.
@column secondary_db_name           Secondary database name.
@column secondary_db_table          Secondary database table.
@column description                 Description.

@see xref
@see unmapped_object
@see protein_align_feature
@see dna_align_feature

*/


CREATE TABLE external_db (

  external_db_id              INT UNSIGNED NOT NULL AUTO_INCREMENT,
  db_name                     VARCHAR(100) NOT NULL,
  db_release                  VARCHAR(255),
  status                      ENUM('KNOWNXREF','KNOWN','XREF','PRED','ORTH', 'PSEUDO') NOT NULL,
  priority                    INT NOT NULL,
  db_display_name             VARCHAR(255),
  type                        ENUM('ARRAY', 'ALT_TRANS', 'ALT_GENE', 'MISC', 'LIT', 'PRIMARY_DB_SYNONYM', 'ENSEMBL'),
  secondary_db_name           VARCHAR(255) DEFAULT NULL,
  secondary_db_table          VARCHAR(255) DEFAULT NULL,
  description                 TEXT,

  PRIMARY KEY (external_db_id),
  UNIQUE KEY db_name_db_release_idx (db_name,db_release)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table biotype
@colour   #1E90FF
@desc Stores data about the biotypes and mappings to Sequence Ontology.

@column biotype_id              Primary key, internal identifier.
@column name                    Ensembl biotype name.
@column object_type             Ensembl object type: 'gene' or 'transcript'.
@column db_type                 Type, e.g. 'cdna', 'core', 'coreexpressionatlas', 'coreexpressionest', 'coreexpressiongnf', 'funcgen', 'otherfeatures', 'rnaseq', 'variation', 'vega', 'presite', 'sangervega'
@column attrib_type_id          Foreign key references to the @link attrib_type table.
@column description             Description.
@column biotype_group           Group, e.g. 'coding', 'pseudogene', 'snoncoding', 'lnoncoding', 'mnoncoding', 'LRG', 'undefined', 'no_group'
@column so_acc                  Sequence Ontology accession of the biotype.

@see attrib_type

*/


CREATE TABLE biotype (
  biotype_id      INT UNSIGNED NOT NULL AUTO_INCREMENT,
  name            VARCHAR(64) NOT NULL,
  object_type     ENUM('gene','transcript') NOT NULL DEFAULT 'gene',
  db_type         set('cdna','core','coreexpressionatlas','coreexpressionest','coreexpressiongnf','funcgen','otherfeatures','rnaseq','variation','vega','presite','sangervega') NOT NULL DEFAULT 'core',
  attrib_type_id  INT DEFAULT NULL,
  description     TEXT,
  biotype_group   ENUM('coding','pseudogene','snoncoding','lnoncoding','mnoncoding','LRG','undefined','no_group') DEFAULT NULL,
  so_acc          VARCHAR(64),
  PRIMARY KEY (biotype_id),
  UNIQUE KEY name_type_idx (name, object_type)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


/**
@table external_synonym
@colour   #1E90FF
@desc Some xref objects can be referred to by more than one name. This table relates names to xref IDs.

@column xref_id           Primary key, internal identifier.
@column synonym           Synonym

@see xref

*/


CREATE TABLE external_synonym (

  xref_id                     INT(10) UNSIGNED NOT NULL,
  synonym                     VARCHAR(100) NOT NULL,

  PRIMARY KEY (xref_id, synonym),
  KEY name_index (synonym)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table identity_xref
@colour   #1E90FF
@desc Describes how well a particular xref object matches the EnsEMBL object.

@column object_xref_id        Primary key, internal identifier. Foreign key references to the @link object_xref table.
@column xref_identity         Percentage identity.
@column ensembl_identity      Percentage identity.
@column xref_start            Xref sequence start.
@column xref_end              Xref sequence end.
@column ensembl_start         Ensembl sequence start.
@column ensembl_end           Ensembl sequence end.
@column cigar_line            Used to encode gapped alignments.
@column score                 Match score.
@column evalue                Match evalue.

@see object_xref

*/


CREATE TABLE identity_xref (

  object_xref_id          INT(10) UNSIGNED NOT NULL,
  xref_identity           INT(5),
  ensembl_identity        INT(5),

  xref_start              INT,
  xref_end                INT,
  ensembl_start           INT,
  ensembl_end             INT,
  cigar_line              TEXT,

  score                   DOUBLE,
  evalue                  DOUBLE,

  PRIMARY KEY (object_xref_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table interpro
@colour   #1E90FF
@desc Allows storage of links to the InterPro database. InterPro is a database of protein families, domains and functional sites in which identifiable features found in known proteins can be applied to unknown protein sequences.
 <a href="http://www.ebi.ac.uk/interpro/">InterPro</a> - The InterPro website

@column interpro_ac               InterPro protein accession number.
@column id                        InterPro protein id.

*/


CREATE TABLE interpro (
  interpro_ac               VARCHAR(40) NOT NULL,
  id                        VARCHAR(40) CHARACTER SET latin1 COLLATE latin1_bin NOT NULL,

  UNIQUE KEY accession_idx (interpro_ac, id),
  KEY id_idx (id)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table object_xref
@colour   #1E90FF
@desc Describes links between EnsEMBL objects and objects held in external databases.
The EnsEMBL object can be one of several types; the type is held in the ensembl_object_type column.
The ID of the particular EnsEMBL gene, translation or whatever is given in the ensembl_id column.
The xref_id points to the entry in the xref table that holds data about the external object.
Each EnsEMBL object can be associated with zero or more xrefs. An xref object can be associated with one or more EnsEMBL objects.

@column object_xref_id            Primary key, internal identifier.
@column ensembl_id                Foreign key references to the @link seq_region, @link transcript, @link gene, @translation tables depending on ensembl_object_type.
@column ensembl_object_type       Ensembl object type: 'RawContig', 'Transcript', 'Gene','Translation'.
@column xref_id                   Foreign key references to the @link xref table.
@column linkage_annotation        Additional annotation on the linkage.
@column analysis_id               Foreign key references to the @link analysis table.

@see xref
@see identity_xref

*/


CREATE TABLE object_xref (

  object_xref_id              INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  ensembl_id                  INT(10) UNSIGNED NOT NULL,
  ensembl_object_type         ENUM('RawContig', 'Transcript', 'Gene',
                                   'Translation', 'Operon', 'OperonTranscript', 'Marker') NOT NULL,
  xref_id                     INT(10) UNSIGNED NOT NULL,
  linkage_annotation          VARCHAR(255) DEFAULT NULL,
  analysis_id                 SMALLINT UNSIGNED,

  PRIMARY KEY (object_xref_id),

  UNIQUE KEY xref_idx (xref_id, ensembl_object_type, ensembl_id, analysis_id),

  KEY ensembl_idx (ensembl_object_type, ensembl_id),
  KEY analysis_idx (analysis_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table ontology_xref
@colour   #1E90FF
@desc This table associates Evidence Tags to the relationship between EnsEMBL objects and ontology accessions (primarily GO accessions).
The relationship to GO that is stored in the database is actually derived through the relationship of EnsEMBL peptides to SwissProt peptides, i.e. the relationship is derived like this:

  ENSP -> SWISSPROT -> GO

And the evidence tag describes the relationship between the SwissProt Peptide and the GO entry.
In reality, however, we store this in the database like this:

  ENSP -> SWISSPROT
  ENSP -> GO

and the evidence tag hangs off of the relationship between the ENSP and the GO identifier.
Some ENSPs are associated with multiple closely related Swissprot entries which may both be associated with the same GO identifier but with different evidence tags.
For this reason a single Ensembl - external db object relationship in the object_xref table can be associated with multiple evidence tags in the ontology_xref table.

@column object_xref_id          Composite key. Foreign key references to the @link object_xref table.
@column source_xref_id          Composite key. Foreign key references to the @link xref table.
@column linkage_type            Composite key. <a href="http://www.geneontology.org/GO.evidence.shtml">Evidence tags</a>

@see object_xref

*/


CREATE TABLE ontology_xref (

  object_xref_id          INT(10) UNSIGNED DEFAULT '0' NOT NULL,
  source_xref_id          INT(10) UNSIGNED DEFAULT NULL,
  linkage_type            VARCHAR(3) DEFAULT NULL,

  KEY source_idx (source_xref_id),
  KEY object_idx (object_xref_id),
  UNIQUE KEY object_source_type_idx (object_xref_id, source_xref_id, linkage_type)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table unmapped_object
@colour   #1E90FF
@desc Describes why a particular external entity was not mapped to an ensembl one.

@column unmapped_object_id         Primary key, internal identifier.
@column type                       Object type: 'xref', 'cDNA', 'Marker'.
@column analysis_id                Foreign key references to the @link analysis table.
@column external_db_id             Foreign key references to the @link external_db table.
@column identifier                 External database identifier.
@column unmapped_reason_id         Foreign key references to the @link unmapped_reason table.
@column query_score                Actual mapping query score.
@column target_score               Target mapping query score.
@column ensembl_id                 Foreign key references to the @link seq_region, @link transcript, @link gene, @translation tables depending on ensembl_object_type.
@column ensembl_object_type        Ensembl object type: 'RawContig', 'Transcript', 'Gene','Translation'.
@column parent                     Foreign key references to the @link dependent_xref table, in case the unmapped object is dependent on a primary external reference which wasn't mapped to an ensembl one.

@see analysis
@see external_db
@see unmapped_reason
*/


CREATE TABLE unmapped_object (

  unmapped_object_id    INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  type                  ENUM('xref', 'cDNA', 'Marker') NOT NULL,
  analysis_id           SMALLINT UNSIGNED NOT NULL,
  external_db_id        INT UNSIGNED,
  identifier            VARCHAR(255) NOT NULL,
  unmapped_reason_id    INT(10) UNSIGNED NOT NULL,
  query_score           DOUBLE,
  target_score          DOUBLE,
  ensembl_id            INT(10) UNSIGNED DEFAULT '0',
  ensembl_object_type   ENUM('RawContig','Transcript','Gene','Translation') DEFAULT 'RawContig',
  parent                VARCHAR(255) DEFAULT NULL,

  PRIMARY KEY (unmapped_object_id),
  UNIQUE KEY unique_unmapped_obj_idx (ensembl_id, ensembl_object_type, identifier, unmapped_reason_id,parent, external_db_id),
  KEY id_idx (identifier(50)),
  KEY anal_exdb_idx (analysis_id, external_db_id),
  KEY ext_db_identifier_idx (external_db_id, identifier)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table unmapped_reason
@colour   #1E90FF
@desc Describes the reason why a mapping failed.

@column unmapped_reason_id           Primary key, internal identifier.
@column summary_description          Summarised description.
@column full_description             Full description.

@see unmapped_object
*/


CREATE TABLE unmapped_reason (

  unmapped_reason_id     INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  summary_description    VARCHAR(255),
  full_description       VARCHAR(255),

  PRIMARY KEY (unmapped_reason_id)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table xref
@colour   #1E90FF
@desc Holds data about objects which are external to EnsEMBL, but need to be associated with EnsEMBL objects.
Information about the database that the external object is stored in is held in the external_db table entry referred to by the external_db column.

@column xref_id                 Primary key, internal identifier.
@column external_db_id          Foreign key references to the @link external_db table.
@column dbprimary_acc           Primary accession number.
@column display_label           Display label for the EnsEMBL web site.
@column version                 Object version.
@column description             Object description.
@column info_type               'PROJECTION', 'MISC', 'DEPENDENT','DIRECT', 'SEQUENCE_MATCH','INFERRED_PAIR', 'PROBE','UNMAPPED', 'COORDINATE_OVERLAP', 'CHECKSUM'.
@column info_text               Text

@see associated_xref
@see dependent_xref
@see gene
@see transcript
@see external_db
@see external_synonym
@see object_xref
@see ontology_xref

*/


CREATE TABLE xref (

   xref_id                    INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
   external_db_id             INT UNSIGNED NOT NULL,
   dbprimary_acc              VARCHAR(512) NOT NULL,
   display_label              VARCHAR(512) NOT NULL,
   version                    VARCHAR(10) DEFAULT NULL,
   description                TEXT,
   info_type                  ENUM( 'NONE', 'PROJECTION', 'MISC', 'DEPENDENT',
                                    'DIRECT', 'SEQUENCE_MATCH',
                                    'INFERRED_PAIR', 'PROBE',
                                    'UNMAPPED', 'COORDINATE_OVERLAP',
                                    'CHECKSUM' ) DEFAULT 'NONE' NOT NULL,
   info_text                  VARCHAR(255) DEFAULT '' NOT NULL,

   PRIMARY KEY (xref_id),
   UNIQUE KEY id_index (dbprimary_acc, external_db_id, info_type, info_text, version),
   KEY display_index (display_label),
   KEY info_type_idx (info_type)

) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@header Misc
@colour #FF8500
*/


/**
@table operon
@colour #FF8500
@desc allows one or more polycistronic transcripts to be grouped together

@column operon_id                   Primary key, internal identifier.
@column analysis_id                 Foreign key references to the @link analysis table.
@column seq_region_id               Foreign key references to the @link seq_region table.
@column seq_region_start            Sequence start position.
@column seq_region_end              Sequence end position.
@column seq_region_strand           Sequence region strand: 1 - forward; -1 - reverse.
@column display_label               Short name for operon
@column stable_id                   Release-independent stable identifier.
@column version                     Stable identifier version number.
@column created_date                Date created.
@column modified_date               Date modified.

@see operon_transcript
@see operon_stable_id
*/


CREATE TABLE operon (
  operon_id                 INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id             INT(10) UNSIGNED NOT NULL,
  seq_region_start          INT(10) UNSIGNED NOT NULL,
  seq_region_end            INT(10) UNSIGNED NOT NULL,
  seq_region_strand         TINYINT(2) NOT NULL,
  display_label             VARCHAR(255) DEFAULT NULL,
  analysis_id               SMALLINT UNSIGNED NOT NULL,
  stable_id                 VARCHAR(128) DEFAULT NULL,
  version                   SMALLINT UNSIGNED DEFAULT NULL,
  created_date              DATETIME DEFAULT NULL,
  modified_date             DATETIME DEFAULT NULL,

  PRIMARY KEY (operon_id),
  KEY seq_region_idx (seq_region_id, seq_region_start),
  KEY name_idx (display_label),
  KEY stable_id_idx (stable_id, version)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table operon_transcript
@colour #FF8500
@desc represents polycistronic transcripts which belong to operons and encode more than one gene

@column operon_transcript_id        Primary key, internal identifier.
@column analysis_id                 Foreign key references to the @link analysis table.
@column seq_region_id               Foreign key references to the @link seq_region table.
@column seq_region_start            Sequence start position.
@column seq_region_end              Sequence end position.
@column seq_region_strand           Sequence region strand: 1 - forward; -1 - reverse.
@column operon_id                   Foreign key references to the @link operon table.
@column display_label               Short name for operon transcript
@column stable_id                   Release-independent stable identifier.
@column version                     Stable identifier version number.
@column created_date                Date created.
@column modified_date               Date modified.

@see operon
@see operon_transcript_stable_id
@see operon_transcript_gene
*/


CREATE TABLE operon_transcript (
  operon_transcript_id      INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  seq_region_id             INT(10) UNSIGNED NOT NULL,
  seq_region_start          INT(10) UNSIGNED NOT NULL,
  seq_region_end            INT(10) UNSIGNED NOT NULL,
  seq_region_strand         TINYINT(2) NOT NULL,
  operon_id                 INT(10) UNSIGNED NOT NULL,
  display_label             VARCHAR(255) DEFAULT NULL,
  analysis_id               SMALLINT UNSIGNED NOT NULL,
  stable_id                 VARCHAR(128) DEFAULT NULL,
  version                   SMALLINT UNSIGNED DEFAULT NULL,
  created_date              DATETIME DEFAULT NULL,
  modified_date             DATETIME DEFAULT NULL,

  PRIMARY KEY (operon_transcript_id),
  KEY operon_idx (operon_id),
  KEY seq_region_idx (seq_region_id, seq_region_start),
  KEY stable_id_idx (stable_id, version)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


/**
@table operon_transcript_gene
@colour #FF8500
@desc allows association of genes with polycistronic transcripts

@column operon_transcript_id      Foreign key references to the @link operon_transcript table.
@column gene_id                   Foreign key references to the @link gene table.

@see operon_transcript
@see gene
*/


CREATE TABLE operon_transcript_gene (
  operon_transcript_id      INT(10) UNSIGNED,
  gene_id                   INT(10) UNSIGNED,

  KEY operon_transcript_gene_idx (operon_transcript_id,gene_id)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


