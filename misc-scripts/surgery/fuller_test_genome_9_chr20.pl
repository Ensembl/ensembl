#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

use Getopt::Std;
use DBI;

use vars qw($opt_d $opt_s $opt_u $opt_h $opt_p $opt_P);

#get command line arguments
getopts('s:d:h:u:p:P:');

my $usage = "Usage: " .
  "test_genome.pl -s srcDB -d destDB -h host -u user [-p pass] [-P port]\n";

my ($user, $pass, $host, $destDB, $srcDB, $port) = 
  ($opt_u, $opt_p, $opt_h, $opt_d, $opt_s, $opt_P);

unless($port) {
  $port = 3306;
}

# If needed command line args are missing print the usage string and quit
$user and $host and $destDB and $srcDB or die $usage;



my $dsn = "DBI:mysql:host=$host;port=$port";

#print "User: $user; Pass: $pass; DSN: $dsn\n";

#
# Connect to the mySQL host
#
my $dbh = DBI->connect( $dsn, $user, $pass, {RaiseError => 1})
  or die "Could not connect to database host : " . DBI->errstr;

print "\nWARNING: If the $destDB database already exists the existing copy \n"
  . "will be destroyed. Proceed (Y/N)? ";

my $key = lc(getc());

unless( $key =~ /y/ ) {
  $dbh->disconnect();
  print "Test Genome Creation Aborted\n";
  exit;
}

print "Proceeding with test genome database $destDB creation\n";  

#
# Create the new database, dropping any existing database
#
$dbh->do("DROP DATABASE $destDB");

$dbh->do( "CREATE DATABASE " . $destDB )
  or die "Could not create database $destDB: " . $dbh->errstr;

#
# Dump the source database table structure (w/o data) and use it to create
# the new database schema
#

# May have to eliminate the -p pass part... not sure

my $rc = 0xffff & system(
  "mysqldump -p$pass -u $user -h $host -P $port --no-data $srcDB | " .
  "mysql -p$pass -u $user -h $host -P port $destDB");

if($rc != 0) {
  $rc >>= 8;
  die "mysqldump and insert failed with return code: $rc";
}  

#
# Create  a temp table to store ids of clones we are interested in
#
$dbh->do("
CREATE TEMPORARY TABLE $destDB.tmp1(
       clone_id INT(10) NOT NULL)
") or die "Could create tmp1 table " . $dbh->errstr;


# Find clones present in a central, exciting 1Mb region of Chromosome 20
$dbh->do("
INSERT INTO $destDB.tmp1
SELECT distinct( c.clone_id )
FROM   $srcDB.contig c, $srcDB.assembly a, $srcDB.chromosome chr
WHERE  a.contig_id = c.contig_id
AND    a.chromosome_id = chr.chromosome_id
AND    chr.name  = '20'
AND    a.chr_end >= 30252000
AND    a.chr_start < 31252001
") or die "Could not do tmp1 chr2 clones insert statement:" . $dbh->errstr;


#Select relevant clones from the source database for the new database
$dbh->do("
INSERT INTO $destDB.clone
SELECT c.*
FROM   $srcDB.clone c, $destDB.tmp1 t
WHERE  c.clone_id = t.clone_id
") or die "Could not do clone insertion statement:" . $dbh->errstr;

#
# Retrieve all contigs on the clones present in the first 1Mbases whether
# or not they are on the golden path
#
$dbh->do("
INSERT INTO $destDB.contig
SELECT c.contig_id, c.name, c.clone_id, c.length, c.embl_offset, c.dna_id
FROM   $srcDB.contig c, $destDB.tmp1 t
WHERE  c.clone_id = t.clone_id
") or die "Could not do contig insertion statement:" . $dbh->errstr;

#
# Create the relevant dna table for the contigs in the test genome
#

$dbh->do("
INSERT INTO $destDB.dna
SELECT d.*
FROM   $srcDB.dna d, $destDB.contig c
WHERE  d.dna_id = c.dna_id
") or die "Could not do dna insertion statement:" . $dbh->errstr;

#
# Copy the entire analysis table (This could be improved I think [mcvicker])
#

$dbh->do("
INSERT INTO $destDB.analysis
SELECT *
FROM $srcDB.analysis
") or die "Could not do analysis insertion statement:" . $dbh->errstr;

#
# Copy the static golden path table
#
$dbh->do("
INSERT INTO $destDB.assembly
SELECT a.*
FROM   $srcDB.assembly a, $srcDB.chromosome chr
WHERE  a.chromosome_id = chr.chromosome_id
AND chr.name = '20'
") or die "Could not do assembly insertion statement for chr2: $dbh->errstr";

#
# Copy overlapping features, repeats, genes, exons, etc
#

# first find gene start end thing

$dbh->do("
CREATE TABLE $destDB.gene_global_start_end
SELECT STRAIGHT_JOIN tr.gene_id,
MIN(IF( a.contig_ori=1,
        (e.contig_start+a.chr_start-a.contig_start),
        (a.chr_start+a.contig_end-e.contig_end)))
   as start,
MAX(IF( a.contig_ori=1,
        (e.contig_end+a.chr_start-a.contig_start),
        (a.chr_start+a.contig_end-e.contig_start)))
   as end,
IF (a.contig_ori=1, e.contig_strand, (-e.contig_strand)) as strand,
chr.name as chromosome

FROM  $srcDB.transcript tr, $srcDB.exon_transcript et,
      $srcDB.exon e, $srcDB.assembly a, $srcDB.chromosome chr
WHERE tr.transcript_id = et.transcript_id
AND   et.exon_id = e.exon_id
AND   e.contig_id = a.contig_id
AND   a.chromosome_id = chr.chromosome_id
GROUP BY tr.gene_id
") or die "Could not do gene_global_start_end table creation statement: $dbh->errstr";

#
# Then make gene list for the area.  Genes have to be completely in
#

$dbh->do("
CREATE TABLE $destDB.gene_list
SELECT gene_id
FROM $destDB.gene_global_start_end g
WHERE g.chromosome = '20'
AND g.end >= 30252000
AND g.start < 31252001
") or die "Could not do gene_list table insertion from chr2: " . $dbh->errstr;

$dbh->do("
ALTER TABLE $destDB.gene_list ADD INDEX gene_id_idx( gene_id )
") or die "Could not prepare genelist table alteration: " . $dbh->errstr;

# Now copy gene

$dbh->do("
INSERT INTO $destDB.gene
SELECT g.* 
FROM $srcDB.gene g, $destDB.gene_list gl
WHERE g.gene_id = gl.gene_id
") or die "Could not do gene insertion statement: " . $dbh->errstr;

# now transcript

$dbh->do("
INSERT INTO $destDB.transcript
SELECT tr.* 
FROM $srcDB.transcript tr, $destDB.gene_list gl
WHERE tr.gene_id = gl.gene_id
") or die "Could not do transcript insertion: " . $dbh->errstr;


# translations

$dbh->do("
INSERT INTO $destDB.translation
SELECT tl.*
FROM $srcDB.translation tl, $destDB.transcript tr
WHERE tr.translation_id = tl.translation_id
") or die "Could not do translation insertion " . $dbh->errstr;


# exon_transcript

$dbh->do("
INSERT INTO $destDB.exon_transcript
SELECT et.*
FROM $srcDB.exon_transcript et, $destDB.transcript tr
WHERE tr.transcript_id = et.transcript_id
") or die "Could not do exon_transcript insertion: " . $dbh->errstr;

# exons are tricky, first make a unique list

$dbh->do ("
CREATE TABLE $destDB.exon_unique
SELECT distinct exon_id
FROM $destDB.exon_transcript
") or die "Could not do exon_unique table create: " . $dbh->errstr;

# then uses this list to copy 

$dbh->do("
INSERT INTO $destDB.exon
SELECT e.*
FROM $srcDB.exon e, $destDB.exon_unique eu
WHERE e.exon_id = eu.exon_id
") or die "Could not do exon insertion: " . $dbh->errstr;


$dbh->do("drop table $destDB.exon_unique")
  or die "Could not drop exon_unique temp table: $dbh->errstr\n";
$dbh->do("drop table $destDB.gene_list")
  or die "Could not drop gene_list temp table: $dbh->errstr\n";
$dbh->do("drop table $destDB.gene_global_start_end")
  or die "Could not drop gene_global_start_end temp table: $dbh->errstr\n";

# meta table

$dbh->do("
INSERT INTO $destDB.meta
SELECT * 
FROM $srcDB.meta
") or die "Could not do meta table insertion: $dbh->errstr\n";

$dbh->do("
INSERT INTO $destDB.chromosome
SELECT * 
FROM $srcDB.chromosome
") or die "Could not do meta table insertion: $dbh->errstr\n";

# object_xref, identity_xref, xref
# gene_description
# external_db

$dbh->do("
INSERT INTO $destDB.gene_description
SELECT gd.* 
FROM $srcDB.gene_description gd, $destDB.gene g
WHERE gd.gene_id = g.gene_id
") or die "Could not do gene description insertion: $dbh->errstr\n";

$dbh->do("
INSERT INTO $destDB.object_xref
SELECT ox.* 
FROM $srcDB.object_xref ox, $destDB.translation tr
WHERE ox.ensembl_id = tr.translation_id
AND   ox.ensembl_object_type = 'Translation'
") or die "Could not do object_xref insertion: $dbh->errstr\n";

$dbh->do("
INSERT INTO $destDB.xref
SELECT x.*
FROM $srcDB.xref x, $destDB.object_xref ox
WHERE x.xref_id = ox.xref_id
") or die "Could not do xref insertion: $dbh->errstr\n";

$dbh->do("
INSERT INTO $destDB.identity_xref 
SELECT ix.* 
FROM $srcDB.identity_xref ix, $destDB.object_xref ox
WHERE ix.object_xref_id = ox.object_xref_id
") or die "Could not do identity insertion: $dbh->errstr\n";

$dbh->do("
INSERT INTO $destDB.external_db
SELECT ed.*
FROM $srcDB.external_db ed
") or die "Could not do external_db insertion: $dbh->errstr\n";

# the interpro we need are those which are in xref


# copy across associated simple features
$dbh->do("
INSERT INTO $destDB.simple_feature
SELECT esf.*
FROM $destDB.contig dc, $srcDB.simple_feature esf
WHERE esf.contig_id = dc.contig_id
") or die "Could not do simple_feature insertion: $dbh->errstr\n";


# copy across associated dna_align_features
$dbh->do("
INSERT INTO $destDB.dna_align_feature
SELECT edna.*
FROM $destDB.contig dc, $srcDB.dna_align_feature edna
WHERE edna.contig_id = dc.contig_id
") or die "Could not do dna_align_feature insertion: $dbh->errstr\n";


# copy across associated protein_align_features
$dbh->do("
INSERT INTO $destDB.protein_align_feature
SELECT eprot.*
FROM $destDB.contig dc, $srcDB.protein_align_feature eprot
WHERE eprot.contig_id = dc.contig_id
") or die "Could not do protein_align_feature insertion: $dbh->errstr\n";


# copy across associated protein_align_features
$dbh->do("
INSERT INTO $destDB.prediction_transcript
SELECT epred.*
FROM $destDB.contig dc, $srcDB.prediction_transcript epred
WHERE epred.contig_id = dc.contig_id
") or die "Could not do prediction_transcript insertion: $dbh->errstr\n";


# copy across associated supporting_features
$dbh->do("
INSERT INTO $destDB.supporting_feature
SELECT esup.*
FROM $destDB.exon dex, $srcDB.supporting_feature esup
WHERE esup.exon_id = dex.exon_id
") or die "Could not do supporting_feature insertion: $dbh->errstr\n";


# copy across associated protein_features
$dbh->do("
INSERT INTO $destDB.protein_feature
SELECT eprot.*
FROM $destDB.translation dt, $srcDB.protein_feature eprot
WHERE eprot.translation_id = dt.translation_id
") or die "Could not do protein_feature insertion: $dbh->errstr\n";


# copy across associated repeat_features
$dbh->do("
INSERT INTO $destDB.repeat_feature
SELECT erep.*
FROM $destDB.contig dc, $srcDB.repeat_feature erep
WHERE erep.contig_id = dc.contig_id
") or die "Could not do protein_align_feature insertion: $dbh->errstr\n";


# copy across associated repeat_consensi
$dbh->do("
INSERT INTO $destDB.repeat_consensus
SELECT econ.*
FROM $destDB.repeat_feature dr, $srcDB.repeat_consensus econ
WHERE econ.repeat_consensus_id = dr.repeat_consensus_id
") or die "Could not do repeat_consensus insertion: $dbh->errstr\n";


# copy across associated karyotype info
$dbh->do("
INSERT INTO $destDB.karyotype
SELECT ek.*
FROM $srcDB.karyotype ek
WHERE ek.chromosome_id = '20'
AND ek.chr_end >= 30252000
AND ek.chr_start < 31252001
") or die "Could not do karyotype insertion: $dbh->errstr\n";


#
# Create a temp table to store superctg_name we are interested in
#
$dbh->do("
CREATE TEMPORARY TABLE $destDB.tmp2(
       superctg_name VARCHAR(20) NOT NULL)
") or die "Could create tmp2 table " . $dbh->errstr;


# grab the unique superctg_names
$dbh->do("
INSERT INTO $destDB.tmp2
SELECT distinct(superctg_name) 
FROM $destDB.assembly
") or die "Could not select the unique superctg_names: $dbh->errstr\n";


# copy across mapfrag
$dbh->do("
INSERT INTO $destDB.mapfrag
SELECT emf.*
FROM $srcDB.mapfrag emf, $destDB.tmp2 sctg
WHERE emf.name = sctg.superctg_name
") or die "Could not do mapfrag insertion: $dbh->errstr\n";


# copy across mapset
$dbh->do("
INSERT INTO $destDB.mapset
SELECT ems.*
FROM $srcDB.mapset ems
") or die "Could not do mapset insertion: $dbh->errstr\n";


# copy across map_density
$dbh->do("
INSERT INTO $destDB.map_density
SELECT emd.*
FROM $srcDB.map_density emd
WHERE emd.chromosome_id = '20'
AND emd.chr_start < 31252001
AND emd.chr_end >= 30252000
") or die "Could not do map_density insertion: $dbh->errstr\n";


# copy across mapannotation
$dbh->do("
INSERT INTO $destDB.mapannotation
SELECT ema.*
FROM $srcDB.mapannotation ema, $destDB.tmp2 sctg
WHERE ema.value = sctg.superctg_name
") or die "Could not do mapannotation insertion: $dbh->errstr\n";


# copy across mapannotationtype
$dbh->do("
INSERT INTO $destDB.mapannotationtype
SELECT emt.*
FROM $srcDB.mapannotationtype emt
") or die "Could not do mapannotationtype insertion: $dbh->errstr\n";


# copy across mapfrag_mapset
$dbh->do("
INSERT INTO $destDB.mapfrag_mapset
SELECT emfms.*
FROM $srcDB.mapfrag_mapset emfms, $destDB.mapfrag dmf
WHERE emfms.mapfrag_id = dmf.mapfrag_id
") or die "Could not do mapfrag_mapset insertion: $dbh->errstr\n";


# copy across dnafrag
$dbh->do("
INSERT INTO $destDB.dnafrag
SELECT edf.*
FROM $srcDB.dnafrag edf
") or die "Could not do dnafrag insertion: $dbh->errstr\n";


# finally, drop the temporary tables
$dbh->do("
DROP TABLE $destDB.tmp1
") or die "Could not drop temporary table tmp1: $dbh->errstr\n";

$dbh->do("
DROP TABLE $destDB.tmp2
") or die "Could not drop temporary table tmp2: $dbh->errstr\n";


#disconnect from the DB
$dbh->disconnect();

print "Test genome database $destDB created\n";

1;
