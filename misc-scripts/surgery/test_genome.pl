# Copyright EMBL-EBI 2001
# Author: Graham McVicker - Based on SQL script by Alistair Rust
# Creation: 04.07.2002
# Last modified: 04.07.2002
#
# File name: test_genome.pl
#
# Given a source and destination database this script will create 
# a ensembl database and perform 2Mbases of test data insertion 
#


use strict;

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


# Find clones present in a central, exciting 1Mb region of Chromosome 2
$dbh->do("
INSERT INTO $destDB.tmp1
SELECT distinct( c.clone_id )
FROM   $srcDB.contig c, $srcDB.assembly a, $srcDB.chromosome chr
WHERE  a.contig_id = c.contig_id
AND    a.chromosome_id = chr.chromosome_id
AND    chr.name  = '2'
AND    chr_start < 224000000
AND    chr_end > 223000001
") or die "Could not do tmp1 chr2 clones insert statement:" . $dbh->errstr;

#Get some sticky exons from chromosome 19
$dbh->do("
INSERT INTO $destDB.tmp1
SELECT distinct(c.clone_id)
FROM   $srcDB.contig c, $srcDB.assembly a, $srcDB.chromosome chr
WHERE  a.contig_id = c.contig_id
AND    chr.name = '19'
AND    chr_end > 57000000
AND    chr_start < 58000000
") or die "Could not do tmp1 chr19 clone insert stmnt:" . $dbh->errstr;

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
SELECT c.contig_id, c.name, c.clone_id, c.length, c.offset, c.corder, c.dna_id,
       c.international_name
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
AND chr.name = '2'
") or die "Could not do assembly insertion statement for chr2: $dbh->errstr";

$dbh->do("
INSERT INTO $destDB.assembly
SELECT a.*
FROM   $srcDB.assembly a, $srcDB.chromosome chr
WHERE  a.chromosome_id = chr.chromosome_id
AND chr.name = '19'
") or die "Could not do assembly insertion statement for chr19: $dbh->errstr";

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
WHERE g.chromosome = '19'
AND g.end < 58000000
AND g.start > 57000000
") or die "Could not do gene_list table create statement: $dbh->errstr";

$dbh->do("
INSERT INTO $destDB.gene_list
SELECT gene_id
FROM $destDB.gene_global_start_end g
WHERE g.chromosome = '2'
AND g.start >= 223000000
AND g.end < 224000000
") or die "Could not do gene_list table insertion from chr2: " . $dbh->errstr;

$dbh->do("
ALTER TABLE $destDB.gene_list ADD INDEX gene_id_idx( gene_id )
") or die "Could not prepare genelist table alteration: " . $dbh->errstr;

# Now copy gene and stable id

$dbh->do("
INSERT INTO $destDB.gene
SELECT g.* 
FROM $srcDB.gene g, $destDB.gene_list gl
WHERE g.gene_id = gl.gene_id
") or die "Could not do gene insertion statement: " . $dbh->errstr;

$dbh->do("
INSERT INTO $destDB.gene_stable_id
SELECT gsi.*
FROM $srcDB.gene_stable_id gsi, $destDB.gene_list gl
WHERE gsi.gene_id = gl.gene_id
") or die "Could not do gene_stable_id: " . $dbh->errstr;

# now transcript, transcript_stable_id

$dbh->do("
INSERT INTO $destDB.transcript
SELECT tr.* 
FROM $srcDB.transcript tr, $destDB.gene_list gl
WHERE tr.gene_id = gl.gene_id
") or die "Could not do transcript insertion: " . $dbh->errstr;

$dbh->do("
INSERT INTO $destDB.transcript_stable_id
SELECT tsi.*
FROM $srcDB.transcript_stable_id tsi, $destDB.transcript tr
WHERE tsi.transcript_id = tr.transcript_id
") or die "Could not do transcript_stable_id insertion: " . $dbh->errstr;

# translations, translation_stable_id

$dbh->do("
INSERT INTO $destDB.translation
SELECT tl.*
FROM $srcDB.translation tl, $destDB.transcript tr
WHERE tr.translation_id = tl.translation_id
") or die "Could not do translation insertion " . $dbh->errstr;

$dbh->do("
INSERT INTO $destDB.translation_stable_id
SELECT tsi.*
FROM $srcDB.translation_stable_id tsi, $destDB.translation tl
WHERE tsi.translation_id = tl.translation_id
") or die "Could not do translation_stable_id insertion" . $dbh->errstr;

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

$dbh->do("
INSERT INTO $destDB.exon_stable_id
SELECT esi.*
FROM $srcDB.exon_stable_id esi, $destDB.exon_unique eu
WHERE esi.exon_id = eu.exon_id
") or die "Could not do exon_stable_id insertion: " . $dbh->errstr;

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







# finally, drop the temporary table
$dbh->do("
DROP TABLE $destDB.tmp1
") or die "Could not drop temporary table tmp1: $dbh->errstr\n";

#disconnect from the DB
$dbh->disconnect();

print "Test genome database $destDB created\n";

1;
