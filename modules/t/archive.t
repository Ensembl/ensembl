## Bioperl Test Harness Script for Modules
##
# CVS Version
# $Id$


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..5\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}
use strict;
use lib 't';
use EnsTestDB;
use Bio::EnsEMBL::DBArchive::Obj;

$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("../sql/archive.sql");
$ens_test->do_sql_file("t/archive.dump");


my $host = $ens_test->host;
my $dbname = $ens_test->dbname;
my $user = $ens_test->user;
my $pass = $ens_test->pass;

my $archive = Bio::EnsEMBL::DBArchive::Obj->new( -host => $host, -dbname => $dbname, -user => $user, -pass => $pass);


print "ok 2\n";

my $seq = Bio::Seq->new( -id => 'silly',-seq => 'ATTCGTTGGGTGGCCCGTGGGTG');

$archive->write_seq($seq,1,'exon','ENSG000000012',1,'AC00012',2);

print "ok 3\n";

my ($seq2) = $archive->get_seq_by_gene_version('ENSG000000012',1,'exon');


if( !defined $seq2 || $seq2->id ne "silly.1" || $seq2->seq ne $seq->seq ) {
  print "not ok 4\n";
  print STDERR "Got ",$seq2->id," with ",$seq2->seq,"\n"
} else {
  print "ok 4\n";
}

my $new_id = $archive->get_new_id_from_old_id('gene','ENSG0000018');

if( $new_id ne 'ENSG0000019' ) {
    print "not ok 5\n";
    print STDERR "Got $new_id for new id\n";
} else {
  print "ok 5\n";
}



