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
BEGIN { $| = 1; print "1..29\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
use lib 't';
use EnsTestDB;
$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/rawcontig.dump");

print STDERR "** This is not testing get_MarkerFeature or get_all_ExternalFeatures\n";
print STDERR "** This is not testing overlap gets\n";
print STDERR "** This test assummes that prediction features are returned in\n";
print STDERR "** get_all_SimilarityFeature call (to be fixed sometime)\n";
    
# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok 2\n";    
$db->static_golden_path_type('UCSC');

my $c = $db->get_Contig('contig1');

die "$0\nError fetching contig1 : $!" unless defined ($c);

print "ok 3\n";

@genes = $c->get_all_Genes();
if( scalar(@genes) != 1 )  {
     print STDERR "Wrong number of genes returned\n";
     print "not ok 4\n";
} else {
     print "ok 4\n";
}

$gene = shift @genes;

if( $gene->id ne 'gene-id-1' ) {
     print STDERR "Wrong gene id\n";
     print "not ok 5\n"
} else {
     print "ok 5\n";
}

@trans = $gene->each_Transcript();
$trans = shift @trans;
@exons = $trans->each_Exon();

if( $trans->id ne 'transcript-1' || $trans->translation->id ne 'trl-1' ||
    scalar(@exons) != 2 ) {
    print STDERR "Something is wrong with the gene get!",$trans->id,":",$trans->translation->id,":",scalar(@exons),"\n";
    print "not ok 6\n";
} else {
    print "ok 6\n";
}

if( $c->has_genes != 1 ) {
    print STDERR "has genes did not return 1\n";
    print "not ok 7\n";
} else {
    print "ok 7\n";
}

$seq = $c->primary_seq();
if( !$seq->isa('Bio::PrimarySeqI') || $seq->length != $c->length ) {
    print STDERR "Got a $seq with length",$seq->length,":",$c->length,"\n";

    print "not ok 8\n";
} else {
    print "ok 8\n";
}

$dbseq = $c->db_primary_seq();
$pseq  = $c->perl_primary_seq();

if( !$dbseq->isa('Bio::PrimarySeqI') || !$pseq->isa('Bio::PrimarySeqI') ||
    $dbseq->seq ne $pseq->seq ) {
    print STDERR "Something up between db and perl primary seq\n";
    print "not ok 9\n";
} else {
    print "ok 9\n";
}

@features = $c->get_all_SeqFeatures();
if( scalar(@features) != 3 ) {
    print STDERR "Did not get the expected 3 features out from all_SeqFeatures\n";
    print "not ok 10\n";
} else {
    print "ok 10\n";
}


@features = $c->get_all_SimilarityFeatures_above_score('swissprot',80);
$f = shift @features;
if( $f->start != 5 || $f->end != 8 ) {
    print STDERR "Did not get the right sequence feature in get_all_SimilarityFeatures_by_score\n";
    print "not ok 11\n";
} else {
    print "ok 11\n";
}

@features = $c->get_all_SimilarityFeatures();
if( scalar(@features) != 2 ) {
    print STDERR "Did not get the expected 2 features out from all_SimilarityFeatures\n";
    print "not ok 12\n";
} else {
    print "ok 12\n";
}


@features = $c->get_all_RepeatFeatures();
$f = $features[0];
if( scalar(@features) != 1 || $f->start != 5 || $f->end != 8 ) {
    print STDERR "Did not get the expected 1 features out from all_RepeatFeatures",$f->start,":",$f->end,":",scalar(@features),"\n";
    print "not ok 13\n";
} else {
    print "ok 13\n";
}


# get_MarkerFeature needs a co-located map database. Test somewhere else.

@features = $c->get_all_PredictionFeatures();
$f = $features[0];
if( scalar(@features) != 1 || $f->start != 2 || $f->end != 8 ) {
    print STDERR "Did not get the expected 1 features out from all_PredictionFeatures\n";
    print "not ok 14\n";
} else {
    print "ok 14\n";
}

if( $c->cloneid ne 'pog' ) {
    print STDERR "Clone is not pog for test contig\n";
    print "not ok 15\n";
} else {
    print "ok 15\n";
}



# chromosome should not be tested here.

if( $c->seq_version != 3 ) {
    print STDERR "contig does not have seq version 3 for test contig\n";
    print "not ok 16\n";
} else {
    print "ok 16\n";
}


if( $c->embl_offset != 500 ) {
    print STDERR "contig does not have embl offset of 500 for test contig\n";
    print "not ok 17\n";
} else {
    print "ok 17\n";
}


# static golden path tests...

if( $c->chromosome ne 'chr2' ) {
   print "not ok 18\n";
} else {
   print "ok 18\n";
}

if( $c->is_static_golden != 1 ) {
   print "not ok 19\n";
} else {
   print "ok 19\n";
}

if( $c->fpc_contig_name ne 'ctg123' ) {
   print "not ok 20\n";
} else {
   print "ok 20\n";
}

if( $c->fpc_contig_start != 3000 ) {
   print "not ok 21\n";
} else {
   print "ok 21\n";
}

if( $c->fpc_contig_end != 3040 ) {
   print "not ok 22\n";
} else {
   print "ok 22\n";
}

if( $c->chr_start != 122300 ) {
   print "not ok 23\n";
} else {
   print "ok 23\n";
}

if( $c->chr_end != 122338 ) {
   print "not ok 24\n";
} else {
   print "ok 24\n";
}

if( $c->static_golden_start != 2  ) {
   print "not ok 25\n";
} else {
   print "ok 25\n";
}

if( $c->static_golden_end != 40  ) {
   print "not ok 26\n";
} else {
   print "ok 26\n";
}

if( $c->static_golden_ori != 1  ) {
   print "not ok 27\n";
} else {
   print "ok 27\n";
}

if( $c->static_golden_type ne 'UCSC'  ) {
   print "not ok 28\n";
} else {
   print "ok 28\n";
}

$c->length();


$seq = $c->get_repeatmasked_seq();

$c->length();


if  (!$seq->isa('Bio::PrimarySeqI') 
        || $seq->length != $c->length 
        || $seq->seq !~ /N/i) 
{
    print STDERR "Got a $seq with length",$seq->length,":",$c->length,"\n";
    print "not ok 29\n";
} 
else 
{
    print "ok 29\n";
}
  
$db = undef;










