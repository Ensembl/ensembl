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
BEGIN { $| = 1; print "1..22\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use lib 't';
use EnsTestDB;
use Bio::EnsEMBL::Archive::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Archive::Seq;
use Bio::EnsEMBL::Archive::VersionedSeq;

$loaded = 1;
print "ok 1\n";    # 1st test passes.

my $hash = {
    'schema_sql'    => ['../sql/archive_new.sql'],
    'module'        => 'Bio::EnsEMBL::Archive::DBSQL::DBAdaptor'
    };

my $ens_test = EnsTestDB->new($hash);
my $db = $ens_test->get_DBSQL_Obj();
print "ok 2\n";

my $seq = Bio::EnsEMBL::Archive::Seq->new(
					  -name => 'test_t',
					  -type => 'transcript',
					  -created => '2001-06-28 13:45:00'
					  );

my $vseq = Bio::EnsEMBL::Archive::VersionedSeq->new(
						    -archive_seq => $seq,
						    -version => 1,
						    -start_clone => 'test_clone',
						    -start => 1,
						    -end_clone => 'test_clone',
						    -end => 9,
						    -sequence => 'ATGCGTGTG',
						    -modified => '2001-06-28 13:45:00',
						    -release_number => 100
						    );


my $seq2 = Bio::EnsEMBL::Archive::Seq->new(
					  -name => 'test_e',
					  -type => 'exon',
					  -created => '2001-06-28 13:45:00'
					  );

my $vseq2 = Bio::EnsEMBL::Archive::VersionedSeq->new(
						    -archive_seq => $seq2,
						    -version => 1,
						    -start_clone => 'test_clone',
						    -start => 1,
						    -end_clone => 'test_clone',
						    -end => 4,
						    -sequence => 'ATGC',
						    -modified => '2001-06-28 13:45:00',
						    -release_number => 100
						    );




$vseq->add_relative($vseq2);

my $vseq3 = Bio::EnsEMBL::Archive::VersionedSeq->new(
						    -archive_seq => $seq,
						    -version => 2,
						    -start_clone => 'test_clone',
						    -start => 1,
						    -end_clone => 'test_clone',
						    -end => 11,
						    -sequence => 'ATGCGTGTGTG',
						    -modified => '2001-06-28 13:45:00',
						    -release_number => 100
						    );

$vseq->add_future_vseq($vseq3);

my $vsda = $db->get_VersionedSeqAdaptor;

my $id = $vsda->store($vseq);
print "ok 3\n";
my $vs = $vsda->fetch_by_dbID($id);
print "ok 4\n";

#Test all Seq methods
if ($vs->archive_seq->name eq 'test_t') {
    print "ok 5\n";
}
else {
    print "not ok 5\n";
}

if ($vs->archive_seq->type eq 'transcript') {
    print "ok 6\n";
}
else {
    print "not ok 6\n";
}

if ($vs->archive_seq->created eq '2001-06-28 13:45:00') {
    print "ok 7\n";
}
else {
    print "not ok 7\n";
}

#Now check all VersionedSeq methods
if ($vs->version == 1) {
    print "ok 8\n";
}
else {
    print "not ok 8\n";
}

if ($vs->start_clone eq 'test_clone') {
    print "ok 9\n";
}
else {
    print "not ok 9\n";
}

if ($vs->start == 1) {
    print "ok 10\n";
}
else {
    print "not ok 10\n";
}

if ($vs->end_clone eq 'test_clone') {
    print "ok 11\n";
}
else {
    print "not ok 11\n";
}

if ($vs->end == 9) {
    print "ok 12\n";
}
else {
    print "not ok 12\n";
}

if ($vs->modified eq '2001-06-28 13:45:00') {
    print "ok 13\n";
}
else {
    print "not ok 13\n";
}

if ($vs->release_number == 100) {
    print "ok 14\n";
}
else {
    print "not ok 14\n";
}

if ($vs->seq eq 'ATGCGTGTG') {
    print "ok 15\n";
}
else {
    print "not ok 15\n";
}

if ($vs->subseq(1,3) eq 'ATG') {
    print "ok 16\n";
}
else {
    print "not ok 16\n";
}

if ($vs->length == 9) {
    print "ok 17\n";
}
else {
    print "not ok 17\n";
}

if ($vs->moltype eq 'dna') {
    print "ok 18\n";
}
else {
    print "not ok 18\n";
}

if ($vs->display_id eq 'test_t') {
    print "ok 19\n";
}
else {
    print "not ok 19\n";
}

if ($vseq->translate->seq eq 'MRV') {
    print "ok 20\n";
}
else {
    print "not ok 20\n";
}

foreach my $r ($vseq->each_relative) {
    if ($r->archive_seq->name eq 'test_e' && $r->seq eq 'ATGC') {
	print "ok 21\n";
    }
    else {
	print "not ok 21\n";
    }
}

foreach my $f ($vseq->each_future_vseq) {
    if ($f->version == 2 && $f->seq eq 'ATGCGTGTGTG') {
	print "ok 22\n";
    }
    else {
	print "not ok 22\n";
    }
}
