use lib 't';
use strict;
use warnings;


BEGIN { $| = 1;  
	use Test;
	plan tests => 26;
}

use MultiTestDB;
use Bio::EnsEMBL::DBSQL::RawContigAdaptor;
use Bio::EnsEMBL::RawContig;
use TestUtils qw(test_getter_setter);
use Bio::Seq;

use Bio::EnsEMBL::Utils::Exception qw(verbose);

######################################################################
# 
# RawContigAdaptor is a deprecated class but needed for backwards 
# compatibility.  These tests ensure that it actually works,
# but verbosity is turned off to avoid all of the deprecated warnings
#
#######################################################################

verbose(-1);


#
#1 slice adaptor compiles
#
ok(1);

my $multi = MultiTestDB->new;
my $db    = $multi->get_DBAdaptor('core');


#
# 2-3 RawContigAdaptor::new
#
my $raw_adaptor = Bio::EnsEMBL::DBSQL::RawContigAdaptor->new($db->_obj);
ok($raw_adaptor->isa('Bio::EnsEMBL::DBSQL::RawContigAdaptor'));
ok($raw_adaptor->db);


#
# 4 fetch_by_name
#
my $contig_name = 'AL031658.11.1.162976';
my $contig = $raw_adaptor->fetch_by_name($contig_name);

ok($contig->name eq $contig_name);


#
# 5 fetch_by_name
#
my $dbID = 368744;
my $related_name = 'AL359765.6.1.13780';
my $contig2 = $raw_adaptor->fetch_by_dbID($dbID);

ok($contig2->name eq $related_name);


#
# 6-11 fetch_all_by_Clone
#

my $clone_name = 'AL359765';
my $clone_adaptor = $db->get_CloneAdaptor;
my $clone = $clone_adaptor->fetch_by_name($clone_name);

ok($clone->isa('Bio::EnsEMBL::Clone'));

my $contigs_from_clone = $raw_adaptor->fetch_all_by_Clone($clone);

ok(scalar(@$contigs_from_clone) == 1);

my $contig3 = $contigs_from_clone->[0];

ok ($contig3->dbID == 368744);
ok ($contig3->name eq "AL359765.6.1.13780");
ok ($contig3->length == 13780);
ok ($contig3->clone->name eq "AL359765.6");


#
# 12-16 fetch_attributes
#
# fetch as an unyet, uncached contig
my $unpop_contig = $raw_adaptor->fetch_by_dbID(469270);

$raw_adaptor->fetch_attributes($unpop_contig);
ok($unpop_contig->name eq 'AL133343.23.1.63609');
ok($unpop_contig->length == 63609);
ok($unpop_contig->dbID == 469270);
ok($unpop_contig->embl_offset == 1);
ok($unpop_contig->clone->name eq "AL133343.23");


#
# 17 fetch_by_all
#

my $retrieved_contigs = $raw_adaptor->fetch_all;
ok(scalar (@$retrieved_contigs) == 12);

my @names = @{$raw_adaptor->fetch_all_names};

ok(scalar(@names) == 12);
#
# 18-25 fetch_filled_by_dbIDs
#

my @contig_list = (368744, 317101);
my $filled_contigs = $raw_adaptor->fetch_filled_by_dbIDs(@contig_list);

# should retrieve just two contigs
ok(scalar(keys %$filled_contigs) == 2);
ok($$filled_contigs{'368744'});
ok($$filled_contigs{'317101'});

# check the contents of one of the retrieved contigs
my $first_contig = @$filled_contigs{'368744'};
ok($first_contig->name eq 'AL359765.6.1.13780');
ok($first_contig->length == 13780);
ok($first_contig->dbID == 368744);
ok($first_contig->embl_offset == 1);
ok($first_contig->clone->name eq "AL359765.6");




##
## 26-32 check out the store
##

## save the original state of the contig table
#$multi->save("core","contig","dna");

#my $dname = 'dummy_contig';
#my $dummy_contig = Bio::EnsEMBL::RawContig->new;
#$dummy_contig->name($dname);
#$dummy_contig->embl_offset(7);
#$dummy_contig->length(24);

#my $seq  = Bio::Seq->new(-seq => 'ATGCAGCTAGCATCGATGACATCG',
#                         -id => 'dummy_contig',
#                         -accession => 'dummy_contig');
#ok($seq);

#$dummy_contig->seq($seq);


#$raw_adaptor->store($dummy_contig, $clone);
#ok($raw_adaptor);

##
## manual check to see whether the contig has gone in
##
#my $sth = $db->prepare("select * from contig");
#$sth->execute;
##print STDERR "Num contigs " . scalar($sth->rows) . "\n";
#ok(scalar($sth->rows) == 13);


##
## and just to check, retrieve the stored items
##

#my $stored_contig = $raw_adaptor->fetch_by_name($dname);
#ok($stored_contig->name eq $dname);
#ok($stored_contig->embl_offset == 7);
#ok($stored_contig->length == 24);
#ok($stored_contig->clone->isa('Bio::EnsEMBL::Clone'));


## restore the contig table
#$multi->restore("core","contig","dna");



##
## 33-40 remove
##

## save the original state of the contig table
#$multi->save("core","contig","dna","repeat_feature","simple_feature",
#             "dna_align_feature","protein_align_feature","prediction_transcript");

## remove the contig AL031658.11.1.162976 from the db
#$raw_adaptor->remove($contig);
#ok($raw_adaptor);

## manual check to see whether the contig and associated
##features have gone in
##

#$sth = $db->prepare("select * from contig");
#$sth->execute;
#ok(scalar($sth->rows) == 11);

#$sth = $db->prepare("select * from dna");
#$sth->execute;
#ok(scalar($sth->rows) == 11);

## contig should have 419 records
#$sth = $db->prepare("select * from repeat_feature");
#$sth->execute;
#ok(scalar($sth->rows) == 1937);

## contig should have 20 records
#$sth = $db->prepare("select * from simple_feature");
#$sth->execute;
#ok(scalar($sth->rows) == 116);

## contig should have 11641 records
#$sth = $db->prepare("select * from dna_align_feature");
#$sth->execute;
#ok(scalar($sth->rows) == 15525);

## contig should have 2507 records
#$sth = $db->prepare("select * from protein_align_feature");
#$sth->execute;
#ok(scalar($sth->rows) == 4727);

## contig should have 30 records
#$sth = $db->prepare("select * from prediction_transcript");
#$sth->execute;
#ok(scalar($sth->rows) == 161);

## restore the original state of the contig table
#$multi->restore("core","contig","dna","repeat_feature","simple_feature",
#             "dna_align_feature","protein_align_feature","prediction_transcript");



verbose(0);