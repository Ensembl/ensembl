use lib 't';
use strict;
use warnings;


BEGIN { $| = 1;  
	use Test;
	plan tests => 20;
}

use MultiTestDB;
use Bio::EnsEMBL::DBSQL::RawContigAdaptor;
use Bio::EnsEMBL::RawContig;
use TestUtils qw(test_getter_setter);
use Bio::Seq;

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
# 6-10 fetch_all_by_Clone
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

#
# 11 fetch_by_all
#

my $retrieved_contigs = $raw_adaptor->fetch_all;
ok(scalar (@$retrieved_contigs) == 12);


#
# 12 fetch_filled_by_dbIDs
#

my @contig_list = (368744, 317101);
my $filled_contigs = $raw_adaptor->fetch_filled_by_dbIDs(@contig_list);

# Still need to process these things.....



#
# methods which may well become deprecated
#
# 13 get_internal_id_by_id
#

my $internal_name = 'AL390298.13.1.49208'; 
my $internal_id = $raw_adaptor->get_internal_id_by_id($internal_name);
ok($internal_id == 376992);


#
# 14 get_id_by_contig_id
#

my $contig_id = 317101;
my $get_name = $raw_adaptor->get_id_by_contig_id($contig_id);
ok($get_name eq 'AL031658.11.1.162976');



#
# XX check out the store
#

my $dname = 'dummy_contig';
my $dummy_contig = Bio::EnsEMBL::RawContig->new;
$dummy_contig->name($dname);
$dummy_contig->embl_offset(7);
$dummy_contig->length(24);

my $seq  = Bio::Seq->new(-seq => 'ATGCAGCTAGCATCGATGACATCG',
                         -id => 'dummy_contig',
                         -accession => 'dummy_contig');
ok($seq);

$dummy_contig->seq($seq);


$raw_adaptor->store($dummy_contig, $clone);
ok($raw_adaptor);

#
# manual check to see whether the contig has gone in
#
my $sth = $db->prepare("select * from contig");
$sth->execute;
#print STDERR "Num contigs " . scalar($sth->rows) . "\n";
ok(scalar($sth->rows) == 13);


#
# and just to check, retrieve the stored items
#

my $stored_contig = $raw_adaptor->fetch_by_name($dname);
ok($stored_contig->name eq $dname);
ok($stored_contig->embl_offset == 7);
ok($stored_contig->length == 24);
ok($stored_contig->clone->isa('Bio::EnsEMBL::Clone'));