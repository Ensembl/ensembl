use lib 't';
use strict;
use warnings;


BEGIN { $| = 1;  
	use Test;
	plan tests => 8;
}

use MultiTestDB;
use Bio::EnsEMBL::DBSQL::RawContigAdaptor;
use TestUtils qw(test_getter_setter);

my ($CHR, $START, $END, $FLANKING) = ("20", 30_252_000, 31_252_001, 1000);

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
# 12-13 fetch_all_by_Clone
#

my $clone_name = 'AL359765';
my $clone_adaptor = $db->get_CloneAdaptor;
my $clone = $clone_adaptor->fetch_by_name($clone_name);

ok($clone);

my $contigs_from_clone = $raw_adaptor->fetch_all_by_Clone($clone);

ok(scalar(@$contigs_from_clone) == 1);


#
# 12-13 fetch_by_all
#

my $retrieved_contigs = $raw_adaptor->fetch_all;
ok(scalar (@$retrieved_contigs) == 12);
