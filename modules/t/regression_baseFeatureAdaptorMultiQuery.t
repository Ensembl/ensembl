use strict;
use warnings;
use Test::More;
use Bio::EnsEMBL::Test::MultiTestDB;
use DBI qw/:sql_types/;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $multi->get_DBAdaptor( "core" );

my $sa = $db->get_SliceAdaptor();
my $ga = $db->get_GeneAdaptor();

my $name  = '20_HAP1';
my $start = undef;
my $end = undef;
my $strand = 1;

my $slice = $sa->fetch_by_region('toplevel', $name, $start, $end, $strand);
$ga->bind_param_generic_fetch('protein_coding', SQL_VARCHAR);
my $genes = $ga->fetch_all_by_Slice_constraint($slice, 'g.biotype =?');

# Logic here is that we divide our slice into 3; prior HAP (Chr20), HAP & post HAP (Chr20)
# and each one needs the bind parameters set for each query. By the time we
# get into BaseFeatureAdaptor's internals we can record and replay them everytime
# we run the same query on a different Slice

is(scalar(@{$genes}), 14, 'Should have 14 genes when we query slices which need normalisation');

done_testing();