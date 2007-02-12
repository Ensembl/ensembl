use strict;

BEGIN { $| = 1;
	use Test ;
	plan tests => 37;
}

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 1; #set to 1 to turn on debug printouts

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $multi->get_DBAdaptor( 'core' );


#
# Test constructor
#

my $csa = $db->get_CoordSystemAdaptor();

ok($csa && $csa->isa('Bio::EnsEMBL::DBSQL::CoordSystemAdaptor'));


#
# Test fetch_by_name()
#

my $cs = $csa->fetch_by_name('chromosome');

ok($cs->name eq 'chromosome');
ok($cs->dbID());
ok($cs->version eq 'NCBI33');
ok(!$cs->is_top_level());
ok(!$cs->is_sequence_level());
ok($cs->is_default());


#
#  Test fetch_all_by_name
#
my @cs_list = @{$csa->fetch_all_by_name('chromosome')};

ok(@cs_list == 1);

ok($cs_list[0]->equals($cs));



#
# Test fetch_by_dbID()
#
$cs = $csa->fetch_by_dbID(3);

ok($cs->name() eq 'clone');
ok($cs->version() eq '');



#
# Test fetch_top_level
#
$cs = $csa->fetch_top_level();

ok($cs->name eq 'toplevel');
ok($cs->is_top_level());
ok($cs->rank == 0);

#
# Test fetch_by_rank
#
$cs = $csa->fetch_by_rank(1);
ok($cs->name() eq 'chromosome' && $cs->rank() == 1);

$cs = $csa->fetch_by_rank(0);
ok($cs->name() eq 'toplevel' && $cs->rank() == 0);


#
# Test fetch_sequence_level
#
$cs = $csa->fetch_sequence_level();

ok($cs->name eq 'contig');
ok($cs->is_sequence_level());


#
# Test fetch_all
#
@cs_list = @{$csa->fetch_all()};
my $prev_cs;

#make sure that they are ordered by rank
foreach my $cs (@cs_list) {
  if($prev_cs) {
    ok($prev_cs->rank < $cs->rank);
  }
  $prev_cs = $cs;
}


#
# Test get_mapping_path
#

my $ctg_cs = $csa->fetch_by_name('contig');
my $chr_cs = $csa->fetch_by_name('chromosome');
my $cln_cs = $csa->fetch_by_name('clone');

my $path = $csa->get_mapping_path($ctg_cs, $chr_cs);

ok(@$path == 2 &&
   $path->[0]->name() eq 'chromosome' &&
   $path->[1]->name() eq 'contig');

$path = $csa->get_mapping_path($chr_cs, $cln_cs);


ok(@$path == 3 &&
   (($path->[0]->name eq 'chromosome' &&  #there are 2 equally valid paths
     $path->[1]->name eq 'contig' &&
     $path->[2]->name eq 'clone')
    ||
    ($path->[0]->name eq 'clone' &&
     $path->[1]->name eq 'contig' &&
     $path->[2]->name eq 'chromosome'))
   );

#
# Test store
#

$multi->save('core', 'coord_system');
$multi->save('core', 'meta');

$cs = Bio::EnsEMBL::CoordSystem->new
  (-NAME            => 'newsystem',
   -VERSION         => 'NCBI35',
   -DEFAULT         => 1,
   -SEQUENCE_LEVEL  => 0,
   -RANK            => 10);

$csa->store($cs);

ok($cs->adaptor == $csa);
ok($cs->dbID());

#now make sure we can retrieve this
$cs = $csa->fetch_by_name('newsystem', 'NCBI35');
ok($cs->name eq 'newsystem');
ok($cs->version eq 'NCBI35');
ok($cs->is_default);
ok(!$cs->is_sequence_level);
ok(!$cs->is_top_level);
ok($cs->rank() == 10);

my $sth = $db->dbc->prepare('SELECT attrib FROM coord_system ' .
                       'WHERE  name = ? and version = ?');
$sth->execute('newsystem', 'NCBI35');

my ($attrib) = $sth->fetchrow_array();
ok($attrib eq 'default_version');
$sth->finish();

#
# Test store_mapping_path
#

my $new_paths = $csa->store_mapping_path( $cs, $cln_cs );
ok( @{$new_paths} == 1 &&
    $new_paths->[0] eq 'clone|newsystem:NCBI35' );

my $new_path = $csa->get_mapping_path( $cs, $cln_cs );
ok( @{$new_path} == 2 &&
    $new_path->[0]->name eq 'clone' &&
    $new_path->[1]->name eq 'newsystem' );

my $new_paths2 = $csa->store_mapping_path( $cs, $cln_cs );
ok( @{$new_paths2} == 0 ); # Should not update second time round

$multi->restore('core', 'coord_system');
$multi->restore('core', 'meta');









