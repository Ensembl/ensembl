use lib 't';
use strict;

BEGIN { $| = 1;
	use Test ;
	plan tests => 9;
}

use MultiTestDB;
use TestUtils qw(debug test_getter_setter);

our $verbose = 0; #set to 1 to turn on debug printouts

my $multi = MultiTestDB->new();
my $db = $multi->get_DBAdaptor( 'core' );

#
# 1 Test AssemblyMapperAdaptor constructor
#
my $asma = $db->get_AssemblyMapperAdaptor();


ok($asma && $asma->isa('Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor'));


#
# 2 Test fetch_by_CoordSystems
#

my $csa = $db->get_CoordSystemAdaptor();
my $ctg_cs = $csa->fetch_by_name('contig');
my $chr_cs = $csa->fetch_by_name('chromosome');
my $cln_cs = $csa->fetch_by_name('clone');

my $asm_mapper = $asma->fetch_by_CoordSystems($ctg_cs, $chr_cs);

ok($asm_mapper && $asm_mapper->isa('Bio::EnsEMBL::AssemblyMapper'));


#
# test db has chr 20  (50KB -> 62MB)
#

#
# 3 Test map
#

my @coords = $asm_mapper->map('20', 500_001, 60_000_000, 1, $chr_cs);
ok(@coords);
print_coords(@coords);


@coords = $asm_mapper->map('AL359765.6.1.13780', 1, 13780, 1, $ctg_cs);
ok(@coords);
print_coords(@coords);

my $cln_mapper = $asma->fetch_by_CoordSystems($cln_cs, $ctg_cs);
@coords = $cln_mapper->map('AL359765.6', 1, 20_000, 1, $cln_cs);
ok(@coords);
print_coords(@coords);


#
# Test list_seq_regions
#

my @seq_regions =
  $asm_mapper->list_seq_regions('20', 500_001, 60_000_000, $chr_cs);
ok(@seq_regions);
my $str = join("\n", "----------", @seq_regions);
debug("$str\n");

@seq_regions =
  $asm_mapper->list_seq_regions('AL359765.6.1.13780', 1, 13780, $ctg_cs);
ok(@seq_regions);
$str = join("\n", "----------", @seq_regions);
debug("$str\n");


#
# Test list_seq_ids
#


my @seq_ids =
  $asm_mapper->list_ids('20', 500_001, 60_000_000, $chr_cs);
ok(@seq_ids);
$str = join("\n", "----------", @seq_ids);
debug("$str\n");

@seq_ids =
  $asm_mapper->list_ids('AL359765.6.1.13780', 1, 13780, $ctg_cs);
ok(@seq_ids);
$str = join("\n", "----------", @seq_ids);
debug("$str\n");



sub print_coords {
  my @coord_list = @_;

  return if(!$verbose);

  foreach my $coord (@coord_list) {
    if($coord->isa('Bio::EnsEMBL::Mapper::Gap')) {
      debug("GAP");
      next;
    }
    debug($coord->id()."\t". $coord->start()."-".$coord->end().
          " (".$coord->strand.")");
  }
}

