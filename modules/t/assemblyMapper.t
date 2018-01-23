# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;

use Test::More;
use Test::Warnings qw(warning);

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0; #set to 1 to turn on debug printouts

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
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
my $chnk_cs = $csa->fetch_by_name('chunk');

my $sa = $db->get_SliceAdaptor();
my $slice = $sa->fetch_by_region('chromosome', '1');

my $asm_mapper =  $asma->fetch_by_CoordSystems($chnk_cs, $chr_cs);
ok( $asm_mapper && $asm_mapper->isa( "Bio::EnsEMBL::ChainedAssemblyMapper" ));

# Test full caching
$asm_mapper->register_all();

#
# Test if the multi mapping works (meta_key=assembly.mapping entry with #)
#

my @coords =  $asm_mapper->map( '1', 1, 50, 1, $chr_cs );

ok( $coords[0]->id() == 965905);
ok( $coords[0]->start() == 10 );
ok( $coords[0]->end() == 59 );
ok( $coords[0]->strand() == 1 );
is( $coords[0]->name(), "multimap_testregion" ); # [ENSCORESW-844]

# Test seq_ids_to_regions with empty cache
$asma->delete_cache();
my @seq_ids = (965905);
$asma->seq_ids_to_regions(\@seq_ids);

@coords = $asm_mapper->map( '1', 1, 50, 1, $chr_cs, $slice );
is(scalar(@coords), 5, "Chained mapping with target slice");

@coords = $asm_mapper->map( "multimap_testregion", 100, 800, 1, $chnk_cs );

ok( $coords[0]->id() == 469271 );  #seq_region_id not name now.
ok( $coords[0]->start() == 91 );
ok( $coords[0]->end() == 200 );
ok( $coords[0]->strand() == 1 );
is( $coords[0]->name(), "1" ); # [ENSCORESW-844]

ok( $coords[1]->isa( "Bio::EnsEMBL::Mapper::Gap" ) );

ok( $coords[2]->id() == 469271);
ok( $coords[2]->start() == 201 );
ok( $coords[2]->end() == 400 );
ok( $coords[2]->strand() == -1 );
is( $coords[2]->name(), "1" ); # [ENSCORESW-844]

ok( $coords[4]->id() == 469282);
ok( $coords[4]->start() == 1 );
ok( $coords[4]->end() == 100 );
ok( $coords[4]->strand() == -1 );
is( $coords[4]->name(), "2" ); # [ENSCORESW-844]

# Test mapping with slice and identical coord system
# will warn about using an implicit mapping
warning { $asm_mapper = $asma->fetch_by_CoordSystems($chr_cs, $chr_cs); };
@coords =  $asm_mapper->map( '1', 1, 50, 1, $chr_cs, 0, $slice );
is( scalar(@coords), 1, "Found 1 mapping when specifying a slice");

$asm_mapper = $asma->fetch_by_CoordSystems($ctg_cs, $chr_cs);
$asm_mapper->register_all();

ok($asm_mapper && $asm_mapper->isa('Bio::EnsEMBL::AssemblyMapper'));


#
# test db has chr 20  (50KB -> 62MB)
#

#
# 3 Test map
#

@coords = $asm_mapper->map('20', 500_001, 60_000_000, 1, $chr_cs);
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

done_testing();
