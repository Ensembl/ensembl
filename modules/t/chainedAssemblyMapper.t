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
use Test::Warnings;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0; # set to 1 to turn on debug printouts

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
my $chr_cs = $csa->fetch_by_name('chromosome');
my $cln_cs = $csa->fetch_by_name('clone');
my $sctg_cs = $csa->fetch_by_name('supercontig');

my $asm_mapper = $asma->fetch_by_CoordSystems($cln_cs, $chr_cs);

ok($asm_mapper && $asm_mapper->isa('Bio::EnsEMBL::ChainedAssemblyMapper'));

my $chr_sctg_mapper = $asma->fetch_by_CoordSystems($chr_cs, $sctg_cs);
ok($chr_sctg_mapper &&
   $chr_sctg_mapper->isa('Bio::EnsEMBL::ChainedAssemblyMapper'));




#
# test db has chr 20  (50KB -> 62MB)
#

#
# 3 Test map
#

my @coords = $asm_mapper->map('20', 500_001, 60_000_000, 1, $chr_cs);
ok(@coords);
debug("MAP 20->clone\n");
print_coords(@coords);
# [ENSCORESW-844]. Test mapped coordinate names
my @coord_names = qw ( AL359765.6 AL031658.11 AL353092.6 AL049539.21 AL121897.32 AL354800.4 AL121583.25 AL034550.31 AL133343.23 AL132653.22 AL035071.17 AL390298.13);
for my $coord (@coords) {
  next if $coord->isa('Bio::EnsEMBL::Mapper::Gap');
  is($coord->name(), shift @coord_names);
}

debug("MAP 'AL359765.6'->chromosome\n");
@coords = $asm_mapper->map('AL359765.6', 1, 13780, 1, $cln_cs);
ok(@coords);
print_coords(@coords);

debug("MAP 20->supercontig\n");
@coords = $chr_sctg_mapper->map('20', 500_001, 60_000_000, 1, $chr_cs);
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
  $asm_mapper->list_seq_regions('AL359765.6', 1, 13780, $cln_cs);
ok(@seq_regions);
$str = join("\n", "----------", @seq_regions);
debug("$str\n");

@seq_regions = 
  $chr_sctg_mapper->list_seq_regions('NT_028392',600_000, 1_000_000, $sctg_cs);
ok(@seq_regions);
$str = join("\n", "----------", @seq_regions);
debug("$str\n");


@seq_regions = 
  $chr_sctg_mapper->list_seq_regions('20', 30_000_000, 31_000_000, $chr_cs);
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
  $asm_mapper->list_ids('AL359765.6', 1, 13780, $cln_cs);
ok(@seq_ids);
$str = join("\n", "----------", @seq_ids);
debug("$str\n");

@seq_ids =
  $asm_mapper->list_ids('AL359765.6', 1, 13780, $cln_cs);
ok(@seq_ids);
$str = join("\n", "----------", @seq_ids);
debug("$str\n");

@seq_ids = 
  $chr_sctg_mapper->list_ids('20', 30_000_000, 31_000_000, $chr_cs);
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
