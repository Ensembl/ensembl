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


my $asma = $db->get_AssemblyMapperAdaptor();

#
# Test fetch_by_CoordSystems
#

my $csa = $db->get_CoordSystemAdaptor();
my $toplevel_cs = $csa->fetch_by_name('toplevel');
my $cln_cs    = $csa->fetch_by_name('clone');
my $superctg_cs = $csa->fetch_by_name('supercontig');

my $cln_toplevel_mapper =
  $asma->fetch_by_CoordSystems($toplevel_cs, $cln_cs);
my $superctg_toplevel_mapper =
  $asma->fetch_by_CoordSystems($toplevel_cs, $superctg_cs);

ok($cln_toplevel_mapper && 
   $cln_toplevel_mapper->isa('Bio::EnsEMBL::TopLevelAssemblyMapper'));


#
# test db has chr 20  (50KB -> 62MB)
#

#
# Test map
#

debug("MAP 'AL359765.6'->toplevel");
my @coords = $cln_toplevel_mapper->map('AL359765.6', 1, 13780, 1, $cln_cs);
print_coords(@coords);
ok(@coords);
# [ENSCORESW-844]. Test mapped coordinate names
is($coords[1]->name(), "20");

debug("MAP NT_028392->toplevel");
@coords = $superctg_toplevel_mapper->map('NT_028392', 600_000, 1e6, 1, 
                                         $superctg_cs);
print_coords(@coords);
ok(@coords);



#
# Test list_seq_regions
#

my @seq_regions =
  $cln_toplevel_mapper->list_seq_regions('AL359765.6', 1, 13780, $cln_cs);
my $str = join("\n", "----------", @seq_regions);
debug("$str\n");
ok(@seq_regions == 1 && $seq_regions[0] eq '20');


@seq_regions = $superctg_toplevel_mapper->list_seq_regions('NT_028392',
                                                           6e5, 1e6,
                                                           $superctg_cs);

$str = join("\n", "----------", @seq_regions);
debug("$str\n");
is(@seq_regions, 1, "Got 1 seq region");
is($seq_regions[0], 20, "Seq region is 20");


#
# Test list_seq_ids
#

my @ids =
  $cln_toplevel_mapper->list_ids('AL359765.6', 1, 13780, $cln_cs);
$str = join("\n", "----------", @ids);
debug("$str\n");
ok(@ids == 1 && $ids[0] == 469283);

@ids = $superctg_toplevel_mapper->list_ids('NT_028392',
                                                   6e5, 1e6,
                                                   $superctg_cs);
$str = join("\n", "----------", @ids);
debug("$str\n");
ok(@ids == 1 && $ids[0] == 469283);

sub print_coords {
  my @coord_list = @_;

  return if(!$verbose);

  foreach my $coord (@coord_list) {
    if($coord->isa('Bio::EnsEMBL::Mapper::Gap')) {
      debug("GAP");
      next;
    }
    my $cs = $coord->coord_system();
    my $cs_str = ($cs) ? $cs->name() . $cs->version() : '';
    debug($coord->id()."\t". $coord->start()."-".$coord->end().
          " (".$coord->strand.") [$cs_str]");
  }
}

done_testing();
