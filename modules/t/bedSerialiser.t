# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

## no critic (RequireFilenameMatchesPackage)
use strict;
use warnings;

use Test::More;
use Test::Warnings;
use Test::Differences;
use IO::String;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Utils::IO::BEDSerializer;
use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Slice;

my $db = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $db->get_DBAdaptor('core');

my $id = 'ENSG00000131044';

my $ga = $dba->get_GeneAdaptor();

{
  my $gene = $ga->fetch_by_stable_id($id);
  delete $gene->{source};
  $gene->{description} = undef; #empty value means don't emit the key

  my $expected = qq{chr20\t30274333\t30300924\tENSG00000131044\t0\t+\n};

  assert_bed($gene, $expected, 'Gene with no source serialises to BED as expected. Source is ensembl');
}

{
  my $cs = $dba->get_CoordSystemAdaptor()->fetch_by_name('chromosome');
  my $feature = Bio::EnsEMBL::Feature->new(
    -SLICE => Bio::EnsEMBL::Slice->new(
      -COORD_SYSTEM => $cs,
      -SEQ => ('A'x10),
      -SEQ_REGION_NAME => 'wibble',
      -START => 1,
      -END => 10
    ),
    -START => 1,
    -END => 10,
    -STRAND => 1,
  );
  my $expected = qq{wibble\t0\t10\t\t0\t+\n};

  assert_bed($feature, $expected, 'Default feature should seralise without attributes but leave a trailing \t');
}
sub assert_bed {
  my ($feature, $expected, $msg) = @_;
  my $fh = IO::String->new();
  my $ser = Bio::EnsEMBL::Utils::IO::BEDSerializer->new($fh);
  $ser->print_feature($feature);
  eq_or_diff(${$fh->string_ref()}, $expected, $msg);
}

done_testing();
