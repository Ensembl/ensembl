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
my $t_id = 'ENST00000310998'; # or maybe ENST00000278995

my $ga = $dba->get_GeneAdaptor();
my $ta = $dba->get_TranscriptAdaptor();

{
  my $gene = $ga->fetch_by_stable_id($id);
  delete $gene->{source};
  $gene->{description} = undef; #empty value means don't emit the key
  my $expected = qq{chr20\t30274333\t30300924\tENSG00000131044\t1000\t+\n};

  assert_bed($gene, $expected, 'Gene with no source serialises to BED as expected. Source is ensembl');
}

# Test transcripts
{
  my $transcript = $ta->fetch_by_stable_id($t_id);
  my $expected = qq{chr20\t30274333\t30298904\tENST00000310998\t1000\t+\t30274333\t30298904\t0,0,0\t6\t92,112,186,69,74,82,\t0,10117,11263,21390,22172,24489,\tC20orf125\tcmpl\tcmpl\t0,2,0,0,0,2,\tprotein_coding\tENSG00000131044\tC20orf125\tprotein_coding\n};
  assert_bed($transcript, $expected, 'Transcript emits as genePred format');
}

{
  my $sa = $dba->get_SliceAdaptor();
  my $slice = $sa->fetch_by_region('chromosome', 12, 1, 10);
  my $feature = Bio::EnsEMBL::Feature->new(
    -SLICE => $slice,
    -START => 1,
    -END => 10,
    -STRAND => 1,
  );
  my $expected = qq{12\t0\t10\t\t1000\t+\n};

  assert_bed($feature, $expected, 'Default feature should seralise without attributes');
}
sub assert_bed {
  my ($feature, $expected, $msg) = @_;
  my $fh = IO::String->new();
  my $ser = Bio::EnsEMBL::Utils::IO::BEDSerializer->new($fh);
  $ser->print_feature($feature);
  eq_or_diff(${$fh->string_ref()}, $expected, $msg);
}

done_testing();
