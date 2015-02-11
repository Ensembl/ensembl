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

use strict;

use Bio::EnsEMBL::Test::TestUtils;

use Test::More;
use Test::Warnings;

use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;

# get a core DBAdaptor
my $dba = $multi->get_DBAdaptor("core");

#
# Test get_ProteinFeatureAdaptor works
#
my $pfa = $dba->get_ProteinFeatureAdaptor();

ok($pfa && ref($pfa) && $pfa->isa('Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor'));

my $pfs = $pfa->fetch_all_by_translation_id(21724);

print_features($pfs);

ok(@$pfs == 15);

sub print_features {
  my $features = shift;
  foreach my $f (@$features) {
	if (defined($f)) {
	  debug($f->start . '-' . $f->end . ' -> ' . $f->hseqname . ':' . $f->hstart . '-' . $f->hend);
	}
  }
}

# test adding and retrieving a new feature
my $start  = 10;
my $end    = 100;
my $hstart = 1;
my $hend   = 90;
my $hstrand = 1;
my $hseqname = 'RF1231';
my $percent_id = 90.8;
my $p_value = '1.52';
my $score   = 50;
my $species = 'Homo_sapiens';
my $hspecies = 'Mus_musculus';
my $hdes = "Hit description";

my $idesc = 'interpro description';
my $interpro_ac = 'interpro accession';

my $analysis = Bio::EnsEMBL::Analysis->new(-LOGIC_NAME => 'test');
$multi->save('core', 'protein_feature', 'meta_coord');


my $f = Bio::EnsEMBL::ProteinFeature->new
  (-START       => $start,
   -END         => $end,
   -ANALYSIS    => $analysis,
   -HSTART      => $hstart,
   -HEND        => $hend,
   -HSEQNAME    => $hseqname,
   -PERCENT_ID  => $percent_id,
   -P_VALUE     => $p_value,
   -SCORE       => $score,
   -SPECIES     => $species,
   -HSPECIES    => $hspecies,
   -HDESCRIPTION=> $hdes,
   -IDESC       => $idesc,
   -INTERPRO_AC => $interpro_ac);
   
$pfa->store($f,21724);

my $pfs = $pfa->fetch_all_by_translation_id(21724);

ok(@$pfs == 16);

my @pfs = grep{$_->hdescription() eq $hdes} @$pfs;

ok(scalar @pfs > 0);

$multi->restore('core', 'protein_feature');

done_testing();
