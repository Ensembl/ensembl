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
use warnings;

use Test::More;
use Test::Warnings;

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor;

# Get a DBAdaptor to from the test system
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;
ok($multi);
my $db = $multi->get_DBAdaptor("core");
ok($db);

# Should get meaningful type back
my $test_adaptor;

$test_adaptor = $db->get_ArchiveStableIdAdaptor();
if (defined $test_adaptor) {
  ok($test_adaptor->isa("Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor"));
}
$test_adaptor = $db->get_ProteinFeatureAdaptor();
ok($test_adaptor->isa("Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor"));
$test_adaptor = $db->get_PredictionTranscriptAdaptor();
ok($test_adaptor->isa("Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor"));
$test_adaptor = $db->get_SequenceAdaptor();
ok($test_adaptor->isa("Bio::EnsEMBL::DBSQL::SequenceAdaptor"));
$test_adaptor = $db->get_GeneAdaptor();
ok($test_adaptor->isa("Bio::EnsEMBL::DBSQL::GeneAdaptor"));
$test_adaptor = $db->get_ExonAdaptor();
ok($test_adaptor->isa("Bio::EnsEMBL::DBSQL::ExonAdaptor"));
$test_adaptor = $db->get_TranscriptAdaptor();
ok($test_adaptor->isa("Bio::EnsEMBL::DBSQL::TranscriptAdaptor"));
$test_adaptor = $db->get_TranslationAdaptor();
ok($test_adaptor->isa("Bio::EnsEMBL::DBSQL::TranslationAdaptor"));
$test_adaptor = $db->get_SliceAdaptor();
ok($test_adaptor->isa("Bio::EnsEMBL::DBSQL::SliceAdaptor"));
$test_adaptor = $db->get_AnalysisAdaptor();
ok($test_adaptor->isa("Bio::EnsEMBL::DBSQL::AnalysisAdaptor"));
$test_adaptor = $db->get_SimpleFeatureAdaptor();
ok($test_adaptor->isa("Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor"));
$test_adaptor = $db->get_RepeatConsensusAdaptor();
ok($test_adaptor->isa("Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor"));
$test_adaptor = $db->get_RepeatFeatureAdaptor();
ok($test_adaptor->isa("Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor"));
$test_adaptor = $db->get_ProteinAlignFeatureAdaptor();
ok($test_adaptor->isa("Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor"));
$test_adaptor = $db->get_DnaAlignFeatureAdaptor();
ok($test_adaptor->isa("Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor"));
$test_adaptor = $db->get_AssemblyMapperAdaptor();
ok($test_adaptor->isa("Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor"));
$test_adaptor = $db->get_DBEntryAdaptor();
ok($test_adaptor->isa("Bio::EnsEMBL::DBSQL::DBEntryAdaptor"));
$test_adaptor = $db->get_KaryotypeBandAdaptor();
ok($test_adaptor->isa("Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor"));
$test_adaptor = $db->get_SupportingFeatureAdaptor();
ok($test_adaptor->isa("Bio::EnsEMBL::DBSQL::SupportingFeatureAdaptor"));
$test_adaptor = $db->get_MarkerFeatureAdaptor();
ok($test_adaptor->isa("Bio::EnsEMBL::Map::DBSQL::MarkerFeatureAdaptor"));
$test_adaptor = $db->get_MarkerAdaptor();
ok($test_adaptor->isa("Bio::EnsEMBL::Map::DBSQL::MarkerAdaptor"));

#ok($test_adaptor->isa("Bio::EnsEMBL::DBSQL::ProxySNPAdaptor"));

# Should get an error if we ask for something non-existent
#eval { $db->get_adaptor("SomeNonExistentType") };
#ok($@);

# Check setting module with good values
ok($db->set_adaptor("ArchiveStableId", "Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor" ));

# Setting an unknown data type should NO LONGER give an error
my $ret = $db->set_adaptor("SomeNonExistentType", "Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor");
ok(defined($ret));


# Setting to a non-subclass of the default should NOT give an error
$ret =  $db->set_adaptor("ArchiveStableId", "Bio::EnsEMBL::DBSQL::GeneAdaptor");
ok(defined($ret));

# Generic adaptors

# Should work OK with subclasses of BaseFeatureAdaptor
my $rfa = Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor->new($db);
my $sfa = Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor->new($db);
ok($db->add_GenericFeatureAdaptor("Repeat", $rfa));
ok($db->add_GenericFeatureAdaptor("Simple", $sfa));


# Check get-ing the above
# by name ...
my %generic_adaptors = %{$db->get_GenericFeatureAdaptors("Simple", "Repeat")};
ok(%generic_adaptors);

# no arg should return all
%generic_adaptors = %{$db->get_GenericFeatureAdaptors()};
my $size = keys(%generic_adaptors);
ok($size == 2);

# requesting one that hasn't been added should throw
eval { my %generic_adpators = $db->get_GenericFeatureAdaptors("Mickey") };
ok($@);

# Slice tests - should these go in slice.t?
my $slice = $db->get_SliceAdaptor()->fetch_by_region('chromosome','X',1,10000);
ok($slice);

my %features = %{$slice->get_generic_features()};
ok(!%features);

done_testing();
