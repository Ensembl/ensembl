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
use Test::Exception;
use Test::Warnings;
use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
ok(1);


#Fetch Transcript by stable_id and ensure it has a slice
my $db = $multi->get_DBAdaptor('core');
my $ta = $db->get_TranscriptAdaptor();
my $tr = $ta->fetch_by_stable_id( "ENST00000222222" );
ok($tr, "Fetched the Transcript by stable_id");
ok($tr->slice(), "Transcript has slice");


#Get TranscriptSupportingFeatureAdaptor
my $tsf_adaptor = $db->get_TranscriptSupportingFeatureAdaptor;


#Fetch all transcript supporting features
my $supporting_features = $tsf_adaptor->fetch_all_by_Transcript($tr);
isa_ok( $supporting_features, 'ARRAY' );
ok(5 == scalar(@$supporting_features), "Fetched all transcript supporting features");


#Fetch all transcript supporting features of type 'protein_align_feature'
my $supporting_features_proteins = $tsf_adaptor->fetch_all_by_Transcript($tr, "protein_align_feature");
isa_ok( $supporting_features_proteins, 'ARRAY' );
ok(3 == scalar(@$supporting_features_proteins), "Fetched all transcript supporting features (protein_align_features) ");


foreach my $supporting_feature(@$supporting_features_proteins){
	ok('Bio::EnsEMBL::DnaPepAlignFeature' eq ref $supporting_feature, "Got back the right ref type: Bio::EnsEMBL::DnaPepAlignFeature");
}

#Fetch all transcript supporting features of type 'dna_align_feature'
my $supporting_features_dnas = $tsf_adaptor->fetch_all_by_Transcript($tr, "dna_align_feature");
isa_ok( $supporting_features_dnas, 'ARRAY' );
ok(2 == scalar(@$supporting_features_dnas), "Fetched all transcript supporting features (dna_align_features) ");


foreach my $supporting_feature(@$supporting_features_dnas){
	ok('Bio::EnsEMBL::DnaDnaAlignFeature' eq ref $supporting_feature, "Got back the right ref type: Bio::EnsEMBL::DnaDnaAlignFeature");
}

#Test exceptions with unknown feature_types
dies_ok { $tsf_adaptor->fetch_all_by_Transcript($tr, "unknown_feature_type") } 'fetch_all_by_Transcript dies with unknown feature_type';
throws_ok {$tsf_adaptor->fetch_all_by_Transcript($tr, "unknown_feature_type") } qr/feature type must be dna_align_feature or protein_align_feature/, 
           'feature type must be dna_align_feature or protein_align_feature';

#Test exception with empty features array
throws_ok { Bio::EnsEMBL::DnaPepAlignFeature->new(-features => [], -align_type => 'ensembl') } qr/features array must not be empty/, 
           'features array must not be empty';

done_testing();
