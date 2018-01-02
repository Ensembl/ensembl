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
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils qw/is_rows/;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $multi->get_DBAdaptor( 'core' );

note 'ISSUE: ENSCORESW-445';
note 'DESC:  Attributes can be duplicated in the schema after adding them to a transformed transcript. Storing the transcript has no affect';

# "Globals"
my $transcript_stable_id = 'ENST00000246203';

# Get adaptors
my $tra = $db->get_TranscriptAdaptor();

# Now going after the transcript and store some new attributes
my $transcript = $tra->fetch_by_stable_id($transcript_stable_id);
$transcript->load();
is_deeply($transcript->get_all_Attributes(), [], 'Transcript has no attributes');
my $basic_attribute = Bio::EnsEMBL::Attribute->new(-CODE => 'note', -VALUE => 'This is a note', -NAME => 'Note', -DESCRIPTION => q{});
$transcript->add_Attributes($basic_attribute);
is_deeply($transcript->get_all_Attributes(), [$basic_attribute], 'Transcript has one attribute') or diag explain $transcript->get_all_Attributes();

# And then we loop through twice creating two new transcripts but we expect the total number of attributes to hit 5
foreach my $i (1..2) {
  note '$transcript->transform("contig") number '.$i;
  #We have an attribute attached; transform() and store the new transcript with an additional attribute
  my $transformed_transcript = $transcript->transform('contig');
  ok(defined $transformed_transcript, 'Transcript projected');
  my $cloned_attr = bless({%{$basic_attribute}}, ref($basic_attribute));
  $cloned_attr->value('Look it is a different value');
  $transformed_transcript->add_Attributes($cloned_attr);
  is_deeply(
    $transformed_transcript->get_all_Attributes(), 
    [$basic_attribute, $cloned_attr], 
    'We should have the old and the cloned new attribute available'
  ) or diag explain $transformed_transcript->get_all_Attributes;
}

done_testing();
