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
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::MiscFeature;
use Bio::EnsEMBL::MiscSet;
use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0; #set to 1 to turn on debug printouts


#test constructor
my $mf = Bio::EnsEMBL::MiscFeature->new(-START => 10,
                                        -END   => 100);

ok($mf->start() == 10 && $mf->end() == 100);



#
# Test add_set, get_set, get_set_codes
#
my $ms1 = Bio::EnsEMBL::MiscSet->new(-CODE => '1mbcloneset',
                                     -NAME => '1mb Clone Set',
                                     -DESCRIPTION => 'This is a 1MB cloneset',
                                     -LONGEST_FEATURE => 1e7);

my $ms2 = Bio::EnsEMBL::MiscSet->new(-CODE => 'tilepath',
                                     -NAME => 'Tiling Path',
                                     -DESCRIPTION => 'NCBI33 Tiling Path',
                                     -LONGEST_FEATURE => 1e6);



$mf->add_MiscSet($ms1);
$mf->add_MiscSet($ms2);


my $ms3 = $mf->get_all_MiscSets($ms1->code)->[0];
my $ms4 = $mf->get_all_MiscSets($ms2->code)->[0];

ok( $ms3 == $ms1);
ok( $ms4 == $ms2);


#
# Test add_attribute, get_attribute_types, get_attribute
#

my $name1 = Bio::EnsEMBL::Attribute->new
  ( -CODE => 'name',
    -VALUE => 'test name'
  );

my $name2 = Bio::EnsEMBL::Attribute->new
  ( -CODE => 'name',
    -VALUE => 'AL4231124.1'
  );


$mf->add_Attribute( $name1 );

ok($mf->display_id eq "test name");

$mf->add_Attribute( $name2 );

my $vers1 = Bio::EnsEMBL::Attribute->new
  ( -CODE => 'version',
    -VALUE => 4
  );

$mf->add_Attribute( $vers1 );


my @attribs = @{$mf->get_all_Attributes('name')};

ok(@attribs == 2);

@attribs = grep {$_ eq $name1 || $_ eq $name2} @attribs;

ok(@attribs == 2);

@attribs = @{$mf->get_all_Attributes('version')};
ok(@attribs == 1 && $attribs[0]->value() eq '4');

@attribs = @{$mf->get_all_Attributes()};
ok( @attribs == 3 );

done_testing();
