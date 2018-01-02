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

my $loaded = 0;
END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0;

$loaded = 1;

ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
ok( $multi );
$multi->hide( "patch", "analysis", "analysis_description", "gene" );

my $db = $multi->get_DBAdaptor( "patch" );
ok($db);


my $analysis_ad = $db->get_AnalysisAdaptor();

ok($analysis_ad);



my $analysis = Bio::EnsEMBL::Analysis->new();

$analysis->logic_name('dummy_analysis');
$analysis->db('dummy');
$analysis->program('dummy');
$analysis->gff_source('dummy');
$analysis->gff_feature('dummy');
$analysis->description( "some funny description" );
$analysis->display_label( "and a label" );
$analysis->displayable( 1 );
$analysis->created( "2005-10-28 10:28:29");
$analysis->web_data({thing => 'blah'});

ok($analysis);

$analysis_ad->store($analysis);

ok(defined $analysis->dbID() );


my $analysis_out = $analysis_ad->fetch_by_logic_name('dummy_analysis');


ok($analysis_out);

is($analysis_out->db, 'dummy', "Db matches");
is($analysis_out->web_data->{thing}, 'blah', "Web data correctly extracted from DB via dirty eval");

ok( check_methods( $analysis_out, "db", "db_file", "dbID", "compare",
		   "logic_name", "parameters", "gff_source", "gff_feature",
		   "module", "module_version", "program_file",
		   "program", "db_version", "adaptor", "display_label", 
		   "displayable", "web_data" ));

is( $analysis_out->description, "some funny description", "Funny description matches" );

# try updating existing description
$analysis->logic_name("new_dummy");
$analysis->description("new description");
$analysis->display_label("new label");
$analysis->displayable(0);
$analysis->web_data("blahblah");
my $dbID = $analysis->dbID();
$analysis_ad->update($analysis);
my $analysis_updated = $analysis_ad->fetch_by_dbID($dbID);
is($analysis_updated->logic_name(), "new_dummy", "Logic name is correct");
is($analysis_updated->description(), "new description", "Description is correct");
is($analysis_updated->display_label(), "new label", "Label is correct");
is($analysis_updated->displayable(), 0, "Displayable is correct");
is($analysis_updated->web_data(), "blahblah", "Web data is correct"); # Note this does not test the eval output, just the set value.

# now try updating analysis that has no existing description
$analysis = Bio::EnsEMBL::Analysis->new();
$analysis->logic_name('dummy_analysis');
$analysis->created( "2005-10-28 10:28:29");
$analysis->displayable(1);
$analysis_ad->store($analysis);
$dbID = $analysis->dbID();
$analysis->description("updated description");
$analysis_ad->update($analysis);
$analysis_updated = $analysis_ad->fetch_by_dbID($dbID);
is($analysis_updated->description(), "updated description", "Description matches the update");
is_rows(2, $db, "analysis_description");

$multi->restore();

# Fetch all analyses for a feature class
my @analyses = @{ $analysis_ad->fetch_all_by_feature_class('Gene') };
is(scalar(@analyses), 6, "Found 6 analyses for genes");

@analyses = @{ $analysis_ad->fetch_all_by_feature_class('DensityFeature') };
is(scalar(@analyses), 6, "Found 6 analyses for density features");



sub check_methods { 
  my $obj = shift;

  my $all_implemented = 1;
  while( my $method = shift ) {
    if( ! $obj->can( $method )) {
      $all_implemented = 0;
      debug( "Analysis doesnt implement $method" );
    }
  }
  return $all_implemented;
}

done_testing();
