
use strict;
use warnings;
use Test::More;

my $loaded = 0;
END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

our $verbose = 0;

$loaded = 1;

ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
ok( $multi );
$multi->hide( "core", "analysis", "analysis_description" );

my $db = $multi->get_DBAdaptor( "core" );
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
$analysis->web_data("blah");

ok($analysis);

$analysis_ad->store($analysis);

ok(defined $analysis->dbID() );


my $analysis_out = $analysis_ad->fetch_by_logic_name('dummy_analysis');


ok($analysis_out);

ok($analysis_out->db eq 'dummy');

ok( check_methods( $analysis_out, "db", "db_file", "dbID", "compare",
		   "logic_name", "parameters", "gff_source", "gff_feature",
		   "module", "module_version", "program_file",
		   "program", "db_version", "adaptor", "display_label", 
		   "displayable", "web_data" ));

ok( $analysis_out->description eq "some funny description" );

# try updating existing description
$analysis->logic_name("new_dummy");
$analysis->description("new description");
$analysis->display_label("new label");
$analysis->displayable(0);
$analysis->web_data("blahblah");
my $dbID = $analysis->dbID();
$analysis_ad->update($analysis);
my $analysis_updated = $analysis_ad->fetch_by_dbID($dbID);
ok($analysis_updated->logic_name() eq "new_dummy");
ok($analysis_updated->description() eq "new description");
ok($analysis_updated->display_label() eq "new label");
ok($analysis_updated->displayable() eq 0);
ok($analysis_updated->web_data() eq "blahblah");

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
ok($analysis_updated->description() eq "updated description");
ok( count_rows( $db, "analysis_description" ) == 2 );




$multi->restore();


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
