
BEGIN { $| = 1;  
	use Test;
	plan tests => 9;
}

my $loaded = 0;
END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 0;

$loaded = 1;

ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
ok( $multi );
$multi->hide( "core", "analysis" );

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


ok($analysis);

$analysis_ad->store($analysis);


ok(defined $analysis->dbID() );


my $analysis_out = $analysis_ad->fetch_by_logic_name('dummy_analysis');


ok($analysis_out);

ok($analysis_out->db eq 'dummy');

ok( check_methods( $analysis_out, "db", "db_file", "dbID", "compare",
		   "logic_name", "parameters", "gff_source", "gff_feature",
		   "module", "module_version", "program_file",
		   "program", "db_version", "adaptor" ));
		   

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


sub debug {
  my $txt = shift;
  if( $verbose ) {
    print STDERR $txt,"\n";
  }
}
