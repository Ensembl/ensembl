use strict;
use warnings;

BEGIN {
	$| = 1;
	use Test;
	plan tests => 12;
}
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Operon;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::DBSQL::OperonAdaptor;
use Bio::EnsEMBL::DBSQL::OperonTranscriptAdaptor;
debug("Startup test");
ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $dba = $multi->get_DBAdaptor("core");

debug("Test database instatiated");
ok($dba);

# get a slice
my $slice = $dba->get_SliceAdaptor()->fetch_by_seq_region_id(469283);

# create operon
my $start         = 31225346;
my $end           = 31225946;
my $strand        = 1;
my $display_label = "accBC";
my $analysis = $dba->get_AnalysisAdaptor->fetch_by_logic_name("Genscan");
my $operon = Bio::EnsEMBL::Operon->new(
	-START         => $start,
	-END           => $end,
	-STRAND        => $strand,
	-SLICE         => $slice,
	-DISPLAY_LABEL => $display_label,
	-ANALYSIS => $analysis);
$operon->add_DBEntry( Bio::EnsEMBL::DBEntry->new( -DBNAME     => 'EMBL',
												  -RELEASE    => 1,
												  -PRIMARY_ID => 'XZY',
												  -DISPLAY_ID => '123',
												  -ANALYSIS   => $analysis ) );

# check it has the correct properties
ok( $display_label, $operon->display_label(),     "Operon name" );
ok( $start,         $operon->seq_region_start(),  "Operon start" );
ok( $end,           $operon->seq_region_end(),    "Operon end" );
ok( $strand,        $operon->seq_region_strand(), "Operon strand" );

my $operon_adaptor = Bio::EnsEMBL::DBSQL::OperonAdaptor->new($dba);

# store operon
$operon_adaptor->store($operon);
ok( defined $operon->dbID() );
# retrieve operon
my $operon2 = $operon_adaptor->fetch_by_dbID( $operon->dbID() );
ok( $operon2->dbID(),             $operon->dbID(),             "Operon ID" );
ok( $operon2->display_label(),    $operon->display_label(),    "Operon name" );
ok( $operon2->seq_region_start(), $operon->seq_region_start(), "Operon start" );
ok( $operon2->seq_region_end(),   $operon->seq_region_end(),   "Operon end" );
ok( $operon2->seq_region_strand(),
	$operon->seq_region_strand(),
	"Operon strand" );

