
## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..5\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}


use lib 't';
use EnsTestDB;
$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new( { 'schema_sql' => ['../sql/table.clean'] } );

print "ok 2\n";

$ens_test->do_sql_file("t/cleangenome.dump");

print "ok 3\n";	

$db = $ens_test->get_DBSQL_Obj;

$sfadp = $db->get_SimpleFeatureAdaptor();

print "ok 4\n";

my $analysis = $db->get_AnalysisAdaptor->fetch_by_dbID(1);

my $sf = Bio::EnsEMBL::SimpleFeature->new();

my $start = 10;
my $end = 20;
my $strand = 1;
my $display_text = "Display text";

$sf->start($start);
$sf->end($end);
$sf->strand($strand);
$sf->display_text($display_text);
$sf->analysis($analysis);

$sfadp->store(1,$sf);

print "ok 5\n";


my @sf = $sfadp->fetch_by_contig_id(1);

if( scalar(@sf) != 1 || $sf[0]->start != $start ||
	$sf[0]->end != $end || $sf[0]->strand != $strand ||
	$sf[0]->display_text ne $display_text ) {
	print "not ok 6\n";
} else {
	print "ok 6\n";
}





