
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

$simplef = $db->get_SimpleFeatureAdaptor;

print "ok 4\n";

