
## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..4\n";
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}


use Bio::EnsEMBL::EMBL_Dump;
use Bio::SeqIO;

use lib 't';
use EnsTestDB;
$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/staticgoldenpath.dump");

my $db = $ens_test->get_DBSQL_Obj;

my $cloneadaptor = $db->get_CloneAdaptor();

my $clone = $cloneadaptor->fetch_by_accession('clone1');

print "ok 2\n";

my $vc = $clone->virtualcontig();

if( !defined $vc ) {
    print "not ok 3\n";
} else {
  print "ok 3\n";
}

&Bio::EnsEMBL::EMBL_Dump::add_ensembl_comments($vc); 
    
$seqout = Bio::SeqIO->new(-format => 'embl', -file => ">t/emblclone.out" );

# adds standard dumping information to the aseqstream to drive off the
# the annseq objects
&Bio::EnsEMBL::EMBL_Dump::ensembl_annseq_output($seqout);


$seqout->write_seq($vc);

print "ok 4\n";


