## Bioperl Test Harness Script for Modules
##
# CVS Version
# $Id$


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..7\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use lib 't';
use EnsTestDB;
$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/diffdump.dump");
    
# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok 2\n";    

my $gene_obj=$db->gene_Obj;

$db->diffdump(1);
$db->diff_fh("t/diff.sql");


$gene = new Bio::EnsEMBL::Gene;
$gene->id('temp_gene');

$tr   = new Bio::EnsEMBL::Transcript;
$tr1   = new Bio::EnsEMBL::Transcript;

$tr->id('temp_transcript_1');
$tr1->id('temp_transcript_2');
$tr->version(1);
$tr1->version(1);

$ex1   = new Bio::EnsEMBL::Exon;
$ex2   = new Bio::EnsEMBL::Exon;
$ex3   = new Bio::EnsEMBL::Exon;
print "ok 3\n";   

$ex1->id("dummy_id_1");
$ex1->version(1);
$ex1->created(time);
$ex1->modified(time);
$ex1->contig_id("AC021078.00006");
$ex1->phase(0);
$ex1->start(8);
$ex1->end(13);
$ex1->strand(1);

$ex2->id("dummy_id_2");
$ex2->created(time);
$ex2->modified(time);
$ex2->version(1);
$ex2->contig_id("AC021078.00006");
$ex2->phase(0);
$ex2->start(18);
$ex2->end(23);
$ex2->strand(1);

$ex3->id("dummy_id_3");
$ex3->created(time);
$ex3->modified(time);
$ex3->version(1);
$ex3->contig_id("AC021078.00006");
$ex3->phase(0);
$ex3->start(26);
$ex3->end(28);
$ex3->strand(1);

$tr->add_Exon($ex1);
$tr->add_Exon($ex2);
$trans = Bio::EnsEMBL::Translation->new();
$trans->id('temp_translation_1');
$trans->version(1);
$trans->start_exon_id('dummy_id_1');
$trans->start(8);
$trans->end_exon_id('dummy_id_2');
$trans->end(23);
$tr->translation($trans);

$tr1->add_Exon($ex1);
$tr1->add_Exon($ex2);
$tr1->add_Exon($ex3);
$trans = Bio::EnsEMBL::Translation->new();

$trans->id('temp_translation_2');
$trans->version(1);
$trans->start_exon_id('dummy_id_1');
$trans->start(8);
$trans->end_exon_id('dummy_id_3');
$trans->end(28);
$tr1->translation($trans);

$gene->add_Transcript($tr);
$gene->add_Transcript($tr1);
print "ok 4\n";

$gene_obj->write($gene);
print "ok 5\n";

$gene_obj->delete($gene->id);
print "ok 6\n";

open (FILE,"<t/diff.sql");
my $ok=0;
while (<FILE>) {
    if (/delete/) {
	$ok++;
    }
    if (/insert/) {
	$ok++;
    }
}

#Check if there were 28 inserts/deletes ;)
if ($ok == 29) {
    print "ok 7\n";
}
else {
    print "not ok 7\n";
    print STDERR "got $ok lines\n";
}

$db = 0;

unlink("t/diff.sql");


