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
BEGIN { $| = 1; print "1..19\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use lib 't';
use EnsTestDB;
$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/newgene.dump");
    
# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok 2\n";    # 1st test passes.

my $gene_obj=$db->gene_Obj;

#Checking get_all_Gene_id;
my @genes=$gene_obj->get_all_Gene_id;
if ((scalar @genes == 1) && ($genes[0] eq 'test_gene')) {
    print "ok 3\n";
}
else {
    print "not ok 3\n";
}

#Checking get method exception
eval{
my $gene = $gene_obj->get('blabla');
};
if ($@) {
    print "ok 4\n";
}
else {
    print "not ok 4\n";
} 

#Checking get method (correct use)
my $gene = $gene_obj->get('test_gene');
if ($gene->isa('Bio::EnsEMBL::Gene')) {
    print "ok 5\n";
}
else {
    print "not ok 5\n";
}

#Check gene methods
if ($gene->id eq 'test_gene') {
    print "ok 6\n";
}
else {
    print "not ok 6\n";
}

if ($gene->version == 1) {
    print "ok 7\n";
}
else {
    print "not ok 7\n";
}

if ($gene->created == 962638206) {
    print "ok 8\n";
}
else {
    print "not ok 8\n";
}  

if ($gene->modified == 962638206) {
    print "ok 9\n";
}
else {
    print "not ok 9\n";
}  
if ($gene->_stored == 962638206) {
    print "ok 10\n";
}
else {
    print "not ok 10\n";
} 

#Not checking clone neighbourhood because they are on their way out of the schema
my $dblink=Bio::Annotation::DBLink->new();
$dblink->database('EMBL_dummy');
$dblink->primary_id('ACC_test');
$dblink->optional_id('ACC_optional');
$dblink->comment('This is a fake dblink object');
$gene->add_DBLink($dblink);
print "ok 11\n";

$ok=0;;
foreach my $link ($gene->each_DBLink) {
    if ($link->database eq 'EMBL_dummy') {
	$ok++;
    }
    if ($link->primary_id eq 'ACC_test') {
	$ok++;
    }
    if ($link->optional_id eq 'ACC_optional') {
	$ok++;
    }
    if ($dblink->comment eq 'This is a fake dblink object') {
	$ok++;
    }
}
if ($ok == 4) {
    print "ok 12\n";
}
else {
    print "not ok 12\n";
}

if ($gene->is_known == 1) {
    print "ok 13\n";
}
else {
    print "not ok 13\n";
}
$ok=0;
foreach my $contig_id ($gene->unique_contig_ids) {
    if ($contig_id =~ /test_contig_1|test_contig_2/){
	$ok++;
    }
}
if ($ok == 2) {
    print "ok 14\n";
}
else {
    print "not ok 14\n";
}

$ok=0;
foreach my $exon ($gene->each_unique_Exon) {
    if ($exon->id =~ /test_exon_1|test_exon_2|test_exon_3/) {
	$ok++;
    }
}
if ($ok == 3) {
    print "ok 15\n";
}
else {
    print "not ok 15\n";
}

$ok=0;
foreach my $exon ($gene->all_Exon_objects) {
    if ($exon->id =~ /test_exon_1|test_exon_2|test_exon_3/) {
	$ok++;
    }
}
if ($ok == 4) {
    print "ok 16\n";
}
else {
    print "not ok 16 $ok\n";
}

#Not testing add_transcript yet
$ok=0;
foreach my $trans ($gene->each_Transcript) {
    if ($trans->id =~ /test_transcript_1|test_transcript_2/) {
	$ok++;
    }
}
if ($ok == 2) {
    print "ok 17\n";
}
else {
    print "not ok 17\n";
}

my $gene=$gene_obj->get_Gene_by_Transcript_id('test_transcript_1');
if ($gene->id eq 'test_gene') {
    print "ok 18\n";
}
else {
    print "not ok 18\n";
}

my $gene=$gene_obj->get_Gene_by_Transcript_id('bla_bla');

if ($gene == undef) {
    print "ok 19\n";
}
else {
    print "not ok 19\n";
}



#Methods to test:
#$gene_obj->delete;
#$gene_obj->delete_Exon;
#$gene_obj->delete_Supporting_Evidence;

#$gene_obj->get_geneids_by_hids;
#$gene_obj->get_array_supporting;
#$gene_obj->get_Gene_by_Transcript_id;
#$gene_obj->get_Gene_by_DBLink; 
#$gene_obj->get_Exon;
#$gene_obj->get_supporting_evidence;
#$gene_obj->get_Transcript;
#$gene_obj->get_Translation;
#$gene_obj->get_Virtual_Contig;
#$gene_obj->get_Transcript_in_VC_coordinates;
#$gene_obj->write;
#$gene_obj->write_Exon;
#$gene_obj->write_supporting_evidence;
#$gene_obj->write_Transcript;
#$gene_obj->write_Translation;
#$gene_obj->get_NewId;
#$gene_obj->get_new_GeneID;
#$gene_obj->get_new_TranscriptID;
#$gene_obj->get_new_ExonID;
#$gene_obj->get_new_TranslationID;

