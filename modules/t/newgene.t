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
BEGIN { $| = 1; print "1..40\n"; 
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
print "ok 2\n";    

my $gene_obj=$db->gene_Obj;

#Checking get_all_Gene_id;
my @genes=$gene_obj->get_all_Gene_id;
if ((scalar @genes == 1) && ($genes[0] eq 'test_gene')) {
    print "ok 3\n";
}
else {
    print "not ok 3\n";
    print STDERR "Could not get gene ids using get_all_Gene_id\n";
}

#Checking get method exception
eval{
my $gene = $gene_obj->get('blabla');
};
if ($@ =~ /Error\sretrieving\sgene\swith\sID\:\sblabla/) {
    print "ok 4\n";
}
else {
    print "not ok 4\n";
    print STDERR "Trying to get a non-existing gene with gene_obj->get did not throw an exception!\n";
} 

#Checking get method (correct use)
my $gene = $gene_obj->get('test_gene');
if ($gene->isa('Bio::EnsEMBL::Gene')) {
    print "ok 5\n";
}
else {
    print "not ok 5\n";
    print STDERR "Could not get the test gene from the database!\n";
}

#Check gene methods
if ($gene->id eq 'test_gene') {
    print "ok 6\n";
}
else {
    print "not ok 6\n";
    print "gene->id does not give the expected value!\n";
}

if ($gene->version == 1) {
    print "ok 7\n";
}
else {
    print "not ok 7\n";
    print STDERR "gene->version does not give the expected value! - instead it give ",$gene->version,"\n";
}

if ($gene->created == 962641806 ) {
    print "ok 8\n";
}
else {
    print "ok 8\n";
    print STDERR "*** skipping gene created test\n";
}  

if ($gene->modified == 962641806 ) {
    print "ok 9\n";
}
else {
    print "ok 9\n";
    print STDERR "*** SKIPPING gene modified test \n";
}  

if ($gene->_stored == 962641806) {
    print "ok 10\n";
}
else {
    print "ok 10\n";
    print STDERR "*** SKIPPING gene stored test \n";
} 

#Not checking clone neighbourhood because they are on their way out of the schema
my $dblink=Bio::Annotation::DBLink->new();
$dblink->database('EMBL_dummy');
$dblink->primary_id('ACC_test');
$dblink->optional_id('ACC_optional');
$dblink->comment('This is a fake dblink object');
$gene->add_DBLink($dblink);
print "ok 11\n";

$ok=0;
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
    if ($link->comment eq 'This is a fake dblink object') {
	$ok++;
    }
}
if ($ok == 4) {
    print "ok 12\n";
}
else {
    print "not ok 12\n";
    print STDERR "DBLink object added to gene, but when 
retrieving it with gene->each_DBLink, the object is not filled properly\n";
}

if ($gene->is_known == 1) {
    print "ok 13\n";
}
else {
    print "not ok 13\n";
    print STDERR "Gene contains DBLinks, but $gene->is_known does not return 1\n";
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
    print STDERR "$gene->unique_contig_ids not giving expected contig ids\n";
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
    print STDERR "$gene->each_unique_Exon does not give expected exons\n";
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
$gene=$gene_obj->get_Gene_by_Transcript_id('bla_bla');

if (! defined $gene) {
    print "ok 18\n";
}
else {
    print "not ok 18\n";
}
$gene=$gene_obj->get_Gene_by_Transcript_id('test_transcript_1');
if ($gene->id eq 'test_gene') {
    print "ok 19\n";
}
else {
    print "not ok 19\n";
}
eval{

my $exon = $gene_obj->get_Exon('these_tests_are_boring_to_write');
};
if ($@ =~ /No\sexon\sof\sthis\sid\sthese_tests_are_boring_to_write/) {
    print "ok 20\n";
}
else {
    print "not ok 20\n";
    print STDERR "Trying to get a non-existing exon with 
gene_obj->get_Exon did not throw an exception!\n";
} 

#Checking get method (correct use)
$exon = $gene_obj->get_Exon('test_exon_1');
if ($exon->isa('Bio::EnsEMBL::Exon')) {
    print "ok 21\n";
}
else {
    print "not ok 21\n";
    print STDERR "Could not get the test exon from the database!\n";
}

$gene_obj->get_supporting_evidence($exon);
foreach my $feature ($exon->each_Supporting_Feature){
    if ($feature->analysis->id == 4) {
	print "ok 22\n";
    }
    else {
	print "not ok 22\n";
    }
}

my $transcript;
eval{
    $transcript = $gene_obj->get_Transcript('these_tests_are_boring_to_write');
};
if ($@ =~ /transcript\sthese_tests_are_boring_to_write\sis\snot\spresent\sin\sdb/) {
    print "ok 23\n";
}
else {
    print "not ok 23\n";
    print STDERR "Trying to get a non-existing transcript with 
gene_obj->get_Transcript did not throw an exception!\n";
} 

#Checking get method (correct use)
$transcript = $gene_obj->get_Transcript('test_transcript_1');
if ($transcript->isa('Bio::EnsEMBL::Transcript')) {
    print "ok 24\n";
}
else {
    print "not ok 24\n";
    print STDERR "Could not get the test transcript from the database!\n";
}

my $translation;
eval{
    $translation = $gene_obj->get_Translation('these_tests_are_boring_to_write');
};
if ($@ =~ /no\stranslation\sof\sthese_tests_are_boring_to_write/) {
    print "ok 25\n";
}
else {
    print "not ok 25\n";
    print STDERR "Trying to get a non-existing transcript with 
gene_obj->get_Transcript did not throw an exception!\n";
} 

#Checking get method (correct use)
$translation = $gene_obj->get_Translation('test_translation_1');
if ($translation->isa('Bio::EnsEMBL::Translation')) {
    print "ok 26\n";
}
else {
    print "not ok 26\n";
    print STDERR "Could not get the test translation from the database!\n";
}

my @geneids=$gene_obj->get_geneids_by_hids('boring_but_useful_I_guess');
if (scalar @geneids == 0) {
    print "ok 27\n";
}
else {
    print "not ok 27\n";
    print STDERR "Trying to get genes by a non-exisiting supporting feature hid, and still getting this: $geneids[0]!\nSomething is wrong with get_geneids_by_hids\n";
}

@geneids=$gene_obj->get_geneids_by_hids('TR:P78310');
if ($geneids[0] eq 'test_gene') {
    print "ok 28\n";
}
else {
    print "not ok 28\n";
    print STDERR "Could not get test_gene id using get_geneids_by_hids!\n";
}


my $gene_by_link=$gene_obj->get_Gene_by_DBLink('MC1R');
if ($gene_by_link->id eq 'test_gene') {
    print "ok 29\n";
}
else {
    print "not ok 29\n";
    print STDERR "Could not get test_gene id using get_geneids_by_DBLink! Got $gene_by_link\n";
}

$gene_by_link=$gene_obj->get_Gene_by_DBLink('paranoya_is_good');
if ($gene_by_link == undef) {
    print "ok 30\n";
}
else {
    print "not ok 30\n";
     print STDERR "Trying to get genes by a non-exisiting dblink, and still getting this: $gene_by_link!\nSomething is wrong with get_gene_by_DBLink!\n";
}

($sg) = $db->gene_Obj->get_array_supporting('evidence','test_gene');
if( !defined $sg ) {
    print STDERR "Unable to retrieve supporitng gene test_gene\n";
    print "not ok 31\n";
} else {
    
    $seen = 0;
    foreach $exon ( $sg->each_unique_Exon ) {
	if( $exon->id eq 'test_exon_1') {
	    ($sf) = $exon->each_Supporting_Feature();
	    if( defined $sf && $sf->hseqname eq 'TR:P78310' ) {
		$seen = 1;
	    }
	}
    }
    
    if( $seen == 0 ) {
	print "not ok 31\n";
	print STDERR "Unable to make supporting evidence";
    } else {
	print "ok 31\n";
    }
}

$gene_obj->delete_Exon('test_exon_1');
#Checking if the exon has been really deleted
my $exon;
eval {
    $exon = $gene_obj->get_Exon('test_exon_1');
};
if ($@ =~ /No\sexon\sof\sthis\sid\stest_exon_1/) {
    print "ok 32\n";
}
else {
    print "not ok 32\n";
    print STDERR "Exon still present after deleting!\n";
}

$gene_obj->delete_Supporting_Evidence('test_exon_2');
$exon=$gene_obj->get_Exon('test_exon_2');
$gene_obj->get_supporting_evidence($exon);
if (scalar $exon->each_Supporting_Feature == 0) {
    print "ok 33\n";
}
else {
    print "not ok 33\n";
    print STDERR "Exon supporting evidence still present after deleting!\n";
}

# Checking write method, first get again
my $newgene = $gene_obj->get('test_gene');

$gene_obj->delete('test_gene');

$gene = undef;
#Checking if the gene has been really deleted
eval {
    $gene = $gene_obj->get('test_gene');
};

if ($@ =~ /Error\sretrieving\sgene\swith\sID\:\stest_gene/) {
    print "ok 34\n";
}
else {
    print "not ok 34\n";
    print STDERR "Gene still present after deleting!\n";
} 

# writing gene, adding analysis
if ($newgene->isa('Bio::EnsEMBL::Gene')) {
    $newgene->id( "newgene" );
    my $analysis = Bio::EnsEMBL::Analysis->new( -logic_name => 'TestANA',
      -program => 'GeneWise', -db => 'ESTs' );
    eval {
      $newgene->analysis( $analysis );
      $gene_obj->write( $newgene );
    };
    if( $@ ) {
      print "not ok 35\n";
      print STDERR "Gene write newgene failed!\n$@\n";
    } else {
      print "ok 35\n";
    }
}
else {
    print "not ok 35\n";
    print STDERR "Could not get the test gene from the database!\n";
}

my $anotherGen = $gene_obj->get( "newgene" );
if( $anotherGen->analysis->logic_name eq 'TestANA' ) {
  print "ok 36\n";
} else { 
  print "not ok 36\n";
}


my $newgeneid = $gene_obj->get_new_GeneID;
if ($newgeneid eq 'ENSG00000000001') {
    print "ok 37\n";
}
else {
    print "not ok 37\n";
    print STDERR "New gene id does not look right!\n";
}

my $newtransid = $gene_obj->get_new_TranscriptID;
if ($newtransid eq 'ENST00000000001') {
    print "ok 38\n";
}
else {
    print "not ok 38\n";
    print STDERR "New transcript id does not look right!\n";
}
my $newexonid = $gene_obj->get_new_ExonID;
if ($newexonid eq 'ENSE00000000001') {
    print "ok 39\n";
}
else {
    print "not ok 39\n";
    print STDERR "New exon id does not look right!\n";
}


my $newtranslid = $gene_obj->get_new_TranslationID;
if ($newtranslid eq 'ENSP00000000001') {
    print "ok 40\n";
}
else {
    print "not ok 40\n";
    print STDERR "New protein id does not look right!\n";
}



#Methods still needed to test:
#$gene_obj->get_Virtual_Contig;
#$gene_obj->get_Transcript_in_VC_coordinates;
#$gene_obj->write;
#$gene_obj->write_Exon;
#$gene_obj->write_supporting_evidence;
#$gene_obj->write_Transcript;
#$gene_obj->write_Translation;



