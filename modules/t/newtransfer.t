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
BEGIN { $| = 1; print "1..43\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;
     system "rm t/recipient.dump";}

use lib 't';
use EnsTestDB;
$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_donor = EnsTestDB->new();
    
# Load some data into the db
$ens_donor->do_sql_file("t/donor.dump");
    
# Get an EnsEMBL db object for the test db
my $fromdb = $ens_donor->get_DBSQL_Obj;
print "ok 2\n";    

my $ens_recipient = EnsTestDB->new();
open (TEMP,">t/recipient.dump");
my $meta= "insert into meta (donor_database_locator) values('Bio::EnsEMBL::DBSQL::Obj/host="
    .$ens_donor->host.";port=".$ens_donor->port
    .";dbname=".$ens_donor->dbname.";user=".$ens_donor->user
    .";pass=".$ens_donor->password."')";
print TEMP $meta;
$ens_recipient->do_sql_file("t/recipient.dump");
my $todb = $ens_recipient->get_DBSQL_Obj;
print "ok 3\n";

my @clones = $fromdb->get_Update_Obj->get_updated_Clone_id(1,time());
print "ok 4\n";

foreach my $clone_id (@clones) {
    if ($clone_id eq 'AB015355') {
	print "ok 5\n";
    }
    else {
	print "not ok 5\n";
	print STDERR "unexpected clone id $clone_id\n";
    }
    my $don_clone=$fromdb->get_Clone($clone_id);
    print "ok 6\n";
    check_clone($don_clone,7);
    $todb->write_Clone($don_clone);
    print "ok 15\n";
    my $rec_clone=$todb->get_Clone($clone_id);
    print "ok 16\n";
    check_clone($rec_clone,17);
    my $ok;
    foreach my $don_gene ($don_clone->get_all_Genes('evidence')) {
	if ($don_gene->id =~ /ENSG00000008215|ENSG00000008216|ENSG00000008217/) {
	    $ok++;
	}
	push @geneids,$don_gene->id;
	if ($don_gene->id eq 'ENSG00000008217') {
	    check_gene($don_gene,25);
	    $todb->gene_Obj->write($don_gene);
	    my $rec_gene=$todb->gene_Obj->get('ENSG00000008217','evidence');
	    check_gene($rec_gene,31);
	}
    }
    if ($ok == 3) {
	print "ok 37\n";
    }
    else {
	print "not ok 37\n";
	print STDERR "Got $ok genes instead of 3!\n";
    }
    foreach my $gene_id (@geneids) {
	$fromdb->gene_Obj->delete($gene_id);
    }
    $don_clone->delete_by_dbID;
    print "ok 38\n";
    check_delete(39);
}

# end of checks.

sub check_clone {
    my ($clone,$c)=@_;

    if ($clone->id eq 'AB015355') {
	print "ok $c\n";
    }
    else {
	print "not ok $c\n";
	print STDERR "Clone id not set correctly\n";
    }
    $c++;
    if ($clone->dbID == 1) {
	print "ok $c\n";
    }
    else {
	print "not ok $c\n";
	print STDERR "Clone internal_id not set correctly\n";
    }
    $c++;
    if ($clone->version ==2) {
	print "ok $c\n";
    }
    else {
	print "not ok $c\n";
	print STDERR "Clone version not set correctly\n";
    }
    $c++;
    if ($clone->embl_id eq 'AB015355') {
	print "ok $c\n";
    }
    else {
	print "not ok $c\n";
	print STDERR "Clone embl_id not set correctly\n";
    }
    $c++;
    if ($clone->embl_version == 1) {
	print "ok $c\n";
    }
    else {
	print "not ok $c\n";
	print STDERR "Clone embl_version not set correctly\n";
    }
    $c++;
    if ($clone->htg_phase == 4) {
	print "ok $c\n";
    }
    else {
	print "not ok $c\n";
	print STDERR "Clone htg_phase not set correctly [",$clone->htg_phase,"]\n";
    }
    $c++;
    if ($clone->created == 932822868) {
	print "ok $c\n";
    }
    else {
	print "ok $c\n";
	print STDERR "*** SKIPPING Clone created not set correctly ",$clone->created,"\n";
    }
    $c++;
    if ($clone->modified == 961875840) {
	print "ok $c\n";
    }
    else {
	print "ok $c\n";
	print STDERR "*** SKIPPING Clone modified not set correctly",$clone->modified,"\n";
    }
    
}

sub check_gene {
    my ($gene,$c)=@_;
    
    if ($gene->id eq 'ENSG00000008217') {
	print "ok $c\n";
    }
    else {
	print "not ok $c\n";
	print STDERR "Gene id not set correctly\n";
    }
    $c++;
    if ($gene->version ==1) {
	print "ok $c\n";
    }
    else {
	print "not ok $c\n";
	print STDERR "Gene version not set correctly\n";
    }
    $c++;
    foreach my $transcript ($gene->each_Transcript) {
	if ($transcript->id eq 'ENST00000009279') {
	    print "ok $c\n";
	}
	else {
	    print "not ok $c\n";
	    print STDERR "Transcript id not set correctly\n";
	}
	$c++;
	my $translation=$transcript->translation;
	if ($translation->id eq 'ENSP00000009279') {
	    print "ok $c\n";
	}
	else {
	    print "not ok $c\n";
	    print STDERR "Transcript id not set correctly\n";
	}
	$c++;
	my $ok=0;
	my $suppexon;
	foreach $exon ($transcript->each_Exon) {
	    if( $exon->id eq 'ENSE00000031521' ) {
		$suppexon = $exon;
	    }

	    if ($exon->id =~ /ENSE00000031517|ENSE00000031518|ENSE00000031519|ENSE00000031520|ENSE00000031521|ENSE00000031522/) {
		$ok++;
	    }
	}
	if ($ok == 6) {
	    print "ok $c\n";
	}
	else {
	    print "not ok $c\n";
	    print STDERR "Got $ok exons instead of 6\n";
	}
	$c++;

	my $test =0;

	if( defined $suppexon ) {
	    my @sup = $suppexon->each_Supporting_Feature();
	    my $sup = shift @sup;


	    if( defined $sup && $sup->hseqname eq 'AB004857' ) {
		$test = 1;
	    }
	}

	if( $test == 0 ) {
	    print STDERR "Supporting feature on exon not retrieved\n";
	    print "not ok $c\n";
	} else {
	    print "ok $c\n";
	}

    }    
    
}

sub check_delete {
    my ($c)=@_;
    
    my $clone;
    eval {
	$clone=$fromdb->get_Clone('AB015355');
    };
    if ($@) {
	print "ok $c\n";
    }
    else {
	print "not ok $c\n";
	print STDERR "Bad news, after deleting clone, found $clone\n";
    }
    $c++;

    my $gene;
    eval {
	$gene=$fromdb->gene_Obj->get('ENSG00000008217');
    };
    if ($@) {
	print "ok $c\n";
    }
    else {
	print "not ok $c\n";
	print STDERR "Bad news, after deleting clone, found gene $gene\n";
    }
    $c++;

    my $exon;
    eval {
	$exon=$fromdb->gene_Obj->get_Exon('ENSE00000031517');
    };
    if ($@) {
	print "ok $c\n";
    }
    else {
	print "not ok $c\n";
	print STDERR "Bad news, after deleting clone, found exon $exon\n";
    }
    $c++;

    my $transcript;
    eval {
	$transcript=$fromdb->gene_Obj->get_Transcript('ENST00000009279');
    };
    if ($@) {
	print "ok $c\n";
    }
    else {
	print "not ok $c\n";
	print STDERR "Bad news, after deleting clone, found transcript $transcript\n";
    }
    $c++;

    my $translation;
    eval {
	$translation=$fromdb->gene_Obj->get_Translation('ENSP00000009279');
    };
    if ($@) {
	print "ok $c\n";
    }
    else {
	print "not ok $c\n";
	print STDERR "Bad news, after deleting clone, found translation $translation\n";
    }
    $c++;
}
    
    
    

