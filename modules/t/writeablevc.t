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
BEGIN { $| = 1; print "1..8\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::DB::ConvertibleVirtualContig;

use lib 't';
use EnsTestDB;
$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/writeablevc.dump");
    
# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;


print "ok 2\n";    


my $contig = $db->get_Contig('contig2');

die "$0\nError fetching contig1 : $!" unless defined ($contig);

print "ok 3\n";

my $wvc     = new Bio::EnsEMBL::DB::ConvertibleVirtualContig(-focuscontig   => $contig,
						 -focusposition => 1,
						 -ori           => 1,
						 -left          => 30,
						 -right         => 30,
                                                 );

die ("$0\nCan't create virtual contig :$!") unless defined ($wvc);

#$wvc->_dump_map(\*STDERR);

$gene = Bio::EnsEMBL::Gene->new();
$gene->id('gene-id-1');
$gene->version(1);
$trans = Bio::EnsEMBL::Transcript->new();
$trans->id('trans-id-1');
$trans->version(1);

$dbl = Bio::Annotation::DBLink->new();
$dbl->database('embl');
$dbl->primary_id('AC000012');
$trans->add_DBLink($dbl);

$dbl = Bio::Annotation::DBLink->new();
$dbl->database('swissprot');
$dbl->primary_id('P000012');
$gene->add_DBLink($dbl);

$trl = Bio::EnsEMBL::Translation->new();
$trl->start_exon_id('exon-1');
$trl->end_exon_id('exon-2');
$trl->start(2);
$trl->end(13);
$trl->id('trl-id');
$trl->version(1);
$trans->translation($trl);
$trans->created(1);
$trans->modified(1);
$gene->add_Transcript($trans);
$gene->created(1);
$gene->modified(1);

$exon = Bio::EnsEMBL::Exon->new();
$exon->id('exon-1');
$exon->start(5);
$exon->end(15);
$exon->strand(1);
$exon->version(1);
$exon->phase(0);
$exon->created(1);
$exon->modified(1);
$exon->contig_id($wvc->id);
$trans->add_Exon($exon);
$exon{'exon-1'} = $exon;
$exon->attach_seq($wvc->primary_seq);

$sf = Bio::EnsEMBL::FeatureFactory->new_feature_pair();

$sf->start(6);
$sf->end(14);
$sf->hstart(100);
$sf->hend(110);
$sf->strand(1);
$sf->seqname($wvc->id);
$sf->hseqname('other');
$sf->score(100);
$sf->primary_tag('similarity');
$sf->source_tag('someone');
$sf->feature2->primary_tag('similarity');
$sf->feature2->source_tag('someone');

$analysis = Bio::EnsEMBL::FeatureFactory->new_analysis();
$analysis->program('program');
$analysis->program_version('version-49');
$analysis->gff_source('source');
$analysis->gff_feature('feature');

$sf->analysis($analysis);
$sf->feature2->analysis($analysis);
$sf->hstrand(1);
$sf->hscore(100);

$exon->add_Supporting_Feature($sf);

$exon = Bio::EnsEMBL::Exon->new();
$exon->id('exon-2');
$exon->start(25);
$exon->end(35);
$exon->strand(-1);
$exon->version(1);
$exon->phase(0);
$exon->created(1);
$exon->modified(1);
$exon->contig_id($wvc->id);
$trans->add_Exon($exon);
$exon{'exon-2'} = $exon;
$exon->attach_seq($wvc->primary_seq);

$newgene = $wvc->convert_Gene_to_raw_contig($gene);
$gene_obj = $db->gene_Obj;
$gene_obj->use_delayed_insert(0);
$gene_obj->write($newgene);
print "ok 4\n";

# let the write complete...
#$sth = $db->prepare("flush tables");
#$sth->execute;

($newgene) = $db->gene_Obj->get_array_supporting('evidence','gene-id-1');

if( !defined $newgene ) {
    print "not ok 5\n";
} else {
    $error = 0;
    foreach $exon ( $newgene->each_unique_Exon ) {
	if( exists $exon{$exon->id} ) {
	    if( $exon->seq->seq ne $exon{$exon->id}->seq->seq ) {
		print STDERR "For exon ",$exon->id," sequences do not agree!\n";
		print STDERR "RC :",$exon->seq->seq,"\nVC :",$exon{$exon->id}->seq->seq,"\n";
		
		$error =1;
	    } else {
		if( $exon->id eq 'exon-1') {
		    @sf = $exon->each_Supporting_Feature();
		    $sf = shift @ sf;
		    if( !defined $sf || !ref $sf || $sf->hseqname() ne 'other' ) {
			print STDERR "did not retrieve supporting evidence on exon-1 [$sf]";
			$error = 1;
		    }
		}
	    }
	} else {
	    print STDERR "Weird exon in output...\n";
	}
    }
    
    $savedgene = $newgene;
    
    if( $error == 1 ) {
	print "not ok 5\n";
    } else {
	print "ok 5\n";
    }
}


$wvc     = new Bio::EnsEMBL::DB::ConvertibleVirtualContig(-focuscontig   => $contig,
						 -focusposition => 1,
						 -ori           => -1,
						 -left          => 30,
						 -right         => 30,
                                                 -gene_obj      => $db->gene_Obj);




$gene = Bio::EnsEMBL::Gene->new();
$gene->id('rgene-id-1');
$gene->version(1);
$trans = Bio::EnsEMBL::Transcript->new();
$trans->id('rtrans-id-1');
$trans->version(1);
$trl = Bio::EnsEMBL::Translation->new();
$trl->start_exon_id('rexon-1');
$trl->end_exon_id('rexon-2');
$trl->start(2);
$trl->end(13);
$trl->id('rtrl-id');
$trl->version(1);
$trans->translation($trl);
$trans->created(1);
$trans->modified(1);
$gene->add_Transcript($trans);
$gene->created(1);
$gene->modified(1);

$exon = Bio::EnsEMBL::Exon->new();
$exon->id('rexon-1');
$exon->start(5);
$exon->end(15);
$exon->strand(1);
$exon->version(1);
$exon->phase(0);
$exon->created(1);
$exon->modified(1);
$exon->contig_id($wvc->id);
$trans->add_Exon($exon);
$exon{'rexon-1'} = $exon;
$exon->attach_seq($wvc->primary_seq);

$exon = Bio::EnsEMBL::Exon->new();
$exon->id('rexon-2');
$exon->start(25);
$exon->end(35);
$exon->strand(-1);
$exon->version(1);
$exon->phase(0);
$exon->created(1);
$exon->modified(1);
$exon->contig_id($wvc->id);
$trans->add_Exon($exon);
$exon{'rexon-2'} = $exon;
$exon->attach_seq($wvc->primary_seq);

$newgene = $wvc->convert_Gene_to_raw_contig($gene);
$db->write_Gene($newgene);



$newgene = $db->get_Gene('rgene-id-1');

$error = 0;
foreach $exon ( $newgene->each_unique_Exon ) {
    if( exists $exon{$exon->id} ) {
	$exon->seq->seq();
	$exon{$exon->id}->seq->seq();
	if( $exon->seq->seq ne $exon{$exon->id}->seq->seq ) {
	    print STDERR "For exon ",$exon->id," sequences do not agree!\n";
	    print STDERR "RC :",$exon->seq->seq,"\nVC :",$exon{$exon->id}->seq->seq,"\n";
	    $error =1;
	} else {
	    print STDERR "Got exon match for",$exon->id,"\n";
	}
    } else {
	print STDERR "Weird exon in output...\n";
    }
}

if( $error == 1 ) {
    print "not ok 6\n";
} else {
    print "ok 6\n";
}

@dblink = $savedgene->each_DBLink();
$dbl = shift @dblink;

if( !defined $dbl || $dbl->database ne 'swissprot' ) {
    print "not ok 7\n";
} else {
  print "ok 7\n";
}

@trans = $savedgene->each_Transcript();
$trans = shift @trans;


@dblink = $trans->each_DBLink();
$dbl = shift @dblink;

if( !defined $dbl || $dbl->database ne 'embl' ) {
    print "not ok 8\n";
}    else {
  print "ok 8\n";
}




