use Test;
use strict;

BEGIN { $| = 1; plan tests => 6;}

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::Seq;
use Bio::EnsEMBL::RawContig;
use Bio::EnsEMBL::Clone;

use Bio::EnsEMBL::Utils::Exception qw(verbose);

######################################################################
# 
# Clone is a deprecated class but needed for backwards 
# compatibility.  These tests ensure that it actually works,
# but verbosity is turned off to avoid all of the deprecated warnings
#
#######################################################################

verbose(-1);

ok(1);

# Database will be dropped when this
# object goes out of scope
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;

ok($multi);

#
## get a core DBAdaptor
##
my $dba = $multi->get_DBAdaptor("core");


#ok($dba);


#my $clone = Bio::EnsEMBL::Clone->new();
#ok($clone);
	

##
## load the clone with some data
##

#$clone->id('dummy_clone');
#$clone->dbID(24);
#$clone->embl_id('dummy_clone');
#$clone->embl_version(1);
#$clone->embl_id(42);
#$clone->version(1);
#$clone->htg_phase(3);
#$clone->created('2002-11-08 12:00:00');
#$clone->modified('2002-11-08 13:00:00');
  
#ok($clone);


##
## create a dummy seq
## 
#my $seq  = Bio::Seq->new(-seq => 'ATGCAGCTAGCATCGATGACATCG',
#                         -id => 'dummy_contig',
#                         -accession => 'dummy_contig');  
#ok($seq);


#my $contig1 = Bio::EnsEMBL::RawContig->new();
 
#my $name =  'dummy_contig';
#$contig1->id($name);
#$contig1->embl_offset(0);
#$contig1->seq($seq);
#ok($contig1);


##
## add the contig to the clone
##
#$clone->add_Contig($contig1);
#ok($clone);


##
## create and add a 2nd contig
##
#my $contig2 = Bio::EnsEMBL::RawContig->new();
#my $name2 =  'dummy_contig2';
#$contig2->id($name2);
#$contig2->embl_offset(0);
#$contig2->seq($seq);

#$clone->add_Contig($contig2);
#ok($clone);


##
## now retrieve the contigs
##
#my $contigs = $clone->get_all_Contigs;
#ok(scalar(@$contigs) == 2);


##
## now let's get a real clone from the test db
#
my $clone_adaptor = $dba->get_CloneAdaptor;
ok($clone_adaptor);


my $real_clone = $clone_adaptor->fetch_by_accession('AL031658');
ok($real_clone);


#
# get all the contigs on this clone
#
my $contigs = $real_clone->get_all_Contigs;
ok(scalar(@$contigs) == 1);


#
# and get all the genes
#
#my $genes = $real_clone->get_all_Genes;
#ok(scalar(@$genes) == 5);


##
## to test the delete method we need to save the tables that we delete from
##

#$multi->save("core","contig","clone","dna","repeat_feature","simple_feature",
#	     "prediction_transcript","protein_align_feature","dna_align_feature");

##
## do the deletion
##
#$clone_adaptor->remove($real_clone);

##
## check the clones
##
#my $sth = $dba->prepare("select * from clone");
#$sth->execute;
##print STDERR "Num clones " . scalar($sth->rows) . "\n";
#ok(scalar($sth->rows) == 11);

##
## check the contigs
##$sth = $dba->prepare("select * from contig");
#$sth->execute;
##print STDERR "Num contigs " . scalar($sth->rows) . "\n";
#ok(scalar($sth->rows) == 11);

##
## check the simple features
##
#$sth = $dba->prepare("select * from simple_feature");
#$sth->execute;
##print STDERR "Num simple_features " . scalar($sth->rows) . "\n";
#ok(scalar($sth->rows) == 116);

##
## check the repeat features
##
#$sth = $dba->prepare("select * from repeat_feature");
#$sth->execute;
##print STDERR "Num repeat_features " . scalar($sth->rows) . "\n";
#ok(scalar($sth->rows) == 1937);

##
## check the protein_align_features
##
#$sth = $dba->prepare("select * from protein_align_feature");
#$sth->execute;
##print STDERR "Num protein_align_features " . scalar($sth->rows) . "\n";
#ok(scalar($sth->rows) == 4727);

##
## check the protein_align_features
##
#$sth = $dba->prepare("select * from dna_align_feature");
#$sth->execute;
##print STDERR "Num dna_align_features " . scalar($sth->rows) . "\n";
#ok(scalar($sth->rows) == 15525);

##
## check the prediction_transcripts
##
#$sth = $dba->prepare("select * from prediction_transcript");
#$sth->execute;
##print STDERR "Num dna_align_features " . scalar($sth->rows) . "\n";
#ok(scalar($sth->rows) == 161);


##
## check the dna
##
#$sth = $dba->prepare("select * from dna");
#$sth->execute;
##print STDERR "Num dna records " . scalar($sth->rows) . "\n";
#ok(scalar($sth->rows) == 11);


## restore the tables for the next test
#$multi->restore("core","contig","clone","dna","repeat_feature","simple_feature",
#		"prediction_transcript","protein_align_feature","dna_align_feature");


##
## just a quick check to see whether the restore has worked
##
#$sth = $dba->prepare("select * from contig");
#$sth->execute;
#ok(scalar($sth->rows) == 12);


#
# This is a snpview method
#
my $raw_contig = undef;
$raw_contig = $real_clone->get_RawContig_by_position(42);
ok($raw_contig);


verbose(0);