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
BEGIN { $| = 1; print "1..10\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::Utils::GTF_handler;
use lib 't';

$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $gtfh=Bio::EnsEMBL::Utils::GTF_handler->new();
my @genes=&parse_genes('genes.gtf',1);

open (DUMP,">t/gene_dump.gtf");
$gtfh->dump_genes(\*DUMP,@genes);
close (DUMP);
print "ok 6\n";
&parse_genes('gene_dump.gtf',6);

sub parse_genes {
    my ($file,$c)=@_;
    
    #Use new parser module, to empty gene array
    my $gtfh=Bio::EnsEMBL::Utils::GTF_handler->new();
    open (PARSE,"t/$file") || die("Could not open $file for Fasta stream reading $!");
    
    @genes=$gtfh->parse_file(\*PARSE);
    close (PARSE);
    $c++;
    print "ok $c\n";

    my $gc;
    my $tc;
    my $ec;
    foreach my $gene (@genes) {
	if ($gene->id =~ /F15G000000000/){
	    $gc++;
	}
	foreach my $trans ($gene->each_Transcript) {
	    if ($trans->id =~ /F15T000000000/){
		$tc++;
	    }
	    foreach my $exon ($gene->each_unique_Exon) {
		if ($exon->id =~ /ENS-F15G0000/) {
		    $ec++;
		}
	    }
	    #print STDERR "   translation start ".$trans->translation->start."\n";
	    #print STDERR "   translation end ".$trans->translation->end."\n";
	}
    }
    $c++;
    if ($gc == 7) {
	print "ok $c\n";
    }
    else {
	print "not ok $c\n";
	print STDERR "Got $gc genes instead of 7!\n";
    }
    $c++;
    if ($tc == 10) {
	print "ok $c\n";
    }
    else {
	print "not ok $c\n";
	print STDERR "Got $tc transcripts instead of 10!\n";
    }
    $c++;
    if ($ec == 21) {
	print "ok $c\n";
    }
    else {
	print "not ok $c\n";
	print STDERR "Got $ec exons instead of 21!\n";
    }
    return @genes;
}




