use strict;

=head1 get_hugo_mapping

=head2 Description

This script reads a set of files produced by HUGO and built a file which will be used to get the DBlink table.The format of the file is the following

SP      P27348  HUGO    YWHAQ 

(known database\t known ac\t hugo or alias\t hugo ac)

=head2 Options

The different options only deal with file names

-nomeid: Hugo file (http://www.gene.ucl.ac.uk/public-files/nomen/nomeids.txt)

-ens1: Hugo file (http://www.gene.ucl.ac.uk/public-files/nomen/ens1.txt)

-ens2: Hugo file (http://www.gene.ucl.ac.uk/public-files/nomen/ens2.txt)

-output: Filename were the output should be written

-dbmap: File which give the corresponding database for each known protein

=cut

use Getopt::Long;

#perl ../../../src/ensembl-live/misc-scripts/protein_match/get_hugo_mapping.pl -ens1 ../secondary/ens1.txt -ens2 ../secondary/ens2.txt -ens4 ../secondary/ens4.txt -ens5 ../secondary/ens5.txt -out hugo.map -dbmap mapdb.map

my ($ens1,$ens2,$ens4,$ens5,$out,$dbmap);

my %map;
my %hugo_sp;
my %hugo_refseq;
my %en2;
my %hugohash;

&GetOptions(
	    'ens1:s'=>\$ens1,
            'ens2:s'=>\$ens2,
	    'ens4:s'=>\$ens4,
	    'ens5:s'=>\$ens5,
	    'dbmap:s'=>\$dbmap,
	    'output:s'=>\$out
            );



open (ENS1,"$ens1") || die "Can't open file $ens1\n";
open (ENS2,"$ens2") || die "Can't open file $ens2\n";
open (ENS4,"$ens4") || die "Can't open file $ens4\n";
open (ENS5,"$ens5") || die "Can't open file $ens5\n";
open (DBMAP,"$dbmap") || die "Can't open file $dbmap\n";
open (OUT,">$out") || die "Can't open output file $out\n";
open (ERROR,">hugo.err") || die "Can't open output file hugo.err\n";

while (<DBMAP>) {
    chomp;
     my ($mapac,$mapdb) = split(/\t/,$_);
     $map{$mapac} = $mapdb;
}

while (<ENS1>) {
    chomp;
    #Get hugo id
    #Get rid of the annoying carriage return!
    $_ =~ s/\r//g;
    my ($hgnc,$sp,$refseq) = split(/\t/,$_);

    if ($sp) {
	print OUT "$map{$sp}\t$sp\tHUGOID\t$hgnc\n";
    }

    if ($refseq) {
	print OUT "$map{refseq}\t$refseq\tHUGOID\t$hgnc\n";
    }

    if ($sp) {
	$hugo_sp{$hgnc} = $sp;
    }
    if ($refseq) {
	$hugo_refseq{$hgnc} = $refseq;
    }
}

while (<ENS2>) {
    chomp;
#Get hugo symbol
    $_ =~ s/\r//g;
    my ($hgnc1,$hugo) = split(/\t/,$_);
    
#    if (!defined $hugo_sp{$hgnc1}) {
#	print ERROR "Can't map back $hugo_sp{$hgnc} (ENS2)\n";
#    }

#    if (!defined $hugo_refseq{$hgnc1}) {
#	print ERROR "Can't map back $hugo_refseq{$hgnc} (ENS2)\n";
#    }

    if ($hugo_sp{$hgnc1}) {
	print OUT "$map{$hugo_sp{$hgnc1}}\t$hugo_sp{$hgnc1}\tHUGOSYMBOL\t$hugo\n";
    }

    if ($hugo_refseq{$hgnc1}) { 
	print OUT "$map{$hugo_refseq{$hgnc1}}\t$hugo_refseq{$hgnc1}\tHUGOSYMBOL\t$hugo\n";
    }

    if (!defined $en2{$hgnc1}) {
	$en2{$hgnc1} = [];
    }
    
    $en2{$hgnc1} = $hugo;
}

while (<ENS4>) {
#Get hugo aliases given a hugo primary accession number. For each primary accession number, the aliases are put in a hash of array
    
    chomp;
    my ($hgnc2, $symbol, $alias, $withdrawn) = split (/\t/,$_);

    if ((defined $hugo_sp{$hgnc2}) && (defined $alias)) {
	my @aliases1 = split (/, /,$alias);
	foreach my $aliase1 (@aliases1) {
	    print OUT "$map{$hugo_sp{$hgnc2}}\t$hugo_sp{$hgnc2}\tHUGOALIAS\t$aliase1\n";
	}
    }
    
    if ((defined $hugo_sp{$hgnc2}) && ($withdrawn =~ /\S+/)) {
	my @withdrawns1 = split (/, /,$withdrawn);
	foreach my $withdrawn1 (@withdrawns1) {
	    print OUT "$map{$hugo_sp{$hgnc2}}\t$hugo_sp{$hgnc2}\tHUGOWITHDRAWN\t$withdrawn1\n";
	}
    }
    
    if ((defined $hugo_refseq{$hgnc2}) && (defined $alias)) {
	my @aliases2 = split (/, /,$alias);
	foreach my $aliase2 (@aliases2) {
	    print OUT "$map{$hugo_sp{$hgnc2}}\t$hugo_sp{$hgnc2}\tHUGOALIAS\t$aliase2\n";
	}
    }
    
    if ((defined $hugo_refseq{$hgnc2}) && ($withdrawn =~ /\S+/)) {
	my @withdrawns2 = split (/, /,$withdrawn);
	foreach my $withdrawn2 (@withdrawns2) {
	    print OUT "$map{$hugo_sp{$hgnc2}}\t$hugo_sp{$hgnc2}\tHUGOWITHDRAWN\t$withdrawn2\n";
	}
    }
}

while (<ENS5>) {
#Use Hugo mapping to get EC numbers    
    chomp;
    my ($hgnc3, $symbol1, $name, $ec, $sp) = split (/\t/,$_);
    
    if ((defined $hugo_sp{$hgnc3}) && (defined $ec)) {
	 print OUT "$map{$hugo_sp{$hgnc3}}\t$hugo_sp{$hgnc3}\tEC\t$ec\n";
     }

    if ((defined $hugo_refseq{$hgnc3}) && (defined $ec)) {
	 print OUT "$map{$hugo_sp{$hgnc3}}\t$hugo_sp{$hgnc3}\tEC\t$ec\n";
     }
}



