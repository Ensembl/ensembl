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

my ($nomeid,$ens1,$ens2,$out,$dbmap);

my %map;
my %en1;
my %en2;
my %hugohash;

&GetOptions(
	    'nomeid:s'=>\$nomeid,
	    'ens1:s'=>\$ens1,
            'ens2:s'=>\$ens2,
	    'output:s'=>\$out,
	    'dbmap:s'=>\$dbmap
            );



open (ENS1,"$ens1") || die "Can't open file $ens1\n";
open (ENS2,"$ens2") || die "Can't open file $ens2\n";
open (NOME,"$nomeid") || die "Can't open file $nomeid\n";
open (DBMAP,"$dbmap") || die "Can't open file $dbmap\n";
open (OUT,">$out") || die "Can't open output file $out\n";


while (<ENS1>) {
    chomp;
    
    #Get rid of the annoying carriage return!
    $_ =~ s/\r//g;
    my ($hgnc,$sp,$refseq) = split(/\t/,$_);
    
    if ($sp) {
	$en1{$sp} = $hgnc;
    }
    if ($refseq) {
	$en1{$refseq} = $hgnc;
    }
}

while (<ENS2>) {
    chomp;
    
    $_ =~ s/\r//g;
    my ($hgnc,$hugo) = split(/\t/,$_);
    
    if (!defined $en2{$hgnc}) {
	$en2{$hgnc} = [];
    }
    
    $en2{$hgnc} = $hugo;
}



while (<NOME>) {
#Get hugo aliases given a hugo primary accession number. For each primary accession number, the aliases are put in a hash of array
    
    chomp;
    my @chunk = split (/\t/,$_);
    
    if ($chunk[8]) {
	$hugohash{$chunk[1]} = $chunk[8];
    }
    
}

while (<DBMAP>) {
    chomp;
    my ($mapac,$mapdb) = split(/\t/,$_);
    my $hugo_ac = $en2{$en1{$mapac}};
    if ($hugo_ac) {

#Print the HUGOs primary accession numbers	
	print OUT "$mapdb\t$mapac\tHUGO\t$hugo_ac\n";
	
	if ($hugohash{$hugo_ac}) {
	    my @syn = split (/, /,$hugohash{$hugo_ac});
	    
	    foreach my $sol (@syn) {
#print the HUGOs aliases
		print OUT "$mapdb\t$mapac\tALIAS\t$sol\n";
	    }
	    
	}
    }
}








