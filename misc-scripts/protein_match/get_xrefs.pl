use strict;

=head1 get_xrefs

=head2 Description
This script take the post processed pmatch output (see process_pmatch.pl) and a file which contains the links of each known gene to other databases (eg: SP to hugo or EMBL) and put them back together in a format suitable for the DBlink tables.

=head2 Options
-mapping:  Name of the file corresponding to postprocessed pmatch
-xrefs: Name of the file linking known genes to other DB
-dbmap: File giving for each known gene its DB
-refseq: If refseq ac is used, file which store for each NP its corresponding NM
-output: Name of the output file 

=head2 Contact
mongin@ebi.ac.uk
birney@ebi.ac.uk

=cut


use Getopt::Long;

my ($mapping,$xrefs,$dbmap,$refseq,$out);

my %map;
my %hash;
my %ref_map;
my %ens2embl;
my %sp2embl;


&GetOptions(
            
            'mapping:s'=>\$mapping,
            'xrefs:s'=>\$xrefs,
	    'dbmap:s'=>\$dbmap,
	    'refseq:s'=>\$refseq,
	    'output:s'=>\$out
            );

open (DBMAP,"$dbmap") || die "Can't open file $dbmap\n"; 
open (XREF,"$xrefs") || die "Can't open file $xrefs\n";
open (MAP,"$mapping") || die "Can't open file $mapping\n";
if ($refseq) {
    open (REFSEQ,"$refseq") || die "Can't open file $refseq\n";
}
open (OUT,">$out") || die "Can't open file $out\n";

while (<DBMAP>) {
    chomp;
#Get put in a hash the corresponding database for an external accession number. Get the infos from a file already processed following this format:
#P31946  SP     
     my ($mapac,$mapdb) = split(/\t/,$_);
     $map{$mapac} = $mapdb;
}

#Read the file by genbank entries (separated by //) 
$/ = "\/\/\n";
while (<REFSEQ>) {
#This subroutine store for each NP (refseq protein accession number) its corresponding NM (DNA accession number)
    my ($prot_ac) = $_ =~ /ACCESSION\s+(\S+)/;
    my ($dna_ac) = $_ =~ /DBSOURCE    REFSEQ: accession\s+(\w+)/;

    $ref_map{$prot_ac} = $dna_ac;
}

#Put back the default (new line) for reading file
$/ = "\n";   
while (<XREF>) {
    chomp;

#SP      P31946  EMBL    X57346 
    my ($xrdb,$xrac,$db,$id) = split (/\t/,$_);
    
    if ($xrdb ne "ENSEMBL") {
	my $both = "$db:$id";
  
	if( !defined $hash{$xrac} ) {
	    $hash{$xrac} = [];
	}
    
	push(@{$hash{$xrac}},$both);
    }
    if ($xrdb eq "ENSEMBL") {
	push(@{$ens2embl{$xrac}},$id);
    }
    if (($xrdb eq "SP") && ($db eq "EMBL")) {
	push(@{$sp2embl{$xrac}},$id);
    }

}

while (<MAP>) {
    chomp;

#P01111  COBP00000000001 100     PRIMARY  
    my ($xr,$ens,$perc,$tag) = split (/\t/,$_);
    if (($tag eq "PRIMARY") || ($tag eq "DUPLICATE")) {

#Its a hack an another solution will have to be found, if the external known gene is a refseq protein accession number get back the equivalent refseq DNA accession number 
	if ($xr =~ /^NP_\d+/) {
	    $xr = $ref_map{$xr};
	}

	if ($sp2embl{$xr}) {
	    my $tot_sp_embl;
	    my $tot_ens_embl;
	    my @sp_embl = @{$sp2embl{$xr}};
	    
	    foreach my $sing1 (@sp_embl) {
		$tot_sp_embl .= $sing1;
	    }
	    

	    my @ens_embl = @{$ens2embl{$xr}};
	    
	    foreach my $sing2 (@sp_embl) {
		$tot_ens_embl .= $sing2;
	    }
	    if ($tot_ens_embl =~ $tot_sp_embl) {
		print OUT "$ens\t$map{$xr}\t$xr\n";
	    }
	}
	
#Print the know gene AC and its database
	print OUT "$ens\t$map{$xr}\t$xr\n";

#Print all of the external database it links to (eg: HUGO)
	foreach my $both (@{$hash{$xr}}){
	    ($a,$b) = split(/:/,$both);
	    print OUT "$ens\t$a\t$b\n";
	}
    }
}


