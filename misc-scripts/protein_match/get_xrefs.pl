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
my %embl_clone;

&GetOptions(
            
            'mapping:s'=>\$mapping,
            'xrefs:s'=>\$xrefs,
	    'dbmap:s'=>\$dbmap,
	    'refseq:s'=>\$refseq,
	    'output:s'=>\$out
            );

#perl ../../../src/ensembl-live/misc-scripts/protein_match/get_xrefs.pl -mapping ../map_outputs/map.total -xrefs ../sec_outputs/xref.total -dbmap ../sec_outputs/mapdb.map -refseq ../primary/hs.gnp -output final.map

open (DBMAP,"$dbmap") || die "Can't open file $dbmap\n"; 
open (XREF,"$xrefs") || die "Can't open file $xrefs\n";
open (MAP,"$mapping") || die "Can't open file $mapping\n";
if ($refseq) {
    open (REFSEQ,"$refseq") || die "Can't open file $refseq\n";
}
open (OUT,">$out") || die "Can't open file $out\n";

#open (CLONE,"clones.txt") || die "Can't open file\n";

#Put in a hash all of the embl clones used by Ensembl
#while (<CLONE>) {
#    chomp;
#    my ($embl_ac,$id) = split(/\t/,$_);
#    print "$embl_ac\n";
#    $embl_clone{$embl_ac}=1;
#}


while (<DBMAP>) {
    chomp;
#Get put in a hash the corresponding database for an external accession number.Get the infos from a file already processed following this format:
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
	my $both = "$db&$id";
  
	if( !defined $hash{$xrac} ) {
	    $hash{$xrac} = [];
	}
    
	push(@{$hash{$xrac}},$both);
    }

#Get the embl clone corresponding for each Ensembl peptides
    #if (($xrdb eq "ENSEMBL")) {
	
	#push(@{$ens2embl{$xrac}},$id);
    #}

#Get the embl ACs for each SP and SPTREMBL proteins
    #if ((($xrdb eq "SP") || ($xrdb eq "SPTREMBL")) && ($db eq "EMBL")) {
	 #print "$id\n";
        #if ($embl_clone{$id}) {
	    
	 #   push(@{$sp2embl{$xrac}},$id);
	#}
    #}
}

while (<MAP>) {
    chomp;
    
#P01111  COBP00000000001 100     PRIMARY  
    my ($xr,$ens,$perc,$tag) = split (/\t/,$_);
    if ($xr =~ /^\w+-\d{2}/) {
	($xr) = $xr =~ /^(\w+)-\d{2}/;
    }
#Hack to be taken away
    #my ($en1,$en2) = $ens =~ /(\w{3})P(\d+)/;
    #my $enst = $en1."T".$en2;
    
#For now take primary or duplicates and only matches which correspond to more than 25% of the external peptide. These criteria will have to be lowered up.
    if ((($tag eq "PRIMARY") || ($tag eq "DUPLICATE")) && ($perc >= 25)) {
	
#Its a hack an another solution will have to be found, if the external known gene is a refseq protein accession number get back the equivalent refseq DNA accession number 
	if ($xr =~ /^NP_\d+/) {
	    $xr = $ref_map{$xr};
	}
	
#If the external peptide correspond to an embl clone, we will take the match only if the Ensembl peptide correspond to the same clone (at least one exon)
	#if ($sp2embl{$xr}) {
	#    print "$xr\t".@{$sp2embl{$xr}}."\n";
	#    my $tot_sp_embl;
	#    my $tot_ens_embl;
	#    my @sp_embl = @{$sp2embl{$xr}};
	
	#    foreach my $sing1 (@sp_embl) {
	#	#print "$sing1\n";
	#	$tot_sp_embl .= $sing1;
	
	#    }
	  	  	    
	#    if ($ens2embl{$enst}) {
	#	my @ens_embl = @{$ens2embl{$enst}};
	    
	#	foreach my $sing2 (@sp_embl) {
	#	    $tot_ens_embl .= $sing2;
	#	}
	#	if ($tot_ens_embl =~ $tot_sp_embl) {
	#	    print  OUT "$ens\t$map{$xr}\t$xr\n";
	#	}
	#	else {
	#	    #print "no\n";
	#	}
	#    }
	#}
	#else {
#Print the know gene AC and its database
	#print OUT "$ens\t$map{$xr}\t$xr\n";
	#}
	
	if (!defined $map{$xr}) {
	    print STDERR "can't find primary $xr\n";
	}

	print OUT "$ens\t$map{$xr}\t$xr\n";
	
	#Print all of the external database it links to (eg: HUGO)
	    foreach my $both (@{$hash{$xr}}){
		if (!defined $hash{$xr}) {
	    print STDERR "can't find Xref $xr\n";
	}
		($a,$b) = split(/&/,$both);
		print OUT "$ens\t$a\t$b\n";
	    }
    }
}


