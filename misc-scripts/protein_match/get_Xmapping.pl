#Contact Emmanuel Mongin (mongin@ebi.ac.uk)
use strict;
use Getopt::Long;
use Bio::SeqIO;

BEGIN {
    my $script_dir = $0;
    $script_dir =~ s/(\S+\/)\S+/$1/;
    unshift (@INC, $script_dir);
    require "mapping_conf.pl";
}

my %conf =  %::mapping_conf; # configuration options

# global vars

#Get general options
my $organism   = $conf{'organism'};
my $sptr_swiss = $conf{'sptr_swiss'};
my $out        = $conf{'x_map_out'};

#Get specific options for human/mouse
my $refseq_gnp = $conf{'refseq_gnp'};

#Get specific options for human
my $ens1       = $conf{'ens1'};
my $ens4       = $conf{'ens4'};
my $refseq_pred = $conf{'refseq_pred_gnp'};
my $go           = $conf{'go'}; 

#Get specific options for the mouse
my $mgi_sp     = $conf{'mgi_sp'};
my $mgi_locus  = $conf{'mgi_locus'};

if ((!defined $organism) || (!defined $sptr_swiss) || (!defined $out)) {
    die "\nSome basic options have not been set up, have a look at mapping_conf\nCurrent set up (required options):\norganism: $organism\nsptr_swiss: $sptr_swiss\nx_map: $out\n";
}

my %refseq_map;
my %sp_db;
my %hugo_id;
my %hugo_syn;

open (OUT,">$out") || die "Can't open $out\n";

#First read the SPTR file in swiss format
print STDERR "Reading SPTR file\n";


my $in  = Bio::SeqIO->new(-file => $sptr_swiss, '-format' =>'swiss');

while ( my $seq = $in->next_seq() ) {

    #For each SP entry, store the information concerning the entry
    my $db;
    my $ac = $seq->accession;
    my $id = $seq->display_id;

    #print STDERR "ID: $id\n";

    my ($displ_id,$tag) = split(/:/,$id);

    
    #Humm not good... will be removed soon
    $tag = $seq->primary_id;

    if ($tag eq "STANDARD") {
	$db = "SWISSPROT";
    }
    elsif ($tag eq "PRELIMINARY") {
	$db = "SPTREMBL";
    }
    else {
	die "Try to load unknown SPTR database\n";
    }
    my $un_ac = "$ac:$db";
    my @secs = $seq->get_secondary_accessions;
    my $syns = join(';',@secs);

#These entries correspond to known genes thus they will be flaged as known
    print OUT "$ac\tSPTR\t$ac\t$db\t$displ_id\t$syns\tKNOWN\n";

    $sp_db{$ac} = $db;

    #Then get info about the Xref mapping
    my @dblink = $seq->annotation->each_DBLink;
    
#All of these entries are flaged as XREF, this mean that they have been retrieved using primary DB (SPTR, Refseq)
    foreach my $link(@dblink) {
	if ($link->database eq "EMBL") {
	    print OUT "$ac\tSPTR\t".$link->primary_id."\t".$link->database."\t".$link->primary_id."\t\tXREF\n";

	    my ($protac) = $link->optional_id =~ /^(\w+).\S+/;
	    if ($protac) {
		print OUT "$ac\tSPTR\t".$protac."\tprotein_id\t$protac\t\tXREF\n";
	    }
	}

	if  ($link->database eq "MIM") {
	    print OUT "$ac\tSPTR\t".$link->primary_id."\t".$link->database."\t".$link->primary_id."\t\tXREF\n";
	}

	if  ($link->database eq "PDB") {
	    print OUT "$ac\tSPTR\t".$link->primary_id."\t".$link->database."\t".$link->primary_id."\t\tXREF\n";
	}

    }



#Get Xref mapping specifically for drosophila

    if ($organism eq "drosophila") {
	#Get the gene name, especially useful to replace Hugo name when these don't exist.
	my @gene_names = $seq->annotation->each_gene_name;
	foreach my $g(@gene_names) {
	    print OUT "$ac\tSPTR\t".$g."\tFlyBase\t".$g."\t\tXREF\n";
	}
    }
    
}

if (($organism eq "human") || ($organism eq "mouse")) {
#Read the refseq file in gnp format
    print STDERR "Reading REFSEQ File\n";
    
    open (REFSEQ,"$refseq_gnp") || die "Can't open $refseq_gnp\n";
    
    $/ = "\/\/\n";
    
    while (<REFSEQ>) {
	my ($prot_ac) = $_ =~ /ACCESSION\s+(\S+)/;
	my ($dna_ac) = $_ =~ /DBSOURCE    REFSEQ: accession\s+(\w+)/;
	
	$refseq_map{$dna_ac} = $prot_ac; 
	
#Its a curated Refseq, flag it as known
	print OUT "$prot_ac\tRefSeq\t$prot_ac\tRefSeq\t$prot_ac\t\tKNOWN\n";
       
	my ($mim) = $_ =~ /\/db_xref=\"MIM:(\d+)/;
	my ($locus) = $_ =~ /\/db_xref=\"LocusID:(\d*)/;

#XREF entries	
	if ($mim) {
	    print OUT "$prot_ac\tRefSeq\t$mim\tMIM\t$mim\t\tXREF\n";
	}
	
	if ($locus) {
	    print OUT "$prot_ac\tRefSeq\t$locus\tLocusLink\t$locus\t\tXREF\n";
	}
    }
    close (REFSEQ);
    
    $/ = "\n";
}
    


#Get Xref mapping specifically for human
if ($organism eq "human") {

#Read the Hugo files
    print STDERR "Reading Hugo files\n";
    
    open (ENS4,"$ens4") || die "Can't open $ens4\n";;
    
    while (<ENS4>) {
	chomp;
	my @array = split(/\t/,$_);
	my $hgnc = $array[0];
	my $id = $array[1];
	#my $syn1 = $array[2];
	#my $syn2 = $array[3];
	
	my $syn1 = join (';',split (/,\s/,$array[2]));
	my $syn2 = join (';',split (/,\s/,$array[3]));
	
	my $syn = "$syn1;$syn2";
	
	$hugo_id{$hgnc} = $id;
	$hugo_syn{$hgnc} = $syn;
	#print $hugo_syn{$hgnc};
    }
    close (ENS4);
    
    open (ENS1,"$ens1") || die "Can't open $ens1\n";
    
    while (<ENS1>) {
	chomp;
	my @array = split(/\t/,$_);
	my $hgnc = $array[0];
	
	if ($array[1]) {
	    #my $db = $sp_db{$array[1]};
	    print OUT "$array[1]\tSPTR\t$hgnc\tHUGO\t$hugo_id{$hgnc}\t$hugo_syn{$hgnc}\tXREF\n";
	}
	
	if ($array[2]) {
	    #my $db = $sp_db{$array[1]};
	    print OUT "$array[2]\tRefSeq\t$hgnc\tHUGO\t$hugo_id{$hgnc}\t$hugo_syn{$hgnc}\tXREF\n";
	}
    }
    close (ENS1);

#Read the file containing the NCBI prediction in gnp format
    open (REFSEQ,"$refseq_pred") || die "Can't open $refseq_pred\n";
    
    print STDERR "Reading Refseq file with NCBI predictions\n";

    $/ = "\/\/\n";
    
    while (<REFSEQ>) {
	my ($prot_ac) = $_ =~ /ACCESSION\s+(\S+)/;
	my ($dna_ac) = $_ =~ /DBSOURCE    REFSEQ: accession\s+(\w+)/;
	
	$refseq_map{$dna_ac} = $prot_ac; 
	
#Its a curated Refseq, flag it as predicted
	print OUT "$prot_ac\tRefSeq_pred\t$prot_ac\tRefSeq_pred\t$prot_ac\t\tPRED\n";
       
	my ($mim) = $_ =~ /\/db_xref=\"MIM:(\d+)/;
	my ($locus) = $_ =~ /\/db_xref=\"LocusID:(\d*)/;
    }
    close (REFSEQ);
    
    $/ = "\n";

    open (GO,"$go") || die "Can't open $go\n";

    print STDERR "Reading GO file\n";

    while (<GO>) {
	chomp;
	my @array = split (/\t/,$_);
	print OUT "$array[1]\tSPTR\t$array[4]\tGO\t$array[4]\t\tXREF\n";	
	}

}

#Get Xref mapping specifically for mouse.
if ($organism eq "mouse") {
    my %mgi2sp;
    print STDERR "Getting Xrefs specifically for mouse\n";
    open (MGISP, "$mgi_sp") || die "Can't open $mgi_sp\n";
    while (<MGISP>) {
	chomp;
	my ($mgi,$rik,$a,$b,$c,$sps) = split (/\t/,$_);
      	
	my @sp = split(/\s/,$sps);
	
#put in hash all of the SP entries which correspond to an MGI (this will be used later)
	$mgi2sp{$mgi} = $sps;
	
	foreach my $s(@sp) {
	    print OUT "$s\tSPTR\t$mgi\tMGI\t$mgi\t\tXREF\n";
	}
    }
    open (MGILOC, "$mgi_locus") || die "Can't open $mgi_locus\n";
    
    while (<MGILOC>) {
#The input file gives us MGI to LOCUS, we want SP to LOCUS, thus we use the hash %mgi2sp
	
	chomp;
	my ($mgi,$locus) = split (/\t/,$_);
	
	if ($mgi2sp{$mgi}) {
#There can be many SPs for one MGI
	    my @swiss = split (/\s/,$mgi2sp{$mgi}); 
	    
	    foreach my $sw(@swiss) {
		print OUT "$sw\tSPTR\t$locus\tLOCUS\t$locus\t\tXREF\n";
	    }
	}
    }
}


print STDERR "The output has been written there: $out\n";






