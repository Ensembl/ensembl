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

my $sptr_swiss = $conf{'sptr_swiss'};
my $refseq_gnp = $conf{'refseq_gnp'};
my $ens1       = $conf{'ens1'};
my $ens4       = $conf{'ens4'};
my $out        = $conf{'x_map'};

my %refseq_map;
my %sp_db;
my %hugo_id;
my %hugo_syn;

open (OUT,">$out") || die "Can't open $out\n";;

#First read the SPTR file in swiss format
print STDERR "Reading SPTR file\n";

my $in  = Bio::SeqIO->new(-file => $sptr_swiss, '-format' =>'swiss');

while ( my $seq = $in->next_seq() ) {

    #For each SP entry, store the information concerning the entry
    my $db;
    my $ac = $seq->accession;
    my $id = $seq->display_id;
    my ($displ_id,$tag) = split(/:/,$id);

    if ($tag eq "STANDARD") {
	$db = "SWISS-PROT";
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
    print OUT "$ac\tSPTR\t$ac\t$db\t$displ_id\t$syns\n";

    $sp_db{$ac} = $db;

    #Then get info about the Xref mapping
    my @dblink = $seq->annotation->each_DBLink;
    
    foreach my $link(@dblink) {
	if ($link->database eq "EMBL") {
	    print OUT "$ac\tSPTR\t".$link->primary_id."\t".$link->database."\t".$link->primary_id."\t\n";

	    my ($protac) = $link->optional_id =~ /^(\w+).\S+/;
	    if ($protac) {
		print OUT "$ac\tSPTR\t".$protac."\tprotein_id\t$protac\t\n";
	    }
	}

	if  ($link->database eq "MIM") {
	    print OUT "$ac\tSPTR\t".$link->primary_id."\t".$link->database."\t".$link->primary_id."\t\n";
	}
    }
}


#Read the refseq file in gnp format
print STDERR "Reading REFSEQ File\n";

open (REFSEQ,"$refseq_gnp") || die "Can't open $refseq_gnp\n";

$/ = "\/\/\n";

while (<REFSEQ>) {
    my ($prot_ac) = $_ =~ /ACCESSION\s+(\S+)/;
    my ($dna_ac) = $_ =~ /DBSOURCE    REFSEQ: accession\s+(\w+)/;

    $refseq_map{$dna_ac} = $prot_ac; 
    
    print OUT "$prot_ac\tRefSeq\t$prot_ac\tRefSeq\t$prot_ac\t\n";
       
    my ($mim) = $_ =~ /\/db_xref=\"MIM:(\d+)/;
    my ($locus) = $_ =~ /\/db_xref=\"LocusID:(\d*)/;

    if ($mim) {
	print OUT "$dna_ac\tRefSeq\t$mim\tMIM\t$mim\t\n";
    }

    if ($locus) {
	print OUT "$dna_ac\tRefSeq\t$locus\tLocusLink\t$locus\t\n";
    }
}
close (REFSEQ);

$/ = "\n";


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
	print OUT "$array[1]\tSPTR\t$hgnc\tHUGO\t$hugo_id{$hgnc}\t$hugo_syn{$hgnc}\n";
    }
    
    if ($array[2]) {
	#my $db = $sp_db{$array[1]};
	print OUT "$array[2]\tRefSeq\t$hgnc\tHUGO\t$hugo_id{$hgnc}\t$hugo_syn{$hgnc}\n";
    }
}
close (ENS1);


