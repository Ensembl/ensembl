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

my $org_list = $conf{'organism_list'};

#Get general options
my $organism   = $conf{'organism'};
my $sptr_swiss = $conf{'sptr_swiss'};
my $out        = $conf{'x_map_out'};

#Get specific options for human/mouse/rat/drosophila
my $refseq_gnp = $conf{'refseq_gnp'};

#Get specific options for human
my $ens1       = $conf{'ens1'};
my $ens4       = $conf{'ens4'};
my $go          = $conf{'go'}; 
my $gkb         = $conf{'gkb'}; 

#Get specific options for the mouse
my $mgi_sp     = $conf{'mgi_sp'};
my $mgi_locus  = $conf{'mgi_locus'};

#Get specific options for anopheles
my $sub_genes  = $conf{'submitted_genes'};

#Get specific options for elegans
my $eleg_nom   = $conf{'eleg_nom'};

#Get specific options for zebrafish
my $zeb_gene    = $conf{'zeb_gene'};
my $zeb_dblink  = $conf{'zeb_dblink'};


my $briggsae_peptides = $conf{'briggsae_hybrid'};

my $help;


&GetOptions(
	    'help' => \$help,
	    );

if ($help) {
    print STDERR $conf{'help'}."\n";
    exit();
}

#Check that the configuration file has been well filled in for each different organism
#Beginning of check

my %check;
my $seenorg = 0;

#Check if the organism is correct
foreach my $or (@{$org_list}) {
    if ($or eq $organism) {
	$seenorg = 1;
    }
}

if ($seenorg == 0) {
    print STDERR "Either the organism name you are using ($organism) is not define or is not allowed\n";
    print STDERR "Here is a list of authorised organisms:\n";
    foreach my $or (@{$org_list}) {
	print STDERR "$or\n";
    }

    exit();
}


#Organism specific checks
if($organism eq "human") {
    $check{'sptr_swiss'} = $conf{'sptr_swiss'};
    $check{'x_map_out'} = $conf{'x_map_out'};
    $check{'refseq_gnp'} = $conf{'refseq_gnp'};
    $check{'ens1'} = $conf{'ens1'};
    $check{'ens4'} = $conf{'ens4'};
    $check{'go'} = $conf{'go'};
    $check{'gkb'} = $conf{'gkb'};
    
    foreach my $k (keys %check) {
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }

}

if ($organism eq "mouse") {
    $check{'sptr_swiss'} = $conf{'sptr_swiss'};
    $check{'x_map_out'} = $conf{'x_map_out'};
    $check{'refseq_gnp'} = $conf{'refseq_gnp'};
    $check{'mgi_sp'} = $conf{'mgi_sp'};
    $check{'mgi_locus'} = $conf{'mgi_locus'};

    foreach my $k (keys %check) {
	print STDERR $check{$k}."\n";
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }

}

if ($organism eq "elegans") {
    $check{'sptr_swiss'} = $conf{'sptr_swiss'};
    $check{'x_map_out'} = $conf{'x_map_out'};
    $check{'elegans_nom'} = $conf{'elegans_nom'};

    foreach my $k (keys %check) {
	print STDERR $check{$k}."\n";
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }

}

if ($organism eq "anopheles") {
    $check{'sptr_swiss'} = $conf{'sptr_swiss'};
    $check{'x_map_out'} = $conf{'x_map_out'};
    $check{'submitted_genes'} = $conf{'submitted_genes'};

    foreach my $k (keys %check) {
	print STDERR $check{$k}."\n";
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }

}

if ($organism eq "drosophila") {
    $check{'sptr_swiss'} = $conf{'sptr_swiss'};
    $check{'x_map_out'} = $conf{'x_map_out'};
    $check{'refseq_gnp'} = $conf{'refseq_gnp'};

    foreach my $k (keys %check) {
	print STDERR $check{$k}."\n";
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }

}

if ($organism eq "rat") {
    $check{'sptr_swiss'} = $conf{'sptr_swiss'};
    $check{'x_map_out'} = $conf{'x_map_out'};
    $check{'refseq_gnp'} = $conf{'refseq_gnp'};

    foreach my $k (keys %check) {
	print STDERR $check{$k}."\n";
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }

}

if ($organism eq "zebrafish") {
    $check{'sptr_swiss'} = $conf{'sptr_swiss'};
    $check{'x_map_out'} = $conf{'x_map_out'};
    $check{'zeb_gene'} = $conf{'zeb_gene'};
    $check{'zeb_dblink'} = $conf{'zeb_dblink'};

    foreach my $k (keys %check) {
	print STDERR $check{$k}."\n";
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }

}

if ($organism eq "fugu") {
  $check{'sptr_swiss'} = $conf{'sptr_swiss'};
  $check{'x_map_out'} = $conf{'x_map_out'};
  
  foreach my $k (keys %check) {
    print STDERR $check{$k}."\n";
    if ($check{$k} !~ /(\S+)/) {
      usage();
    }
  } 
}


#End of check

if ((!defined $organism) || (!defined $sptr_swiss) || (!defined $out)) {
    die "\nSome basic options have not been set up, have a look at mapping_conf\nCurrent set up (required options):\norganism: $organism\nsptr_swiss: $sptr_swiss\nx_map: $out\n";
}

my %refseq_map;
my %sp_db;
my %hugo_id;
my %hugo_syn;

open (OUT,">$out") || die "Can't open OUTFILE $out\n";

#First read the SPTR file in swiss format
print STDERR "Reading SPTR file\n";

my ($swiss, $ac, $id) = &parse_sp_file($sptr_swiss);
&process_parsed_sp($swiss, $ac, $id, \*OUT);


if (($organism eq "human") || ($organism eq "mouse") || ($organism eq "rat") || ($organism eq "zebrafish") || ($organism eq "drosophila")) {
  #Read the refseq file in gnp format
  print STDERR "Reading REFSEQ File\n";
  
  open (REFSEQ,"$refseq_gnp") || die "Can't open Refseq gnp $refseq_gnp\n";
  
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
    
    open (ENS4,"$ens4") || die "Can't open hugo ens4 $ens4\n";;
    
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
    
    open (ENS1,"$ens1") || die "Can't open hugo ens1 $ens1\n";
    
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

    open (GO,"$go") || die "Can't open Go file $go\n";

    print STDERR "Reading GO file\n";

    while (<GO>) {
	chomp;
	my @array = split (/\t/,$_);
	$array[9] =~ s/\'/\\\'/g;
	print OUT "$array[1]\tSPTR\t$array[4]\tGO\t$array[6]\t$array[9]\tXREF\n";	
	}
    
    close (GO);

#Get the GKB (Genome Knowledge Database) mapping
#Contacts: Imre Vastrik <vastrik@ebi.ac.uk>, Hester <eschmidt@ebi.ac.uk>

    open (GKB,"$gkb") || die "Can't open GKB data: $gkb";

    while (<GKB>) {
	chomp;
	my ($sp,$id) = split;
	print OUT "$sp\tSPTR\t$id\tGKB\t$id\t\tXREF\n";
    }
}

#Get Xref mapping specifically for mouse.
if ($organism eq "mouse") {
    my %mgi2sp;
    # 20/02/04:th MGI changed files for MGI to SP mappings.  File format also changed
    # moving $sps from 6th to 7th column.
    open (MGISP, "$mgi_sp") || die "Can't open $mgi_sp\n";
    while (<MGISP>) {
      chomp;
      my ($mgi,$rik,$a,$b,$c,$d,$sps) = split (/\t/,$_);
      my @sp = split(/\s/,$sps);
      #put in hash all of the SP entries which correspond to an MGI (this will be used later)
      $mgi2sp{$mgi} = $sps;
    }
    close MGISP;
    open (MGILOC, "$mgi_locus") || die "Can't open $mgi_locus\n";
    my %mgi_got;
    my %mgi_syns;
    while (<MGILOC>) {
      #The input file gives us MGI to LOCUS, we want SP to LOCUS, thus we use the hash %mgi2sp
      chomp;
      my ($mgi,$locus) = split (/\t/,$_);
      if($mgi_got{$mgi}){
	if(!$mgi_syns{$mgi}){
	  $mgi_syns{$mgi} = [];
	  push(@{$mgi_syns{$mgi}}, $locus);
	}else{
	  push(@{$mgi_syns{$mgi}}, $locus);
	}	    
      }else{
	$mgi_got{$mgi} = $locus;
      }      
    }
    close MGILOC;
	
    my @mgi_ids = keys(%mgi_got);
    foreach my $mgi(@mgi_ids){
      my @syns;
      my $syns;
      if($mgi_syns{$mgi}){
	@syns = @{$mgi_syns{$mgi}};
	$syns = join(';',@syns);
      }
      my $locus = $mgi_got{$mgi};
      if ($mgi2sp{$mgi}) {
	#There can be many SPs for one MGI
	my @swiss = split (/\s/,$mgi2sp{$mgi}); 
	
	foreach my $sw(@swiss) {
	  print OUT "$sw\tSPTR\t$mgi\tMarkerSymbol\t$locus\t$syns\tXREF\n";
	}
      }
    }
}

#Get specific xmapping for anopheles.
if($organism eq "anopheles") {
    print STDERR "Getting Xref specifically for anopheles\n";
    open (SUB,"$sub_genes") || die "Can't open $sub_genes\n";

    while(<SUB>) {
#This input uses fasta files with the header as follow:
#>AC\tGene name
#This files can be dumped using the following utility scripts in /ensembl-genename/scripts
	chomp;
	if ($_ =~ />/) {
	    my ($ac,$name) = $_ =~ />(\S+)\t(\S+)/;
	    print OUT "$ac\tAnopheles_symbol\t$ac\tAnopheles_symbol\t$ac\t\tKNOWN\n";
	}
    }
}



#Get specific xmapping for zebrafishq
if($organism eq "zebrafish") {
    print STDERR "Getting Xref specifically for zebrafish";
    open (ZEBGENE,"$zeb_gene") || die "Can't open $zeb_gene";
    open (ZEBLINK,"$zeb_dblink") || die "Can't open $zeb_dblink";
    
    my %map;
    
     while (<ZEBLINK>) {
	chomp;
	my ($ac,$db,$ext_ac) = split /\t/;
        my (@accs) = split /,/, $ext_ac if ($ext_ac =~ /,/);
        @accs = ($ext_ac) unless ($ext_ac =~ /,/);

        if (($db eq "SWISS-PROT") || ($db eq "RefSeq") || ($db eq "LocusLink") || ($db eq "Genbank") || ($db eq "UniGene")) {
	    foreach my $newac (@accs) {
                push @{$map{$ac}}, "$newac:$db";
            }
        }
    }
    
    while (<ZEBGENE>) {
	chomp;
        next if /ZDB\S+\s+\w\w:/;
        next if /cnf/;
        my ($ac,$desc,$id,$num) = split /\t/;
	if (exists $map{$ac}) {
            my @accs = @{$map{$ac}} ;

            for (my $k = 0; $k < scalar @accs; $k++) {
                my $pri = $accs[$k];
                print STDERR "$_\n" unless ($pri);
                next unless ($pri);
	        my ($displ_id,$tag) = split(/:/,$pri);
	        if ($tag eq "SWISS-PROT") {
	            $tag = "SPTR";
	        }
	        print OUT "$displ_id\tSPTR\t$ac\tZFIN_ID\t$id\t\tXREF\n";
	        print OUT "$displ_id\tSPTR\t$ac\tZFIN_AC\t$ac\t\tXREF\n";
            }

        }
    }
    close (ZEBGENE);
}


if($organism eq 'briggsae'){
  my $in  = Bio::SeqIO->new(-file => $briggsae_peptides, '-format' =>'fasta');
 BRIGGSAE: while(my $seq = $in->next_seq){
    
    #ID CBG11531 desc CBP02734 (cb25.fpc2454.en7794a/cb25.fpc2454.tw352/cb25.fpc2454.gc383) W = 100.00; B = 99.7 
    #print STDERR "ID ".$seq->id." desc ".$seq->desc."\n";
    my $id = $seq->id;
    my @values = split /\s+/, $seq->desc;
    my $display_id = $values[0];
    my $syns = $values[1];
    $syns =~ s/[\(\)]//g;
    my @secs = split /\//, $syns;
    my $syn = join(';',@secs);
    print OUT $display_id."\tSPTR\t".$display_id."\tBRIGGSAE_HYBRID\t".$id."\t".$syn."\tXREF\n";
  }
}

close(OUT);
print STDERR "The output has been written there: $out\n";


sub parse_sp_file{
  my ($file) = @_;
  #print STDERR "opening file ".$file."\n";
  open(FH, $file) or die("couldn't open file ".$file." $!");
  my $counter = 0;
  my %swiss;
  my %ac;
  my %id;
  my $id = 0;
  my $ac = 0;
  
 SWISS: while(<FH>){
    chomp;
    if($_ =~ /^ID/ && $id == 0){
      my @values = split;
      if($id{$counter}){
	print STDERR "something odd going on ".$counter." element already defined ".$id{$counter}."\n";
	exit;
      }else{
	$id{$counter} = $values[1];
      }
      $id = 1;
    }
    if($_ =~ /^AC/ && $ac == 0){
      my @values = split;
      if($ac{$counter}){
	print STDERR "something odd going on ".$counter." element already defined ".$ac{$counter}."\n";
	exit;
      }else{
	$ac{$counter} = $values[1];
      }
      $ac = 1;
    }
    if(($_ =~ /^ID/) || ($_ =~ /^AC/) || ($_ =~ /^DR/) || ($_ =~ /^GN/)){
      if(!$swiss{$counter}){
	$swiss{$counter} = [];
      }
      push @{$swiss{$counter}}, $_;
    }
    if($_ =~ /http/){
      next SWISS;
    }
    if($_ =~ /ftp/){
      next SWISS;
    }
    if($_ =~ /\/\//){
      $counter++;
      $id = 0;
      $ac = 0;
    }
  }
  close(FH);
  return \%swiss, \%ac, \%id;
}


sub process_parsed_sp{
  my ($swiss, $ac, $id, $out) = @_;

  my %swiss = %{$swiss};
  my %ac = %{$ac};
  my %id = %{$id};
  foreach my $entry(keys(%swiss)){
    my @id_syns;
    my @ac_syns;
    my @embl_ids;
    my @protein_ids;
    my @mims;
    my @pdbs;
    my @genenames;
    my $id = $id{$entry};
    my $ac = $ac{$entry};
    my $tag;
    my $db;
    $ac =~ s/\;//;
    my @lines = @{$swiss{$entry}};
#    print STDERR "Entry ".$entry." has ".@lines." lines ".$ac.":".$id."\n";
  SP_ENTRY:foreach my $line(@lines){

      if(!$line){
	next SP_ENTRY;
      } 
      if($line =~ /^ID/){
	#print STDERR $line."\n";
	my @values = split /\s+/, $line;
#	print STDERR "have ".@values." values\n";
	if($values[1] ne $id){
	  push(@id_syns, $id);
	}
	if($values[2]){
	  $tag = $values[2];
	}
      }
      if($line =~ /^AC/){
	#AC   Q04456;
	my @values = split /\s+/, $line;
	shift @values;
	foreach my $v(@values){
	  $v =~ s/\;//;
	  if($v ne $ac){
	    push(@ac_syns, $v);
	  }
	}
    }
      
     
      if($line =~ /^GN/){
	  my @values = split /\s+/, $line;
	  my $v = $values[1];
	  $v =~ s/\.//;
	  if (($v =~ /EBIG/) || ($v =~ /AGCG/) || ($v =~ /ENSANGG/)) {
	      $tag = "Prediction-SPTREMBL";
	  }
      }
      

      if($line =~ /^DR/){
	if($line =~ /EMBL/){
	  #DR   EMBL; M96144; AAA28056.1; -. 
	  #print STDERR $line."\n";
	  my($foo, $bar, $embl, $protein) = split /\s+/, $line;
	  $embl =~ s/\;//;
	  $protein =~ s/\;//;
	  push(@embl_ids, $embl);
	  #print STDERR "have protein ".$protein."\n";
	  push(@protein_ids, $protein) unless($protein eq '-');
	}
	if($line =~ /MIM/){
	  #DR   MIM; 115200; -.
	  my @values = split /\s+/, $line;
	  my $mim = $values[2];
	  $mim =~ s/\;//;
	  push(@mims, $mim);
	}
	if($line =~ /PDB/){
	  #DR   PDB; 1ABQ; 15-OCT-95.
	  my @values = split /\s+/, $line;
	  my $mim = $values[2];
	  $mim =~ s/\;//;
	  push(@pdbs, $mim);
	}
      }
      if($organism eq 'drosophila'){
	if($line =~ /^GN/){
	  $line =~ s/GN//;
	  my @values = split /\s+OR\s+/, $line;
	  foreach my $g(@values){
	    $g =~ s/\.//;
	    $g =~ s/\s+//g;
	    push(@genenames, $g);
	  }
	}
      }
    }
    
    if(!$tag){
      die "you have no data tag can't decide if data is from swissprot or trembl for id ".$id." ac ".$ac." entry ".$entry."\n"
    }
    if($tag =~ /STANDARD/){
      $db = 'SWISSPROT';
    }elsif($tag =~ /PRELIMINARY/){
      $db = 'SPTREMBL';
    }elsif ($tag =~ /Prediction-SPTREMBL/){
	$db = 'prediction-SPTREMBL';
    }else {
      die("can't deal with tag ".$tag."\n");
    }
    my $ac_syns = join(';',@ac_syns);
    print $out $ac."\tSPTR\t".$ac."\t".$db."\t".$id."\t".$ac_syns."\tKNOWN\n";
    foreach my $embl_acc(@embl_ids){
      print $out $ac."\tSPTR\t".$embl_acc."\tEMBL\t".$embl_acc."\t\tXREF\n";
    }
    foreach my $protein(@protein_ids){
      print $out $ac."\tSPTR\t".$protein."\tprotein_id\t".$protein."\t\tXREF\n";
    }
    foreach my $mim(@mims){
      print $out $ac."\tSPTR\t".$mim."\tMIM\t".$mim."\t\tXREF\n";
    }
    foreach my $mim(@pdbs){
      print $out $ac."\tSPTR\t".$mim."\tPDB\t".$mim."\t\tXREF\n";
    }
    foreach my $mim(@genenames){
      print $out $ac."\tSPTR\t".$mim."\tFlyBase\t".$mim."\t\tXREF\n";
    }
  }
}

sub usage {
    
  print STDERR <<HELP

Usage: get_Xmapping.pl 
One of the element of the configuration file has not been properly loaded
for the organism $organism
Please fill in properly your configuration file

Here is your set up:
HELP
;

 foreach my $k (keys %check) {
	print STDERR "$k:\t$check{$k}\n";
    }



  exit();
}

