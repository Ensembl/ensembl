#!/usr/local/bin/perl -w
use strict;
use Getopt::Long;
=head1 NAME

  prepare_proteome.pl

=head1 SYNOPSIS
 
  prepare_proteome.pl

=head1 DESCRIPTION
    
  This script is actually a copy of Val Curven script used for the Genebuilt, its function is the same

  prepare_proteome.pl prepares a fasta file of protein sequences from swissprot and refseq 
  input files (also in fasta format). This file is needed for pmatch comparisons and its 
  creation is the first part of the GeneBuild.

  The file has a description line consisting solely of the accession number after the leading >
  All U are replaced by X to prevent pmatch complaining.

  The final part of the script does a tiny pmatch test run to reveal any problems persisting in 
  the file that would prevent later full scale pmatches from running smoothly.

=head1 OPTIONS
  
  Options are to be set in GB_conf.pl
  The important ones for this script are:
     refseq      location of refseq file in fasta format
     sptr        location of swissprot file in fasta format
     pfasta      where to write the clean fasta file
     pmatch      location of the pmatch executable

     eg.
	    'refseq'      => '/work2/vac/GeneBuild/rf.fa',
	    'sptr'        => '/work2/vac/GeneBuild/sptr.fa',
	    'pfasta'      => '/work2/vac/GeneBuild/human_proteome.fa',
	    'pmatch'      => '/work2/vac/rd-utils/pmatch',
  
=cut


BEGIN {
  # oooh this is not nice
  my $script_dir = $0;
  $script_dir =~ s/(\S+\/)\S+/$1/;
  unshift (@INC, $script_dir);
  require "mapping_conf.pl";
}


#Global var
my %conf     = %::mapping_conf;

my $org_list = $conf{'organism_list'};
my $refseq   = $conf{'refseq_fa'};
my $sptr     = $conf{'sptr_fa'};
my $protfile = $conf{'total_known_fa'};
my $pmatch   = $conf{'pmatch'};
my $organism = $conf{'organism'};
my $refseq_pred = $conf{'refseq_pred_fa'};
my $sub_genes = $conf{'submitted_genes'};
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
    $check{'refseq_gnp'} = $conf{'refseq_gnp'};
    $check{'total_known_fa'} = $conf{'total_known_fa'};
    $check {'pmatch'} = $conf{'pmatch'};
    
    foreach my $k (keys %check) {
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }
    
}

if ($organism eq "mouse") {
    $check{'sptr_swiss'} = $conf{'sptr_swiss'};
    $check{'refseq_gnp'} = $conf{'refseq_gnp'};
    $check{'total_known_fa'} = $conf{'total_known_fa'};
    $check {'pmatch'} = $conf{'pmatch'};

    foreach my $k (keys %check) {
	print STDERR $check{$k}."\n";
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }

}

if ($organism eq "elegans") {
    $check{'sptr_swiss'} = $conf{'sptr_swiss'};
    $check{'total_known_fa'} = $conf{'total_known_fa'};
    $check {'pmatch'} = $conf{'pmatch'};

    foreach my $k (keys %check) {
	print STDERR $check{$k}."\n";
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }

}

if ($organism eq "anopheles") {
    $check{'sptr_swiss'} = $conf{'sptr_swiss'};
    $check{'submitted_genes'} = $conf{'submitted_genes'};
    $check{'total_known_fa'} = $conf{'total_known_fa'};
    $check {'pmatch'} = $conf{'pmatch'};

    foreach my $k (keys %check) {
	print STDERR $check{$k}."\n";
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }
}

if ($organism eq "drosophila") {
    $check{'sptr_swiss'} = $conf{'sptr_swiss'};
    $check{'refseq_gnp'} = $conf{'refseq_gnp'};
    $check{'total_known_fa'} = $conf{'total_known_fa'};
    $check {'pmatch'} = $conf{'pmatch'};
    
    foreach my $k (keys %check) {
	print STDERR $check{$k}."\n";
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }

}

if ($organism eq "rat") {
    $check{'sptr_swiss'} = $conf{'sptr_swiss'};
    $check{'refseq_gnp'} = $conf{'refseq_gnp'};
    $check{'total_known_fa'} = $conf{'total_known_fa'};
    $check {'pmatch'} = $conf{'pmatch'};

    foreach my $k (keys %check) {
	print STDERR $check{$k}."\n";
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }

}

if ($organism eq "zebrafish") {
    $check{'sptr_swiss'} = $conf{'sptr_swiss'};
    $check{'total_known_fa'} = $conf{'total_known_fa'};
    $check {'pmatch'} = $conf{'pmatch'};

    foreach my $k (keys %check) {
	print STDERR $check{$k}."\n";
	if ($check{$k} !~ /(\S+)/) {
	    usage();
	}
    }

}


#End of checks


print STDERR "Prepare proteome using following files: SPTR ($sptr)\nREFSEQ ($refseq)\nPROTFILE ($protfile)\n";

if (($organism eq "human") || ($organism eq "mouse") || ($organism eq "rat") || ($organism eq "drosophila")) {
    &parse_refseq;
}

if ($organism eq "anopheles") {
    &parse_sub;
}


if($organism eq 'briggsae'){
   &parse_briggsae_sptr;
   &parse_blessed;
   &test_protfile;
   exit(0);
}

&parse_sptr;
&test_protfile;


### END MAIN

sub parse_refseq {

  open (IN, "<$refseq") or die "Can't open $refseq\n";
  open (OUT, ">>$protfile") or die "Can't open $protfile\n";

  while(<IN>){
    # eg >gi|4501893|ref|NP_001094.1| actinin, alpha 2 [Homo sapiens]
    if(/^>/){
      if(/^>\w+\|\w+\|\w+\|(\S+)\.\d+\|/){
        print OUT ">$1\n";
      }
      else {
        print OUT $_;
      }
    }
    else {
      # sequence - sub U by X
      s/U/X/g;
      print OUT $_;
    }
  }
  close IN;
  close OUT;

}

sub parse_sptr {

  open (IN, "<$sptr") or die "Can't open $sptr\n";
  open (OUT, ">>$protfile") or die "Can't open $protfile\n";
  print STDERR "parsing sptr\n";
  while(<IN>){
    #>ACE1_CAEBR Q27459 AAB41269.1 Desc: Acetylcholinesterase 1 precursor (EC 3.1.1.7) (AChE 1).
    # eg >143G_HUMAN (Q9UN99) 14-3-3 protein gamma
    if(/^>\S+\s+\((\S+)\)/){
     

      if($1 eq 'P17013'){
	die("DYING: $sptr still contains P17013. \nThis will probably cause problems with pmatch.\nYou should REMOVE IT AND RERUN prepare_proteome!\n");
      }
	
      if($1 eq 'Q99784'){
	die("DYING: $sptr still contains Q99784. \nThis will probably cause problems with pmatch.\nYou should REMOVE IT AND RERUN prepare_proteome!\n");
      }
      print OUT ">$1\n";
    }
    else {
      print OUT $_;
    }
  }
  
  close IN;
  close OUT;

}

sub parse_briggsae_sptr {

  open (IN, "<$sptr") or die "Can't open $sptr $!";
  open (OUT, ">>$protfile") or die "Can't open $protfile $!";
  print STDERR "parsing sptr\n";
  while(<IN>){
    #>ACE1_CAEBR Q27459 AAB41269.1 Desc: Acetylcholinesterase 1 precursor (EC 3.1.1.7) (AChE 1).
    # eg >143G_HUMAN (Q9UN99) 14-3-3 protein gamma
    if(/^>\S+\s+(\S+)/){
     

      if($1 eq 'P17013'){
	die("DYING: $sptr still contains P17013. \nThis will probably cause problems with pmatch.\nYou should REMOVE IT AND RERUN prepare_proteome!\n");
      }
	
      if($1 eq 'Q99784'){
	die("DYING: $sptr still contains Q99784. \nThis will probably cause problems with pmatch.\nYou should REMOVE IT AND RERUN prepare_proteome!\n");
      }
      print OUT ">$1\n";
    }
    else {
      print OUT $_;
    }
  }
  
  close IN;
  close OUT;

}


sub parse_sub {
    #Simply add the file to the main file which is going to be used for the mapping
    open (IN, "$sub_genes") || die "Can't open $sub_genes\n";
    open (OUT, ">>$protfile") || die "Can't open $protfile\n";

    while(<IN>) {
	print OUT $_;
    }
    
    close IN;
    close OUT;

}

sub parse_blessed {

  open (IN, "<$refseq") or die "Can't open $refseq\n";
  open (OUT, ">>$protfile") or die "Can't open $protfile\n";

  while(<IN>){
    #>CBG00032 CBP00001 (cb25.fpc0002.en5152a) W = 89.6; B = 95.9
    if(/^>/){
      if(/^>\w+\s+(\w+)/){
	print OUT ">$1\n";
      }
      else {
	print OUT $_;
      }
    }
    else {
      # sequence - sub U by X
      s/U/X/g;
      print OUT $_;
    }
  }
  close IN;
  close OUT;

}

sub parse_refseq_pred {
    
  open (IN, "<$refseq_pred") or die "Can't open $refseq_pred\n";
  open (OUT, ">>$protfile") or die "Can't open $protfile\n";

  while(<IN>){
    # eg >gi|4501893|ref|NP_001094.1| actinin, alpha 2 [Homo sapiens]
    if(/^>/){
      if(/^>\w+\|\w+\|\w+\|(\S+)\.\d+\|/){
	print OUT ">$1\n";
      }
      else {
	print OUT $_;
      }
    }
    else {
      # sequence - sub U by X
      s/U/X/g;
      print OUT $_;
    }
  }
  close IN;
  close OUT;

}


sub parse_hybrid{
  open (IN, "<$refseq_pred") or die "Can't open $refseq_pred\n";
  open (OUT, ">>$protfile") or die "Can't open $protfile\n";

  while(<IN>){
    # eg >gi|4501893|ref|NP_001094.1| actinin, alpha 2 [Homo sapiens]
    if(/^>/){
      if(/^>\w+\|\w+\|\w+\|(\S+)\.\d+\|/){
	print OUT ">$1\n";
      }
      else {
	print OUT $_;
      }
    }
    else {
      # sequence - sub U by X
      s/U/X/g;
      print OUT $_;
    }
  }
  close IN;
  close OUT;

}


sub test_protfile {

  # set up a temporary file
  my $time = time;
  chomp ($time);
  my $tmpfile = "cup.$$.$time.fa";
  open (SEQ, ">$tmpfile") or die "can't open $tmpfile\n";
  print SEQ ">test_seq\n";
  print SEQ 'cctgggctgcctggggaagcacccagggccagggagtgtgaccctgcaggctccacacaggactgccagaggcacac';
  close SEQ;

  # do a pmatch test run
  print "starting pmatch test ... \n";
  open(PM, "$pmatch -D $protfile $tmpfile | ") or die "Can't run $pmatch\n";
  while(<PM>) {
    print $_;
  }
  close PM;

  # tidy up
  unlink $tmpfile;

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
