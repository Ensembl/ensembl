#!/usr/local/bin/perl -w
# $Header$

# Create a file dbname\taccno\tdescription from SWISSPROT/SPTrEMBL/RefSeq flat file.
# Used by gene-descriptions.pl

# Usage: sp-descriptions.pl < swissprot-flatfile [RefSeq-flatfile] > sp-descriptions.dat

use strict;

my $swissprot_db_name = "SWISSPROT";
my $sptrembl_db_name = "SPTREMBL";
my $refseq_db_name = "RefSeq";

my $db;
my $acc; 
my $desc = "";
my %sp_desc;

while (defined (my $line = <>)) {

  if ($line =~ /^ID\s{3}.*$/o) {
    if ($line =~ /^ID\s{3}.*\s+PRELIMINARY;\s+.*$/o) {
      $db = $sptrembl_db_name;
    } elsif ($line =~ /^ID\s{3}.*\s+STANDARD;\s+.*$/o) { 
      $db = $swissprot_db_name;
    } else {
      chomp $line;
      warn "\nCould not recognize line : \"$line\"\nCheck the input file format.\n\n";
      exit 1;
    }
  }

  if ($line =~ /^LOCUS\s+\S+\s+.*$/o) {
    $db = $refseq_db_name;
  }
  
  if ($db eq $swissprot_db_name ||
      $db eq $sptrembl_db_name) {
    
    if ($line =~ /^AC\s{3}(\S+);.*$/o &&
       ! defined $acc)  {
      $acc = $1;
    }
    
    if ($line =~ /^DE\s{3}(\S.*\S)$/o) {
      $desc .= " " unless ($desc eq "");
      $desc .= uc $1;
    }
    
    if ($line =~ /^\/\/$/o) {
      $desc =~ s/\{.*\}//g; 
      $sp_desc{"$db\t$acc"} = $desc; 
      $db = undef;
      $acc = undef; 
      $desc = ""; 
    }
    
  } elsif ($db eq $refseq_db_name) {

    if ($line =~ /^DEFINITION\s+(\S.*\S)$/o) {
      $desc .= uc $1;
      while (defined ($line = <>)) {
	last if ($line =~ /^ACCESSION\s+(\S+)\s*.*$/o);
	if ($line =~ /^\s+(\S.*\S)$/o) {
	  $desc .= " ".uc $1;
	}
      }
    }

    if ($line =~ /^ACCESSION\s+(\S+)\s*.*$/o)  {
      $acc = $1;
    }

    if ($line =~ /^\/\/$/o) {
      $desc =~ s/\s*\[.*\]//g; 
      $sp_desc{"$db\t$acc"} = $desc; 
      $db = undef;
      $acc = undef; 
      $desc = ""; 
    }
  }

  
}

foreach my $db_acc (keys %sp_desc)  {
  my $desc = $sp_desc{$db_acc};
  print "$db_acc\t$desc\n";
}

exit 0;
