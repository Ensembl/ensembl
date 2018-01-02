#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

my $dir = shift;
my $submit = shift;

opendir(DIR, $dir);
my @allfiles = readdir DIR;
closedir DIR;

#xref_0_dna.3100129-307.out
#ExonerateGappedBest1_dna_307.map
#xref_0_dna.3100129-307.err

my %hash;
my %max;
my %fasta;
my $method;

$max{'pep'} = 0;
$max{'dna'} = 0;

foreach my $file(@allfiles){
  if($file eq '.' || $file eq '..'){
    next;
  }elsif(-d $file){
    next;
  }else{
    my $type;
    if($file =~ /dna.fasta/){
      $fasta{'dna'} = $file;
      next;
    }
    elsif($file =~ /protein.fasta/){
      $fasta{'pep'} = $file;
      next;
    }
    elsif($file =~ /dna/){
      $type = 'dna';
      if(!defined($method)){
	if($file =~ /(.*)_dna_\d+.map/){
	  $method= $1;
	}
      }
    }
    elsif($file =~ /pep/){
      $type = 'pep';
    }
    else{
      next;
    }

    if($file =~ /_(\d+).map/){
      $hash{$type}{'map'}{$1} = $1;
    }
    elsif($file =~ /-(\d+).err/){
      if($max{$type} < $1){
	$max{$type} = $1
      }
      $hash{$type}{'err'}{$1} = $1;
    }
    elsif($file =~ /-(\d+).out/){
      $hash{$type}{'out'}{$1} = $1;
    }
  }
}

print "Method $method\n";

if(defined($submit)){
  open(OUT,">SUBMIT")|| die "could not open SUBMIT output file\n";
}


foreach my $type ('pep', 'dna'){
  print "max ".$type."-->".$max{$type}."\n";
  for (my $i=1; $i <= $max{$type}; $i++){
    
    if(!defined($hash{$type}{'map'}{$i})){
      print $i,".map\n";
      if(defined($submit)){
	print OUT "/usr/local/ensembl/bin/exonerate-0.9.0 xref/xref_0_dna.fasta ";
	print OUT $dir."/homo_sapiens_dna.fasta --querychunkid ".$i." --querychunktotal ";
	print OUT $max{$type}." --showvulgar false --showalignment FALSE --ryo ";
	print OUT '"xref:%qi:%ti:%ei:%ql:%tl:%qab:%qae:%tab:%tae:%C:%s\n"';
	print OUT ' --model affine:local --subopt no --bestn 1 | grep "^xref" > ';
	print OUT $dir."/".$method."_".$type."_".$i.".map\n";
      }
    }

    if(!defined($hash{$type}{'err'}{$i})){
      print $i,".err\n";
    }
    if(!defined($hash{$type}{'out'}{$i})){
      print $i,".out\n";
    }
  }
}


