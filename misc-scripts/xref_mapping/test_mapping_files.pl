use strict;
use warnings;

my $dir = shift;

opendir(DIR, $dir);
my @allfiles = readdir DIR;
closedir DIR;

#xref_0_dna.3100129-307.out
#ExonerateGappedBest1_dna_307.map
#xref_0_dna.3100129-307.err

my %hash;
my %max;

foreach my $file(@allfiles){
  if($file eq '.' || $file eq '..'){
    next;
  }elsif(-d $file){
    next;
  }else{
    my $type;
    if($file =~ /dna/){
      $type = 'dna';
    }
    elsif($file =~ /pep/){
      $type = 'pep';
    }
    elsif($file =~ /fasta/){
      next;
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


foreach my $type ('pep', 'dna'){
  print "max ".$type."-->".$max{$type}."\n";
  for (my $i=1; $i <= $max{$type}; $i++){
    if(!defined($hash{$type}{'map'}{$i})){
      print $i,".map\n";
    }
    if(!defined($hash{$type}{'err'}{$i})){
      print $i,".err\n";
    }
    if(!defined($hash{$type}{'out'}{$i})){
      print $i,".out\n";
    }
  }
}
