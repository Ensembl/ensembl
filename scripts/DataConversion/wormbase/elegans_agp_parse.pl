#!/usr/local/ensembl/bin/perl -w

# mouse agp -> ensembl sgp

$| = 1;

#print STDERR "have entered script\n";
my $agp = shift or die;
my $raw = shift or die;
my $chromosome = shift or die;

my @sgp;
my %chr_id;
my %ctg;

open RAW, "< $raw" or die "couldn't open file ".$raw." @!";
open AGP, "< $agp" or die;
open CHR, "< $chromosome" or die;

while (<RAW>) { 
  chomp;
  #print;
  #print "\n";
  my ($contig_id, $iid, $length) = split;
  #print "have ".$contig_id." ".$iid." ".$length."\n";
  #AL031185.2.1.4059
  #my ($id) = $contig_id =~ /(\w+\.\d+)\.\d+\.\d+/;   # strip off leading 'sc<date>_'
  $ctg{$contig_id} = [ $iid, $length ];
}

while(<CHR>){
  chomp;
  my ($id, $name) = split;
  #print "have ".$id." name ".$name."\n";
  $chr_id{$name} = $id; 
}

while (<AGP>) {
  
  chomp;
  #I	47490	107680	3	F	AC024796.1	1	60191	+
  #print;
  #print "\n";
  my ($chr, $chr_start, $chr_end, $gap,  $contig, $raw_start, $raw_end, $raw_ori) =
    (split)[0, 1, 2, 4, 5, 6, 7, 8];
  #next if $raw_start eq 'fragment';
  #next if $raw_start eq 'clone';
  #next if $raw_start eq 'contig';
  if($gap eq 'N'){
    next;
  }
  if($contig eq '.'){
    next;
  }
  #$chr =~ s/chr//;
  
  #next if $chr eq 'NA_contam';
  #print $chr." start ".$chr_start." end ".$chr_end." contig ".$contig." start ".$raw_start." end ".$raw_end." ori ".$raw_ori."\n";
  if ($raw_ori eq '+') {
    $raw_ori = 1;
  }
  elsif ($raw_ori eq '-') {
    $raw_ori = -1;
  }
  else {
    $raw_ori = 1;
    #print "$chr Contig  $contig  $chr_start \n";
    #print "Warning assumed default orientation for $contig\n";
  }
  my $chromosome_id = $chr_id{$chr};
  my $ctg_id    = $ctg{$contig}->[0];
  die "Don't have raw contig for $contig\n" unless $ctg_id;
  #print $contig."\n" if(!$ctg_id);
  #my $ctg_start = $ctg{$contig}->[1];
  my $ctg_end   = $ctg{$contig}->[1];
  if($raw_end > $ctg_end){
    next;
    #warn "contig ".$contig." length ".$ctg_end." is less than agp length ".$raw_end." won't work\n";
  }
  #if ($raw_start != $ctg_start || $raw_start != $ctg_start) {
  # raw coords from AGP assumed to be from 1 to length of contig
  #	print "Warning: my assumption about contig coords for $contig is wrong\n";
  #    }
  
  # $sgp[0] = qq{};
  $sgp[0] = $chromosome_id;   # fpc = contig id
  $sgp[1] = $chr_start;
  $sgp[2] = $chr_end;
  $sgp[3] = $chr;
  $sgp[4] = $chr_start;
  $sgp[5] = $chr_end;   # fpc start
  $sgp[6] = 1;   # fpc end
  $sgp[7] = $ctg_id;          # fpc ori
  $sgp[8] = $raw_start;   # raw start
  $sgp[9] = $raw_end;   # raw end
  $sgp[10] = $raw_ori;
  $sgp[11] = qq{elegans_90};
  
  print join("\t", @sgp), "\n";
}

close RAW;
close AGP;
close CHR;
