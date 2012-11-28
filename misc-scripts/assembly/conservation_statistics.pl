#!/usr/bin/env perl

## Script used in conjunction with exon_conservation_check.pl. Requires
## Tie::IxHash for key ordering.

## RUN:
##     conservation_statistics.pl my_conservation.log

use strict;
use warnings;
use Tie::IxHash;

tie my %exon_states, 'Tie::IxHash', (
  '!!' => 'Protein Coding MisMatch',
  '%%' => 'Non-Coding MisMatch',
  '??' => 'Missed mapping',
  '¤¤' => 'Eval error'
);
tie my %transcript_states, 'Tie::IxHash', (
  '@@' => 'Protein Coding Protein MisMatch',
  '££' => 'Protein Coding cDNA MisMatch',
  '**' => 'Non-Coding MisMatch',
  'XX' => 'Missed mapping',
  '±±' => 'Eval error'
);

my $file = $ARGV[0];
die "Cannot find $file" if ! -f $file;
my %states = map { $_, 0 } (keys %exon_states, keys %transcript_states);
my %totals;
open my $fh, '<', $file or die "Cannot open file '$file': $!";
while(my $line = <$fh>) {
  chomp $line;
  if($line =~ /(^[!%?¤@£*X±]{2})/xms) {
    $states{$1}++;
    next;
  }
  if($line =~ /^Total ((?:exons|transcripts)) : (\d+)$/) {
    $totals{$1} = $2;
  }
}
close $fh;

my $format = "%s : %d (%.2f%%)\n";

print "Exon summary\n";
foreach my $key (keys %exon_states) {
  my $count = $states{$key};
  my $percentage = ($count / $totals{exons})*100;
  printf $format, $exon_states{$key}, $count, $percentage;
}

print "\nTranscript summary\n";
foreach my $key (keys %transcript_states) {
  my $count = $states{$key};
  my $percentage = ($count / $totals{transcripts})*100;
  printf $format, $transcript_states{$key}, $count, $percentage;
}