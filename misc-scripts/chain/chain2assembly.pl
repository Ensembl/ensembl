use strict;
use warnings;

use Getopt::Long;



my ($seq_region_in, $chain_in);


GetOptions('seq_region_in=s'  => \$seq_region_in,
           'chain_in=s'       => \$chain_in);

$seq_region_in || usage();
$chain_in || usage();

open(SEQ_REGION_IN, $seq_region_in) or die("Could not open file $seq_region_in");
open(CHAIN_IN,      $chain_in)      or die("Could not open file $chain_in");

my %seq_region_map;

my @row;
while(<SEQ_REGION_IN>) {
  chomp;
  @row = split("\t", $_);
  $seq_region_map{$row[1]} = $row[0];
}

close SEQ_REGION_IN;


my $is_header_line = 1;
my ($asm_sr_id, $cmp_sr_id, $asm_start, $cmp_start, $asm_strand, $cmp_strand, $ori);
my ($match, $insert, $delete, $asm_begin, $asm_end, $cmp_begin, $cmp_end, $align_id);

LINE:
while(<CHAIN_IN>) {
  chomp;

  if($is_header_line) {
    my @header = split;

    my $asm_sr_name = $header[2];
    my $cmp_sr_name = $header[7];
    $asm_sr_name =~ s/chr//; #take the 'chr' off of chromosome names
    $cmp_sr_name =~ s/chr//;
    $asm_sr_id = $seq_region_map{$asm_sr_name};
    $cmp_sr_id = $seq_region_map{$cmp_sr_name};

    if(!$asm_sr_id) {
      warn("Unknown asm_seq_region $asm_sr_name");
      $is_header_line = 0;
      next LINE;
    }
    if(!$cmp_sr_id) {
      warn("Unknown cmp_seq_region $cmp_sr_name");
      $is_header_line = 0;
      next LINE;
    }

    $asm_strand = ($header[4] eq '+') ? 1 : -1;
    $cmp_strand = ($header[9] eq '+') ? 1 : -1;
    $ori        = $asm_strand * $cmp_strand;

    #start is end if we are on negative strand
    #start coordinates in file start at 0 not 1
    $asm_start = ($asm_strand == 1) ? $header[5] +1 : $header[6];
    $cmp_start = ($cmp_strand == 1) ? $header[10]+1 : $header[11]; 

    $is_header_line = 0;

    $align_id = $header[12];
    
  } else {
    if($_) {
      next LINE if(/^\#/); #skip comments

      next LINE if(!$asm_sr_id || !$cmp_sr_id);

      ($match,$insert,$delete) = split;

      if($asm_strand == 1) {
        $asm_begin  = $asm_start;
        $asm_end    = $asm_start + $match - 1;
        $asm_start += $match;
        $asm_start += $insert if($insert);
      } else {
        $asm_begin  = $asm_start - $match + 1; 
        $asm_end    = $asm_start;
        $asm_start -= $match;
        $asm_start -= $insert if($insert);
      }
      if($cmp_strand == 1) {
        $cmp_begin  = $cmp_start;
        $cmp_end    = $cmp_start + $match - 1;
        $cmp_start += $match;
        $cmp_start += $delete if($delete);
      } else {
        $cmp_begin  = $cmp_start - $match + 1;
        $cmp_end    = $cmp_start;
        $cmp_start -= $match;
        $cmp_start -= $delete if($delete);
      }

      print join("\t", 
                 $asm_sr_id, $cmp_sr_id, 
                 $asm_begin, $asm_end, 
                 $cmp_begin, $cmp_end,
                 $ori, "\n"); 

    } else {
      #empty line, next line is header line
      $is_header_line = 1;
    }
  }
}

close(CHAIN_IN);

sub usage {
  print "\nusage: perl chain2assembly.pl -seq_region_in <filename> -chain_in <filename>\n";
  print "\n  seq_region_in must be a tab delimited text file dumped from an ensembl core".
        "\n  database (with revision > 19) that contains all of the seq_region names and ".
        "\n  identifiers used in the chain_in file.\n" .
        "\n  chain_in must be a chain format alignment file.\n" .
        "\n  assembly table data in tab delimited text format is printed to STDOUT\n";
  exit 0;
}
