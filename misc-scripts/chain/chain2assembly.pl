use strict;
use warnings;

use Getopt::Long;



my ($seq_region_in_1, $seq_region_in_2, $chain_in);


GetOptions('seq_region_in_1=s'  => \$seq_region_in_1,
           'seq_region_in_2=s'  => \$seq_region_in_2,
           'chain_in=s'       => \$chain_in);

$seq_region_in_1 || usage();
$seq_region_in_2 || usage();
$chain_in || usage();

open(SEQ_REGION_IN_1, $seq_region_in_1)
  or die("Could not open file $seq_region_in_1");
open(SEQ_REGION_IN_2, $seq_region_in_2)
  or die("Could not open file $seq_region_in_2");
open(CHAIN_IN,      $chain_in)      or die("Could not open file $chain_in");

my %seq_region_map_1;
my %seq_region_map_2;

my @row;

while(<SEQ_REGION_IN_1>) {
  chomp;
  @row = split("\t", $_);
  $seq_region_map_1{$row[1]} = $row[0];
}

close SEQ_REGION_IN_1;

while(<SEQ_REGION_IN_2>) {
  chomp;
  @row = split("\t", $_);
  $seq_region_map_2{$row[1]} = $row[0];
}

close SEQ_REGION_IN_2;


my $is_header_line = 1;
my ($asm_sr_id, $cmp_sr_id, $asm_start, $cmp_start, $asm_strand, $cmp_strand,
    $ori, $match, $insert, $delete, $asm_begin, $asm_end, $cmp_begin,
    $cmp_end, $align_id, $asm_sr_len, $cmp_sr_len);


my %seen_coordinates;
my $line_num = 0;

my ($asm_sr_name,$cmp_sr_name);

LINE:
while(<CHAIN_IN>) {
  $line_num++;
  chomp;

  if($is_header_line) {
    my @header = split;

    $asm_sr_name = $header[2];
    $cmp_sr_name = $header[7];
    $asm_sr_name =~ s/chr//; #take the 'chr' off of chromosome names
    $cmp_sr_name =~ s/chr//;
    $asm_sr_id = $seq_region_map_1{$asm_sr_name};
    $cmp_sr_id = $seq_region_map_2{$cmp_sr_name};

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

    $asm_sr_len = $header[3];
    $cmp_sr_len = $header[8];

    $asm_strand = ($header[4] eq '+') ? 1 : -1;
    $cmp_strand = ($header[9] eq '+') ? 1 : -1;
    $ori        = $asm_strand * $cmp_strand;

    #start is (length-start) if we are on negative strand
    #start coordinates in file start at 0 not 1
    $asm_start = ($asm_strand == 1) ? $header[5]+1  : $asm_sr_len-$header[5];
    $cmp_start = ($cmp_strand == 1) ? $header[10]+1 : $cmp_sr_len-$header[10];

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

      if($seen_coordinates{"$asm_sr_id:$asm_begin"}) {
        my $first_line = $seen_coordinates{"$asm_sr_id:$asm_begin"};
         die("Duplicate coordinates [$asm_sr_name $asm_begin]. " .
             "At line=$line_num and line=$first_line\n");
      }

      if($seen_coordinates{"$cmp_sr_id:$cmp_begin"}) {
        my $first_line = $seen_coordinates{"$cmp_sr_id:$cmp_begin"};
         die("Duplicate coordinates [$cmp_sr_name $cmp_begin]. " .
             "At line=$line_num and line=$first_line\n");
      }

      if($seen_coordinates{"$asm_sr_id:$asm_end"}) {
        my $first_line = $seen_coordinates{"$asm_sr_id:$asm_end"};
         die("Duplicate coordinates [$asm_sr_name $asm_end]. " .
             "At line=$line_num and line=$first_line\n");
      }

      if($seen_coordinates{"$cmp_sr_id:$cmp_end"}) {
        my $first_line = $seen_coordinates{"$cmp_sr_id:$cmp_end"};
         die("Duplicate coordinates [$cmp_sr_name $cmp_end]. " .
             "At line=$line_num  and line=$first_line\n");
      }

      print join("\t",
                 $asm_sr_id, $cmp_sr_id,
                 $asm_begin, $asm_end,
                 $cmp_begin, $cmp_end,
                 $ori, "\n");

      $seen_coordinates{"$asm_sr_id:$asm_begin"} = $line_num;
      $seen_coordinates{"$cmp_sr_id:$cmp_begin"} = $line_num;
      $seen_coordinates{"$asm_sr_id:$asm_end"}   = $line_num;
      $seen_coordinates{"$cmp_sr_id:$cmp_end"}   = $line_num;

    } else {
      #empty line, next line is header line
      $is_header_line = 1;
    }
  }
}

close(CHAIN_IN);

sub usage {
  print "\nusage: perl chain2assembly.pl -seq_region_in_1 <filename> -seq_region_in_2 <filename> -chain_in <filename>\n";
  print "\n  seq_region_in must be a tab delimited text file dumped from an ensembl core".
        "\n  database (with revision > 19) that contains all of the seq_region names and ".
        "\n  identifiers used in the chain_in file.\n" .
        "\n  chain_in must be a chain format alignment file.\n" .
        "\n  assembly table data in tab delimited text format is printed to STDOUT\n";
  exit 0;
}
