use strict;
use warnings;

use Getopt::Long;

my ($seq_region_in_1, $seq_region_in_2, $axt_in);


GetOptions('seq_region_in_1=s'  => \$seq_region_in_1,
           'seq_region_in_2=s'  => \$seq_region_in_2);


$seq_region_in_1 || usage();
$seq_region_in_2 || usage();

open(SEQ_REGION_IN_1, $seq_region_in_1)
  or die("Could not open file $seq_region_in_1");
open(SEQ_REGION_IN_2, $seq_region_in_2)
  or die("Could not open file $seq_region_in_2");

my %seq_region_map_1;
my %seq_region_map_2;

my @row;

while(<SEQ_REGION_IN_1>) {
  chomp;
  @row = split("\t", $_);
  $seq_region_map_1{$row[1]} = [$row[0],$row[2]];
}

close SEQ_REGION_IN_1;

while(<SEQ_REGION_IN_2>) {
  chomp;
  @row = split("\t", $_);
  $seq_region_map_2{$row[1]} = [$row[0],$row[2]];
}

close SEQ_REGION_IN_2;

my ($asm_sr_name, $cmp_sr_name, $asm_sr_id, $cmp_sr_id, $asm_start,
    $cmp_start, $asm_strand, $cmp_strand,
    $ori, $match, $insert, $delete, $asm_begin, $asm_end, $cmp_begin,
    $cmp_end, $align_id, $asm_sr_len, $cmp_sr_len);

my $line_num = 0;

my %seen_coords;

while(<>) {
  next if((++$line_num-1) % 4); # read only every fourth line

  chomp;
  my @line = split;

  $asm_sr_name = $line[1];
  $cmp_sr_name = $line[4];
  $asm_sr_name =~ s/chr//;
  $cmp_sr_name =~ s/chr//;

  if(!$seq_region_map_1{$asm_sr_name}) {
    warn("Unknown seq_region $asm_sr_name");
    next;
  }

  ($asm_sr_id, $asm_sr_len) = @{$seq_region_map_1{$asm_sr_name}};

  if(!$seq_region_map_2{$cmp_sr_name}) {
    warn("Unknown seq_region $cmp_sr_name");
    next;
  }

  ($cmp_sr_id, $cmp_sr_len) = @{$seq_region_map_2{$cmp_sr_name}};

  $ori = ($line[7] eq '+') ? 1 : -1;

  $asm_begin = $line[2]+1;
  $asm_end   = $line[3]+1;
  $cmp_begin = ($ori == 1) ? $line[5]+1 : $cmp_sr_len - $line[6];
  $cmp_end   = ($ori == 1) ? $line[6]+1 : $cmp_sr_len - $line[5];

  if($seen_coords{"$asm_sr_id:$asm_begin"}) {
    my $first_line = $seen_coords{"$asm_sr_id:$asm_begin"};
    die("Duplicate coordinates [$asm_sr_name $asm_begin] " .
        "At line=$line_num  and line=$first_line");
  }

  if($seen_coords{"$cmp_sr_id:$cmp_begin"}) {
    my $first_line = $seen_coords{"$cmp_sr_id:$cmp_begin"};
    die("Duplicate coordinates [$cmp_sr_name $cmp_begin] " .
        "At line=$line_num and line=$first_line");
  }

  if($seen_coords{"$asm_sr_id:$asm_end"}) {
    my $first_line = $seen_coords{"$asm_sr_id:$asm_end"};
    die("Duplicate coordinates [$asm_sr_name $asm_end] " .
        "At line=$line_num  and line=$first_line");
  }

  if($seen_coords{"$cmp_sr_id:$cmp_end"}) {
    my $first_line = $seen_coords{"$cmp_sr_id:$cmp_end"};
    die("Duplicate coordinates [$cmp_sr_name $cmp_end]. " .
        "At line=$line_num and line=$first_line");
  }

  print join("\t",
             $asm_sr_id, $cmp_sr_id,
             $asm_begin, $asm_end,
             $cmp_begin, $cmp_end,
             $ori, "\n");

  $seen_coords{"$asm_sr_id:$asm_begin"} = $line_num;
  $seen_coords{"$cmp_sr_id:$cmp_begin"} = $line_num;
  $seen_coords{"$asm_sr_id:$asm_end"}   = $line_num;
  $seen_coords{"$cmp_sr_id:$cmp_end"}   = $line_num;
}


sub usage {
  print STDERR "gzcat axtfile.gz | perl axt2assembly.pl -seq_region_1 <file> -seq_region_2 <file>\n";
  exit(1);
}
