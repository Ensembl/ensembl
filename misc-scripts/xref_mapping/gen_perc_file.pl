
use strict;
use warnings;
  
use Getopt::Long;
use Cwd;
use XrefMapper::db;
  
use vars qw(@INC);
  

my $acc_file; # xref 2 acc list.
my $outfile;  # where to dump output to.
my $map_dir;  # directory map files are in.

GetOptions ('acc_file=s'              => \$acc_file,
	    'out_file=s'              => \$outfile,
	    'map_dir=s'               => \$map_dir);

open(FILE,"<".$acc_file) 
  || die "Could not open $acc_file";

<FILE>;

my %xref2acc;
while(<FILE>){
  chomp;
  my @arr = split;
  $xref2acc{$arr[0]} = $arr[1];
}

close FILE;


open(TEMP,">".$outfile) || die "COuld not open $outfile for writing\n";

my $count=0;
my %seen;
foreach my $file (glob("$map_dir/*.map")) {
    
  open(FILE, $file);
  while (<FILE>) {
#    $total_lines++;
    chomp();
    my ($label, $query_id, $target_id, $identity, $query_length, $target_length,
	$query_start, $query_end, $target_start, $target_end, $cigar_line, $score)
      = split(/:/, $_);
    
    if(defined($xref2acc{$query_id})){
      
      
      my $query_identity = int (100 * $identity / $query_length);
      my $target_identity = int (100 * $identity / $target_length);
      
      if(!defined($seen{$xref2acc{$query_id}})){
	print TEMP $query_identity."\t".$target_identity."\t".$xref2acc{$query_id}."\n";
      }
      $seen{$xref2acc{$query_id}} = 1;
    }
    else{
      $count++;
    }
  }
}


if($count){
  print STDERR "$count do not have accessions ??\n";
}
$count=0;
foreach my $acc (values %xref2acc){
  if(!defined($seen{$acc})){
    print TEMP "0\t0\t$acc\n";
    $count++;
  }
}
print STDERR "$count do not have any hits at all\n";
