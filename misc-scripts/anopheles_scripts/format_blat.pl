use strict;

my ($input) = @ARGV;

print STDERR "$input\n";

open (IN,"$input") || die;

my %score;
my %load;

while(<IN>) {
    chomp;

    my ($match,$a,$b,$c,$f,$g,$h,$k,$strand,$hid,$r,$hstart,$hend,$id,$m,$start,$end,$blockcount,$blocksize,$w,$tstart) = split;

    if ($hid =~ /^AG/) {

 #   print STDERR "HID: $hid\n";

	if ($match > $score{$hid}) {
#	print STDERR "$hid\n";
	    
	    $load{$hid} = $_;
	    $score{$hid} = $match;
	    
	}
    }
}


foreach my $k (keys %load) {
 #   print "$k\n";

    my ($match,$a,$b,$c,$f,$g,$h,$l,$strand,$hid,$r,$hstart,$hend,$id,$m,$start,$end,$blockcount,$blocksize,$w,$tstart) = split (/\t/,$load{$k});

    my $idt = $match * 100 / $r;

    my @blocks = split(/,/,$blocksize);
    my @targets = split(/,/,$tstart);
    my @queries = split(/,/,$w);

    my $count = 0;

    if ($strand eq "+") {
	$strand = 1;
    }
    else {
	$strand = -1;
    }

    if ($idt >= 70) {
	
	$start = $start + 1;
	$hstart = $hstart + 1;


#	foreach my $bl(@blocks) {
#	    my $start = $targets[$count] + 1;
#	    my $startq = $queries[$count] + 1;
#	    my $end = $start + $bl;
#	    my $endq = $startq + $bl;
	    print "\\N\t$id\t$start\t$end\t$strand\t$k\t50\t$idt\n";
#	    $count++;
#	}
#    }
#    print "SCORE $k\t$match\n"; 
	}
}











