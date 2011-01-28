=head1 NAME

  format_blat.pl

=head1 SYNOPSIS
 
  format_blat.pl blat.out > blat_formated.txt

=head1 DESCRIPTION

  Formats a blat output into a format ready to go to mapfag. This is currently used to  produce the markers mapping

=head1 CONTACT

  dev@ensembl.org
  mongin@ebi.ac.uk

=cut

use strict;

my ($input) = @ARGV;

print STDERR "Blat output $input\n";

open (IN,"$input") || die;

my %score;
my %load;

while(<IN>) {
    chomp;

    my ($match,$a,$b,$c,$f,$g,$h,$k,$strand,$hid,$r,$hstart,$hend,$id,$m,$start,$end,$blockcount,$blocksize,$w,$tstart) = split;

    if ($hid =~ /^AG/) {
	if ($match > $score{$hid}) {
	    $load{$hid} = $_;
	    $score{$hid} = $match;
	}
    }
}


foreach my $k (keys %load) {

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
	    print "\\N\t$id\t$start\t$end\t$strand\t$k\t50\t$idt\n";
	}
}











