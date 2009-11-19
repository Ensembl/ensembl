#!/usr/local/ensembl/bin/perl
while(<>) {
    chomp;
    ($f,$c,$rid,$cstart,$cend,$fstart,$fend,$rstart,$rend,$ori,$type) = split;

    my $cdiff=$cend-$cstart;
    my $fdiff=$fend-$fstart;
    my $rdiff=$rend-$rstart;

    if (($cdiff != $fdiff) || ($fdiff != $rdiff)) {
	print STDERR "Error in line: $_\n got cdiff:$cdiff, fdiff:$fdiff, rdiff:$rdiff\n";
    }
}
