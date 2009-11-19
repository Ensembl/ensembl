#!/usr/local/ensembl/bin/perl

open(ERROR,">nonunique.txt");

while(<>) {
    chomp;
    ($f,$c,$rid,$cstart,$cend,$fstart,$fend,$rstart,$rend,$ori,$type) = split;

    if( defined $h{$rid} ) {
	print ERROR "$_\n";
	next;
    }
    print "$_\n";
    $h{$rid} = 1;
}
