#!/usr/local/bin/perl

while (my $line = <>) {
    if ($line =~ /Got gene with id (\S+),/) {
	$gene = $1;
	print "$gene\n";
    }
}
   



