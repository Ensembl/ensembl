#!/usr/local/bin/perl
my @clones;
my $oldclone="";
while (my $line = <>) {
    if ($line =~ /Got two exons, \S+\[(\S+)\./) {
	$clone = $1;
	if ($clone ne $oldclone) {
	    push @clones,$clone;
	}
	$oldclone=$clone;
    }
    elsif ($line =~ /translation has stop codons. \(in clone (\S+)\)/) {
	$clone = $1;
	if ($clone ne $oldclone) {
	    push @clones,$clone;
	}
	$oldclone=$clone;
    }
    if ($clone ne $oldclone) {
	push @clones,$clone;
    }
    $oldclone=$clone;
}   
foreach $clone (@clones) {
    print "$clone\n";
}


