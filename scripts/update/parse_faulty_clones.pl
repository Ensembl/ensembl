#!/usr/local/bin/perl
my @clones;
my $oldclone="";
while (my $line = <>) {
    if ($line =~ /Found bug in clone (\S+)/) {
	$clone = $1;
	if ($clone ne $oldclone) {
	    push @clones,$clone;
	}
	$oldclone=$clone;
    }
}   
foreach $clone (@clones) {
    print "$clone\n";
}


