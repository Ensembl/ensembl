#!/usr/bin/env perl
# $Id$
# quick script to remap fasta files with COBP id's. Use as 
#
#   megamap.pl peptide.map gene.map  < input > output

use strict; 

my %map;
foreach my $f (@ARGV) {
    open (FILE,"<$f") || die "$f:$!";
    while (<FILE>) {
	chomp;
	$_ =~ /(\S+)\s+(\S+)/;
	$map{$1}=$2;
    }
    close (FILE) || die "$f:$!";
}

my $word_delim = '[-> \t:.\n;,_]';
my $ensid_regexp = '^([A-Z]{3})([PG])(0\d{10})$'; # '; # fool emacs

my $remapped=0;
while (<STDIN>) {
    if (/^>/) {
        my @words = split /($word_delim)/;
        foreach my $w ( @words ) {
            if (length ($w) > 1 && $map{$w} ) {
                s/$w/$map{$w}/g;
                $remapped++;
            }
        }
    }
    print;
}
warn "mapped $remapped identifiers\n";
