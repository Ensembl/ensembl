#!/usr/bin/env perl
# $Id$
# quick script to remap fasta files with COBP id's. Use as
#
#   fasta-remap.pl peptide.map gene.map  < input > output

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
my $ensid_regexp = '^[A-Z]{3})([PG])(0\d{10}$'; # '; # fool emacs

my $remapped=0;
my @notmapped=();
while (<STDIN>) {
    if (/^>/) {
        my @words = split /($word_delim)/;
        foreach my $w ( @words ) {
            if (length ($w) > 1 && $w =~ /$ensid_regexp/  ) {
                if (defined $map{$w}) {
                    s/$w/$map{$w}/g;
                    $remapped++;
                } else {
                    push(@notmapped, $w);
                }
            }
        }
    }
    print;
}
warn "Mapped $remapped identifiers\n";
if (@notmapped) {
    warn "Not mapped:\n" , join("\n", @notmapped) ,
    "\n" , int(@notmapped), " were not mapped\n";
}
