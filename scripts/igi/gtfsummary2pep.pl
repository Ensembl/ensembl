#!/usr/local/bin/perl

# $Id$ 

# This script takes a gtf summary file and peptide files, and produces an
# ipi file. It is derived from gtfsummary2gtf.pl, and could possibly do
# with some refactoring.

# Written by Philip lijnzaad@ebi.ac.uk
# Copyright EnsEMBL http://www.ensembl.org

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Utils::igi_utils;

my $usage = "$0 options gtf.summary  peptidefiles > peptidefile\n";
my $help;
my @argv_copy = @ARGV;

&GetOptions( 'h|help' => \$help
	     );
die $usage if $help;

(my $summary = shift) || die $usage;

open(SUMMARY, "< $summary") || die "Summary file $summary: $!";

my $all_igis = Bio::EnsEMBL::Utils::igi_utils::read_igis_from_summary(\*SUMMARY); # all igi's in one big hash
# (a hash ref of igi to native_ids)

my $blurp = blurp();
print $blurp;
print <<MOREBLURP
### IPI file {based on summary/containing peptides} from @ARGV.
MOREBLURP
  ;

foreach $peptides ( @ARGV) {
    open(SOURCE, "<$peptides") || die "$peptides:$!";
    read_and_print(\*SOURCE, $all_igis);
    close(SOURCE);
}

### simple log message
sub blurp {
    my (@stuff)  = ("on", `date`, "by", `whoami`, "@",  `hostname`);
    foreach (@stuff) {chomp};
    my $s =  '### Produced by $Id$  ' . "\n" 
      . "### run " . join(' ',@stuff) .  "\n"
      . "### for EnsEMBL (http://www.ensembl.org)\n"
      . "### argument(s): ". join(' ', @argv_copy) . "\n";
    foreach (@argv_copy) {   $s .= "### ", `ls -l $_`; }
    $s .= "###\n";
    $s;
}                                       # blurp

sub read_and_print {
    my ($IN, $all_igis) = @_;

    local( $/ ) = "\n>";
    while( $_=<> ) {
        s/>//g;
        my $n =(index $_, "\n");
        my $firstline =  substr($_, 0, $n);
        my $rest = substr($_, $n+1);

        my ($native_id, $rest) = ( />([^ \t])+(.*)$/);
        die "no native_id for $." unless $native_id;

        print "first line: '$firstline'\n";
        print "rest: '$rest'\n";
    }
}                                       # read_and_print
