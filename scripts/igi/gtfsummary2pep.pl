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

# put all igi's in one big hash(ref) (a table of of igi to lists of
# native_ids):
my $all_igis = 
  Bio::EnsEMBL::Utils::igi_utils::read_igis_from_summary(\*SUMMARY); 

my $igis_of_natives = igi_of_native($all_igis);
my $all_peptides={};

foreach my $file ( @ARGV) {
    warn "reading $file ...\n";
    read_peptide_file($all_peptides, $file, $igis_of_natives);
    warn " done\n";
}

# my $blurp = blurp();
# print $blurp;
# print <<MOREBLURP
# ### IPI file {based on: see arguments line)
# MOREBLURP
#   ;

print_peptides($all_igis, $all_peptides);

# ### simple log message
# sub blurp {
#     my (@stuff)  = ("on", `date`, "by", `whoami`, "@",  `hostname`);
#     foreach (@stuff) {chomp};
#     my $s =  '### Produced by $Id$  ' . "\n" 
#       . "### run " . join(' ',@stuff) .  "\n"
#       . "### for EnsEMBL (http://www.ensembl.org)\n"
#       . "### argument(s): ". join(' ', @argv_copy) . "\n";
#     foreach (@argv_copy) {   $s .= "### ", `ls -l $_`; }
#     $s .= "###\n";
#     $s;
# }                                       # blurp

# invert the igis hash so we can lookup the igi given a native id:
sub igi_of_native  {
    my ($all_igis) = @_;
    my %igis_of_natives = undef;
    foreach my $igi (keys %$all_igis) {
        my @natives = @{$all_igis->{$igi}};
        foreach my $nat ( @natives ) {  # $nat includes the source prefix,
                                        # e.g. FGENH:S.C16000125
            $igis_of_natives{$nat} = $igi;
        }
    }
    return \%igis_of_natives;
}

# add peptides all into one big ugly hash
sub read_peptide_file {
    my ($peptides, $file, $igis_of_natives) = @_;
    my ($kept, $ignored) = (0,0);

    open(IN, $file) || die "$file: $!";

    local( $/ ) = "\n>";
  PEPTIDE:
    while( <IN> ) {
        s/>//g;
        my $n =(index $_, "\n");
        my $firstline =  substr($_, 0, $n);
        my $sequence = substr($_, $n+1);
        my $seqlen = ($sequence =~ tr /a-zA-Z//);

        my ($native_id, $rest) = ( $firstline =~ /(\S+)(.*)$/);
        unless ($native_id) {
            die "no native_id for $." 
        }

        my $igi = $igis_of_natives->{$native_id};
        
        if (! defined($igi)) {; # didn't make it into the IGI set
            $ignored++;
            next PEPTIDE ;
        }
        $kept++;
        push @{$peptides->{$igi}}, [$seqlen, $native_id, $rest, $sequence ];
    }
    warn "From $file, kept $kept , ignored $ignored peptides\n";
    close(IN);
    return;
}                                       # read_peptide_file


sub print_peptides { 
    my ($igis,$peptides) = @_; 
    my ($ignored, $kept) = (0,0);

  IGI:
    foreach my $igi (sort { substr($a,5) <=> substr($b,5) } keys %$igis) {
        # (this sorts on the numerical part of the "igi3_" identifier)
        my $peps= $peptides->{$igi};
        if (!defined $peps) {
            warn "no  peptides for $igi\n";
            next IGI;
        }
        my (@peps) = @$peps;
        
        my $max = -1;
        my $maxpep = undef;
        # get the longest (or ENS if equal;-) peptide:
        my @all_natives; # 
        foreach my $pep (@peps) {
            my $len = $ {$pep}[0];
            my $native = $ {$pep}[1];
            push(@all_natives, $native);
            if ($len > $max  || ( $len == $max && $native =~ /^ENS/ )  ) {
                $max=$len;
                $maxpep = $pep;
            }
        }
        my ($len, $native_id, $rest, $sequence) = @$maxpep;
        my @others = grep ( $_ ne $native_id, @all_natives);
        my $others = join(',', @others);
        my $newheader = ">$igi $native_id; length: $len; others: $others; original: $rest";
#        print "$newheader\n$sequence";
    }
}                                       # print_peptides
