#!/usr/local/bin/perl

# $Id$ 

# This script takes a gtf summary file and a merged gft file, and from
# this, filters out the features that don't appear in the summary
# file. This is what's called the final gtf file.

# Written by Philip lijnzaad@ebi.ac.uk
# Copyright EnsEMBL http://www.ensembl.org

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Utils::igi_utils;

my $usage = "$0 options < merged-file.gtf gtf.summary  > gtfile\n";
my $help;
my @argv_copy = @ARGV;

&GetOptions( 'h|help' => \$help
	     );
die $usage if $help;

(my $summary = shift) || die $usage;

open(SUMMARY, "< $summary") || die "$summary: $!";

my $all_igis = Bio::EnsEMBL::Utils::igi_utils::read_igis_from_summary(\*SUMMARY); # all igi's in one big hash

my $blurp = blurp();
print $blurp;
print <<MOREBLURP
### Merged GTF file based on summary file $summary; see the header of that
### file to know how many sources had to have predicted this gene in order
### to 'make' it into this file
MOREBLURP
  ;

GTF_LINE:
while (<>) {
    # taken largely from Bio::EnsEMBL::Utils::GTF_handler:
    next GTF_LINE if /^                 #/;
      next GTF_LINE if /^\s*$/;
    chomp;
    
    my @fields = split "\t", $_;
    my ($seq_name, $source, $feature,
        $start,  $end,    $score,
        $strand, $phase,  $group_field)  = @fields;
    $feature = lc $feature;

    unless ($group_field) {
        warn("line $.: no group field: skipping : '$_'");
        next GTF_LINE ;
    }

    # Extract the extra information from the final field of the GTF line.
    my ($igi, $gene_name, $native_ids, $transcript_id, $exon_num, $exon_id) =
      Bio::EnsEMBL::Utils::igi_utils::parse_group_field($group_field);
    if ( $all_igis->{$igi} ) {          
        # this is a feature we want; rearrange things a bit so the order is
        # more consistent:
        my @fields = ($seq_name, $source, $feature,
                      $start,  $end,    $score,
                      $strand, $phase);

        my $native_id; 
        unless ($native_ids) {         
            die("line has no gene_id: '$_'\n");
        }
        
        if (int(@$native_ids) > 1 ) {   # this would be bizarre, but never
                                        # mind
            warn("Line with several gene_ids (taking first one): '$_'\n");
            next GTF_LINE;
        }
        $native_id =  ${$native_ids}[0];
        
        my $rest =  "igi_id \"$igi\"; gene_id \"$native_id\"; ";
        $rest .= "gene_name \"$gene_name\"; " if $gene_name;
        $rest .= "transcript_id \"$transcript_id\"; " if $transcript_id;
        $rest .= "transcript_id \"$transcript_id\"; " if $transcript_id;
        $rest .= "exon_id \"$exon_id\"; " if $exon_id;
        $rest .= "exon_number $exon_num; " if defined($exon_num);
        push @fields, $rest;
        print join("\t", @fields), "\n";
    }
}                                       # while(<>)


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

