#!/usr/local/bin/perl

# $Id$ 

# This script takes a gtf summary file and a merged gft file, and  from this,
# filters out the features that don't appear in the summary file.

# Written by Philip lijnzaad@ebi.ac.uk
# Copyright EnsEMBL http://www.ensembl.org

use strict;
use Getopt::Long;



my $usage = "$0 options < merged-file.gtf gtf.summary  > gtfile\n";
my $help;
my @argv_copy = @ARGV;

&GetOptions( 'h|help' => \$help
	     );
die $usage if $help;

(my $summary = shift) || die $usage;

open(SUMMARY, "< $summary") || die "$summary: $!";

my %all_igis = undef;
read_igis(\*SUMMARY);                   # all igi's in one big hash

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
    my ($igi, $gene_name, $native_id, $transcript_id, $exon_num, $exon_id) =
      parse_group_field($group_field);
    if ( $all_igis{$igi} ) {
        # rearrange things a bit so the order is more consistent:
        my @fields = ($seq_name, $source, $feature,
                      $start,  $end,    $score,
                      $strand, $phase);
        my $rest =  "igi_id \"$igi\"; gene_id \"$native_id\"; ";
        $rest .= "gene_name \"$gene_name\"; " if $gene_name;
        $rest .= "transcript_id \"$transcript_id\"; " if $transcript_id;
        $rest .= "transcript_id \"$transcript_id\"; " if $transcript_id;
        $rest .= "exon_id \"$exon_id\"; " if $exon_id;
        $rest .= "exon_number $exon_num; " if $exon_num;
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

### put all the igi's found in the summary into one big hash
sub  read_igis {
    my ($IN)  = @_; 

    SUMMARY_LINE:
    while (<$IN>) {
        next SUMMARY_LINE if /^#/;
          next SUMMARY_LINE if /^\s*$/;
        chomp;
    
        my @fields = split "\t", $_;
        my ($seq_name, $source, $feature,
            $start,  $end,    $score,
            $strand, $phase,  $group_field)  = @fields;
        $feature = lc $feature;
        
        unless ($group_field) {
            warn("no group field: skipping : '$_'\n");
            next SUMMARY_LINE ;
        }
        
        # Extract the extra information from the final field of the GTF line.
        my ($igi, $gene_name, $native_id, $transcript_id, $exon_num, $exon_id) =
          parse_group_field($group_field);
        
        $all_igis{$igi}++;                  # c'est tout
    }                                   # while <$IN>
}                                       # read_summary

# following will have to be factored out into a igi-utils.pm at some
# point, since also used by stats-from-merge-files.pl
sub parse_group_field {
    my( $group_field ) = @_;
    
    my ($igi, $gene_name, $native_id, $transcript_id, $exon_num, $exon_id);

    # Parse the group field
    foreach my $tag_val (split /;/, $group_field) {

        # Trim trailing and leading spaces
        $tag_val =~ s/^\s+|\s+$//g;

        my($tag, $value) = split /\s+/, $tag_val, 2;

        # Remove quotes from the value
        $value =~ s/^"|"$//g;
        $tag = lc $tag;

        if ($tag eq 'igi_id') {
            $igi = $value;
        }
        elsif ($tag eq 'gene_name') {
            $gene_name = $value;
        }
        elsif ($tag eq 'gene_id') {
            $native_id = $value;
        }
        elsif ($tag eq 'transcript_id') {
            $transcript_id = $value;
        }
        elsif ($tag eq 'exon_number') {
            $exon_num = $value;
        }
        elsif ($tag eq 'exon_id') {
            $exon_id = $value;
        }
        else {
            #warn "Ignoring group field element: '$tag_val'\n";
        }
    }
    return($igi, $gene_name, $native_id, $transcript_id, $exon_num, $exon_id);
}                                       # parse_group_field
