#!/usr/bin/perl -w

use strict;
use warnings;
use diagnostics;

use Data::Dumper;
use Getopt::Std;

# $Id$
#
# Author:  Andreas Kahari, <andreas.kahari@ebi.ac.uk>
#
# This is a wrapper around exonerate that produces a comma
# separated file of alignment information, including cigar lines
# and separate percentage of identity for query and target
# sequences.  Run with -h or -? for usage information.
#
#----------------------------------------------------------------------
#
# Format of output:
#
# Col	Value				Type
# 1	Query sequence ID		string
# 2	Query percentage of identity	%g-type float <= 100
# 3	Alignment start in query	integer
# 4	Alignment end in query		integer
# 5	Target sequence ID		string
# 6	Target percentage of identity	%g-type float <= 100
# 7	Alignment start in target	integer
# 8	Alignment end in target		integer
# 9	Score				integer
# 10	Cigar line			string
#
#----------------------------------------------------------------------
#
# Variables:
#
# $q_fa:   Query FastA file (-q flag).
# $t_fa:   Target FastA file (-t flag).
# $e_cmd:  Which exonerate executable to pick up (-e flag).
# $e_ryo:  Exonerate output format specification.
# $[qts]_min: Minimum values on query, target, and score percentages
#             (min score is as a percentage of best score per query)
#             (-Q, -T, and -S flags).
#
# $fin:    File that is passed to the perl script rather than running
#          exonerate (-f flag).
# $finter: File where the intermediate exonerate output should be
#          stored (this file is the suitable for use with the -f flag
#          in later runs) (-i flag).
# $fout:   Where the output goeth (usually stdout) (-o flag).
#
#----------------------------------------------------------------------

my $q_fa='/acari/work4/mongin/anopheles_mai/mapping/Primary/peptides_new.fa';
my $t_fa='/acari/work4/mongin/anopheles_mai/mapping/Primary/total.fa';
my $e_cmd='/acari/work2/gs2/gs2/local/OSF1/bin/exonerate';
my $q_min=25;
my $t_min=25;
my $s_min=99;

sub usage
{
	my $padding = ' ' x length($0);

	print STDERR <<EOT;
Usage:	$0 [-h?]
	$0 [-v] [-e path] [-q path] [-t path] [-i path]
	$padding [-o path] [-Q val] [-T val] [-S val] [-- opts]	
	$0 [-v] [-f path] [-o path] [-Q val] [-T val] [-S val]
	$0 [-v] -d

-h, -?   Show usage information (this text).

-e path  Explicit path to exonerate executable.

-q path  Path to FastA file containing query sequences.

-t path  Path to FastA file containing target sequences.

-f path  Don't run exonerate, read from this properly formatted
	 file instead.

-i path  Save intermediate exonerate output to the specified
	 file.  This file is properly formatted for use with
	 the -f flag in later runs.

-o path  Save the final output to the specified file instead of dumping
         it on stdout.

-d       Read configuration from "mapping_config.pl".  The
	 configuration in that file overruns any other configuration
	 option specified on the command line.

-v       Be verbose.

-Q val   Don't consider alignments which have less than this
	 query percentage of identity.

-T val   Don't consider alignments which have less than this
	 target percentage of identity.

-S val   Don't consider alignments which have a score less than
	 this many percent of the best score for each query
	 sequence.

opts     Extra options to pass to exonerate, e.g. -M 512 -s 200
	 (must not affect format of output!).

Default values:
  -e $e_cmd
  -q $q_fa
  -t $t_fa
  -Q $q_min
  -T $t_min
  -S $s_min
EOT
}

my %opts;
if (!getopts('h?e:q:t:f:o:i:vQ:T:S:d', \%opts) || defined $opts{h}) {
    usage();
    exit(1);
}

my $fin    = $opts{f} if (defined $opts{f});
my $finter = $opts{i} if (defined $opts{i});
my $fout   = $opts{o} if (defined $opts{o});

$e_cmd = $opts{e} if (defined $opts{e});
$q_fa  = $opts{q} if (defined $opts{q});
$t_fa  = $opts{t} if (defined $opts{t});
$q_min = $opts{Q} if (defined $opts{Q});
$t_min = $opts{T} if (defined $opts{T});
$s_min = $opts{S} if (defined $opts{S});

if (defined $opts{d} && $opts{d} == 1) {
    our %mapping_conf;
    require 'mapping_conf.pl';

    $e_cmd = $mapping_conf{exonerate}  if (exists $mapping_conf{exonerate});
    $q_fa  = $mapping_conf{ensembl_predictions}      if (exists $mapping_conf{ensembl_predictions});
    $t_fa  = $mapping_conf{total_known_fa}
	if (exists $mapping_conf{total_known_fa});
    $fout  = $mapping_conf{mapping_out} if (exists $mapping_conf{mapping_out});
    $q_min = $mapping_conf{min_known_idt}  if (exists $mapping_conf{min_known_idt});
    $t_min = $mapping_conf{min_ensembl_idt} if (exists $mapping_conf{min_ensembl_idt});
}

if (defined($fin) && length($fin) != 0 &&
    defined($finter) && length($finter) != 0) {
    die "Can't specify -f and -i at the same time\n";
} elsif (! -R $q_fa) {
    die "Can't find or read file '$q_fa'\n";
} elsif (! -R $t_fa) {
    die "Can't find or read file '$t_fa'\n";
} elsif (! -X $e_cmd ) {
    die "Can't find or execute command '$e_cmd'\n";
}

# Set default options for exonerate using environment variables (may be
# overridden by command line options).
$ENV{EXONERATE_EXONERATE_FSMMEMORY}      = 512;
$ENV{EXONERATE_EXONERATE_HSPDROPOFF}     = 5;
$ENV{EXONERATE_EXONERATE_HSPTHRESHOLD}   = 30;
$ENV{EXONERATE_EXONERATE_MODEL}          = 'affine:local';
$ENV{EXONERATE_EXONERATE_PROTEINWORDLEN} = 6;
$ENV{EXONERATE_EXONERATE_WORDTHRESHOLD}  = 3;

my $e_ryo =
    '%qi\\\t%qal\\\t%ql\\\t%qab\\\t%qae\\\t' .
    '%ti\\\t%tal\\\t%tl\\\t%tab\\\t%tae\\\t' .
    '%p\\\t%s\\\t%C\\\n';

my $e_opt =
    "--showalignment no --showvulgar no -q $q_fa -t $t_fa --ryo $e_ryo " .
    "@ARGV";

my $cmd;
if (defined($fin) && length($fin) != 0) {
    $cmd = "cat $fin";
} elsif (defined($finter) && length($finter) != 0) {
    $cmd = "$e_cmd $e_opt | tee $finter";
} else {
    $cmd = "$e_cmd $e_opt";
}

if (!defined($fout) || length($fout) == 0) {
   $fout = '/dev/tty';
}

if (defined($opts{v}) && $opts{v} == 1) {
    warn "Will execute and parse the output from " .
	"the following command:\n$cmd\n";
    warn "Output goes to '$fout'\n";
}

open(IN,  "$cmd |") or die "Pipe failed: $!";
open(OUT, ">$fout") or die "Can not open '$fout' for writing: $!";

$s_min /= 100.0;
my %r;

while (defined(my $line = <IN>)) {
    # Perl script to calculate the percentage of identity
    # and reformat cigar lines.  Takes tab-delimited list in
    # specific format from exonerate as input on stdin and
    # writes comma-separated output to stdout.

    next if (($line =~ /^Message:/) || ($line =~ /^--/));
    chomp($line);

    # Pick out the individual fields (variable names correspond to
    # exonerate format parameters).
    my ($qi, $qal, $ql, $qab, $qae,
    	$ti, $tal, $tl, $tab, $tae,
    	$p,  $s,   $C) = split /\t/, $line;

    # Calculate the percent of identity (and skip to the next alignment
    # if it is not good enough).
    my $qp = $qal * $p / $ql; next if ($qp < $q_min);
    my $tp = $tal * $p / $tl; next if ($tp < $t_min);

    # Only store the incoming result if its score is larger than 95% (or
    # whatever is in $s_min) of the best score found so far for this
    # query ID.  If the score is larger than the best score, or if there
    # is no best score yet, then update the best score.

    if (!defined $r{$qi}{BS}) {
	# This is the first score.
	$r{$qi}{BS} = $s;
    }
    if ($s >= $s_min * $r{$qi}{BS}) {
	# This is a good enough score, as far as we know now.  The
	# alignment will not be output at the end if it is found to have
	# a too low score after all.

	# If this query-target pair has been seen before.  Replace
	# the previous pair only if the new score is better.  If this
	# query-target pair has not been seen before, store the data.

	if ((!exists $r{$qi}{ALGN}{$ti}) ||
	    (exists $r{$qi}{ALGN}{$ti} && $r{$qi}{ALGN}{$ti}{s} < $s)) {
	    $r{$qi}{ALGN}{$ti} = {
		qi => $qi, qp => $qp, qab => $qab, qae => $qae,
		ti => $ti, tp => $tp, tab => $tab, tae => $tae,
		s  => $s,  C  => $C };
	}

	if ($s > $r{$qi}{BS}) {
	    # This score was better than the best score so far for this
	    # query.  Do not try to delete bad alignments now (those that
	    # does not fulfill the new requirements), that would be too
	    # slow.
	    $r{$qi}{BS} = $s;
	}
    }
}

close(IN);

# Output

foreach my $q (values %r) {
    foreach my $t (values %{ $q->{ALGN} }) {
	next if ($t->{s} < $s_min * $q->{BS});

	# Reformat the cigar line.
	$t->{C} =~ s/([MDI]+) ([0-9]+) ?/$2$1/g;  # flip
	$t->{C} =~ s/([MDI])1([MDI])/$1$2/g;      # no lone 1
	$t->{C} =~ s/^1([MDI])/$1/;               # no lone 1 at start

	printf OUT "%s,%g,%d,%d,%s,%g,%d,%d,%d,%s\n",
	    $t->{qi}, $t->{qp}, $t->{qab}, $t->{qae},
	    $t->{ti}, $t->{tp}, $t->{tab}, $t->{tae},
	    $t->{s}, $t->{C};
    }
}

close(OUT);
