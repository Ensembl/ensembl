#!/usr/local/bin/perl -w

# $Id$
#
# Author: Andreas Kahari <andreas.kahari@ebi.ac.uk>
#
# This is a  wrapper around Richard Durbin's
# pmatch code (fast protein matcher).

use strict;
use warnings;

use Data::Dumper;	# For debugging output

use Getopt::Std;
use File::Temp qw(tempfile);

use Bio::EnsEMBL::Mapper;
use Bio::Range;

my $pmatch_cmd	= '/nfs/disk5/ms2/bin/pmatch';
my $pmatch_opt	= '-T 14';
my ($unused_fh, $pmatch_out) = tempfile("pmatch_XXXXX",
	DIR => '/tmp', UNLINK => 0);

my $datadir	= '/acari/work4/mongin/final_build/release_mapping/Primary';
my $target	= $datadir . '/final.fa';
my $query	= $datadir . '/sptr_ano_gambiae_19_11_02_formated.fa';
my $q_thr;
my $t_thr;
my $output;

# Set defaults
my %opts = (
	'c'	=> $pmatch_cmd,
	'k'	=> '0',
	'd'     => '0', 
	'p'	=> '2',
	'q'	=> $query,
	't'	=> $target
);

if (!getopts('c:kdp:q:t:', \%opts)) {
	print STDERR <<EOT;
Usage: $0 [-c path] [-k] [-d] [-p num] [-q path] [-t path]

-c path	Use the pmatch executable located at 'path' rather than at
]	'$pmatch_cmd'.

-p num	Report the targets that are 'num' percent within the best
	matching target.  Default is 2%.

-k	Keep the pmatch output file.
	The name of the file will be printed on stderr.

-q path	Take query FastA file from 'path' rather than from
	'$query'.

-t path	Take target FastA file from 'path' rather than from
	'$target'.

-d      Dump an output file which will be used for the known gene mapping

EOT
	die;
}

# Override defaults
$pmatch_cmd	= $opts{'c'};
$query		= $opts{'q'};
$target		= $opts{'t'};

# In that case the know gene mapping is used thus the option
# from the conf file ovveride the ones from the command line
if ($opts{'d'} == 1) {
	print STDERR "Warning: Any values given to the command line\n" .
		"option will be overideen by the values stored\n" .
		"in mapping_conf.pm\n";    

	BEGIN {
		my $script_dir = $0;
		$script_dir =~ s/(\S+\/)\S+/$1/;
		unshift (@INC, $script_dir);
		require "mapping_conf.pl";
	}

	my %conf =  %::mapping_conf;

	$query = $conf{'query'};
	$target = $conf{'pmatch_input_fa'};
	$output = $conf{'pmatch_out'};

	$t_thr = $conf{'target_idt'};
	$q_thr = $conf{'query_idt'};

	$pmatch_cmd = $conf{'pmatch'};
}


if (system("$pmatch_cmd $pmatch_opt $target $query >$pmatch_out") != 0) {
	# Failed to run pmatch command
	die($!);
}

open(PMATCH, $pmatch_out) or die($!);

my %hits;

# Populate the %hits hash.
while (defined(my $line = <PMATCH>)) {
	my ($length,
		$qid, $qstart, $qend, $qperc, $qlen,
		$tid, $tstart, $tend, $tperc, $tlen) = split(/\s+/, $line);

	if (!exists($hits{$qid}{$tid})) {
		$hits{$qid}{$tid} = {
			QID	=> $qid,
			TID	=> $tid,
			QLEN	=> $qlen,
			TLEN 	=> $tlen,
			HITS	=> [ ]
		};
	}

	push(@{ $hits{$qid}{$tid}{HITS} }, {
		QSTART	=> $qstart,
		QEND	=> $qend,
		TSTART	=> $tstart,
		TEND	=> $tend });
}

close(PMATCH);

if (!$opts{'k'}) {
	unlink($pmatch_out);	# Get rid of pmatch output file
} else {
	print(STDERR "$pmatch_out\n");
}


my $r1 = new Bio::Range();  # Outside loop to avoid unnecessary object creation.
my $r2 = new Bio::Range();

foreach my $query (values(%hits)) {
	foreach my $target (values(%{ $query })) {

		foreach my $c ('Q', 'T') {

			my $overlap = 0;	# Total query overlap length
			my $totlen = 0;		# Total hit length

			my @pair;
			foreach my $hit (
				sort { $a->{$c . 'START'} <=>
				       $b->{$c . 'START'} }
				@{ $target->{HITS} }) {

				$totlen += $hit->{$c . 'END'} -
					   $hit->{$c . 'START'} + 1;

				shift(@pair) if (scalar(@pair) == 2);
				push(@pair, $hit);
				next if (scalar(@pair) != 2);

				$r1->start($pair[0]{$c . 'START'});
				$r1->end($pair[0]{$c . 'END'});

				$r2->start($pair[1]{$c . 'START'});
				$r2->end($pair[1]{$c . 'END'});

				$overlap += $r1->intersection($r2);
			}

			# Calculate the query and target identities
			$target->{$c . 'IDENT'} =
				100*($totlen - $overlap)/$target->{$c . 'LEN'};
		}
	}
}

my %goodhits;
foreach my $query (values(%hits)) {
	my $best;
	my $priority = 0;
	foreach my $target (
		sort { $b->{QIDENT} <=> $a->{QIDENT} } values %{ $query }) {

		$best = $target->{QIDENT} if (!defined($best));

		last if ($target->{QIDENT} < $best - $opts{'p'});

		$goodhits{$target->{QID}}{$target->{TID}} = $target;

		foreach my $hit (@{ $target->{HITS} }) {

		    if ($opts{'d'} != 1) {
			printf("%s\t%d\t%s\t%d\t%d\t%d\t%d\n", 
				$target->{QID}, $priority, $target->{TID},
				$hit->{QSTART}, $hit->{QEND},
				$hit->{TSTART}, $hit->{TEND});
		    }
		}
		++$priority;

		# One mapping might look like this in the output:
		#
		# Q8WR21  0       10395   1       114     1       114
		# Q8WR21  0       10395   116     147     116     147
		# Q8WR21  1       10394   1       114     10      123
		# Q8WR21  1       10394   116     146     125     155
		#
		#
		# Columns are:
		#
		# 1. Query ID
		# 2. Priority (lower means better query sequence identity)
		# 3. Target ID
		# 4/5. Query start/stop
		# 6/7. Target start/stop

		#push(@{ $maps{$target->{QID}} }, $map);
	}
}

if ($opts{'d'} == 1) {
    open (OUT,">$output") || die "Can't open output file: $output\n";
	if (!defined($q_thr) || !defined($t_thr) || !defined($output)) {
		die "You need to define:\n" .
		    "query perc idt: $q_thr\n" .
		    "target perc idt: $t_thr\n" .
		    "Output file: $output\n";
	}

	foreach my $query (values(%goodhits)) {
		foreach my $target (values(%{ $query })) {

			if ($target->{'QIDENT'} >= $q_thr &&
			    $target->{'TIDENT'} >= $t_thr) {

				printf OUT "%s\t%s\t%.1f\t%.1f\n",
					$target->{'QID'},
					$target->{'TID'},
					$target->{'QIDENT'},
					$target->{'TIDENT'};
			}

		}
	}
}

# Example use of map.  Map [1,1000] on 'O61479' through whatever the
# best map ('[0]') maps to:
#print Dumper(
	#$maps{O61479}[0]->map_coordinates('O61479', 1, 1000, 1, 'query')
#);
