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

sub overlap
{
	# Returns the length of the overlap of the two ranges
	# passed as argument.  A range is a two element array.

	my $first  = shift;
	my $second = shift;

	# Order them so that $first starts first.
	if ($first->[0] > $second->[0]) {
		($first, $second) = ($second, $first);
	}

	# No overlap
	return 0 if ($first->[1] < $second->[0]);

	# Partial overlap
	return ($first->[1] - $second->[0] + 1) if ($first->[1] < $second->[1]);

	# Full overlap
	return ($first->[1] - $first->[0] + 1);
}

my $pmatch_cmd	= '/nfs/disk5/ms2/bin/pmatch';
my $pmatch_opt	= '-T 14';
my ($unused_fh, $pmatch_out) = tempfile("pmatch_XXXXX",
	DIR => '/tmp', UNLINK => 0);

my $datadir	= '/acari/work4/mongin/final_build/release_mapping/Primary';
my $target	= $datadir . '/final.fa';
my $query	= $datadir . '/sptr_ano_gambiae_19_11_02_formated.fa';

# Set defaults
my %opts = (
	'b'	=> '0',
	'c'	=> $pmatch_cmd,
	'k'	=> '0',
	'p'	=> '2',
	'q'	=> $query,
	't'	=> $target );

if (!getopts('bc:kp:q:t:', \%opts)) {
	print STDERR <<EOT;
Usage: $0 [-b] [-c path] [-k] [-p num] [-q path] [-t path]

-b	Display header for each column.

-c path	Use the pmatch executable located at 'path' rather than at
	'$pmatch_cmd'.

-p num	Report the targets that are 'num' percent within the best
	matching target.  Default is 2%.

-k	Keep the pmatch output file.
	The name of the file will be printed on stderr.

-q path	Take query FastA file from 'path' rather than from
	'$query'.

-t path	Take target FastA file from 'path' rather than from
	'$target'.

EOT
	die;
}

# Override defaults
$pmatch_cmd	= $opts{'c'};
$query		= $opts{'q'};
$target		= $opts{'t'};

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

				my $o = overlap([$pair[0]{$c . 'START'},
						 $pair[0]{$c . 'END'}],
						[$pair[1]{$c . 'START'},
						 $pair[1]{$c . 'END'}]);
				$overlap += $o;
			}

			# Calculate the query and target identities
			$target->{$c . 'IDENT'} =
				100*($totlen - $overlap)/$target->{$c . 'LEN'};
		}
	}
}

if ($opts{'b'}) {	# Display header
	printf("%8s%8s%8s%8s%8s%8s\n",
		'QID', 'QLEN', 'QIDENT',
		'TID', 'TLEN', 'TIDENT');
}

my %maps;

foreach my $query (values(%hits)) {
	my $best;
	foreach my $target (
		sort { $b->{QIDENT} <=> $a->{QIDENT} } values %{ $query }) {

		$best = $target->{QIDENT} if (!defined($best));

		last if ($target->{QIDENT} < $best - $opts{'p'});

		#printf("%8s%8d%8.3f%8s%8d%8.3f\n",
			#$target->{QID}, $target->{QLEN},
			#$target->{QIDENT},
			#$target->{TID}, $target->{TLEN},
			#$target->{TIDENT});

		my $map = new Bio::EnsEMBL::Mapper('query', 'target');

		foreach my $hit (@{ $target->{HITS} }) {

			$map->add_map_coordinates(
				$target->{QID},
				$hit->{QSTART}, $hit->{QEND}, 1,
				$target->{TID},
				$hit->{TSTART}, $hit->{TEND});
		}

		push(@{ $maps{$target->{QID}} }, $map);
	}
}

# Example use of map.  Map [1,1000] on 'Q27250' through whatever the
# best map ('[0]') maps to:
print Dumper(
	$maps{Q27250}[0]->map_coordinates('Q27250', 1, 1000, 1, 'query')
);
