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

my $pmatch_cmd	= '/nfs/disk5/ms2/bin/pmatch';
my $pmatch_opt	= '-T 14';
my $pmatch_out	= '/tmp/pmatch_out.' . $$;	# Will be unlinked

my $datadir	= '/acari/work4/mongin/final_build/release_mapping/Primary';
my $target	= $datadir . '/final.fa';
my $query	= $datadir . '/sptr_ano_gambiae_19_11_02_formated.fa';

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

if (system("$pmatch_cmd $pmatch_opt $target $query >$pmatch_out") == -1) {
	# Failed to run pmatch
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
			QLEN	=> $qlen,
			TLEN 	=> $tlen,
			HITS	=> [ ]
		};
	}

	push(@{ $hits{$qid}{$tid}{HITS} }, {
		LENGTH	=> $length,
		QSTART	=> $qstart,
		QEND	=> $qend,
		TSTART	=> $tstart,
		TEND	=> $tend });
}

close(PMATCH);

unlink($pmatch_out);	# Get rid of pmatch output file

# Sort the %hits hash on QSTART, then on TSTART.
# Also calculate the lengths of any overlaps in the query or target.
foreach my $query (values(%hits)) {
	foreach my $target (values(%{ $query })) {
		@{ $target->{HITS} } =
			sort { $a->{QSTART} <=> $b->{QSTART}  ||
			       $a->{TSTART} <=> $b->{TSTART} }
					@{ $target->{HITS} };

		# Figure out how long the query and target overlaps are.
		# The first hit for any $qid<->$tid combination will never
		# have QOVERLAP or TOVERLAP keys (since it's not preceeded
		# by another hit).  This loop also calculates the total
		# length of the hits and the overlaps.

		my $qoverlap = 0;	# Total query overlap length
		my $toverlap = 0;	# Total target overlap length
		my $totlen = 0;		# Total hit length

		my @pair;
		foreach my $hit (@{ $target->{HITS} }) {
			$totlen += $hit->{LENGTH};

			shift(@pair) if (scalar(@pair) == 2);
			push(@pair, $hit);
			next if (scalar(@pair) != 2);

			$qoverlap +=
				overlap([$pair[0]{QSTART}, $pair[0]{QEND}],
					[$pair[1]{QSTART}, $pair[1]{QEND}]);

			$toverlap +=
				overlap([$pair[0]{TSTART}, $pair[0]{TEND}],
					[$pair[1]{TSTART}, $pair[1]{TEND}]);
		}

		# Calculate the query and target identities
		$target->{QIDENT} = 100*($totlen - $qoverlap)/$target->{QLEN};
		$target->{TIDENT} = 100*($totlen - $toverlap)/$target->{TLEN};
	}
}

print Dumper(\%hits);	# Produce debugging output
