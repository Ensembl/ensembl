#!/usr/local/bin/perl -w

# Author: Andreas Kahari <andreas.kahari@ebi.ac.uk>
#
# This is a  wrapper around Richard Durbin's
# pmatch code (fast protein matcher).

use strict;
use warnings;

use Data::Dumper;	# For debugging output

my $pmatch_cmd	= '/nfs/disk5/ms2/bin/pmatch';
my $pmatch_opt	= '-T 14';

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



open(PMATCH, "$pmatch_cmd $pmatch_opt $target $query |") or die($!);

my %hits;

# Populate the %hits hash.
while (defined(my $line = <PMATCH>)) {
	my ($length,
		$qid, $qstart, $qend, $qperc, $qlen,
		$tid, $tstart, $tend, $tperc, $tlen) = split(/\s+/, $line);

	$hits{$qid}{$tid} = [ ] if (!exists($hits{$qid}{$tid}));

	
	push(@{ $hits{$qid}{$tid} }, {
		LENGTH	=> $length,
		#not needed# QID	=> $qid,
		QSTART	=> $qstart,
		QEND	=> $qend,
		QPERC	=> $qperc,
		QLEN	=> $qlen,
		#not needed# TID	=> $tid,
		TSTART	=> $tstart,
		TEND	=> $tend,
		TPERC	=> $tperc,
		TLEN	=> $tlen });
}

close(PMATCH);

# Sort the %hits hash on QSTART, then on TSTART.
# Also calculate the lengths of any overlaps in the query or target.
foreach my $query (values(%hits)) {
	foreach my $target (values(%{ $query })) {
		@{ $target } =
			sort { $a->{QSTART} <=> $b->{QSTART}  ||
			       $a->{TSTART} <=> $b->{TSTART} } @{ $target };

		# Figure out how long the query and target overlaps are.
		# The first hit for any $qid<->$tid combination will never
		# have QOVERLAP or TOVERLAP keys (since it's not preceeded
		# by another hit).
		my @pair;
		foreach my $hit (@{ $target }) {
			shift(@pair) if (scalar(@pair) == 2);
			push(@pair, $hit);
			next if (scalar(@pair) != 2);

			$hit->{QOVERLAP} =
				overlap([$pair[0]{QSTART}, $pair[0]{QEND}],
					[$pair[1]{QSTART}, $pair[1]{QEND}]);
			$hit->{TOVERLAP} =
				overlap([$pair[0]{TSTART}, $pair[0]{TEND}],
					[$pair[1]{TSTART}, $pair[1]{TEND}]);
		}
	}
}

print Dumper(\%hits);	# Produce debugging output
