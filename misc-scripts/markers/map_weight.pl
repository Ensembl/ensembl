use strict;

my $marker_feature_in = shift;
my $marker_hits_in = shift;

open MARKER_FEATURES, $marker_feature_in or die;
open MARKER_HITS, $marker_hits_in or die;

my %marker_hit_hash;
while(<MARKER_HITS>) {
	my ($id, $count) = split;
	$marker_hit_hash{$id} = $count;
}

close MARKER_HITS;

while(<MARKER_FEATURES>) {
	chomp;
	my @row = split;
	$row[-1] = $marker_hit_hash{$row[1]};
	print join("\t", @row) . "\n";
}

close MARKER_FEATURES;
