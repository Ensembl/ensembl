#!/usr/bin/env perl

package Dummy;

sub start { my $self = shift; return $self->{start}; }

my $array = [ bless({ start => 1, end => 3 }, 'Dummy'), bless({ start => 10, end => 3 }, 'Dummy'), bless({ start => 5, end => 3 }, 'Dummy'), bless({ start => 4, end => 3 }, 'Dummy')];
	      
my $s = [ sort { $a->start <=> $b->start } @{$array} ];
# my $s = [ map { $_->[1] } sort { $a->[0] <=> $b->[0] } map { [ $_->start, $_ ] } @{$array} ];

use Data::Dumper; print Dumper $s, "\n";

my $start = 10.1;
my $end = 15.5;

for my $i ($start .. $end+1) {
  print $i, " ";
}
print "\n";

my $results = [ bless({ start => 1, end => 3 }, 'Dummy'), bless({ start => 10, end => 3 }, 'Dummy') ];
push @{$results}, @{$array};
print Dumper $results;

my $duplicates = [ $array->[0], $array->[1], $array->[0], $array->[2], $array->[0], $array->[3] ];
print Dumper $duplicates, "\n";
my $u = uniq($duplicates);
print Dumper $u, "\n";

sub uniq {
  my $intervals = shift;

  use Tie::RefHash;
  tie my %seen, 'Tie::RefHash';
  
  return [ grep { ! $seen{ $_ }++ } @{$intervals} ];
}

1;
