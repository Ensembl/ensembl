#!/usr/bin/perl -w

# $Id$
#
# Don't use this program unless you know exactly what you and it
# are doing!
#
# Use like this:
#
#   planner.pl <databases.txt >do_it_all.sh
#
# Please take a look at the whole of "do_it_all.sh" before
# running it!  It assumes that the current directory contains
# certain directories and files, it removes files, and it
# contains no safty catches whatsoever.
#
# See "databases.txt.example" for an example of "databases.txt".
#
# Author: Andreas Kahari, <andreas.kahari@ebi.ac.uk>
#

use strict;
use warnings;

my %thing;

while (defined(my $line = <>)) {
    chomp $line;

    $line =~ /^(.*)_([1-9][0-9]*_[0-9]*[a-z]?)$/;
    if (!defined($1) || !defined($2) || $1 eq "" || $2 eq "") {
	printf STDERR "Problems: '%s'\n", $line;
	die;	# FIXME:  Doesn't catch everything...
    }
    my $db = $1;
    my $v  = $2;

    push @{ $thing{$db} }, $v;
}

foreach my $db (keys %thing) {
    my @pair;
    print <<EOT;
# DATABASE: $db
EOT
    foreach my $v (@{ $thing{$db} }) {
	shift(@pair) if (scalar @pair == 2);
	push(@pair, [ $db, $v ]);
	next if (scalar @pair != 2);

	my $p0 = $pair[0][0] . '_' . $pair[0][1];
	my $p1 = $pair[1][0] . '_' . $pair[1][1];
	my $d  = $p0 . '_delta_' . $pair[1][1];
	print <<EOT;
# DELTA: $pair[0][1] -> $pair[1][1]
if [ ! -d databases/$p0 -a ! -f databases/$p0.done ]; then
  scp -c none -r ecs3:/mysqla/current/var/$p0 databases/
fi
if [ ! -d databases/$p1 -a ! -f databases/$p1.done ]; then
  scp -c none -r ecs3:/mysqla/current/var/$p1 databases/
fi
if [ ! -f deltas/$d.txt ]; then
  /usr/bin/time ./build.pl -c ./xdelta.osf -s databases -d deltas \\
    $pair[0][0] $pair[0][1] $pair[1][1] 2>&1 | \\
    tee deltas/$d.txt
  rm -rf databases/$p0
  touch databases/$p0.done
fi
EOT
    }
    if (defined $pair[1]) {
	print <<EOT;
if [ -d databases/$pair[1][0]_$pair[1][1] ]; then
  rm -rf databases/$pair[1][0]_$pair[1][1]
  touch  databases/$pair[1][0]_$pair[1][1].done
fi
EOT
    }
}
