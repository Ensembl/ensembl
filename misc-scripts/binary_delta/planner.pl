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
    foreach my $v (@{ $thing{$db} }) {
	print <<EOT;
scp -c none -r ecs3:/mysqla/current/var/${db}_$v databases/
EOT
    }

    for (my $i = 0; $i < scalar @{ $thing{$db} } - 1; ++$i) {
	print <<EOT;
/usr/bin/time perl ./build.pl -c ./xdelta.osf -s databases -d deltas \\
   $db $thing{$db}[$i] $thing{$db}[$i + 1] 2>&1 | \\
   tee deltas/${db}_$thing{$db}[$i]_delta_$thing{$db}[$i + 1].txt
EOT
    }

    foreach my $v (@{ $thing{$db} }) {
	print <<EOT;
rm -rf databases/${db}_$v
EOT
    }
}
