#!/usr/bin/perl -w

# $Id$
#
# Don't use this program unless you know exactly what you and it
# are doing!
#
# Use like this:
#
#   planner.pl [options] <databases.txt >do_it_all.sh
#   sh -x ./do_it_all.sh
#
# [options] are
#   -a	Generate code that only runs apply.pl (for verification of
#	generated deltas).
#   -b	Generate code that only runs build.pl
#
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

use Getopt::Std;

my %opts;
if (!getopts('ab', \%opts)) {
    die "Options are -a and -b.  See head of script.\n";
}
if (!$opts{'a'} && !$opts{'b'}) {
    # If neither -a nor -b is specified, do both.
    $opts{'a'} = 1;
    $opts{'b'} = 1;
}

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

if ($opts{'a'}) {
    print "APPLY=YES\n";
} else {
    print "APPLY=NO\n";
}
if ($opts{'b'}) {
    print "BUILD=YES\n";
} else {
    print "BUILD=NO\n";
}

foreach my $db (keys %thing) {
    my @pair;
    print <<EOT;
# DATABASE: $db
EOT
    foreach my $v (@{ $thing{$db} }) {
	shift(@pair) if (scalar @pair == 2);
	push(@pair, $v);
	next if (scalar @pair != 2);

	my $v0 = $pair[0]; my $db0 = $db . '_' . $v0;
	my $v1 = $pair[1]; my $db1 = $db . '_' . $v1;
	my $d  = $db . '_' . $v0 . '_delta_' . $v1;
	print <<EOT;
# DELTA: $v0 -> $v1
if [ \\( \\( "x\$BUILD" = "xYES" -a ! -f deltas/${d}_build.txt \\) -o \\
         \\( "x\$APPLY" = "xYES" -a ! -f deltas/${d}_apply.txt \\) \\) -a \\
       ! -d databases/$db0 ]; then
  # Get older revision (needed for build and apply)
  scp -c none -r ecs3:/mysqla/current/var/$db0 databases/
fi

if [ "x\$BUILD" = "xYES" -a ! -f deltas/${d}_build.txt -a \\
     ! -d databases/$db1 ]; then
  # Get newer revision (only needed for build)
  scp -c none -r ecs3:/mysqla/current/var/$db1 databases/
fi

if [ "x\$BUILD" = "xYES" -a ! -f deltas/${d}_build.txt ]; then
  # Compute delta
  /usr/bin/time perl -w ./build.pl -c ./xdelta.osf -s databases -d deltas \\
    $db $v0 $v1 2>&1 | \\
    tee deltas/${d}_build.txt
fi

if [ "x\$APPLY" = "xYES" -a ! -f deltas/${d}_apply.txt ]; then
  # Apply the delta as a test
  /usr/bin/time perl -w ./apply.pl -c ./xdelta.osf -s deltas -d databases \\
    $db $v0 $v1 2>&1 | \\
     tee deltas/${d}_apply.txt
fi

# Remove older revision and new revision built by apply.pl
rm -rf databases/$db0
if [ "x\$APPLY" = "xYES" ]; then
  rm -rf databases/${db1}.????
fi

EOT
    }
    if (defined $pair[1]) {
	my $db1 = $db . '_' . $pair[1];
	print <<EOT;
if [ -d databases/$db1 ]; then
  # Remove newer (latest) revision
  rm -rf databases/$db1
fi

EOT
    }
}
