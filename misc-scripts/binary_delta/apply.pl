#!/usr/bin/perl -w

# $Id$
#
# apply.pl
#
# A program that uses a previously created set of binary delta
# files to produce a new revision of an ensembl database out of
# an older revision of the same database.  The delta files must
# have been built with the build.pl Perl program.
#
# See also apply.README
#
# Author: Andreas Kahari, <andreas.kahari@ebi.ac.uk>
#

use strict;
use warnings;

use File::Basename;
use File::Copy;

use Getopt::Std;

use aux qw(:default :apply);

my %opts;
my $xdelta_cmd	= $opts{'c'} = 'xdelta';
my $src_prefix	= $opts{'s'} = '.';
my $dst_prefix	= $opts{'d'} = '.';

if (!getopts('c:s:d:', \%opts)) {
    usage_apply(\%opts);
    die;
}

$xdelta_cmd = $opts{'c'};
$src_prefix = $opts{'s'};
$dst_prefix = $opts{'d'};

if ($#ARGV != 2) {
    usage_apply(\%opts);
    die;
}

my $db = $ARGV[0];
my $v1 = $ARGV[1]; my $v1_dir = sprintf "%s/%s_%s", $dst_prefix, $db, $v1;
my $v2 = $ARGV[2]; my $v2_dir = sprintf "%s/%s_%s", $dst_prefix, $db, $v2;

my $delta_dir = sprintf "%s/%s_%s_delta_%s", $src_prefix, $db, $v1, $v2;

die "$v1_dir: $!" if (! -d $v1_dir);
die "$delta_dir: $!" if (! -d $delta_dir);

while (-d $v2_dir) {
    $v2_dir = sprintf "%s.%04d", $v2_dir, int(rand(10000));
}

printf STDERR "Creating the directory '%s'\n", $v2_dir;
mkdir($v2_dir) or die $!;

my $v1_all_size = 0;
my $v2_all_size = 0;
my $delta_all_size = 0;

foreach my $info_file (glob($delta_dir . '/*.info')) {
    my $base_name  = basename($info_file);

    $base_name     =~ s/\.info$//;

    my $v1_file    = sprintf "%s/%s", $v1_dir, $base_name;
    my $v2_file    = sprintf "%s/%s", $v2_dir, $base_name;
    my $delta_file = sprintf "%s/%s", $delta_dir, $base_name;

    printf "Processing '%s'\n", $base_name;

    open(INFO, $info_file) or die $!;

    my $patch_command = <INFO>; chomp $patch_command;

    my $v1_line = <INFO>; chomp $v1_line;
    my ($v1_sum, $v1_size) = split /\s+/, $v1_line;

    my $v2_line = <INFO>; chomp $v2_line;
    my ($v2_sum, $v2_size) = split /\s+/, $v2_line;

    my $delta_line = <INFO>; chomp $delta_line;
    my ($delta_sum, $delta_size) = split /\s+/, $delta_line;

    close INFO;

    if ($v1_sum ne '(none)' && $v1_sum ne make_checksum($v1_file)) {
	print "\tChecksum mismatch for old file\n";
	print "\tCan not continue\n";
	die;
    } elsif ($v1_sum ne '(none)' && $v1_size != (stat $v1_file)[7]) {
	print "\tSize mismatch for old file\n";
	print "\tCan not continue\n";
	die;
    } else {
	print "\tChecksum and size ok for old file\n";
    }

    if ($delta_sum ne '(none)' && $delta_sum ne make_checksum($delta_file)) {
	print "\tChecksum mismatch for delta file\n";
	print "\tCan not continue\n";
	die;
    } elsif ($delta_sum ne '(none)' && $delta_size != (stat $delta_file)[7]) {
	print "\tSize mismatch for delta file\n";
	print "\tCan not continue\n";
	die;
    } else {
	print "\tChecksum and size ok for delta file\n";
    }

    if ($patch_command eq 'PATCH') {
	print "\tPatching file\n";
	system($xdelta_cmd, 'patch', $delta_file, $v1_file, $v2_file);
    } elsif ($patch_command eq 'COPY') {
	print "\tCopying old file\n";
	copy($v1_file, $v2_file);
    } elsif ($patch_command eq 'ZIP') {
	print "\tDecompressing compressed file\n";
	do_decompress($delta_file, $v2_file);
    } else {
	warn "\tStrange patch command: $patch_command\n";
    }

    if ($v2_sum ne '(none)' && $v2_sum ne make_checksum($v2_file)) {
	print "\tChecksum mismatch for new file\n";
	print "\tCan not continue\n";
	die;
    } elsif ($v2_sum ne '(none)' && $v2_size != (stat $v2_file)[7]) {
	print "\tSize mismatch for new file\n";
	print "\tCan not continue\n";
	die;
    } else {
	print "\tChecksum and size ok for new file\n";
    }
}
