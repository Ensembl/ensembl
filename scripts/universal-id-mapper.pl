#!/usr/bin/env perl
# $Id$
# 
# quick script for id-mapping pretty much anything to anything. 
#
# The idea is that this script finds, in STDIN, *all* and *only*, the ID's
# given in the file.map argument, and replaces all of them on STDOUT.
#
# To make this process safer (i.e., if things on STDIN look like ID's but
# really should be left alone), make the regexp args (see Usage or -h) more
# restrictive. This will be especially needed for id's that are just
# numbers. 
#
# To see what happens, use the verbose option, which writes all
# unmatcheable mappable-ids-or-at-least-they-looked-like-it to stderr. 
#
# To make the things faster, use the -o option (or its opposite, -n)
#
# Since the script uses STDIN and STDOUT, you can set up a pipe of
# different mappers (operating on either the same or different id's). 

use strict; 
use Getopt::Std;

my $opts='hi:d:vo:n:';
use vars qw($opt_h $opt_i $opt_d $opt_v $opt_o $opt_n);

my $dflt_i = "[A-Z]{3,6}[PGET]\\d{11}";
my $dflt_d = "[-\\\\> \t:.;,/]+";

my ($prog) = ( $0 =~ m|/([^/]*)$| );

my $Usage=<<END_USAGE;

Usage:

  $prog [ options ]  *.map  < file-with-old-ids  > file-with-new-ids

  Options:
   -h       : this message
   -i REGEXP: use REGEXP to regcognize IDs to be mapped (default: $dflt_i)
   -d REGEXP: use REGEXP as word delimiter (default: $dflt_d)
   -v       : verbose mode (prints un-mapped IDs to stderr)
   -o REGEXP: only map lines matching REGEXP (amongst others for speed)
   -n REGEXP: map all lines apart from those matching REGEXP

The *.map files should contain pairs of whitespace delimited words. 

END_USAGE
#'; # pacify emacs

if ( !getopts($opts) || $opt_h ||  @ARGV==0 ) {
    die $Usage; 
}

die "combining options -o and -n doesn't make sense\n$Usage" 
    if ($opt_o && $opt_n);

my $word_delim = ($opt_d || $dflt_d);
my $not_allowed = '[a-zA-Z0-9_]';
die "word delim $word_delim may not contain any of $not_allowed (to avoid problems with '\\b' word boundaries"
 if $word_delim =~ /$not_allowed/;

my $id_regexp =  ($opt_i || $dflt_i);

# read all the map files:
my %map;
foreach my $f (@ARGV) {
    open (FILE,"<$f") || die "$f:$!";
    while (<FILE>) {
	chomp;
	my ($old, $new) = ( /(\S+)\s+(\S+)/ );
        die "$0: didn't find old:$old -> new:$new, file $f, line $."
          unless $old && $new;
	$map{$old}=$new;
    }
    close (FILE) || die "$f:$!";
}

my $remapped=0;
my $notmapped=0;
LINE:
while (<STDIN>) {

    if ($opt_n && /$opt_n/ || ($opt_o && !/$opt_o/) ) {
        print;
        next LINE;
    }
    chomp;
    if ( /$id_regexp/ ) { # ... but only if anything looks like an ID
        my @words = split /$word_delim/;
        foreach my $w ( @words ) {
            if ( $w =~ /$id_regexp/) {
                if (defined $map{$w}) {
                    s/\b$w\b/$map{$w}/g; # works on $_
                    $remapped++;
                } else {
                    warn $w, "\n" if $opt_v;
                    $notmapped++;
                }
            }
        }
    }
    print; 
    print "\n";
}
my $total=$remapped+$notmapped;
warn "Mapped $remapped out of $total identifiers\n";
exit ($remapped != $total);             # i.e, 0 if OK, 1 otherwise
