#! /usr/local/bin/perl
#GTF merging script, takes two arguments, input gtf file and prefix string

use Getopt::Std;
use Bio::EnsEMBL::Utils::GTF_Merge('gtf_merge');

my $def_prefix='foobar';
my $Usage = "Usage: $0 [ -p prefix (default=$def_prefix] inputfile(s) > outputfile\n";
my $opts = 'hp:';

getopts($opts) || die $Usage;

die $Usage if $opt_h;
# die "Need -p prefix!\n$Usage" unless $opt_p;

die $usage unless  @ARGV >= 1;

my $prefix = ($opt_p || $def_prefix);


my  $cmd = "sort -k1,1 -k7,7 -k4,4n " ;
   # i.e., first on fpccontigid, then on strand, then by start position

$pipe = $cmd . join(' ', @ARGV) . " | ";

open (IN,$pipe) || die "can't popen: $pipe";

&gtf_merge(\*IN,\*STDOUT,$prefix);

