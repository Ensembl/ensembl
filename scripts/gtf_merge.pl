#/usr/local/bin/perl
#GTF merging script, takes two arguments, input gtf file and prefix string

use Bio::EnsEMBL::Utils::GTF_Merge('gtf_merge');

my $inputfile1= shift(@ARGV);
my $prefix = shift(@ARGV);
my $output_file="t/merged.gtf";
open (FILE,"<$inputfile1");
open (OUT,">$output_file");
&gtf_merge(\*FILE,\*OUT,'TEST_MERGE');

