use strict;

=head1 get_all_external

=head2 Description

This script will call all of the other script involved in the mapping of cross references.

=head2 Options

-input: Full pathname for the directory where the input files are stored.
-output: Full pathname for the directory where the output files should be stored

=head2 Filenames

The filenames of the input files are hardcoded on purpose. The input files should have the same names than the ones hardcoded in this script.

=head2 Contact

mongin@ebi.ac.uk
birney@ebi.ac.uk

=cut

use Getopt::Long;

my ($input,$output);

&GetOptions(
	    'input:s'=>\$input,
	    'output:s'=>\$output
	    );


##########################################
#First should run the script to get DBMAP#
##########################################

#Think to add $input and $output

print STDERR "Getting external DBs (EMBL and MIM) from Swiss-Prot\n";
my $get_embl_mim = "perl get_embl_mim_mapping.pl -sp sp_trembl.swiss -dbmap dbmap.extmap -output embl_mim.extmap";
system($get_embl_mim) == 0 or die "$0\Error running '$get_embl_mim'";

print STDERR "Getting external DBs (MIM and LOCUS) from Refseq\n";
my $get_mim_locus = "perl get_embl_mim_mapping.pl -refseq refseq.genbank -dbmap dbmap.extmap -output mim_locus.extmap";
system($get_mim_locus) == 0 or die "$0\Error running '$get_mim_locus'";

print STDERR "Getting HUGOs using HUGO mapping with SP and Refseq\n";
my $get_hugo = "perl get_hugo_mapping.pl -nomeid nomeids.txt -ens1 ens1.txt -ens2 ens2.txt -dbmap dbmap.extmap -output hugo.extmap";
system($get_hugo) == 0 or die "$0\Error running '$get_hugo'";






