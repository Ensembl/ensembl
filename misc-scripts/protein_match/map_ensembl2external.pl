use strict;

=head1 map_ensembl2external.pl

=head2 Description

This script runs pmatch and postprocess it and put the results of pmatch together with with the external mapping to produce a file ready to go into TranscriptDBlink 

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

#Think to add $input and $output

print STDERR "Running script process_pmatch.pl\n";

#The output has to be directed in process_pmatch

my $process_pmatch = "perl process_pmatch.pl -ens ensembl.fas -sp sp_trembl.fas -refseq refseq.fas";
system($process_pmatch) == 0 or die "$0\Error running '$process_pmatch'";

print STDERR "Running script get_xrefs.pl\n";

my $get_xrefs "perl get_xrefs.pl -mapping ens_mapping.map -xrefs ext_mapping.extmap -dbmap dbmap.extmap -refseq ???? -output transcriptdblink.txt";
system($get_xrefs) == 0 or die "$0\Error running '$get_xrefs'";  

