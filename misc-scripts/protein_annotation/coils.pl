#! /usr/local/bin/perl

#
# Perl script make_other_rdb
#
# Cared for by Mhairi Marshall <mm1@sanger.ac.uk>
#
# Copyright Kevin Howe & Mhairi Marshall
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

make_other_rdb - DESCRIPTION 

This script populates the pfamseq table in the pfamdev relation datbase.

=head1 SYNOPSIS

make_other_rdb.pl <pfamseq fasta file>

Options:

-help : Show this help message
-files : Use the ultra-fast (but mysql specific) file import mechanism to load data


=head1 DESCRIPTION
    
This script populates the other_region table in the pfam relational database.
It should be run once at each major release (i.e. each time the underlying 
database changes).

Note: This script creates temporary files in the cwd

The following region types are currently supported:

low_complexity
coiled_coil !!!
trans_membrane 
signal_peptide

				    
=head1 CONTACT
						    
Mail pfam@sanger.ac.uk with any queries


=cut
    

## Note: This code splits pfamseq into about 5 chunks
## before feeding each file to seg. This is because
## seg segmentation faults with large input files


# Let the code begin...

    
use strict;
use Getopt::Long;
use FileHandle;




my ( $help,
     $path,
     $use_files,
     $read_fh, 
     $write_fh,
     $fname,
     $total, 
     $rdb,
     $filecount,
     $entries,
     @pfamseqlist,
     @datfiles,
     @list);



### GIVE IT A FILENAME ####

&GetOptions( 'files=s' => \$use_files,
	     'help' => \$help);

my $pfamseqfasta = shift;

(my $pf_file) = $pfamseqfasta =~ /chunk\.(\d+)/;

exec('perldoc', $0) if $help or not defined $pfamseqfasta;




$@ and die "Fatal error: failed to empty pfamseq table before re-inserting [$@]";

			
my $coilsdir='/usr/local/ensembl/data/coils';
$ENV{'COILSDIR'}=$coilsdir;
open (_NCOILS, "/usr/local/ensembl/bin/ncoils -f < $pfamseqfasta  |" ) or die "Fatal: could not open $pfamseqfasta + NCOILS";


my($whole_sequence, $coil_name);

my($data_file) = "coils_upload.dat";

my($data_do) = 0;

### CLOSE STDERR
#open SAVEERR, ">&STDERR";
#open STDERR, "/work1/birney/mongin/src/ensembl-live/misc-scripts/protein_annotation";
open DATABASE, ">outputs/coils.$pf_file";

while(<_NCOILS>) {



	if ($_ =~ /^\>/) {
		&find_coils($coil_name, $whole_sequence ) if($whole_sequence =~ /x/);

		$coil_name = $1 if ($_ =~ /^>(\S+)/);

		$whole_sequence = "";



	} else {

		

		chop($_);
		$whole_sequence .= $_;

	}







} #/ end WHILE

## Add the last one

&find_coils($coil_name, $whole_sequence ) if($whole_sequence =~ /x/);
close(DATABASE);


### Re-open standard error
open STDERR, ">&SAVEERR";
exit(0);           



### FIND THE COIL LOCATION
#
#
#


sub find_coils {

	my( %x_position);

	my($seq_id) = shift;
	my($one_sequence) = shift;




	
	my($first) = 1;
	my($num) = 1;
	my($prev);



	while ($one_sequence =~ m/x/g) {

		my($x_pos);

		$x_pos =  pos $one_sequence;


		## FIRST ONE

		if($first) {
	
			$first = 0;
			$x_position{$num} = $x_pos;
			$prev = $x_pos;

		} else {

			### AMEND TO PREVIOUS
			if ( ($prev + 1) eq ($x_pos)) {
				$prev = $x_pos;
				if ($num eq 8) {
	
				}
				
			} else {

			### WHOLE NEW REGION

				$x_position{$num} .= "~" . $prev;
				$num++;
				$x_position{$num} = $x_pos;
				$prev = $x_pos;


			} #/ end IF

		} #/ end IF


	

	} #/ end WHILE	


	## ADD THE LAST ONE
	$x_position{$num} .= "~" . $prev;


	foreach my $x_arr (sort keys %x_position) {

		if ($x_position{$x_arr} =~ /~/) {

			my($from, $to) = split(/~/, $x_position{$x_arr});

#Print in the following format to allow an easy loading into Ensembl database
#FORMAT: translation\tseq_start\tseq_end\tanalysis\thstart\thend\thid\tscore\tevalue\tperc_id\n
			print DATABASE "NULL\t$seq_id\t$from\t$to\t14\t0\t0\tcoils\t0\t0\t0\n";


			#print DATABASE "\\N\t$seq_id\t$from\t$to\tcoiled_coil\tncoils\t \n";
		} #/ end IF


	
	} #/ end FOR
	


	

} #/ end SUB













