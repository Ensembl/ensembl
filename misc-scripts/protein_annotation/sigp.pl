#! /usr/local/bin/perl

#
# Perl script **** sigp.pl ***********
#
# Cared for 
#
# Copyright
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
coiled_coil
trans_membrane 
signal_peptide ########## this one !!!!!!!!!

				    
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

&GetOptions( 'files=s' => \$use_files,
	     'help' => \$help);

my $pfamseqfasta = shift;

exec('perldoc', $0) if $help or not defined $pfamseqfasta;



### THIS IS THE MINIMAL SCRIPT NEEDED TO SPLIT THE FILES FOR SIGNALP TO USE !!!

my($truncate);
my($file_num) = 0;

my($data_do) = 0;


open(_PFAMSEQ, $pfamseqfasta) or die "Fatal: Could not open $pfamseqfasta";
while(<_PFAMSEQ>) {


  if ($_ =~ /^\>/) {

    if (not ($entries % 1500)) {
      close(DATABASE);
  

      $data_do = 1;
      $file_num++;
     
      my($file) = "chunks/sigp_split." . $file_num;
      open(DATABASE, ">$file");
    }

    $entries++;
    $truncate = 1;

    my ($ac) = $_  =~ /^\>ENSP(\d+)/;
    
		
    print DATABASE $_;

  }  else {

    if($truncate) {
			my($sub_str);

			## Get only what we need !!
			$sub_str  = substr $_, 0, 50;	
	
			print DATABASE $sub_str . "\n";

			### Stops reading in rest of sequence as only need first 50 aa's
			$truncate = 0;
	

		      }
  } #/ end IF
}

                 

close(_PFAMSEQ);

close(DATABASE);
 
