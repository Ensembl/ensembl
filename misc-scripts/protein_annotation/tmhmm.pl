#! /usr/local/bin/perl

#
# Perl script tmhmm.pl
#
# Cared for by Mhairi Marshall
#
# Copyright 
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

make_other_rdb - DESCRIPTION 

This script populates the pfamseq table in the pfamdev relation datbase.

=head1 SYNOPSIS

transmem.pl <pfamseq fasta file>

Options:

-help : Show this help message
-files : Use the ultra-fast (but mysql specific) file import mechanism to load data

#### WHAT IS THIS !!!



=head1 DESCRIPTION
    
This script populates the other_region table in the pfam relational database.
It should be run once at each major release (i.e. each time the underlying 
database changes).

Note: This script creates temporary files in the cwd

The following region types are currently supported:

low_complexity
coiled_coil
trans_membrane  ****** <<<--- THIS ONE --->>> ********
signal_peptide

				    
=head1 CONTACT
						    
Mail pfam@sanger.ac.uk with any queries


=cut
    

    
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

  
#open SAVEERR, ">&STDERR";
#open STDERR, "/dev/null";

my $tmhmmfasta = shift;
(my $pf_file) = $tmhmmfasta =~ /chunk\.(\d+)/;

## Where output is written to
open(_TMFILE, ">outputs/tmhmm.$pf_file") or die "canna open tmhhm file to write $!\n";
close(_TMFILE);

open(_TMFILE, ">>outputs/tmhmm.$pf_file") or die "canna open tmhhm file to write $!\n";



### *** Emmanuel you need to put in the path to decodeanhmm, options & model otherwise it wont work.

open (_TMHMM, "cat $tmhmmfasta | /usr/local/ensembl/bin/decodeanhmm -f /usr/local/ensembl/lib/TMHMM2.0.options -modelfile /usr/local/ensembl/lib/TMHMM2.0.model | ") or die "canna open the file $!\n";


#open STDERR, ">&SAVEERR";
my($id, %transmem );
my $last_membrane = 0;

while (<_TMHMM>) {


   if( $_ =~ /^\>(\S+)/) { 
   
     $id = $1;

   } elsif ($_ =~ /%pred/){
   		   
	my($junk, $values) = split(/:/, $_);
	my(@trans) = split(/,/, $values);

	foreach (@trans) {

	  $_ = substr($_, 1);

	  my($orien, $start, $end) = split(/ /, $_);

	   $orien = uc($orien);
	  if($orien =~ /M/i) {
	 
	  $transmem{$orien} = $start . " " . $end;
	  } else {

	  $transmem{$orien} = $start;
	  }

	 
	  if ($orien =~ /M/i) {
	    $last_membrane = 1;
	  } elsif ($last_membrane) {
	    my($start, $end) = split(/ /, $transmem{'M'});
	    #print _TMFILE "\\N\t$id\t$start\t$end\ttransmembrane\ttmhmm\t\t";

	    print _TMFILE "NULL\t$id\t$start\t$end\t18\t0\t0\ttransmembrane\t0\t0\t0\n";

	     #if($transmem{'I'} < $transmem{'O'}) {
	      # print _TMFILE  " I-O\n";
	     #} else {
	      # print  _TMFILE " O-I\n";
	     #}
	   
	    $last_membrane =0;
	  }

	}	
   }
      
}

#open SAVEERR, ">&STDERR";
 # open STDERR, "/dev/null";

close(_TMHMM);


close(_TMFILE);

    






