#! /usr/local/bin/perl

#
# Perl script **** sigp.pl ***********
#
# Cared for by Mhairi Marshall (mm1@sanger.ac.uk)
#
# Copyright Mhairi Marshall
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

sigp_readsplit.pl - DESCRIPTION 

This script reads in part of the pfamseq data and prints the sig_p data to a file .

=head1 SYNOPSIS



Options:

-help : Show this help message
-files : Use the ultra-fast (but mysql specific) file import mechanism to load data


=head1 DESCRIPTION
    
This script populates the other_region table in the pfam relational database.


Note: This script creates temporary files in the cwd

The following region types are currently supported:

low_complexity
coiled_coil
trans_membrane 
signal_peptide ########## this one !!!!!!!!!

				    
=head1 CONTACT
						    
Mail pfam@sanger.ac.uk with any queries


=cut
    

## Note: This code reads in a file that contains part of the pfamseq data. Then runs
## the sig_p program and prints to another file the sigp data.


# Let the code begin...


#use lib '/nfs/disk45/klh/pfam/scripts/pfamrdb';


### USE MINE !!!
#use lib '/work1/birney/mongin/src/modules/pfamrdb';

#use lib '/nfs/disk65/mm1/pfam/scripts/Modules';

#use lib '/nfs/disk92/PerlSource/Bioperl/Releases/bioperl-0.05/';

#use lib '/nfs/disk100/pubseq/Pfam/scripts/Modules';
#use lib '/nfs/disk100/pubseq/Pfam/scripts/pfamrdb/';
    
use strict;
use Getopt::Long;
use FileHandle;

#use Bio::Pfam;

#use UpdateRDB;


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


### Which file is to be used 
my $file_num = shift;

exec('perldoc', $0) if $help or not defined $file_num;



my($data_do) = 0;


### THIS FILE HAS A LOT OF WRITING TO IT - 
 my($tmp_file) = "sigp_tmp_dat" . $file_num;



#my($read_file) = "sigp_split" . $file_num;
my($read_file) = "chunks/".$file_num;

open(_PFAMSEQ, "<$read_file") or die "Fatal: Could not open the file $read_file !!";

my($file) = "outputs/sigp_dat." . $file_num;
open(DATABASE, ">$file");
close(DATABASE);

my $entries = 1;
my @store_sigp_vals;
while(<_PFAMSEQ>) {



  if ($_ =~ /^\>/) {

    ## every 1000 entries print sigp output to file
	if (not ($entries % 1000)) {

		open(DATABASE, ">>$file");
		print "print data to file: $file !! entries: $entries\n";
		foreach  (@store_sigp_vals) {

		  print DATABASE $_;

		}
		close(DATABASE);

		@store_sigp_vals = ();

	}

	$entries++;
     
    open(TMP, ">$tmp_file");
    print TMP $_ ;
	

  }  else {
   
	
    print TMP $_;
    close(TMP);
    my $sigp_val = find_sigp($tmp_file);
    push @store_sigp_vals, $sigp_val;
		    
  }
 
}                 

close(_PFAMSEQ);
 
 open(DATABASE, ">>$file") or die "canna open the file \n";


### Print last of the sigp data to file
foreach  (@store_sigp_vals) {
  print DATABASE $_;

}

close(DATABASE);


## Find the sigp regions and print then to a file

sub find_sigp {

  my($tmp_file) = @_;
  
  open(_SIGP, "/usr/local/ensembl/bin/signalp -t euk $tmp_file | ") or die "Could not open the file \n";
  
  my($name);
  
  my (%line_store, $line_counter);
  $line_counter = 1;
  
  my $sigp_data;
  while(<_SIGP>) {
    
      #print STDERR "$_\n";
    
    $name = "ENSP".$1 if( $_ =~ /^>(\S+)/);
     
      
  #$name = $1 if( $_ =~ /^>(\S+)/);
  
    print STDERR "$name\n";
        
    $line_store{$line_counter} = $_ if ($_ =~ /max./i);
    if ($_ =~ /Most likely/i) {
      
      
      my($numbers) = m/(\d+)/;
      
      ### SEE if is a valid sigp region  
      if ( ( $line_store{1} =~ /YES/i  ) || ($line_store{2} =~ /YES/i ) || ($line_store{3} =~ /YES/i ) ) {
	#$sigp_data = "\\N\t$name\t1\t$numbers\tsig_p\tsignalp\t \t \t \n";
	$sigp_data = "NULL\t$name\t1\t$numbers\t19\t0\t0\tsignalp\t0\t0\t0\n";
      }
    } #/ end if Actual output
    
    if ($line_counter =~ /3/) {
      $line_counter = 1;
    } else {
      $line_counter++;
  }

  } #/ end WHILE
  close(_SIGP);

  return $sigp_data;

} #/ end find signal peptide 

