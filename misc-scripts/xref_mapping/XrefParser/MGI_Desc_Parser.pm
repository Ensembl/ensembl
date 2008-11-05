package XrefParser::MGI_Desc_Parser;

use strict;
use File::Basename;

use base qw( XrefParser::BaseParser );

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


#my $dbi2;

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print STDERR "\nUsage: MGI_Desc_Parser.pm file <source_id> <species_id>\n\n";
    exit(1);
  }

  run(@ARGV);
}

sub run {

  my $self = shift if (defined(caller(1)));

  my $source_id = shift;
  my $species_id = shift;
  my $files       = shift;
  my $release_file   = shift;
  my $verbose       = shift;

  my $file = @{$files}[0];
  my $syn_file = @{$files}[1];

  if(!defined($source_id)){
    $source_id = XrefParser::BaseParser->get_source_id_for_filename($file);
  }
  if(!defined($species_id)){
    $species_id = XrefParser::BaseParser->get_species_id_for_filename($file);
  }



  my $mgi_io = $self->get_filehandle($file);

  if ( !defined $mgi_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }

  my $xref_count =0;
  my $syn_count =0;

  my %acc_to_xref;
#MGI Marker associations to Sequence (GenBank or RefSeq) information (tab-delimited)
#MGI Marker Accession ID	Marker Symbol	Status	Marker Type	Marker Name	cM Position	Chromosome 	GenBank Accession IDs
#(space-delimited)	Unigene ID
#(if any)	RefSeq ID
#(if any)  
  while ( $_ = $mgi_io->getline() ) {
    chomp;
    if($_ =~ /^ MGI:/){
      my ($junk, $acc, $chr, $pos, $label, $status, @part_desc) = split(/\s+/,$_);
    
      my $desc= join(" ",@part_desc);
      $acc_to_xref{$acc} = $self->add_xref($acc,"",$label,$desc,$source_id,$species_id);
      if($verbose and $desc eq ""){
	print "$acc has no description\n";
      }
      $xref_count++;
    }
  }
  
  $mgi_io->close();
  
  print $xref_count." MGI Description Xrefs added\n" if($verbose);

  #
  # Now process the synonyms
  #
  my $mgi_io = $self->get_filehandle($syn_file);

  if ( !defined $mgi_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }


  my $syn_count = 0;
  while ( $_ = $mgi_io->getline() ) {
    chomp;
    if($_ =~ /^ MGI:/){
     
      my ($junk, $acc, $chr, $pos, $symbol, @part_synonym) = split(/\s+/,$_);
      my $syn = join(" ",@part_synonym);
    
      if(defined($acc_to_xref{$acc})){
        $self->add_synonym($acc_to_xref{$acc}, $syn);
        $syn_count++;
      }
#      Lots of withdrawn entrys.
#      else{
#        print "Could not find xref for $acc to add synonym $syn\n" if($verbose);
#      } 
    }
  }
  $mgi_io->close();
  print $syn_count." synonyms added\n" if($verbose);

  return 0; #successful
}
	



1;

