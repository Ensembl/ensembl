package XrefParser::MGI_Desc_Parser;

use strict;
use warnings;
use Carp;
use File::Basename;

use base qw( XrefParser::BaseParser );

use Bio::EnsEMBL::DBSQL::DBAdaptor;


sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];
  my $syn_file = @{$files}[1];

  my $mgi_io = $self->get_filehandle($file);

  if ( !defined $mgi_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }

  my $xref_count =0;

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
    
      my $type = pop @part_desc; # last array element is the type.
      my $desc= join(" ",@part_desc);
      $acc_to_xref{$acc} = $self->add_xref($acc,"",$label,$desc,$source_id,$species_id,"MISC");
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
  my $mgi_syn_io = $self->get_filehandle($syn_file);

  if ( !defined $mgi_syn_io ) {
    print STDERR "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }


  my $syn_count = 0;
  while ( $_ = $mgi_syn_io->getline() ) {
    chomp;
    if($_ =~ /^ MGI:/){
      my ($junk, $acc, $chr, $pos, $symbol, @part_synonym) = split(/\s+/,$_);
      my $syn = join(" ",@part_synonym);

      if(defined($acc_to_xref{$acc})){
        $self->add_synonym($acc_to_xref{$acc}, $syn);
        $syn_count++;
      }
    }
  }
  $mgi_syn_io->close();
  print $syn_count." synonyms added\n" if($verbose);

  return 0; #successful
}
	



1;

