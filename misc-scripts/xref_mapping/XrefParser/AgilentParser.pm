package XrefParser::AgilentParser;

use strict;
use File::Basename;
use Carp;
use base qw( XrefParser::BaseParser );

# OParser for FASTA-format probe mappings from Agilent
# >A_23_P253586
# CTGTCAGGATTCTAGAACTTCTAAAATTAAAGTTTGGGGAAATCAGTAGCTCTGATGAGA
# >A_23_P217507
# AGAAAGACGTTTTCCAACATGTAGAACTGCTTTTTAACTGGAGGAAAAATACTTCAGGAG

sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $release_file = $ref_arg->{rel_file};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) or (!defined $release_file)){
    croak "Need to pass source_id, species_id, files and rel_file as pairs";
  }
  $verbose |=0;

  my $file = @{$files}[0];

  my @xrefs;

  my $ag_io = $self->get_filehandle($file);

  if ( !defined $ag_io ) {
      print STDERR "Could not open $file\n";
      return 1;
  }

  my $probe;
  while ( $_ = $ag_io->getline() ) {

    chomp;

    my $xref;

    # strip ^M at end of line
    $_ =~ s/\015//g;

    if(/^>(.+)/){
      $probe = $1;
    }
    else{
      my $sequence = $_;

      $sequence =~ s/\n//g;

      # build the xref object and store it
      $xref->{ACCESSION}     = $probe;
      $xref->{LABEL}         = $probe;
      $xref->{SEQUENCE}      = $sequence;
      $xref->{SOURCE_ID}     = $source_id;
      $xref->{SPECIES_ID}    = $species_id;
      $xref->{SEQUENCE_TYPE} = 'dna';
      $xref->{STATUS}        = 'experimental';
      
      push @xrefs, $xref;
    }
  }

  $ag_io->close();


  XrefParser::BaseParser->upload_xref_object_graphs(\@xrefs);

  print scalar(@xrefs) . " Agilent xrefs succesfully parsed\n" if($verbose);

  return 0;
}

1;
