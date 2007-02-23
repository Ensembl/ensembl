package XrefParser::AgilentParser;

use strict;
use File::Basename;

use base qw( XrefParser::BaseParser );

# OParser for FASTA-format probe mappings from Agilent
# >A_23_P253586
# CTGTCAGGATTCTAGAACTTCTAAAATTAAAGTTTGGGGAAATCAGTAGCTCTGATGAGA
# >A_23_P217507
# AGAAAGACGTTTTCCAACATGTAGAACTGCTTTTTAACTGGAGGAAAAATACTTCAGGAG

sub run {

  my ($self, $file, $source_id, $species_id) = @_;

  my @xrefs;

#  local $/ = "\n>";

  my $ag_io = $self->get_filehandle($file);

  if ( !defined $ag_io ) {
      print "Could not open $file\n";
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

  print scalar(@xrefs) . " Agilent xrefs succesfully parsed\n";

  XrefParser::BaseParser->upload_xref_object_graphs(\@xrefs);

  print "Done\n";
  return 0;
}

1;
