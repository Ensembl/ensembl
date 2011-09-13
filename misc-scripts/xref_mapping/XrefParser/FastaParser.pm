package XrefParser::FastaParser;

use strict;
use Carp;
use Bio::SeqIO;
use File::Basename;

use base qw( XrefParser::BaseParser );

# Fasta file format, e.g.
# >foo peptide sequence for the foo gene
# MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAG
# PTRTVDTKQAHELAKSYGIPFIETSAKTRQGVEDAFYTLVREIRQYRMKKLNSSDDGTQG
# CMGLPCVVM

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

  
  my $sio = Bio::SeqIO->new(-format=>'fasta' , -file=>$file );
  my %species_tax_id = %{$self->get_taxonomy_from_species_id($species_id)};
  
  my @xrefs;
  while( my $seq = $sio->next_seq ) {

    # Test species if available
    if( my $sp = $seq->species ){
      if( my $tax_id = $sp->ncbi_taxid ){
        next if (!defined $species_tax_id{$tax_id});
      }
    }

    # build the xref object and store it
    my $xref;
    $xref->{ACCESSION}     = $seq->display_name;
    $xref->{LABEL}         = $seq->display_name;
    $xref->{DESCRIPTION}   = $seq->description;
    $xref->{SEQUENCE}      = $seq->seq;
    $xref->{SOURCE_ID}     = $source_id;
    $xref->{SPECIES_ID}    = $species_id;
    $xref->{SEQUENCE_TYPE} = $seq->alphabet eq 'protein' ? 'peptide' : 'dna';
    $xref->{STATUS}        = 'experimental';
    if( my $v = $seq->version ){ $xref->{VERSION} = $v };
    push @xrefs, $xref;

  }

  print scalar(@xrefs) . " Fasta xrefs succesfully parsed\n" if($verbose);

  $self->upload_xref_object_graphs(\@xrefs);

  print scalar(@xrefs) . " Fasta xrefs succesfully loaded\n" if($verbose);


  return 0; #successful
}

1;
