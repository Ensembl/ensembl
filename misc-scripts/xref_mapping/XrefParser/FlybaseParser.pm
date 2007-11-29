# $Id$

package XrefParser::FlybaseParser;

use strict;
use warnings;

use Carp;

use base qw( XrefParser::BaseParser );

# The object types we'd like to parse.
our %object_types = ( gene       => 1,
                      mRNA       => 1,
                      miRNA      => 1,
                      ncRNA      => 1,
                      protein    => 1,
                      pseudogene => 1,
                      rRNA       => 1,
                      snRNA      => 1,
                      snoRNA     => 1,
                      tRNA       => 1 );

# This is some statistics from the 5.4 file 'dmel-all-r5.4.gff.gz',
# looking at the FlyBase (2nd column) object types (3rd column) gene,
# mRNA, miRNA, ncRNA, protein, pseudogene, rRNA, snRNA, snoRNA, and tRNA
# only:
#
# cnt   Dbxref name             source_name in Xref database
# 57884 FlyBase_Annotation_IDs  (special, see below)
# 33441 GB_protein              protein_id
# 14324 GB                      EMBL
# 13951 flygrid                 FlyGrid
# 13265 FlyBase                 flybase_annotation_id
# 12768 dedb                    dedb
# 11745 UniProt/TrEMBL          Uniprot/SPTREMBL
# 10076 INTERPRO                Interpro
# 8077  orthologs               SKIPPED
# 2089  UniProt/Swiss-Prot      Uniprot/SWISSPROT
# 1596  bdgpinsituexpr          bdgpinsituexpr
# 1207  hybrigenics             SKIPPED
# 787   if                      SKIPPED
# 290   MITODROME               SKIPPED
# 153   TF                      SKIPPED
# 82    EPD                     SKIPPED
# 80    MIR                     SKIPPED
# 61    PDB                     SKIPPED
# 56    MEROPS                  SKIPPED
# 17    GCR                     SKIPPED
# 15    Rfam                    SKIPPED
# 12    NRL_3D                  SKIPPED
# 11    GO                      SKIPPED
#
# The Dbxref name 'FlyBase_Annotation_IDs' will be associated with the
# source_names FlyBaseCGID_{gene,transcript,translation} depending on
# the type of 'ID' of the line.
#
# Likewise, the source_names FlyBaseName_{gene,transcript,translation}
# will be associated with the 'Name' of each entry depending on the type
# of 'ID'.
#
# ... and the source_names flybase_{gene,transcript,translation}_id will
# be associated with the 'ID' of each entry depending on the type of
# 'ID'.

# This hash will translate the Dbxref names in the data file into source
# names known by the Xref system.
our %source_name_map = ( 'FlyBase'    => 'flybase_annotation_id',
                         'GB'         => 'EMBL',
                         'GB_protein' => 'protein_id',
                         'INTERPRO'   => 'Interpro',
                         'UniProt/Swiss-Prot' => 'Uniprot/SWISSPROT',
                         'UniProt/TrEMBL'     => 'Uniprot/SPTREMBL',
                         'bdgpinsituexpr'     => 'bdgpinsituexpr',
                         'dedb'               => 'dedb',
                         'flygrid'            => 'FlyGrid' );

# This is for source_ids that depend on the type of 'ID' of the line.
our %special_source_name_map = (
                                'gene' => {
                                         'Dbxref' => 'FlyBaseCGID_gene',
                                         'Name'   => 'FlyBaseName_gene',
                                         'ID'     => 'flybase_gene_id'
                                },
                                'transcript' => {
                                   'Dbxref' => 'FlyBaseCGID_transcript',
                                   'Name'   => 'FlyBaseName_transcript',
                                   'ID'     => 'flybase_transcript_id'
                                },
                                'translation' => {
                                  'Dbxref' => 'FlyBaseCGID_translation',
                                  'Name'   => 'FlyBaseName_translation',
                                  'ID'     => 'flybase_translation_id'
                                } );

# This hash will eventually be populated with the source_id for the
# sources above.
our %source_id;

sub get_source_id_for_source_name {
  my $self = shift;
  my ($source_name) = @_;

  if ( !defined( $source_id{$source_name} ) ) {
    $source_id{$source_name} =
      $self->SUPER::get_source_id_for_source_name(@_);
  }

  if ( !defined( $source_id{$source_name} )
       || $source_id{$source_name} < 0 )
  {
    croak(
       sprintf( "Can not find source_id for source '%s'", $source_name )
    );
  }

  return $source_id{$source_name};
}

sub run {
  my $self = shift;
  my ( $source_id, $species_id, $data_file, $release_file ) = @_;

  # Fetch a hash of the already stored Uniprot accessions.
  my %uniprot_xref_ids =
    %{ $self->get_valid_code( 'uniprot', $species_id ) };

  my $data_io = $self->get_filehandle($data_file);

  while ( defined( my $line = $data_io->getline() ) ) {
    # Skip comment lines at the start of the file.
    if ( substr( $line, 0, 1 ) eq '#' ) { next }

    chomp($line);

    # Split each line into fields.
    my @fields = split( /\t/, $line );

    # Only pick out the interesting lines.
    if (
         !(    defined( $fields[1] )
            && defined( $fields[2] )
            && $fields[1] eq 'FlyBase'
            && exists( $object_types{ $fields[2] } ) ) )
    {
      next;
    }

    # Go though each attribute (from the 9th field), split them up into
    # key-value pairs and store them.
    my %attributes;
    foreach my $attribute ( split( /;/, $fields[8] ) ) {
      my ( $key, $value ) = split( /=/, $attribute );
      if ( $key ne '' && $value ne '' ) {
        $attributes{$key} = $value;
      }
    }

    my $id = $attributes{'ID'};

    my $type;
    if    ( substr( $id, 0, 4 ) eq 'FBgn' ) { $type = 'gene' }
    elsif ( substr( $id, 0, 4 ) eq 'FBtr' ) { $type = 'transcript' }
    elsif ( substr( $id, 0, 4 ) eq 'FBpp' ) { $type = 'translation' }
    else                                    { $type = 'unknown' }

    # For the 'Dbxref' and 'Ontology_term' attributes, split them up on
    # commas, divide into key-value pairs, and store them.
    foreach my $attribute_key ( 'Dbxref', 'Ontology_term' ) {
      if ( exists( $attributes{$attribute_key} ) ) {
        my %tmphash;
        foreach
          my $subattribute ( split( /,/, $attributes{$attribute_key} ) )
        {
          my ( $key, $value ) = split( /:/, $subattribute, 2 );
          push( @{ $tmphash{$key} }, $value );
        }

        # Replace the attribute entry with the hash.
        $attributes{$attribute_key} = \%tmphash;
      }
    }

    my $dbxref = $attributes{'Dbxref'};

    #-------------------------------------------------------------------
    # Store Xrefs and Direct Xrefs for all the interesting Dbxref
    # entries.
    #-------------------------------------------------------------------
    foreach my $dbxref_name ( keys( %{$dbxref} ) ) {
      my $source_id;
      if ( exists( $source_name_map{$dbxref_name} ) ) {
        $source_id =
          $self->get_source_id_for_source_name(
                                       $source_name_map{$dbxref_name} );

        # Treat Uniprot differently.
        if ( substr( $dbxref_name, 0, 7 ) eq 'UniProt' ) {
          foreach my $accession ( @{ $dbxref->{$dbxref_name} } ) {
            if ( exists( $uniprot_xref_ids{$accession} ) ) {
              $self->add_direct_xref( $uniprot_xref_ids{$accession},
                                      $id, $type, '' );
            }
          }
        } else {
          foreach my $accession ( @{ $dbxref->{$dbxref_name} } ) {
            my $xref_id =
              $self->add_xref( $accession, '', $accession, '',
                               $source_id, $species_id );
            $self->add_direct_xref( $xref_id, $id, $type, '' );
          }
        }
      }
    }

    #-------------------------------------------------------------------
    # Store Xrefs and Direct Xrefs for the GO 'Ontology_term' entries.
    #-------------------------------------------------------------------
    if ( exists( $attributes{'Ontology_term'}{'GO'} ) ) {
      my $source_id = $self->get_source_id_for_source_name('GO');

      foreach my $accession ( @{ $attributes{'Ontology_term'}{'GO'} } )
      {
        my $xref_id =
          $self->add_xref( $accession, '', $accession, '', $source_id,
                           $species_id );
        $self->add_direct_xref( $xref_id, $id, $type, '' );
      }
    }

    #-------------------------------------------------------------------
    # Store Xrefs and Direct Xrefs for the 'FlyBase_Annotation_IDs'
    # Dbxref entry (depends on type of 'ID').
    #-------------------------------------------------------------------
    if ( exists( $dbxref->{'FlyBase_Annotation_IDs'} ) ) {
      my $source_id =
        $self->get_source_id_for_source_name(
                            $special_source_name_map{$type}{'Dbxref'} );

      foreach my $accession ( @{ $dbxref->{'FlyBase_Annotation_IDs'} } )
      {
        my $xref_id =
          $self->add_xref( $accession, '', $accession, '', $source_id,
                           $species_id );
        $self->add_direct_xref( $xref_id, $id, $type, '' );
      }

    }

    #-------------------------------------------------------------------
    # Store Xref and Direct Xref for the 'Name' (depends on type of
    # 'ID').
    #-------------------------------------------------------------------
    {
      my $source_id =
        $self->get_source_id_for_source_name(
                              $special_source_name_map{$type}{'Name'} );
      my $xref_id =
        $self->add_xref( $attributes{'Name'}, '', $attributes{'Name'},
                         $source_id, $species_id );
      $self->add_direct_xref( $xref_id, $id, $type, '' );
    }

    #-------------------------------------------------------------------
    # Store Xref and Direct Xref for the 'ID' (depends on type of 'ID').
    #-------------------------------------------------------------------
    {
      my $source_id =
        $self->get_source_id_for_source_name(
                                $special_source_name_map{$type}{'ID'} );
      my $xref_id =
        $self->add_xref( $id, '', $id, $source_id, $species_id );
      $self->add_direct_xref( $xref_id, $id, $type, '' );
    }

  } ## end while ( defined( my $line...
  $data_io->close();

} ## end sub run

1;
