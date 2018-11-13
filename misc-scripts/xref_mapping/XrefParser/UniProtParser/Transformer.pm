=head1 LICENSE

See the NOTICE file distributed with this work for additional
information regarding copyright ownership.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut



package XrefParser::UniProtParser::Transformer;

use strict;
use warnings;

use Carp;
use Readonly;


# FIXME: this belongs in BaseParser
Readonly my $ERR_SOURCE_ID_NOT_FOUND => -1;

Readonly my $PROTEIN_ID_SOURCE_NAME => 'protein_id';
Readonly my $UNIPROT_GN_SOURCE_NAME => 'Uniprot_gn';

# FIXME: this should probably be combined with
# Extractor::%supported_taxon_database_qualifiers to make sure
# database qualifiers stay in sync
Readonly my %taxonomy_ids_from_taxdb_codes
  => {
      # NCBI taxon codes and Ensembl taxonomy IDs are identical
      'NCBI_TaxID' => sub { return $_[0]; },
    };

Readonly my %whitelisted_crossreference_sources
  => (
      'ChEMBL'                => 1,
      'EMBL'                  => 1,
      'Ensembl'               => 1,
      'MEROPS'                => 1,
      'PDB'                   => 1,
      $PROTEIN_ID_SOURCE_NAME => 1,
      $UNIPROT_GN_SOURCE_NAME => 1,
    );

Readonly my $MAX_TREMBL_EVIDENCE_LEVEL_FOR_STANDARD => 2;
Readonly my %source_selection_criteria_for_status
  => (
      'Reviewed'   => [ 'Uniprot/SWISSPROT',
                        sub {
                          return 'sequence_mapped';
                        }, ],
      'Unreviewed' => [ 'Uniprot/SPTREMBL', sub {
                          my ( $level ) = @_;
                          return ( $level <= $MAX_TREMBL_EVIDENCE_LEVEL_FOR_STANDARD ) ?
                            'sequence_mapped' :
                            "protein_evidence_gt_$MAX_TREMBL_EVIDENCE_LEVEL_FOR_STANDARD";
                        }, ],
    );

Readonly my %protein_id_extraction_recipe_for_database
  => (
      'ChEMBL' => \&_get_protein_id_xref_from_embldb_xref,
      'EMBL'   => \&_get_protein_id_xref_from_embldb_xref,
    );
sub _get_protein_id_xref_from_embldb_xref {
  my ( $protein_id, $linkage_source_id, $source_id ) = @_;

  # Strip the version number, if any, from the protein ID. At the same
  # time, filter out entries with no ID - in which case the ID is a
  # lone hyphen.
  # FIXME:
  #  - are versioned primary IDs still a thing? There are no such
  #    entries in the Swiss-Prot file
  #  - ditto primary ID being absent
  my ( $unversioned_protein_id )
    = ( $protein_id =~ m{
                          \A
                          # Allow hyphens if they are not
                          # the first character
                          ( [^-.] [^.]+ )
                      }msx );

  if ( ! defined $unversioned_protein_id ) {
    return;
  }

  my $xref_link = {
                   'ACCESSION'         => $unversioned_protein_id,
                   'LABEL'             => $protein_id,
                   'LINKAGE_SOURCE_ID' => $linkage_source_id,
                   'SOURCE_ID'         => $source_id,
              };
  return $xref_link;
}



sub new {
  my ( $proto, $arg_ref ) = @_;

  my $self = {
              'dbh'        => $arg_ref->{'dbh'},
              'species_id' => $arg_ref->{'species_id'},
              'maps'       => {},
            };
  my $class = ref $proto || $proto;
  bless $self, $class;

  $self->_load_maps( $arg_ref->{'baseParser'} );

  return $self;
}


# Transforms extracted record into form that can be consumed by
# BaseParser::upload_xref_object_graphs().
sub transform {
  my ( $self, $extracted_record ) = @_;

  $self->{'extracted_record'} = $extracted_record;

  # Only proceed if at least one taxon code in the entry maps
  # to a valid Ensembl species ID
  my $xref_multiplicity = $self->_recognised_taxon_ids();
  if ( ! $xref_multiplicity ) {
    return;
  }

  # Skip unreviewed entries
  if ( $self->_entry_is_unreviewed() ) {
    return;
  }

  my ( $accession, @synonyms )
    = @{ $extracted_record->{'accession_numbers'} };
  my $source_id = $self->_get_source_id();

  my $xref_graph_node
    = {
       'ACCESSION'     => $accession,
       'DESCRIPTION'   => $extracted_record->{'description'},
       'INFO_TYPE'     => 'SEQUENCE_MATCH',
       'LABEL'         => $accession,
       'SEQUENCE'      => $extracted_record->{'sequence'},
       'SEQUENCE_TYPE' => 'peptide',
       'SOURCE_ID'     => $source_id,
       'SPECIES_ID'    => $self->{'species_id'},
       'STATUS'        => 'experimental',
       'SYNONYMS'      => \@synonyms,
       '_multiplicity' => $xref_multiplicity, # hint for Loader
     };

  # UniProt Gene Names links come from the 'gene_names' fields
  my $genename_dependent_xrefs
    = $self->_make_links_from_gene_names( $accession, $source_id );
  # Do not assign an empty array to DEPENDENT_XREFS, current insertion code
  # doesn't like them.
  if ( scalar @{ $genename_dependent_xrefs } > 0 ) {
    push @{ $xref_graph_node->{'DEPENDENT_XREFS'} }, @{ $genename_dependent_xrefs };
  }
  # FIXME: if this outputs multiple dependent xrefs, all their links
  # will point to the first xref object created. This is due to a bug
  # in BaseParser.

  # All other xref links come from crossreferences
  my ( $direct_xrefs, $dependent_xrefs )
    = $self->_make_links_from_crossreferences( $accession, $source_id );
  # Do not assign empty arrays to FOO_XREFS, current insertion code
  # doesn't like them.
  if ( scalar @{ $direct_xrefs } > 0 ) {
    push @{ $xref_graph_node->{'DIRECT_XREFS'} }, @{ $direct_xrefs };
  }
  if ( scalar @{ $dependent_xrefs } > 0 ) {
    push @{ $xref_graph_node->{'DEPENDENT_XREFS'} }, @{ $dependent_xrefs };
  }

  # FIXME: fix BaseParser to make dep-xref insertion via the graph replay-safe!!!

  return $xref_graph_node;
}


# FIXME: description
sub get_source_id_map {
  my ( $self ) = @_;

  # Just in case, even though we presently call _load_maps() in the
  # constructor so it shouldn't be possible to call this method before
  # maps have been loaded
  if ( ! exists $self->{'maps'}->{'named_source_ids'} ) {
    croak 'Source-ID map is missing';
  }

  return $self->{'maps'}->{'named_source_ids'};
}


sub _load_maps {
  my ( $self, $baseParserInstance ) = @_;

  my $taxonomy_ids_for_species
    = $baseParserInstance->get_taxonomy_from_species_id( $self->{'species_id'},
                                                         $self->{'dbh'} );
  # If the map is empty, something is wrong
  if ( scalar keys %{ $taxonomy_ids_for_species } == 0 ) {
    croak "Got zero taxonomy_ids for species_id '"
      . $self->{'species_id'} . q{'};
  }
  $self->{'maps'}->{'taxonomy_ids_for_species'}
    = $taxonomy_ids_for_species;

  my $source_id_map
    = {
       'Uniprot/SWISSPROT'
       => {
           'direct'          => undef,
           'sequence_mapped' => undef,
         },
       'Uniprot/SPTREMBL'
       => {
           'direct'            => undef,
           "protein_evidence_gt_$MAX_TREMBL_EVIDENCE_LEVEL_FOR_STANDARD" => undef,
           'sequence_mapped'   => undef,
         },
     };
  while ( my ( $source_name, $pri_ref ) = each %{ $source_id_map } ) {
    foreach my $priority ( keys %{ $pri_ref } ) {
      $pri_ref->{$priority}
        = $baseParserInstance->get_source_id_for_source_name( $source_name,
                                                              $priority,
                                                              $self->{'dbh'} );
      if ( $pri_ref->{$priority} == $ERR_SOURCE_ID_NOT_FOUND ) {
        croak "No source ID found for source $source_name, priority $priority";
      }
    }
  }
  $self->{'maps'}->{'named_source_ids'}
    = $source_id_map;

  my %dependent_source_map
    = $baseParserInstance->get_xref_sources( $self->{'dbh'} );
  $self->{'maps'}->{'dependent_sources'}
    = \%dependent_source_map;

  return;
}


# Returns true if the current record describes an entry derived from
# Ensembl, false otherwise.
sub _entry_is_from_ensembl {
  my ( $self ) = @_;

  # As of end of October 2018, the old parser's way of identifying
  # proteins derived from Ensembl by searching for a comment topic
  # "CAUTION" stating "The sequence shown here is derived from an
  # Ensembl" no longer works because there seem to be no entries
  # containing this string any more. Therefore, for the time being
  # this test always returns false.

  return 0;
}


# Returns true if the current record describes an entry tagged as
# unreviewed, false otherwise.
sub _entry_is_unreviewed {
  my ( $self ) = @_;

  # This is the way the old UniProtParser identified unreviewed
  # entries. FIXME: is this still a thing? As of October 2018 there
  # are NO such entries in either the first ~1000 lines of the TrEMBL
  # file or anywhere in the SwissProt one.
  my $accession_numbers = $self->{'extracted_record'}->{'accession_numbers'};
  if ( lc( $accession_numbers->[0] ) eq 'unreviewed' ) {
    return 1;
  }

  return 0;
}


# Translate quality of the extracted entry into the matching Ensembl
# source_id, optionally with an override of priority.
sub _get_source_id {
  my ( $self, $priority_override ) = @_;

  my $source_id_map = $self->{'maps'}->{'named_source_ids'};

  my $entry_quality = $self->{'extracted_record'}->{'quality'};
  my $criteria = $source_selection_criteria_for_status{ $entry_quality->{'status'} };
  my $priority_mapper = $criteria->[1];

  my $source_name = $criteria->[0];
  my $priority = $priority_override
    // $priority_mapper->( $entry_quality->{'evidence_level'} );

  return $source_id_map->{$source_name}->{$priority};
}


# Make xrefs from 'crossreferences' entries in the extracted record,
# in a form suitable to attaching to the main xref's graph node as
# consumed by upload_xref_object_graphs(). Ensembl crossreferences
# become direct xrefs, everything else - dependent ones. If requested
# we additionally generate protein_id dependent xrefs from appropriate
# sources, i.e. EMBL and ChEMBL at present.
sub _make_links_from_crossreferences {
  my ( $self, $xref_accession, $xref_source_id ) = @_;

  my $crossreferences = $self->{'extracted_record'}->{'crossreferences'};
  my $dependent_sources = $self->{'maps'}->{'dependent_sources'};

  my @direct_xrefs;
  my @dependent_xrefs;

 REF_SOURCE:
  while ( my ( $source, $entries ) = each %{ $crossreferences } ) {

    if ( ! $whitelisted_crossreference_sources{ $source } ) {
      next REF_SOURCE;
    }

    if ( $source eq 'Ensembl' ) {

    DIRECT_XREF:
      foreach my $direct_ref ( @{ $entries } ) {
        my $xref_link
          = {
             'STABLE_ID'    => $direct_ref->{'id'},
             'ENSEMBL_TYPE' => 'Translation',
             'LINKAGE_TYPE' => 'DIRECT',
             'SOURCE_ID'    => $self->_get_source_id( 'direct' ),
           };
        push @direct_xrefs, $xref_link;
      }

    }
    elsif ( exists $dependent_sources->{$source} ) {

    DEPENDENT_XREF:
      foreach my $dependent_ref ( @{ $entries } ) {
        my $xref_link
          = {
             'ACCESSION'         => $xref_accession,
             'LINKAGE_SOURCE_ID' => $xref_source_id,
             'SOURCE_ID'         => $dependent_sources->{$source},
           };
        push @dependent_xrefs, $xref_link;

        my $protein_id_xref_maker
          = $protein_id_extraction_recipe_for_database{ $source };
        if ( $whitelisted_crossreference_sources{ $PROTEIN_ID_SOURCE_NAME }
             && ( defined $protein_id_xref_maker ) ) {

          # Entries for the source 'protein_id' are constructed from
          # crossreferences to other databases
          my $protein_id_xref
            = $protein_id_xref_maker->( $dependent_ref->{'id'},
                                        $xref_source_id,
                                        $dependent_sources->{$PROTEIN_ID_SOURCE_NAME}
                                     );
          if ( defined $protein_id_xref ) {
            push @dependent_xrefs, $protein_id_xref;
          }

        }

      }

    }
  }

  return ( \@direct_xrefs, \@dependent_xrefs );
}


# Make Uniprot_gn dependent xrefs from 'gene_names' entries in the
# extracted record, in a form suitable to attaching to the main xref's
# graph node as consumed by upload_xref_object_graphs().
sub _make_links_from_gene_names {
  my ( $self, $xref_accession, $xref_source_id ) = @_;

  my @genename_xrefs;

  # Are we supposed to process this xref source to begin with?
  if ( ! $whitelisted_crossreference_sources{ 'Uniprot_gn' } ) {
    return [];
  }

  # UniProt Gene Name xrefs are dependent so in order to avoid
  # circular dependencies, do not generate them for proteins derived
  # from Ensembl.
  if ( $self->_entry_is_from_ensembl() ) {
    return [];
  }

  my $gene_names = $self->{'extracted_record'}->{'gene_names'};
  my $dependent_sources = $self->{'maps'}->{'dependent_sources'};

 GN_ENTRY:
  foreach my $gn_entry ( @{ $gene_names } ) {
    if ( ! exists $gn_entry->{'Name'} ) {
      next GN_ENTRY;
    }

    my $name = $gn_entry->{'Name'};
    my $xref = {
                'ACCESSION'         => $xref_accession,
                'LABEL'             => $name,
                'SOURCE_ID'         => $dependent_sources->{'Uniprot_gn'},
                'LINKAGE_SOURCE_ID' => $xref_source_id,
              };

    my $synonyms = $gn_entry->{'Synonyms'};
    if ( defined $synonyms ) {
      push @{ $xref->{'SYNONYMS'} }, @{ $synonyms };
    }

    push @genename_xrefs, $xref;
  }

  return \@genename_xrefs;
}


# Translate extracted taxon codes into Ensembl taxonomy IDs, then
# return the number of taxons matching the species ID under
# consideration.
sub _recognised_taxon_ids {
  my ( $self ) = @_;

  my $taxon_codes = $self->{'extracted_record'}->{'taxon_codes'};
  my $tid4s_map = $self->{'maps'}->{'taxonomy_ids_for_species'};

  my @taxonomy_ids;
  foreach my $taxon ( @{ $taxon_codes } ) {
    my $code_mapper
      = $taxonomy_ids_from_taxdb_codes{ $taxon->{'db_qualifier'} };
    push @taxonomy_ids, $code_mapper->( $taxon->{'taxon_code'} );
  }

  my $recognised_taxonomy_ids = 0;
  foreach my $taxonomy_id ( @taxonomy_ids ) {
    $recognised_taxonomy_ids += ( exists $tid4s_map->{$taxonomy_id} );
  }

  return $recognised_taxonomy_ids;
}


1;
