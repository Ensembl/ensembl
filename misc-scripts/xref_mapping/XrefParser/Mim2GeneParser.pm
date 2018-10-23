=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

package XrefParser::Mim2GeneParser;

# For non-destructive substitutions in regexps (/r flag)
require 5.014_000;

use strict;
use warnings;

use Carp;
use File::Basename;
use List::Util;
use POSIX qw(strftime);
use Readonly;
use Text::CSV;

use parent qw( XrefParser::BaseParser );


# FIXME: this belongs in BaseParser
Readonly my $ERR_SOURCE_ID_NOT_FOUND => -1;

Readonly my $EXPECTED_NUMBER_OF_COLUMNS => 5;


sub run {

  my ( $self, $ref_arg ) = @_;
  my $general_source_id = $ref_arg->{source_id};
  my $species_id        = $ref_arg->{species_id};
  my $files             = $ref_arg->{files};
  my $verbose           = $ref_arg->{verbose};
  my $dbi               = $ref_arg->{dbi};

  if ( ( !defined $general_source_id ) or
       ( !defined $species_id ) or
       ( !defined $files ) )
  {
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $dbi //= $self->dbi;
  $verbose //= 0;


  my $csv = Text::CSV->new({
                            sep_char => "\t",
                          })
    || croak 'Failed to initialise CSV parser: ' . Text::CSV->error_diag();

  my $filename = @{$files}[0];

  my $eg_io = $self->get_filehandle($filename);
  if ( !defined $eg_io ) {
    croak "Could not open file '${filename}'";
  }

  my $mim_gene_source_id =
    $self->get_source_id_for_source_name( 'MIM_GENE', undef, $dbi );
  my $mim_morbid_source_id =
    $self->get_source_id_for_source_name( 'MIM_MORBID', undef, $dbi );
  my $entrez_source_id =
    $self->get_source_id_for_source_name( 'EntrezGene', undef, $dbi );
  if ( ( $mim_gene_source_id == $ERR_SOURCE_ID_NOT_FOUND )
       || ( $mim_morbid_source_id == $ERR_SOURCE_ID_NOT_FOUND )
       || ( $entrez_source_id == $ERR_SOURCE_ID_NOT_FOUND ) ) {
    croak 'Failed to retrieve all source IDs';
  }

  # This will be used to prevent insertion of duplicates
  $self->get_dependent_mappings( $mim_gene_source_id, $dbi );
  $self->get_dependent_mappings( $mim_morbid_source_id, $dbi );

  # FIXME: should we abort if any of these comes back empty?
  my (%mim_gene) =
    %{ $self->get_valid_codes( "MIM_GENE", $species_id, $dbi ) };
  my (%mim_morbid) =
    %{ $self->get_valid_codes( "MIM_MORBID", $species_id, $dbi ) };
  my (%entrez) =
    %{ $self->get_valid_codes( "EntrezGene", $species_id, $dbi ) };

  # Initialise all counters to 0 so that we needn't handle possible undefs
  # while printing the summary
  my %counters = (
                  'all_entries'                          => 0,
                  'dependent_on_entrez'                  => 0,
                  'direct_ensembl'                       => 0,
                  'expected_in_mim_gene_but_not_there'   => 0,
                  'expected_in_mim_morbid_but_not_there' => 0,
                  'missed_master'                        => 0,
                  'missed_omim'                          => 0,
                );

 RECORD:
  while ( my $line = $csv->getline( $eg_io ) ) {

    my ( $is_comment, $is_header )
      = ( $line->[0] =~ m{
                           \A
                           ([#])?
                           \s*
                           (MIM[ ]Number)?  # FIXME: this is an assumption regarding header contents.
                                            # See if $line has split to the right number of columns instead?
                       }msx );
    if ( $is_comment ) {
      if ( ( scalar @{ $line } == $EXPECTED_NUMBER_OF_COLUMNS )
           && ( ! is_header_file_valid( $line ) ) ) {
        croak "Malformed or unexpected header in Mim2Gene file '${filename}'";
      }
      next RECORD;
    }

    if ( scalar @{ $line } != $EXPECTED_NUMBER_OF_COLUMNS ) {
      croak ' Line ' . $csv->record_number()
        . " of input file '${filename}' has an incorrect number of columns";
    }

    # Do not modify the contents of @{$line}, only the output - hence the /r.
    my ( $omim_id, $type, $entrez_id, $hgnc_symbol, $ensembl_id )
      = map { s{\s+\z}{}rmsx } @{ $line };

    $counters{'all_entries'}++;

    # No point in doing anything if we have no matching MIM xref...
    if ( ( !defined $mim_gene{$omim_id} ) and
         ( !defined $mim_morbid{$omim_id} ) )
    {
      $counters{'missed_omim'}++;
      next RECORD;
    }

    # An unknown type might indicate the change of input format,
    # therefore make sure the user notices it. That said, do not
    # bother we do not have an xref this entry would operate on anyway
    # - which is why we only check this AFTER the missed_omim block
    # above.
    if ( ( $type ne 'gene')
         && ( $type ne 'gene/phenotype' )
         && ( $type ne 'predominantly phenotypes' )
         && ( $type ne 'phenotype' ) ) {
      croak "Unknown type $type";
    }

    if ( $ensembl_id ) {
      $counters{'direct_ensembl'}++;

      # FIXME: *lots* of duplication wrt dependent xrefs

      if ( $type eq 'phenotype' ) {
        # The only type with a different source order

        if ( defined( $mim_morbid{$omim_id} ) ) {
          foreach my $mim_id ( @{ $mim_morbid{$omim_id} } ) {
            $self->add_to_direct_xrefs({
                                        'stable_id'  => $ensembl_id,
                                        'type'       => 'gene',
                                        'acc'        => $mim_id,
                                        'source_id'  => $mim_morbid_source_id,
                                        'species_id' => $species_id,
                                        'dbi'        => $dbi,
                                      });
          }
        }
        else {
          $counters{'expected_in_mim_morbid_but_not_there'}++;
          foreach my $mim_id ( @{ $mim_gene{$omim_id} } ) {
            $self->add_to_direct_xrefs({
                                        'stable_id'  => $ensembl_id,
                                        'type'       => 'gene',
                                        'acc'        => $mim_id,
                                        'source_id'  => $mim_gene_source_id,
                                        'species_id' => $species_id,
                                        'dbi'        => $dbi,
                                      });
          }
        }
      }
      else {
        if ( defined( $mim_gene{$omim_id} ) ) {
          foreach my $mim_id ( @{ $mim_gene{$omim_id} } ) {
            $self->add_to_direct_xrefs({
                                        'stable_id'  => $ensembl_id,
                                        'type'       => 'gene',
                                        'acc'        => $mim_id,
                                        'source_id'  => $mim_gene_source_id,
                                        'species_id' => $species_id,
                                        'dbi'        => $dbi,
                                      });
          }
        }
        else {
          $counters{'expected_in_mim_gene_but_not_there'}++;
          foreach my $mim_id ( @{ $mim_morbid{$omim_id} } ) {
            $self->add_to_direct_xrefs({
                                        'stable_id'  => $ensembl_id,
                                        'type'       => 'gene',
                                        'acc'        => $mim_id,
                                        'source_id'  => $mim_morbid_source_id,
                                        'species_id' => $species_id,
                                        'dbi'        => $dbi,
                                      });
          }
        }
      }

    }
    elsif ( defined $entrez{$entrez_id} ) {
      $counters{'dependent_on_entrez'}++;

      # FIXME: *lots* of duplication wrt direct xrefs

      if ( $type eq 'phenotype' ) {
        # The only type with a different source order

        if ( defined( $mim_morbid{$omim_id} ) ) {
          foreach my $ent_id ( @{ $entrez{$entrez_id} } ) {
            foreach my $mim_id ( @{ $mim_morbid{$omim_id} } ) {
              $self->add_dependent_xref({
                                         'master_xref_id' => $ent_id,
                                         'acc'            => $mim_id,
                                         'source_id'      => $mim_morbid_source_id,
                                         'species_id'     => $species_id,
                                         'linkage'        => $entrez_source_id,
                                         'dbi'            => $dbi,
                                       });
            }
          }
        }
        else {
          $counters{'expected_in_mim_morbid_but_not_there'}++;
          foreach my $ent_id ( @{ $entrez{$entrez_id} } ) {
            foreach my $mim_id ( @{ $mim_gene{$omim_id} } ) {
              $self->add_dependent_xref({
                                         'master_xref_id' => $ent_id,
                                         'acc'            => $mim_id,
                                         'source_id'      => $mim_gene_source_id,
                                         'species_id'     => $species_id,
                                         'linkage'        => $entrez_source_id,
                                         'dbi'            => $dbi,
                                       });
            }
          }
        }
      }
      else {
        if ( defined( $mim_gene{$omim_id} ) ) {
          foreach my $ent_id ( @{ $entrez{$entrez_id} } ) {
            foreach my $mim_id ( @{ $mim_gene{$omim_id} } ) {
              $self->add_dependent_xref({
                                         'master_xref_id' => $ent_id,
                                         'acc'            => $mim_id,
                                         'source_id'      => $mim_gene_source_id,
                                         'species_id'     => $species_id,
                                         'linkage'        => $entrez_source_id,
                                         'dbi'            => $dbi,
                                       });
            }
          }
        }
        else {
          $counters{'expected_in_mim_gene_but_not_there'}++;
          foreach my $ent_id ( @{ $entrez{$entrez_id} } ) {
            foreach my $mim_id ( @{ $mim_morbid{$omim_id} } ) {
              $self->add_dependent_xref({
                                         'master_xref_id' => $ent_id,
                                         'acc'            => $mim_id,
                                         'source_id'      => $mim_morbid_source_id,
                                         'species_id'     => $species_id,
                                         'linkage'        => $entrez_source_id,
                                         'dbi'            => $dbi,
                                       });
            }
          }
        }
      }

    }
    else {
      $counters{'missed_master'}++;
    }

  } ## end record loop

  $csv->eof || croak 'Error parsing CSV: ' . $csv->error_diag();
  $eg_io->close();

  if ( $verbose ) {
    print 'Processed ' . $counters{'all_entries'} . " entries. Out of those\n"
      . "\t" . $counters{'missed_omim'} . " had missing OMIM entries,\n"
      . "\t" . $counters{'direct_ensembl'} . " were direct gene xrefs,\n"
      . "\t" . $counters{'dependent_on_entrez'} . " were dependent EntrezGene xrefs,\n"
      . "\t" . $counters{'missed_master'} . " had missing master entries.\n"
      . "\t * * *\n"
      . "\t" . $counters{'expected_in_mim_gene_but_not_there'}
      . " were expected in MIM_GENE but were not found there,\n"
      . "\t" . $counters{'expected_in_mim_morbid_but_not_there'}
      . " were expected in MIM_MORBID but were not found there.\n";
  }

  return 0;
} ## end sub run


=head2 is_file_header_valid
  Arg [1]    : String file header line
  Example    : if (!is_file_header_valid($header_line)) {
                 croak 'Bad header';
               }
  Description: Verifies if the header of a Mim2Gene file follows expected
               syntax.
               We do not check the number of columns because that is what
               we use to *detect* the header in the first place.
  Return type: boolean
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub is_header_file_valid {
  my ( $header ) = @_;

  my @fields_ok;

  # FIXME: we should probably use a loop + a lookup list for the code
  # below to avoid using hard-coded field indices

  # This one will likely have a hash prepended to it
  my $mim_number_ok = ( $header->[0] =~ m{ \A [#]? \s* MIM[ ]Number }msx );
  push @fields_ok, $mim_number_ok;

  my $mim_type_ok = ( $header->[1] =~ m{ MIM[ ]Entry[ ]Type }msx );
  push @fields_ok, $mim_type_ok;

  my $entrez_id_ok = ( $header->[2] =~ m{ Entrez[ ]Gene[ ]ID }msx );
  push @fields_ok, $entrez_id_ok;

  my $hgnc_id_ok = ( $header->[3] =~ m{ Approved[ ]Gene[ ]Symbol }msx );
  push @fields_ok, $hgnc_id_ok;

  my $ensembl_id_ok = ( $header->[4] =~ m{ Ensembl[ ]Gene[ ]ID }msx );
  push @fields_ok, $ensembl_id_ok;

  # All fields must be in order
  return List::Util::all { $_ } @fields_ok;
}

1;
