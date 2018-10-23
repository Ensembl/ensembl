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

    # ...or no Ensembl ID or EntrezGene xref to match it to
    # FIXME: this number might be underestimated because it doesn't
    # check if Ensembl IDs from the input file actually exist
    if ( ( ! $ensembl_id )
         && ( ( ! $entrez_id ) || ( ! defined $entrez{$entrez_id} ) ) ) {
      $counters{'missed_master'}++;
      next RECORD;
    }

    # An unknown type might indicate the change of input format,
    # therefore make sure the user notices it. That said, do not
    # bother we do not have an xref this entry would operate on anyway
    # - which is why we only check this after the preceding two
    # presence checks.
    if ( ( $type ne 'gene')
         && ( $type ne 'gene/phenotype' )
         && ( $type ne 'predominantly phenotypes' )
         && ( $type ne 'phenotype' ) ) {
      croak "Unknown type $type";
    }

    if ( $type eq 'phenotype' ) {
      # The only type with a different source order

      if ( defined( $mim_morbid{$omim_id} ) ) {
        foreach my $mim_id ( @{ $mim_morbid{$omim_id} } ) {
          $self->process_xref_entry({
            'mim_id'           => $mim_id,
            'mim_source_id'    => $mim_morbid_source_id,
            'species_id'       => $species_id,
            'ensembl_id'       => $ensembl_id,
            'entrez_xrefs'     => $entrez{$entrez_id},
            'entrez_source_id' => $entrez_source_id,
            'dbi'              => $dbi,
            'counters'         => \%counters
          });
        }
      }
      else {
        $counters{'expected_in_mim_morbid_but_not_there'}++;
        foreach my $mim_id ( @{ $mim_gene{$omim_id} } ) {
          $self->process_xref_entry({
            'mim_id'           => $mim_id,
            'mim_source_id'    => $mim_gene_source_id,
            'species_id'       => $species_id,
            'ensembl_id'       => $ensembl_id,
            'entrez_xrefs'     => $entrez{$entrez_id},
            'entrez_source_id' => $entrez_source_id,
            'dbi'              => $dbi,
            'counters'         => \%counters
          });
        }
      }
    }
    else {
      if ( defined( $mim_gene{$omim_id} ) ) {
        foreach my $mim_id ( @{ $mim_gene{$omim_id} } ) {
          $self->process_xref_entry({
            'mim_id'           => $mim_id,
            'mim_source_id'    => $mim_gene_source_id,
            'species_id'       => $species_id,
            'ensembl_id'       => $ensembl_id,
            'entrez_xrefs'     => $entrez{$entrez_id},
            'entrez_source_id' => $entrez_source_id,
            'dbi'              => $dbi,
            'counters'         => \%counters
          });
        }
      }
      else {
        $counters{'expected_in_mim_gene_but_not_there'}++;
        foreach my $mim_id ( @{ $mim_morbid{$omim_id} } ) {
          $self->process_xref_entry({
            'mim_id'           => $mim_id,
            'mim_source_id'    => $mim_morbid_source_id,
            'species_id'       => $species_id,
            'ensembl_id'       => $ensembl_id,
            'entrez_xrefs'     => $entrez{$entrez_id},
            'entrez_source_id' => $entrez_source_id,
            'dbi'              => $dbi,
            'counters'         => \%counters
          });
        }
      }
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

  Readonly my @field_patterns
    => (
        qr{ \A [#]? \s* MIM[ ]Number }msx,
        qr{ MIM[ ]Entry[ ]Type }msx,
        qr{ Entrez[ ]Gene[ ]ID }msx,
        qr{ Approved[ ]Gene[ ]Symbol }msx,
        qr{ Ensembl[ ]Gene[ ]ID }msx,
      );

  my $header_field;
  foreach my $pattern (@field_patterns) {
    $header_field = shift @{ $header };
    # Make sure we run the regex match in scalar context
    push @fields_ok, scalar ( $header_field =~ m{ $pattern }msx );
  }

  # All fields must have matched
  return List::Util::all { $_ } @fields_ok;
}


=head2 process_xref_entry

  Arg [1]    : HashRef list of named arguments: FIXME
  Example    : $self->process_xref_entry({...});
  Description: Wrapper around the most frequently repeated bit of
               run(): if $ensembl_id is defined insert a direct MIM
               xref, otherwise loop over the list of matching
               EntrezGene xrefs and insert dependent MIM xrefs.
               In either case increment the correct counter.
  Return type: none
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub process_xref_entry {
  my ( $self, $arg_ref ) = @_;

  if ( $arg_ref->{'ensembl_id'} ) {
    $arg_ref->{'counters'}->{'direct_ensembl'}++;
    $self->add_to_direct_xrefs({
      'stable_id'  => $arg_ref->{'ensembl_id'},
      'type'       => 'gene',
      'acc'        => $arg_ref->{'mim_id'},
      'source_id'  => $arg_ref->{'mim_source_id'},
      'species_id' => $arg_ref->{'species_id'},
      'dbi'        => $arg_ref->{'dbi'},
    });
  }
  else {
    foreach my $ent_id ( @{ $arg_ref->{'entrez_xrefs'} } ) {
      $arg_ref->{'counters'}->{'dependent_on_entrez'}++;
      $self->add_dependent_xref({
        'master_xref_id' => $ent_id,
        'acc'            => $arg_ref->{'mim_id'},
        'source_id'      => $arg_ref->{'mim_source_id'},
        'species_id'     => $arg_ref->{'species_id'},
        'linkage'        => $arg_ref->{'entrez_source_id'},
        'dbi'            => $arg_ref->{'dbi'},
      });
    }
  }

  return;
}


1;
