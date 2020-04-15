=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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
use List::Util;
use Text::CSV;

use parent qw( XrefParser::BaseParser );


my $EXPECTED_NUMBER_OF_COLUMNS = 6;



=head2 run

  Arg [1]    : HashRef standard list of arguments from ParseSource
  Example    : $m2g_parser->run({ ... });
  Description: Extract mappings between OMIM genes and other gene
               identifiers from a tab-delimited file downloaded from
               the DBASS Web site, then insert corresponding links
               into the xref database:
                - for entries mapped to Ensembl genes, we create
                  gene_direct_xref links;
                - otherwise, if an entry is mapped to an EntrezGene ID
                  that exists in the xref database we creare a
                  dependent_xref link.
               In either case we update info_type of OMIM xrefs
               accordingly.

               DEPENDENCIES: This parser must be run after:
                - MIMParser - without existing OMIM entries this
                  parser does nothing;
                - EntrezGeneParser - otherwise there will be no
                  dependent-xref links.

               mim2gene.txt begins with several lines of comments
               which start with a hash; the last of these comment
               lines contains a tab-separated list of column names.

               The rest of the file are the following columns:
                1) OMIM number
                2) OMIM entry type
                3) EntrezGene ID
                4) HGNC gene symbol
                5) Ensembl gene ID
               The former two are mandatory, the latter can be empty
               strings.

  Return type: none
  Exceptions : throws on all processing errors
  Caller     : ParseSource in the xref pipeline
  Status     : Stable

=cut

sub run {

  my ( $self, $ref_arg ) = @_;
  my $general_source_id = $ref_arg->{source_id};
  my $species_id        = $ref_arg->{species_id};
  my $files             = $ref_arg->{files};
  my $verbose           = $ref_arg->{verbose} // 0;
  my $dbi               = $ref_arg->{dbi} // $self->dbi;

  if ( ( !defined $general_source_id ) or
       ( !defined $species_id ) or
       ( !defined $files ) )
  {
    confess "Need to pass source_id, species_id and files as pairs";
  }

  my $csv = Text::CSV->new({
                            sep_char => "\t",
                          })
    || confess 'Failed to initialise CSV parser: ' . Text::CSV->error_diag();

  my $filename = @{$files}[0];

  my $m2g_io = $self->get_filehandle($filename);
  if ( !defined $m2g_io ) {
    confess "Could not open file '${filename}'";
  }

  my $mim_gene_source_id =
    $self->get_source_id_for_source_name( 'MIM_GENE', undef, $dbi );
  my $mim_morbid_source_id =
    $self->get_source_id_for_source_name( 'MIM_MORBID', undef, $dbi );
  my $entrez_source_id =
    $self->get_source_id_for_source_name( 'EntrezGene', undef, $dbi );

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
                  'missed_master'                        => 0,
                  'missed_omim'                          => 0,
                );

 RECORD:
  while ( my $line = $csv->getline( $m2g_io ) ) {

    my ( $is_comment )
      = ( $line->[0] =~ m{
                           \A
                           ([#])?
                       }msx );
    if ( $is_comment ) {
      # At present we identify the header line among other comments by
      # checking if it has the expected number of tab-delimited
      # columns, which of course means we cannot identify header lines
      # with too few or too many column names. However, this should be
      # mostly harmless - something would have to be very, very wrong
      # with the input file for the header to have the wrong number of
      # column names without a change in the number of actual columns
      # in data rows.
      if ( ( scalar @{ $line } == $EXPECTED_NUMBER_OF_COLUMNS )
           && ( ! is_file_header_valid( @{ $line } ) ) ) {
        confess "Malformed or unexpected header in Mim2Gene file '${filename}'";
      }
      next RECORD;
    }

    if ( scalar @{ $line } != $EXPECTED_NUMBER_OF_COLUMNS ) {
      confess ' Line ' . $csv->record_number()
        . " of input file '${filename}' has an incorrect number of columns";
    }

    # Do not modify the contents of @{$line}, only the output - hence the /r.
    my ( $omim_acc, $entrez_id, $type, $source, $medgen, $comment )
      = map { s{\s+\z}{}rmsx } @{ $line };

    $counters{'all_entries'}++;

    # No point in doing anything if we have no matching MIM xref...
    if ( ( !defined $mim_gene{$omim_acc} ) &&
         ( !defined $mim_morbid{$omim_acc} ) )
    {
      $counters{'missed_omim'}++;
      next RECORD;
    }

    # ...or no EntrezGene xref to match it to
    if ( ( ( ! $entrez_id ) || ( ! defined $entrez{$entrez_id} ) ) ) {
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
      confess "Unknown type $type for MIM Number '${omim_acc}' "
        . "(${filename}:" . $csv->record_number() . ")";
    }

    # With all the checks taken care of, insert the mappings. We check
    # both MIM_GENE and MIM_MORBID every time because some MIM entries
    # can appear in both.
    foreach my $mim_xref_id ( @{ $mim_gene{$omim_acc} } ) {
      $self->process_xref_entry({
        'mim_xref_id'      => $mim_xref_id,
        'mim_source_id'    => $mim_gene_source_id,
        'entrez_xrefs'     => $entrez{$entrez_id},
        'entrez_source_id' => $entrez_source_id,
        'counters'         => \%counters,
        'dbi'              => $dbi,
      });
    }
    foreach my $mim_xref_id ( @{ $mim_morbid{$omim_acc} } ) {
      $self->process_xref_entry({
        'mim_xref_id'      => $mim_xref_id,
        'mim_source_id'    => $mim_morbid_source_id,
        'entrez_xrefs'     => $entrez{$entrez_id},
        'entrez_source_id' => $entrez_source_id,
        'counters'         => \%counters,
        'dbi'              => $dbi,
      });
    }

  } ## end record loop

  $csv->eof || confess 'Error parsing CSV: ' . $csv->error_diag();
  $m2g_io->close();

  if ( $verbose ) {
    print 'Processed ' . $counters{'all_entries'} . " entries. Out of those\n"
      . "\t" . $counters{'missed_omim'} . " had missing OMIM entries,\n"
      . "\t" . $counters{'direct_ensembl'} . " were direct gene xrefs,\n"
      . "\t" . $counters{'dependent_on_entrez'} . " were dependent EntrezGene xrefs,\n"
      . "\t" . $counters{'missed_master'} . " had missing master entries.\n";
  }

  return 0;
} ## end sub run


=head2 is_file_header_valid

  Arg [1..N] : list of column names provided by Text::CSV::getline()
  Example    : if ( ! is_file_header_valid( $csv->getline( $fh ) ) {
                 confess 'Bad header';
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

sub is_file_header_valid {
  my ( @header ) = @_;

  my @field_patterns
    = (
        qr{ \A [#]? \s* MIM[ ]number }msx,
        qr{ GeneID }msx,
        qr{ type }msx,
        qr{ Source }msx,
        qr{ MedGenCUI }msx,
        qr{ Comment }msx,
      );

  my $header_field;
  foreach my $pattern (@field_patterns) {
    $header_field = shift @header;
    # Make sure we run the regex match in scalar context
    return 0 unless scalar ( $header_field =~ m{ $pattern }msx );
  }

  # If we have made it this far, all should be in order
  return 1;
}


=head2 process_xref_entry

  Arg [1]    : HashRef list of named arguments: FIXME
  Example    : $self->process_xref_entry({...});
  Description: Wrapper around the most frequently repeated bit of
               run(): loop over the list of matching
               EntrezGene xrefs and insert dependent MIM xrefs.
  Return type: none
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub process_xref_entry {
  my ( $self, $arg_ref ) = @_;

  foreach my $ent_id ( @{ $arg_ref->{'entrez_xrefs'} } ) {
    $arg_ref->{'counters'}->{'dependent_on_entrez'}++;
    $self->add_dependent_xref_maponly( $arg_ref->{'mim_xref_id'},
                                       $arg_ref->{'mim_source_id'},
                                       $ent_id,
                                       $arg_ref->{'entrez_source_id'},
                                       $arg_ref->{'dbi'},
                                       1
                                    );
  }

  return;
}


1;
