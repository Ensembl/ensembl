=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

package XrefParser::MIMParser;

use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
use File::Basename;

use parent qw( XrefParser::BaseParser );

# This parser will read xrefs from a record file downloaded from the
# OMIM Web site. They should be assigned to two different xref
# sources: MIM_GENE and MIM_MORBID. MIM xrefs are linked to EntrezGene
# entries so the parser does not match them to Ensembl; this will be
# taken care of when EntrezGene antries are matched.
#
# OMIM records are multiline. Each record begins with a specific tag
# line and consists of a number of fields. Each field starts with its
# own start-tag line (i.e. the data proper only appears after a
# newline) and continues until the beginning of either the next field
# in the same record, the next record, or the end-of-input tag. The
# fields are expected to appear IN FIXED, SPECIFIC ORDER. The overall
# structure looks as follows:
#
#   *RECORD*
#   *FIELD* NO
#   *FIELD* TI
#   *FIELD* TX
#   *FIELD* RF
#   *FIELD* CD
#   *FIELD* ED
#   *RECORD*
#   *FIELD* NO
#   *FIELD* TI
#   ...
#   *FIELD* CD
#   *FIELD* ED
#   *THEEND*
#
# All the data relevant to the parser can be found in the TI field.


sub run {

  my ( $self, $ref_arg ) = @_;
  my $general_source_id = $ref_arg->{source_id};
  my $species_id        = $ref_arg->{species_id};
  my $files             = $ref_arg->{files};
  my $verbose           = $ref_arg->{verbose};
  my $dbi               = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

  if ( ( !defined $general_source_id ) or
       ( !defined $species_id ) or
       ( !defined $files ) )
  {
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose |= 0;

  my $file = @{$files}[0];

  my %old_to_new;
  my %removed;
  my @sources;

  push @sources, $general_source_id;

  my $gene_source_id =
    $self->get_source_id_for_source_name( "MIM_GENE", undef, $dbi );
  push @sources, $gene_source_id;
  my $morbid_source_id =
    $self->get_source_id_for_source_name( "MIM_MORBID", undef, $dbi );
  push @sources, $morbid_source_id;
  if ( ( $gene_source_id == -1 ) || ( $morbid_source_id == -1 ) ) {
    croak 'Failed to retrieve MIM source IDs';
  }

  print "sources are:- " . join( ", ", @sources ) . "\n" if ($verbose);

  IO::Handle->input_record_separator('*RECORD*');

  my $mim_io = $self->get_filehandle($file);

  if ( !defined $mim_io ) {
    print {*STDERR} "ERROR: Could not open $file\n";
    return 1;    # 1 is an error
  }

  my $gene          = 0;
  my $phenotype     = 0;
  my $removed_count = 0;

  $mim_io->getline();    # first record is empty with *RECORD* as the
                         # record seperator

  while ( my $input_record = $mim_io->getline() ) {

    my ( $number )
      = ( $input_record =~ m{
                              [*]FIELD[*]\s+NO\n
                              (\d+)
                          }msx );
    if ( defined $number ) {

      my ( $ti )
        = ( $input_record =~ m{
                                [*]FIELD[*]\sTI\n
                                (.+)\n             # Grab the whole TI field
                                [*]FIELD[*]\sTX
                            }msx );
      if ( defined $ti ) {
        # Remove line breaks. FIXME: in some places it will result in words being concatenated
        $ti =~ s{\n}{}gmsx;

        # Extract the 'type' and the whole description
        my ( $type, $long_desc ) =
          ( $ti =~ m{
                      ([#%+*^]*)  # type of entry
                      \d+         # number (should be the same as from NO)
                      \s+         # normally just one space
                      (.+)        # description of entry
                  }msx );
        if ( !defined( $type ) ) {
          print {*STDERR} 'Failed to extract record type and description from TI field';
          return 1;
        }

        # Use the first block of text as description
        my @fields = split( qr{;;}msx, $long_desc );
        my $label = $fields[0] . " [" . $type . $number . "]";

        if ( $type eq q{*} ) {     # gene only
          $gene++;
          $self->add_xref(
                           { acc        => $number,
                             label      => $label,
                             desc       => $long_desc,
                             source_id  => $gene_source_id,
                             species_id => $species_id,
                             dbi        => $dbi,
                             info_type  => "DEPENDENT" } );
        }
        elsif ( ( !defined $type ) or
                ( $type eq q{} )  or
                ( $type eq q{#} ) or
                ( $type eq q{%} ) )
        {    #phenotype only
          $phenotype++;
          $self->add_xref(
                           { acc        => $number,
                             label      => $label,
                             desc       => $long_desc,
                             source_id  => $morbid_source_id,
                             species_id => $species_id,
                             dbi        => $dbi,
                             info_type  => "DEPENDENT" } );
        }
        elsif ( $type eq q{+} ) {    # both
          $gene++;
          $phenotype++;
          $self->add_xref(
                           { acc        => $number,
                             label      => $label,
                             desc       => $long_desc,
                             source_id  => $gene_source_id,
                             species_id => $species_id,
                             dbi        => $dbi,
                             info_type  => "DEPENDENT" } );

          $self->add_xref(
                           { acc        => $number,
                             label      => $label,
                             desc       => $long_desc,
                             source_id  => $morbid_source_id,
                             species_id => $species_id,
                             dbi        => $dbi,
                             info_type  => "DEPENDENT" } );
        }
        elsif ( $type eq q{^} ) {
          my ( $new_number ) = ( $long_desc =~ m{
                                                  MOVED\sTO\s
                                                  (\d+)
                                              }msx );
          if ( defined $new_number ) {
            if ( $new_number ne $number ) {
              $old_to_new{$number} = $new_number;
            }
          }
          # Both leading and trailing whitespace has been removed
          # so don't bother with another regex match, just compare.
          elsif ( $long_desc eq 'REMOVED FROM DATABASE' ) {
            $removed{$number} = 1;
            $removed_count++;
          }
          else {
            print {*STDERR} "Unsupported type of a '^' record: '${long_desc}'\n";
            return 1;
          }

        }
      } ## end if (defined $ti)
    } ## end if (defined $number)
  } ## record loop

  $mim_io->close();

  my $syn_count = 0;
  foreach my $mim ( keys %old_to_new ) {
    my $old = $mim;
    my $new = $old_to_new{$old};
    while ( defined( $old_to_new{$new} ) ) {
      $new = $old_to_new{$new};
    }
    if ( !defined( $removed{$new} ) ) {
      $self->add_to_syn_for_mult_sources( $new, \@sources, $old,
                                          $species_id, $dbi );
      $syn_count++;
    }
  }
  print "$gene genemap and $phenotype phenotype MIM xrefs added\n"
    if ($verbose);
  print "added $syn_count synonyms (defined by MOVED TO)\n"
    if ($verbose);
  return 0;    #successful
} ## end sub run


1;
