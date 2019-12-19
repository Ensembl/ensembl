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

use parent qw( XrefParser::BaseParser );


my $QR_TI_FIELD_TERMINATORS
  = qr{
         (?:                # The TI field spans from *FIELD* TI until:
           [*]FIELD[*]      #  - the next field in same record, or
         | [*]RECORD[*]     #  - the end of current record, or
         | [*]THEEND[*]     #  - the end of input file
         )
     }msx;



=head2 run

  Arg [1]    : HashRef standard list of arguments from ParseSource
  Example    : $omim_parser->run({ ... });
  Description: Extract Online Mendelian Inheritance in Man entries
               from a text file downloaded from the OMIM Web site,
               then insert corresponding xrefs into the xref
               database. Note that all the xrefs produced by this
               parser are unmapped and tagged as such; links of
               appropriate type will be inserted by Mim2GeneParser.

               OMIM entries can either represent a unique locus,
               describe a disorder, or both. In Ensembl these are
               assigned, respectively, to: the source MIM_GENE,
               the source MIM_MORBID, or independently into both.

               OMIM records are multiline. Each record begins with a
               specific tag line and consists of a number of
               fields. Each field starts with its own start-tag line
               (i.e. the data proper only appears after a newline) and
               continues until the beginning of either the next field
               in the same record, the next record, or the
               end-of-input tag. The overall structure looks as
               follows:

                 *RECORD*
                 *FIELD* NO
                 *FIELD* TI
                 *FIELD* TX
                 ...
                 *RECORD*
                 *FIELD* NO
                 *FIELD* TI
                 ...
                 *RECORD*
                 *FIELD* NO
                 ...
                 *FIELD* CD
                 *FIELD* ED
                 *THEEND*

               All the data relevant to the parser can be found in the
               TI field.

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

  my $filename = @{$files}[0];

  my %old_to_new;
  my %removed;
  my %counters;
  my @sources;

  push @sources, $general_source_id;

  my $gene_source_id =
    $self->get_source_id_for_source_name( "MIM_GENE", undef, $dbi );
  push @sources, $gene_source_id;
  my $morbid_source_id =
    $self->get_source_id_for_source_name( "MIM_MORBID", undef, $dbi );
  push @sources, $morbid_source_id;

  my %TYPE_SINGLE_SOURCES = (
    q{*} => $gene_source_id,
    q{}  => $morbid_source_id,
    q{#} => $morbid_source_id,
    q{%} => $morbid_source_id,
  );

  if ($verbose) {
    print "sources are: " . join( ", ", @sources ) . "\n";
  }

  IO::Handle->input_record_separator('*RECORD*');

  my $mim_io = $self->get_filehandle($filename);
  if ( !defined $mim_io ) {
    confess "Failed to acquire a file handle for '${filename}'";
  }

  $mim_io->getline();    # first record is empty with *RECORD* as the
                         # record seperator

 RECORD:
  while ( my $input_record = $mim_io->getline() ) {

    my $ti = extract_ti( $input_record );
    if ( ! defined $ti ) {
      confess 'Failed to extract TI field from record';
    }

    my ( $type, $number, $long_desc ) = parse_ti( $ti );
    if ( ! defined $type ) {
      confess 'Failed to extract record type and description from TI field';
    }

    # Use the first block of text as description
    my @fields = split( qr{;;}msx, $long_desc );
    my $label = $fields[0] . " [" . $type . $number . "]";

    my $xref_object = {
      acc        => $number,
      label      => $label,
      desc       => $long_desc,
      species_id => $species_id,
      dbi        => $dbi,
      info_type  => 'UNMAPPED',
    };

    if ( exists $TYPE_SINGLE_SOURCES{$type} ) {
      my $type_source = $TYPE_SINGLE_SOURCES{$type};

      $xref_object->{'source_id'} = $type_source;
      $counters{ $type_source }++;
      $self->add_xref($xref_object);

    }
    elsif ( $type eq q{+} ) {    # both gene and phenotype

      $xref_object->{'source_id'} = $gene_source_id;
      $counters{ $gene_source_id }++;
      $self->add_xref($xref_object);

      $xref_object->{'source_id'} = $morbid_source_id;
      $counters{ $morbid_source_id }++;
      $self->add_xref($xref_object);

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
        $counters{ 'removed' }++;
      }
      else {
        confess "Unsupported type of a '^' record: '${long_desc}'\n";
      }

    }

  } ## record loop

  $mim_io->close();

  # Generate synonyms from "MOVED TO" entries
  foreach my $mim ( keys %old_to_new ) {
    my $old = $mim;
    my $new = $old_to_new{$old};

    # Some entries in the MIM database have been moved multiple times,
    # and we want each of the synonyms created this way to point to
    # the *current* accession instead of one another. Keep traversing
    # the chain of renames until we have reached the end, i.e. until
    # $new is no longer a valid key in %old_to_new.
    # FIXME: this is not entirely efficient, especially for long
    # rename chains, because the foreach loop processes every single
    # key of %old_to_new (i.e. every single "MOVED TO" entry) from
    # scratch - even though some of them might have already been
    # encountered in the process of traversing the change chains of
    # previously encountered keys. Some sort of a cache pointing each
    # of previously encountered keys to their respective final values,
    # might be in order here.
    # FIXME: If we do implement such a cache, compare performance for
    # retrieving original keys in random order vs in descending
    # numerical order. On the one hand starting with high accessions
    # will likely allow us to process rename chains from shorter to
    # longer ones, thus, maximising the use of the cache; on the other
    # there is the O(n log n) cost of sorting to take into account.
    while ( defined( $old_to_new{$new} ) ) {
      $new = $old_to_new{$new};
    }

    # With the latest value of $new no longer pointing to anything in
    # %old_to_new, we have got two options: either we have finally
    # reached an up-to-date entry number or the entry has ultimately
    # been removed from the database. See if we have logged the
    # removal, if we haven't add the synonyms - letting Ensembl figure
    # out by itself to which of the three (two???) sources the
    # relevant xrefs belong.
    if ( !defined( $removed{$new} ) ) {
      $self->add_to_syn_for_mult_sources( $new, \@sources, $old,
                                          $species_id, $dbi );
      $counters{ 'synonyms' }++;
    }
  }

  if ($verbose) {
    print $counters{ $gene_source_id } . ' genemap and '
      . $counters{ $morbid_source_id } . " phenotype MIM xrefs added\n"
      . $counters{ 'synonyms' } . " synonyms (defined by MOVED TO) added\n"
      . $counters{ 'removed' } . " entries removed\n";
  }

  return 0;
} ## end sub run


=head2 extract_ti

  Arg [1]    : String $input_record A single OMIM record
  Example    : my $ti_string = extract_ti( $omim_record );
  Description: Scan the provided record for the TI field and extract
               its contents, regardless of where in the record that
               field appears or the position of the record in the
               file.
  Return type: String
  Exceptions : none
  Caller     : MIMParser::run()
  Status     : Stable

=cut

sub extract_ti {
  my ( $input_record ) = @_;

  my ( $ti )
    = ( $input_record =~ m{
                            [*]FIELD[*]\sTI\n
                            (.+?)              # (important: NON-greedy match)
                            \n?
                            $QR_TI_FIELD_TERMINATORS
                        }msx );

  return $ti;
}



=head2 parse_ti

  Arg [1]    : String $ti Contents of a single TI field
  Example    : my ( $type_symbol, $omim_number, $description )
                 = parse_ti( $ti_string );
  Description: Extract the type symbol, the entry number and the
               description from the contents of an OMIM record's TI
               field. The description is *not* split into possible
               individual components, we do however remove line breaks
               from multiline entries.
  Return type: Array
  Exceptions : none
  Caller     : MIMParser::run()
  Status     : Stable

=cut

sub parse_ti {
  my ( $ti ) = @_;

  # Remove line breaks, making sure we do not accidentally concatenate words
  $ti =~ s{
            (?:
              ;;\n
            | \n;;
            )
        }{;;}gmsx;
  $ti =~ s{\n}{ }gmsx;

  # Extract the 'type' and the whole description
  my @captures = ( $ti =~ m{
                             \A
                             ([#%+*^]*)  # type of entry
                             (\d+)       # accession number, same as in NO
                             \s+         # normally just one space
                             (.+)        # description of entry
                         }msx );

  return @captures;
}



1;
