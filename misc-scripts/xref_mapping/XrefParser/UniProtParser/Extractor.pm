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



package XrefParser::UniProtParser::Extractor;

use strict;
use warnings;

# For non-destructive substitutions in regexps (/r flag)
require 5.014_000;

use Carp;
use List::Util;
use Readonly;
use charnames ':full';


# Note that care must be taken when adding new prefixes to this list
# because some of them - for instance the Rx family of fields,
# describing publications - are not compatible with the current way of
# processing.
# Syntax: 1 is mandatory field, 0 - an optional one
Readonly my %prefixes_of_interest
  => (
      'ID'  => 1,
      'AC'  => 1,
      'DE'  => 1,
      'GN'  => 0,
      'OX'  => 1,
      'CC'  => 0,
      'DR'  => 1,
      'PE'  => 1,
      'SQ'  => 1,
      q{  } => 1,
    );

# FIXME: this should probably be combined with
# Transformer::%taxonomy_ids_from_taxdb_codes to make sure
# database qualifiers stay in sync
Readonly my %supported_taxon_database_qualifiers
  => (
      'NCBI_TaxID' => 1,
    );



# FIXME: at the moment the extractor only handles one input file at a
# time for backwards compatibility with the old UniProtParser, that
# said there is no reason for it not to be able to handle multiple
# files. Should we want to support this, we should stop opening the
# filehandle in the constructor because it would no longer be a fixed
# property of the object.
sub new {
  my ( $proto, $arg_ref ) = @_;
  my $file_names = $arg_ref->{'file_names'};
  my $baseParserInstance = $arg_ref->{'baseParser'};

  my $filename = $file_names->[0];
  my $filehandle = $baseParserInstance->get_filehandle( $filename );
  if ( !defined $filehandle ) {
    croak "Failed to acquire a file handle for '${filename}'";
  }

  # Keep the file name for possible debugging purposes, unless we can
  # somehow retrieve it from _io_handle
  my $self = {
              'input_name' => $filename,
              '_io_handle' => $filehandle,
            };
  my $class = ref $proto || $proto;
  bless $self, $class;
  return $self;
}


sub close_input {
  my ( $self ) = @_;

  $self->{'_io_handle'}->close();

  return;
}


sub extract {
  my ( $self ) = @_;

  $self->_record_has_all_needed_fields();

  my $entry_object
    = {
       'accession_numbers' => $self->_get_accession_numbers(),
       'comments'          => $self->_get_comments(),
       'crossreferences'   => $self->_get_database_crossreferences(),
       'description'       => $self->_get_description(),
       'gene_names'        => $self->_get_gene_names(),
       'quality'           => $self->_get_quality(),
       'sequence'          => $self->_get_sequence(),
       'taxon_codes'       => $self->_get_taxon_codes(),
     };

  return $entry_object;
}


sub get_uniprot_record {
  my ( $self ) = @_;

  my $io = $self->{'_io_handle'};
  my $uniprot_record = {};

 INPUT_LINE:
  while ( my $file_line = $io->getline() ) {
    chomp $file_line;

    my ( $prefix, $content )
      = ( $file_line =~ m{ \A
                           ([A-Z /]{2})  # prefix
                           (?:
                             \s{3}       # leading spaces are important in e.g. DE, FT
                             (.+)        # content
                           )?            # end-of-record line will not have any of this
                           \z
                       }msx );
    if ( ! defined $prefix ) {
      croak 'Malformed prefix';
    }

    if ( $prefix eq q{//} ) {
      # End of record, return what we have got so far
      $self->{'record'} = $uniprot_record;
      return 1;
    }

    # Do not waste time and memory on fields we do not need
    if ( ! exists $prefixes_of_interest{$prefix} ) {
      next INPUT_LINE;
    }

    if ( ! exists $uniprot_record->{$prefix} ) {
      $uniprot_record->{$prefix} = [];
    }
    push @{ $uniprot_record->{$prefix} }, $content;

  }

  # If we began parsing fields but have never reached the //,
  # something is very wrong
  if ( scalar keys %{ $uniprot_record } > 0 ) {
    croak "Incomplete input record";
  }

  # EOF
  return 0;
}



# Parse the AC fields of the current record and produce a list of
# UniProt accession numbers. The list will reflect the order in which
# accession numbers appeared in the record, which as of October 2018
# is: primary accession number first, then all the secondary ones in
# alphanumerical order.
sub _get_accession_numbers {
  my ( $self ) = @_;

  my $ac_fields = $self->{'record'}->{'AC'};
  my @numbers
    = split( qr{\s* ; \s*}msx, join( q{}, @{ $ac_fields } ) );
  # FIXME: we should probably make this persist until a new record has
  # been loaded

  return \@numbers;
}


# Parse the CC fields of the current record and produce a hash mapping
# comments to their respective topics. Optionally, the returned data
# can be restricted to specific topics.
sub _get_comments {
  my ( $self, @topics_to_return ) = @_;

  my $cc_fields = $self->{'record'}->{'CC'};

  # CC is an optional field
  if ( ! defined $cc_fields ) {
    return {};
  }

  # FIXME: we should probably make this persist until a new record has
  # been loaded
  my %comments_by_topic;

  # FIXME: get rid of extra whitespace from the beginning of
  # continuation lines
  my @topic_lines = split( qr{\s* -!- \s*}msx,
                           join( q{ }, @{ $cc_fields } )
                        );
  foreach my $line ( @topic_lines ) {

    my ($topic, $content)
      = ( $line =~ m{
                      \A

                      # Topic name: one or more words in ALL CAPS
                      ( [A-Z ]+ )

                      # Separator
                      \s*
                      :
                      \s*

                      # Everything else is content
                      ( .+ )
                  }msx );

    # This will bypass the first empty line of @topic_lines
    if ( defined $topic ) {
      # FIXME: can a record have multiple entries with the same topic???
      $comments_by_topic{$topic} = $content;
    }
  }

  # If we haven't been asked for specific topics, return all comments.
  if ( scalar @topics_to_return == 0 ) {
    return \%comments_by_topic;
  }

  my %comments_of_interest;
  foreach my $topic ( @topics_to_return ) {
    $comments_of_interest{$topic} = $comments_by_topic{$topic};
  }
  return \%comments_of_interest;
}


# Parse the DR fields of the current record, break them into
# constituent parts and produce a list (or to be precise, a hash of
# arrayrefs) of database cross-references grouped by reference.
sub _get_database_crossreferences {
  my ( $self ) = @_;

  # Use named sequences for square brackets to avoid excessive
  # escaping as well as for better readability.
  Readonly my $isoform_field_pattern
    => qr{
           \s*
           \N{LEFT SQUARE BRACKET}
           \s*
           ( [^\N{RIGHT SQUARE BRACKET}]+ )
           \s*
           \N{RIGHT SQUARE BRACKET}
           \s*
       }msx;

  my $dr_fields = $self->{'record'}->{'DR'};

  # FIXME: we should probably make this persist until a new record has
  # been loaded
  my $crossreferences = {};

  foreach my $dr_line ( @{ $dr_fields } ) {
    my ( $res_abbrev, $res_id, @opts ) = split( qr{ ;\s* }msx, $dr_line);

    my ( $last_opt, $isoform )
      = ( $opts[-1] =~ m{
                          ( .+ )  # will grab all dots but the last one
                          [.]
                          (?:
                            $isoform_field_pattern
                          )?
                          \z
                      }msx );
    if ( ! defined $last_opt ) {
      croak "Mailformed final-option match in:\n\t$dr_line";
    }

    # At the very least, strips the trailing dot
    $opts[-1] = $last_opt;

    my $crossref
      = {
         'id'            => $res_id,
         'optional_info' => \@opts,
       };
    if ( defined $isoform ) {
      $crossref->{'target_isoform'} = $isoform;
    }

    # There CAN be multiple cross-references to a database
    push @{ $crossreferences->{$res_abbrev} }, $crossref;
  }

  return $crossreferences;
}


# Parse the DE fields of the current record and produce a description
# string compatible with the output of the old UniProtParser, namely:
#  - the description begins with a semicolon-separated list of
#    top-level names (i.e. ones not belonging to a Contains or
#    Includes section), in the order they appear in the record;
#  - this list is followed by a space and a space-separated list of
#    names from Contains and Includes sections, again in the order
#    they appear in the record;
#  - we process both RecNames and SubNames, and both types are
#    considered of equal priority (i.e. ultimately it is their order
#    that matters);
#  - in either case we only consider full names;
#  - evidence codes, PubMed references etc. are discarded.
# Note that unlike most field parsers implemented so far, this one
# does NOT attempt to fully process the syntax of DE fields; this is
# in order to avoid messing with whitespace-defined context only to
# discard most of the field data anyway. Keep this in mind should you
# want to extend the parsing to e.g. extract EC numbers, in which case
# you will likely have to implement a full parser.
sub _get_description {
  my ( $self ) = @_;

  my $de_fields = $self->{'record'}->{'DE'};

  Readonly my $description_name_value
    => qr{
           ( [^;\N{LEFT CURLY BRACKET}]+ )
           # FIXME: the match will fail if there is no
           # whitespace before the left curly bracket
           (?: ; | \s+\N{LEFT CURLY BRACKET} )
       }msx;

  my @names;
  my @subdescs;

  foreach my $line ( @{ $de_fields } ) {
    my ( $indent, $content )
      = ( $line =~ m{
                      \A
                      ( \s* )  # FIXME: explain
                      (?:
                        RecName | SubName
                      )
                      :
                      \s*
                      Full=
                      $description_name_value
                  }msx );
    if ( defined $indent ) {
      if ( $indent eq q{} ) {
        push @names, $content;
      }
      else {
        push @subdescs, $content;
      }
    }
  }

  my $description
    = join( q{ }, (
                   join( q{;}, @names ),
                   @subdescs
                 ) );
  # Make sure we do not return an empty string
  return ( length( $description ) > 0 ) ? $description : undef;
}


# Parse the GN fields of the current record and produce a structured
# list of names of genes coding the protein sequence in question.
sub _get_gene_names {
  my ( $self ) = @_;

  my $gn_fields = $self->{'record'}->{'GN'};

  # GN is an optional field
  if ( ! defined $gn_fields ) {
    return [];
  }

  my @gn_text_entries;

  # Gene-name entries can span multiple GN names but we cannot simply
  # concatenate them all because there can be multiple gene names per
  # record. In order to merge data correctly we must look for the
  # name-separator line.
  my $current_entry = q{};
  foreach my $line ( @{ $gn_fields } ) {
    # This is what the separator line looks like
    if ( $line eq 'and' ) {
      push @gn_text_entries, $current_entry;
      $current_entry = q{};
    }
    else {
      # No need for extra spaces here
      $current_entry .= $line;
    }
  }
  # Make sure the last entry makes it in as well
  push @gn_text_entries, $current_entry;

  my $gene_names = [];
  foreach my $entry ( @gn_text_entries ) {
    my $parsed_entry = {};

    my @entry_captures = ( $entry =~ m{
                                        \s*
                                        ( [^=]+ )
                                        \s*
                                        =
                                        \s*
                                        ( [^;]+ )
                                        \s*
                                        ;
                                    }gmsx );

    while ( my ( $key, $value ) = splice( @entry_captures, 0, 2 ) ) {
      my @split_value = split( qr{ \s*,\s* }msx, $value );

      $parsed_entry->{$key}
        = ( $key eq 'Name' ) ? $value : \@split_value;
    }

    # UniProt-KB User Manual states a "Synonyms" token can only be
    # present if there is a "Name" token.
    if ( ( exists $parsed_entry->{'Synonyms'} )
         && ( ! exists $parsed_entry->{'Name'} ) ) {
      croak "Malformed input: found 'Synonyms' but no 'Name' in:\n\t$entry";
    }

    push @{ $gene_names }, $parsed_entry;
  }

  return $gene_names;
}


# Obtain quality information for the current record. This consists of
# two parts: status (i.e. whether the entry has been reviewed or not)
# from the ID line and evidence level from the PE line.
sub _get_quality {
  my ( $self ) = @_;

  Readonly my $id_status_field
    => qr{
           (?: Unreviewed )
         | (?: Reviewed )
       }msx;

  # These is only one ID line
  my $id_line = $self->{'record'}->{'ID'}->[0];
  my ( $entry_status )
    = ( $id_line =~ m{
                       \A

                       # UniProt name
                       [0-9A-Z_]+

                       \s+

                       ( $id_status_field )
                       \s*
                       ;
                   }msx );
  if ( ! defined $entry_status ) {
    croak "Invalid entry status in:\n\t$id_line";
  }

  # Likewise, these is only one PE line
  my $pe_line = $self->{'record'}->{'PE'}->[0];
  my ( $evidence_level )
    = ( $pe_line =~ m{
                       \A

                       ( [1-4] )
                       \s*
                       :
                   }msx );
  if ( ! defined $evidence_level ) {
    croak "Invalid protein evidence level in:\n\t$pe_line";
  }

  return {
          'status'         => $entry_status,
          'evidence_level' => $evidence_level,
        };
}


# Parse the sequence ('  ') fields of the current record and produce
# its polypeptide sequence. as a continuous string i.e. without
# decorative whitespace.
sub _get_sequence {
  my ( $self ) = @_;

  my $sequence_fields = $self->{'record'}->{ q{  } };

  # Concatenate the sequence into a single continuous string. We do
  # not expect to see more than whitespace at a time so instead of
  # trying to match as long a string of them as possible in order to
  # minimise the number of independent substitutions, we always match
  # on one in order to avoid unnecessary backtracking.
  # Note that we use non-destructive substitution.
  my $sequence = ( join( q{}, @{ $sequence_fields } ) =~ s{ \s }{}grmsx );

  # We could in principle directly return substitution result but then
  # we would have to make sure we always call get_sequence() in scalar
  # context. Safer to simply return a scalar instead.
  return $sequence;
}



# Parse the OX field of the current record and produce a list of
# database_qualifier/taxon_code pairs.
sub _get_taxon_codes {
  my ( $self ) = @_;

  # There is only one OX line
  my $ox_line = $self->{'record'}->{'OX'}->[0];

  # On the one hand, according to UniProt-KB User Manual from October
  # 2018 there should only be a single taxon code per OX line and the
  # current SwissProt data file doesn't contain any records which
  # violate this. On the other hand we haven't checked this in the
  # TrEMBL file because of its size and the old UniProtParser did have
  # support for synonyms (albeit using different syntax -
  # "db=id1, id2, ...;" rather than "db1=id1; db2=id2; ...").
  # The code below assumes there might be multiple entries present, if
  # you want to force it to only ever look for one simply drop the /g
  # modifier from the regex match.

  Readonly my $taxon_db_entry
    => qr{
           # Database qualifier. Chances are the list of
           # allowed characters will change should DBs
           # other than NCBI ever become supported here.
           ( [A-Za-z_]+ )

           \s*  # just in case
           =
           \s*  # same

           # Taxon ID. This is almost certainly NCBI-specific.
           ( [0-9]+ )
       }msx;
  Readonly my $evidence_code_list
    => qr{
           # As of October 2018, this syntax is not declared in
           # UniProt-KB User Manual yet frequently encountered in data
           # files.  Use named sequences for curly brackets to avoid
           # excessive escaping as well as for better readability.
           \N{LEFT CURLY BRACKET}
           \s*
           [^\N{RIGHT CURLY BRACKET}]+
           \s*
           \N{RIGHT CURLY BRACKET}
           \s*
       }msx;
  my @ox_captures
    = ( $ox_line =~ m{
                       $taxon_db_entry
                       \s*

                       # Optional things (e.g. evidence codes) we do
                       # not need now but must parse on the off chance
                       # someone decides e.g. that curly braces
                       # constitute quotes so it is okay to have a
                       # semicolon between them
                       (?:
                         $evidence_code_list
                         | [^;]+
                       )?

                       # End of record
                       ;
                   }gmsx );

  my @extracted_taxon_codes;

  # Reminder: if you change the number of capture groups above
  # remember to adapt both the variable list *and* the third argument
  # to splice().
 TAXON_ENTRY:
  while ( my ( $db_qualifier, $taxon_code ) = splice( @ox_captures, 0, 2 ) ) {

    if ( ( ! defined $db_qualifier )
         || ( ! exists $supported_taxon_database_qualifiers{$db_qualifier} ) ) {
      # Abort on malformed or new database qualifiers
      croak "Cannot use taxon-DB qualifier '${db_qualifier}'";
    }
    elsif ( ! $supported_taxon_database_qualifiers{$db_qualifier} ) {
      # Known but of no interest. Ignore it.
      next TAXON_ENTRY;
    }

    if ( ! defined $taxon_code ) {
      croak "Failed to extract taxon code from:\n\t${ox_line}";
    }

    # FIXME: further processing?

    push @extracted_taxon_codes, {
                                  'db_qualifier' => $db_qualifier,
                                  'taxon_code'     => $taxon_code,
                                };

  }

  return \@extracted_taxon_codes;
}


sub _record_has_all_needed_fields {
  my ( $self ) = @_;

  # Only check mandatory fields
  my @needed_fields = grep { $prefixes_of_interest{$_} } %prefixes_of_interest;

  my $has_all
    = List::Util::all { exists $self->{'record'}->{$_} } @needed_fields;
  if ( ! $has_all ) {
    croak 'One or more required fields missing in record';
 }

  return;
}


1;
