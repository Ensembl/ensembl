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

=head1 NAME

XrefParser::JGI_ProteinParser

=head1 DESCRIPTION

Parser for JGI-1.0 protein files with gene description, FASTA format.

WARNING: this is an extremely simplistic implementation of a FASTA
parser, for instance it does not treat strings beginning with ; as
comments. As of September 2019 it (still) works for JGI data, though.

=head1 SYNOPSIS

  my $parser = XrefParser::JGI_ProteinParser->new($db->dbh);
  $parser->run({
    source_id  => 70,
    species_id => 7719,
    files      => [ "ciona.prot.fasta.gz" ],
  });

=cut

package XrefParser::JGI_ProteinParser;

# For non-destructive substitutions in regexps (/r flag)
require 5.014_000;

use strict;
use warnings;

use Carp;

use parent qw( XrefParser::BaseParser );


=head2 run

  Arg []     : HashRef standard list of arguments from ParseSource
  Example    : $jgi_parser->run({ ... });
  Description: Parse FASTA input file containing JGI-1.0 protein data,
               extract seq xrefs and add them to the xref DB
  Return type: Int; 0 upon success
  Exceptions : throws on all processing errors
  Caller     : ParseSource in the xref pipeline

=cut

sub run {
  my ( $self, $ref_arg ) = @_;

  my $source_id  = $ref_arg->{source_id};
  my $species_id = $ref_arg->{species_id};
  my $files      = $ref_arg->{files};
  my $verbose    = $ref_arg->{verbose} // 0;
  my $dbi        = $ref_arg->{dbi} // $self->dbi;

  if ( ( !defined $source_id ) or
       ( !defined $species_id ) or
       ( !defined $files ) )
  {
    confess 'Need to pass source_id, species_id and files as pairs';
  }

  my $file = @{$files}[0];

  my $file_io = $self->get_filehandle($file);
  if ( !defined $file_io ) {
    confess "Could not open $file\n";
  }
  IO::Handle->input_record_separator("\n>");

  my @xrefs;

 RECORD:
  while ( my $input_data = $file_io->getline() ) {

    my ( $accession, $sequence )
      = ( $input_data =~ m{
                            # Header line. The first record will
                            # have a > but since we use "\n>" as
                            # record separator, further ones will not
                            # contain it.
                            \A >? \s* ci0100 ( \w+? ) \n

                            # Sequence data. Can span multiple
                            # lines. Err on the side of caution and
                            # assume there CAN be records with no
                            # sequence data at all (hence the *), such
                            # records would be useless for xref
                            # generation but at least they shoudn't
                            # trigger parsing errors. By specifying
                            # "not >" as our character class we avoid
                            # having to chomp the input record.
                            ( [^>]* )
                        }msx );

    if ( !defined $accession ) {
      # Is it the file header? If so, just skip it
      if ( $input_data =~ m{ \A File: }msx ) {
        next RECORD;
      }
      # Otherwise, alert the user of parsing problems
      else {
        confess "Can't parse FASTA entry: $input_data";
      }
    }

    # Build an xref object (getting rid of whitespace from the
    # sequence in the process) and store it
    push @xrefs,
      { ACCESSION     => $accession,
        SEQUENCE      => ( $sequence =~ s{ \s }{}grmsx ),
        SOURCE_ID     => $source_id,
        SPECIES_ID    => $species_id,
        SEQUENCE_TYPE => 'peptide',
      };

  } ## end while ( my $input_data = $file_io...)

  $file_io->close();

  $self->upload_xref_object_graphs( \@xrefs, $dbi );

  if ( $verbose ) {
    print scalar(@xrefs) . " JGI_ xrefs succesfully parsed\n";
  }

  return 0;
} ## end sub run


1;
