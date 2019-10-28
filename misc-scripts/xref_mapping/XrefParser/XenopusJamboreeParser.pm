=head1 LICENSE

See the NOTICE file distributed with this work for additional information
regarding copyright ownership.
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

=head1 NAME

XrefParser::XenopusJamboreeParser

=head1 DESCRIPTION

A parser class to parse the Xenbase source file.

-species = xenopus_tropicalis
-species_id = 8364
-data_uri = ftp://ftp.xenbase.org/pub/GenePageReports/GenePageEnsemblModelMapping.txt
-file_format = TSV
-columns = [acc label desc stable_id]


=head1 SYNOPSIS

  my $parser = XrefParser::XenopusJamboreeParser->new($db->dbh);
  $parser->run({
    source_id  => 150,
    species_id => 8364,
    files      => ["xenopusjamboree.txt"],
  });

=cut

package XrefParser::XenopusJamboreeParser;

use strict;
use warnings;

use Carp;
use Text::CSV;

use parent qw( XrefParser::BaseParser );


=head2 run
  Description: Runs the XenopusJamboreeParser
  Return type: N/A
  Caller     : internal

=cut

sub run {
  my ( $self, $ref_arg ) = @_;
  my $source_id  = $ref_arg->{source_id};
  my $species_id = $ref_arg->{species_id};
  my $files      = $ref_arg->{files};
  my $verbose    = $ref_arg->{verbose} // 0;
  my $dbi        = $ref_arg->{dbi} // $self->dbi;

  if ( ( !defined $source_id )
    or ( !defined $species_id )
    or ( !defined $files ) )
  {
    confess 'Need to pass source_id, species_id and files as pairs';
  }

  my $file = @{$files}[0];

  my $file_io = $self->get_filehandle($file);
  if ( !defined $file_io ) {
    confess "Could not open $file\n";
  }

  my $input_file = Text::CSV->new({
    sep_char       => "\t",
    empty_is_undef => 1,
  }) || confess "Cannot use file $file: " . Text::CSV->error_diag();

  my $count = 0;
  while ( my $data = $input_file->getline($file_io) ) {

    my ( $accession, $label, $desc, $stable_id ) = @{$data};

    # If there is a description, trim it a bit
    if ( defined $desc ) {
      $desc = parse_description( $desc );
    }

    if ( $label eq 'unnamed' ) {
      $label = $accession;
    }

    $self->add_to_direct_xrefs({
      stable_id  => $stable_id,
      type       => 'gene',
      acc        => $accession,
      label      => $label,
      desc       => $desc,
      dbi        => $dbi,
      source_id  => $source_id,
      species_id => $species_id,
    });
    $count++;
  }

  $input_file->eof
    || confess "Error parsing file $file: " . $input_file->error_diag();
  $file_io->close();

  if ($verbose) {
    print $count . " XenopusJamboreeParser xrefs succesfully parsed\n";
  }

  return 0;
} ## end sub run


=head2 parse_description
  Description: Extract description information from
               Xenopus downloaded file
  Return type: N/A
  Caller     : internal

=cut

sub parse_description {
  my ( $desc ) = @_;

  # Remove some provenance information encoded in the description
  $desc =~ s{ \s* \[ .* \] }{}msx;

  # Remove labels of type 5 of 14 from the description
  $desc =~ s{ , \s+\d+\s+ of \s+\d+ }{}msx;

  return $desc;
}


1;
