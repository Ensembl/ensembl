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

=head1 NAME

XrefParser::MGIParser

=head1 DESCRIPTION

A parser class to parse the MGI (official) source,
creating a DIRECT xref between MGI accession and ensembl mouse gene stable id ENSMUSG*

-species = mus_musculus
-species_id = 10090
-data_uri = http://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt
-file_format = TSV
-columns = [accession symbol name position chrom ens_gene_stableid] ##ignore other columns

=head1 SYNOPSIS

  my $parser = XrefParser::MGIParser->new($db->dbh);
  $parser->run({
    source_id  => 55,
    species_id => 10090,
    files      => ["MRK_ENSEMBL.rpt"],
  });
=cut

package XrefParser::MGIParser;

use strict;
use warnings;
use Carp;
use Text::CSV;

use parent qw( XrefParser::BaseParser );

=head2 run
  Arg [1]    : HashRef standard list of arguments from ParseSource
  Example    : $mgi_parser->run({ ... });
  Description: Runs the MGIParser
  Return type: 0 on success
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

  #synonyms; move this to SynonymAdaptor?!
  my $syn_hash = $self->get_ext_synonyms( 'MGI', $dbi );

  #Init input file
  my $input_file = Text::CSV->new(
    {
      sep_char           => "\t",
      empty_is_undef     => 1,
      strict             => 1,
      allow_loose_quotes => 1,
    }
  ) or confess "Cannot use file $file: " . Text::CSV->error_diag();

  my $count     = 0;
  my $syn_count = 0;

  while ( my $data = $input_file->getline($file_io) ) {
    my $acc     = $data->[0];
    my $ensid   = $data->[5];

    my $xref_id = $self->add_xref(
      {
        acc        => $acc,
        version    => 0,
        label      => $data->[1],
        desc       => $data->[2],
        source_id  => $source_id,
        species_id => $species_id,
        info_type  => 'DIRECT',
        dbi        => $dbi,
      }
    );

    $self->add_direct_xref( $xref_id, $ensid, 'Gene', undef, $dbi );
    if ( exists $syn_hash->{$acc} ) {
      foreach my $syn ( @{ $syn_hash->{$acc} } ) {
        $self->add_to_syn( $acc, $source_id, $syn, $species_id, $dbi );
        $syn_count += 1;
      }
    }
    $count += 1;

  }
  $input_file->eof
    || confess "Error parsing file $file: " . $input_file->error_diag();
  $file_io->close();

  if ($verbose) {
    print "$count direct MGI xrefs added\n";
    print $syn_count. " synonyms added\n";
  }
  return 0;

}

1;
