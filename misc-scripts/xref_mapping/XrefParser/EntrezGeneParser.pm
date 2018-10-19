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

package XrefParser::EntrezGeneParser;

use strict;

use Carp;
use POSIX qw(strftime);
use File::Basename;
use Text::CSV;

use parent qw( XrefParser::BaseParser );

sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $species_name = $ref_arg->{species};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose} || 0;
  my $dbi          = $ref_arg->{dbi} || $self->dbi;

  croak "Need to pass source_id, species_id, files and rel_file as pairs"
    unless defined $source_id and defined $species_id and defined $files;

  my $file = @{$files}[0];

  my $wiki_source_id = $self->get_source_id_for_source_name("WikiGene", undef, $dbi);

  my %species_tax_id = %{$self->get_taxonomy_from_species_id($species_id, $dbi)};
  $species_tax_id{$species_id} = $species_id if defined $species_name;
  

  my $eg_io = $self->get_filehandle($file);
  croak "ERROR: Could not open $file\n" unless defined $eg_io;

  my $input_file = Text::CSV->new({ sep_char           => "\t",
				    empty_is_undef     => 1,
				    allow_loose_quotes => 1 # turns out file has unescaped quotes
				  }) or croak "Cannot use file $file: ".Text::CSV->error_diag ();

  # process header
  $input_file->column_names( @{ $input_file->getline( $eg_io ) } );

  # read data and load xrefs
  my $xref_count = 0;
  my $syn_count  = 0;
  my %seen; # record already processed xrefs

  while ( my $data = $input_file->getline_hr( $eg_io ) ) {
    next unless exists $species_tax_id{$data->{'#tax_id'}};

    my $acc = $data->{'GeneID'};
    next if $seen{$acc};

    my $symbol = $data->{'Symbol'};
    my $desc   = $data->{'description'};

    $self->add_xref({ acc        => $acc,
		      label      => $symbol,
		      desc       => $desc,
		      source_id  => $source_id,
		      species_id => $species_id,
                      dbi        => $dbi,
		      info_type  =>"DEPENDENT"} );

    $self->add_xref({ acc        => $acc,
		      label      => $symbol,
		      desc       => $desc,
		      source_id  => $wiki_source_id,
		      species_id => $species_id,
                      dbi        => $dbi,
		      info_type  => "DEPENDENT" } ); #,"From EntrezGene $acc");
    $xref_count++;

    my (@syn) = split(/\|/, $data->{'Synonyms'});
    foreach my $synonym (@syn){
      if($synonym ne "-"){
	$self->add_to_syn($acc, $source_id, $synonym, $species_id, $dbi);
	$syn_count++;
      }
    }

    $seen{$acc} = 1;
  }

  $input_file->eof or croak "Error parsing file $file: " . $input_file->error_diag();
  $eg_io->close();

  print $xref_count." EntrezGene Xrefs added with $syn_count synonyms\n" if $verbose;

  return 0; # success
}

 
1;
