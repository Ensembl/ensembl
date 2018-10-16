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

package XrefParser::EntrezGeneParser;

use strict;
use POSIX qw(strftime);
use File::Basename;
use Carp;
use base qw( XrefParser::BaseParser );

sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $species_name = $ref_arg->{species};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};
  my $dbi          = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

  croak "Need to pass source_id, species_id, files and rel_file as pairs"
    unless defined $source_id and defined $species_id and defined $files;

  $verbose |=0;

  my $file = @{$files}[0];

  my $wiki_source_id = $self->get_source_id_for_source_name("WikiGene", undef, $dbi);

  my %species_tax_id = %{$self->get_taxonomy_from_species_id($species_id, $dbi)};
  $species_tax_id{$species_id} = $species_id if defined $species_name;
  

  my $eg_io = $self->get_filehandle($file);
  croak "ERROR: Could not open $file\n" unless defined $eg_io;

  my $head = $eg_io->getline(); # first record are the headers
  chomp $head;

  # process header to figure out correct indexes to use, in case they change
  my (@arr) = split(/\s+/, $head);

  my $gene_id_index = -2;
  my $gene_symbol_index = -2;
  my $gene_desc_index = -2;
  my $gene_tax_id_index = -2;
  my $gene_synonyms_index = -2;

  foreach (my $i=0; $i<= $#arr; $i++){
    # Format change by Entrez, first header element is #tax_id
    if($arr[$i] eq "#tax_id") {
      $gene_tax_id_index = $i;
    } elsif($arr[$i] eq "GeneID") {
      $gene_id_index = $i;
    } elsif($arr[$i] eq "Symbol") {
      $gene_symbol_index = $i;
    } elsif($arr[$i] eq "description") {
      $gene_desc_index = $i;
    } elsif($arr[$i] eq "Synonyms") {
      $gene_synonyms_index = $i;
    }
  }
  if( $gene_id_index       == -2 ||
      $gene_symbol_index   == -2 ||
      $gene_desc_index     == -2 ||
      $gene_synonyms_index == -2 ||
      $gene_tax_id_index == -2) {
    croak "HEADER\n$head\n\n" .
	"Unable to get all the indexes needed\n" .
	"gene_id = $gene_id_index\n" .
	"tax_id = $gene_tax_id_index\n" .
	"symbol = $gene_symbol_index\n" .
	"desc = $gene_desc_index\n" .
	"synonyms = $gene_synonyms_index";
  }

  my $xref_count = 0;
  my $syn_count  = 0;
  my %seen; # record already processed xrefs

  while ( $_ = $eg_io->getline() ) {
    chomp;
    my (@arr) = split(/\t/, $_);
    next unless defined $species_tax_id{$arr[$gene_tax_id_index]};

    my $acc = $arr[$gene_id_index];
    next if $seen{$acc};

    my $symbol = $arr[$gene_symbol_index];
    my $desc   = $arr[$gene_desc_index];

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

    my (@syn) = split(/\|/ ,$arr[$gene_synonyms_index]);
    foreach my $synonym (@syn){
      if($synonym ne "-"){
	$self->add_to_syn($acc, $source_id, $synonym, $species_id, $dbi);
	$syn_count++;
      }
    }

    $seen{$acc} = 1;
  }

  $eg_io->close();

  print $xref_count." EntrezGene Xrefs added with $syn_count synonyms\n" if $verbose;

  return 0; #successful
}

 
1;
