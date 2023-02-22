=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

XrefParser::WormBaseCuratedSynonymsParser

=head1 DESCRIPTION

This parser will read and create dependent xrefs from a simple
tab delimited file.

=head1 SYNOPSIS

  my $parser = XrefParser::WormBaseCuratedSynonymsParser->new($db->dbh);
  $parser->run({
    source_id  => 11,
    species_id => 9606,
    files      => [ "/nfs/production/flicek/wormbase/parasite/data/synonyms/strongyloides_stercoralis_prjeb528/curated_synonyms.tsv" ],
  });

=cut

package XrefParser::WormBaseCuratedSynonymsParser;

use strict;
use warnings;
use Carp;
use File::Basename;

use XrefParser::BaseParser;

use base qw( XrefParser::BaseParser );

sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files)){
    croak "Need to pass source_id, species_id and files as pairs";
  }

  my %src_ids;
  my $sth = $self->dbi()->prepare("SELECT xref_id FROM xref WHERE accession=? AND source_id=? AND species_id=$species_id");
  $src_ids{curated_gene_synonyms} = $self->get_source_id_for_source_name("curated_gene_synonyms");
  my $data = $self->get_data(@$files);
  for my $gene_id (keys %$data){
    for my $synonym (@{$data->{$gene_id}->{synonyms}}) {
        $self->add_xref_and_direct_xref(
            $sth, $species_id, "gene", $src_ids{curated_gene_synonyms},
            $gene_id, $synonym->{synonym_id}
        );
    }
  }
}

sub get_data {
  my ($self, $file) = @_;
  my $pep_io = $self->get_filehandle($file) or croak "Could not open: $file";

  my $data = {};

  while ( $_ = $pep_io->getline() ) {
    chomp;
    next if /^\/\//;
    my ($stabid, $synonym) = split(/\t/, $_);
    chomp($synonym);
    chomp($stabid);
    $data->{$stabid}->{synonyms} //=[];
    push @{$data->{$stabid}->{synonyms}}, {
        synonym_id => $synonym
    }
  }
  $pep_io->close();
  return $data;
}

sub add_xref_and_direct_xref {
  my ($self, $sth, $species_id, $object_type, $source_id,  $object_id, $label, $primary_id) = @_;
  $primary_id //= $object_id;
  return unless $label;
  $sth->execute($primary_id, $source_id);
  $self->add_direct_xref(
      ($sth->fetchrow_array())[0]
      || $self->add_xref({
           acc => $object_id,
           label => $label,
           source_id => $source_id,
           species_id => $species_id,
           info_type  => "DIRECT"
      })
  , $primary_id, $object_type, "");
}
1;
