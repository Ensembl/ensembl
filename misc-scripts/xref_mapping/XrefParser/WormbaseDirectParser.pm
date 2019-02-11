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

package XrefParser::WormbaseDirectParser;

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

  my $file = @{$files}[0];
  my @fields = qw/wormbase_gene wormbase_gseqname wormbase_locus wormbase_transcript wormbase_cds wormpep_id protein_id/;
  my %src_ids;
  my $sth = $self->dbi()->prepare("SELECT xref_id FROM xref WHERE accession=? AND source_id=? AND species_id=$species_id");
  for my $field (@fields){
    $src_ids{$field} = $self->get_source_id_for_source_name($field); 
  }
  my $data = $self->get_data(@$files);
  for my $gene_id (keys %$data){
    $self->add_xref_and_direct_xref(
      $sth, $species_id, "gene", $src_ids{wormbase_gene},
      $gene_id,  $gene_id
    );
    $self->add_xref_and_direct_xref(
      $sth, $species_id, "gene", $src_ids{wormbase_gseqname},
      $gene_id, $data->{$gene_id}->{wormbase_gseqname}
    );
    $self->add_xref_and_direct_xref(
      $sth, $species_id, "gene", $src_ids{wormbase_locus}, 
      $gene_id, $data->{$gene_id}->{wormbase_locus}
    );
    for my $transcript (@{$data->{$gene_id}->{transcripts}}){
      $self->add_xref_and_direct_xref(
        $sth, $species_id, "transcript", $src_ids{wormbase_transcript}, 
        $transcript->{transcript_id}, $transcript->{transcript_id}
      );
      $self->add_xref_and_direct_xref(
        $sth, $species_id, "transcript", $src_ids{wormbase_cds},
        $transcript->{wormbase_cds}, $transcript->{wormbase_cds}, $transcript->{transcript_id}
      );
      $self->add_xref_and_direct_xref(
        $sth, $species_id, "translation", $src_ids{wormpep_id}, 
        $transcript->{wormpep_id}, $transcript->{wormpep_id}, $transcript->{transcript_id}
      );
      $self->add_xref_and_direct_xref(
        $sth, $species_id, "translation", $src_ids{protein_id}, 
        $transcript->{protein_id}, $transcript->{protein_id}, $transcript->{transcript_id}
      );
    } 
  }
}
sub get_data {
  my ($self, $file) = @_;
  my $pep_io = $self->get_filehandle($file) or croak "Could not open: $file";

  my $data = {};

  while ( $_ = $pep_io->getline() ) {
    next if /^\/\//;
    my ($gseqid, $wbgeneid, $locus, $wbtranscript, $wormpep, $insdc_parent, $insdc_locus_tag, $protein_id, $uniprot_id) = split(/\t/, $_);

    $data->{$wbgeneid}->{transcripts} //=[];
    push @{$data->{$wbgeneid}->{transcripts}}, {
      transcript_id => $wbtranscript,
      ($wormpep ne '.' && $wbtranscript =~ /^(.*?)(\.\d+)?$/ ? (wormbase_cds => $1 ) : ()),
      ($wormpep ne '.' ? (wormpep_id => $wormpep) : ()),
      ($protein_id ne '.' ? (protein_id => $protein_id) : ()),
    };
    $data->{$wbgeneid}->{wormbase_gseqname} = $gseqid;
    $data->{$wbgeneid}->{wormbase_locus} =  $locus if $locus ne '.'; 
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
