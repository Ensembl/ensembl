=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

# Parse RefSeq data from central database to create species specific xrefs.

package XrefParser::RefSeqGPFFParser;

use strict;
use warnings;
use Carp;
use File::Basename;

use base qw( XrefParser::BaseParser );

sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $dbi          = $ref_arg->{dbi};
  my $xref_source  = $ref_arg->{xref_source};

  if((!defined $source_id) or (!defined $species_id) or (!defined $xref_source)){
    croak "Need to pass source_id, species_id and xref_source";
  }

  my @source_ids;
  my $peptide_source_id =
    $self->get_source_id_for_source_name('RefSeq_peptide', undef, $dbi);
  push @source_ids, $peptide_source_id;
  my $mrna_source_id =
    $self->get_source_id_for_source_name('RefSeq_mRNA','refseq', $dbi);
  push @source_ids, $mrna_source_id;
  my $ncrna_source_id =
    $self->get_source_id_for_source_name('RefSeq_ncRNA', undef, $dbi);
  push @source_ids, $ncrna_source_id;

  my $pred_peptide_source_id =
    $self->get_source_id_for_source_name('RefSeq_peptide_predicted', undef, $dbi);
  push @source_ids, $pred_peptide_source_id;
  my $pred_mrna_source_id =
    $self->get_source_id_for_source_name('RefSeq_mRNA_predicted','refseq', $dbi);
  push @source_ids, $pred_mrna_source_id;
  my $pred_ncrna_source_id =
    $self->get_source_id_for_source_name('RefSeq_ncRNA_predicted', undef, $dbi);
  push @source_ids, $pred_ncrna_source_id;

  my $entrez_source_id = $self->get_source_id_for_source_name('EntrezGene', undef, $dbi);
  push @source_ids, $entrez_source_id;
  my $wiki_source_id = $self->get_source_id_for_source_name('WikiGene', undef, $dbi);
  push @source_ids, $wiki_source_id;

  # Retrieve existing NCBIGene xrefs
  my (%entrez)     = %{$self->get_acc_to_label("EntrezGene",$species_id, undef, $dbi)};

  my $get_source_version_sql = "SELECT source_release FROM source WHERE source_id = ?";
  my $get_source_version_sth = $xref_source->prepare($get_source_version_sql);
  my ($release_version);
  foreach my $id (@source_ids) {
    $get_source_version_sth->execute($id);
    $get_source_version_sth->bind_columns(\$release_version);
    $self->set_release( $id, $release_version, $dbi);
  }

  my $get_xref_sql = "SELECT x.accession, x.version, x.label, x.description, x.info_type, ".
  "px.sequence, px.sequence_type, p.accession2, ".
  "x2.accession, x2.species_id, dx.linkage_source_id ".
  "FROM primary_xref px, pairs p, xref x2, xref x, dependent_xref dx ".
  "WHERE dx.master_xref_id = x.xref_id AND x.xref_id = px.xref_id AND x.accession = p.accession1 AND dx.dependent_xref_id = x2.xref_id ".
  "AND x.species_id = ? AND x.source_id = ?";
  my $get_xref_sth = $xref_source->prepare($get_xref_sql);
  my ($accession, $version, $label, $description, $xref_source_id, $info_type, $parsed_seq, $type, $refseq_pair, $dep_accession, $dep_species_id, $linkage_source_id);

  my @xrefs;
  my $count = 0;
  foreach my $xref_source_id (@source_ids) {
    $get_xref_sth->execute($species_id, $xref_source_id);
    $get_xref_sth->bind_columns(\$accession, \$version, \$label, \$description, \$info_type, \$parsed_seq, \$type, \$refseq_pair, \$dep_accession, \$dep_species_id, \$linkage_source_id);
    while ($get_xref_sth->fetch()) {
      my $xref = {};
      $count++;
  
      $xref->{ACCESSION} = $accession;
      $xref->{VERSION} = $version;
      $xref->{LABEL} = $label;
      $xref->{DESCRIPTION} = $description;
      $xref->{SOURCE_ID} = $xref_source_id;
      $xref->{SEQUENCE} = $parsed_seq;
      $xref->{SEQUENCE_TYPE} = $type;
      $xref->{SPECIES_ID} = $species_id;
      $xref->{INFO_TYPE} = $info_type;
      $xref->{PAIR} = $refseq_pair;
  
      my %dep;
      my %dep2;
      if (defined $entrez{$dep_accession}) {
        if ($dep_species_id != $species_id) { next; }
        $dep{SOURCE_ID} = $entrez_source_id;
        $dep{LINKAGE_SOURCE_ID} = $linkage_source_id;
        $dep{ACCESSION} = $dep_accession;
        $dep{LABEL} = $entrez{$dep_accession};
        push @{$xref->{DEPENDENT_XREFS}}, \%dep;
  
        my %dep2;
        $dep2{SOURCE_ID} = $wiki_source_id;
        $dep2{LINKAGE_SOURCE_ID} = $linkage_source_id;
        $dep2{ACCESSION} = $dep_accession;
        $dep2{LABEL} = $entrez{$dep_accession};
        push @{$xref->{DEPENDENT_XREFS}}, \%dep2;
      }
      push @xrefs, $xref;
  
      if ($count > 1000) {
        $self->upload_xref_object_graphs( \@xrefs, $dbi );
        $count = 0;
        undef @xrefs;
      }
    }
  }
  $get_xref_sth->finish();

  $self->upload_xref_object_graphs(\@xrefs, $dbi) if scalar(@xrefs) > 0;

  return 0; # successful

}

1;
