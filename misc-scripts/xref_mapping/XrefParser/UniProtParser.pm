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

# Parse UniProt (SwissProt & SPTrEMBL) files to create xrefs.
#
# Files actually contain both types of xref, distinguished by ID line;
#
# ID   CYC_PIG                 Reviewed;         104 AA.  Swissprot
# ID   Q3ASY8_CHLCH            Unreviewed;     36805 AA.  SPTrEMBL



package XrefParser::UniProtParser;

use strict;
use warnings;
use Carp;
use POSIX qw(strftime);
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
  my $sp_source_id = $self->get_source_id_for_source_name('Uniprot/SWISSPROT','sequence_mapped', $dbi);
  push @source_ids, $sp_source_id;
  my $sptr_source_id = $self->get_source_id_for_source_name('Uniprot/SPTREMBL', 'sequence_mapped', $dbi);
  push @source_ids, $sptr_source_id;
  my $sptr_non_display_source_id = $self->get_source_id_for_source_name('Uniprot/SPTREMBL', 'protein_evidence_gt_2', $dbi);
  push @source_ids, $sptr_non_display_source_id;
  my $sp_direct_source_id = $self->get_source_id_for_source_name('Uniprot/SWISSPROT', 'direct', $dbi);
  push @source_ids, $sp_direct_source_id;
  my $sptr_direct_source_id = $self->get_source_id_for_source_name('Uniprot/SPTREMBL', 'direct', $dbi);
  push @source_ids, $sptr_direct_source_id;
  my $isoform_source_id = $self->get_source_id_for_source_name('Uniprot_isoform');
  push @source_ids, $isoform_source_id;

  my $get_source_version_sql = "SELECT source_release FROM source WHERE source_id = ?";
  my $get_source_version_sth = $xref_source->prepare($get_source_version_sql);
  my ($release_version);
  foreach my $id (@source_ids) {
    $get_source_version_sth->execute($id);
    $get_source_version_sth->bind_columns(\$release_version);
    $self->set_release( $id, $release_version, $dbi);
  }

  my $get_xref_sql = "SELECT xref_id, accession, version, label, description, info_type ".
  "FROM xref WHERE species_id = ? AND source_id = ?";
  my $get_xref_sth = $xref_source->prepare($get_xref_sql);
  my $get_dependent_sql = "SELECT x.xref_id, x.accession, x.source_id, x.species_id, dx.linkage_source_id FROM xref x, dependent_xref dx ".
  "WHERE dx.dependent_xref_id = x.xref_id and dx.master_xref_id = ?";
  my $get_dependent_sth = $xref_source->prepare($get_dependent_sql);
  my $get_sequence_sql = "SELECT sequence, sequence_type, status FROM primary_xref WHERE xref_id = ?";
  my $get_sequence_sth = $xref_source->prepare($get_sequence_sql);
  my $get_synonym_sql = "SELECT synonym FROM synonym WHERE xref_id = ?";
  my $get_synonym_sth = $xref_source->prepare($get_synonym_sql);
  my $get_direct_sql = "SELECT ensembl_stable_id, linkage_xref FROM translation_direct_xref WHERE general_xref_id = ?";
  my $get_direct_sth = $xref_source->prepare($get_direct_sql);
  my ($xref_id, $accession, $version, $label, $description, $info_type, $parsed_seq, $type, $status, $dep_xref_id, $dep_accession, $dep_source_id, $dep_species_id, $linkage_source_id, $synonym, $stable_id, $linkage_xref);

  my @xrefs;
  my $count = 0;

  foreach my $xref_source_id (@source_ids) {
    $get_xref_sth->execute($species_id, $xref_source_id);
    $get_xref_sth->bind_columns(\$xref_id, \$accession, \$version, \$label, \$description, \$info_type);
    while ($get_xref_sth->fetch) {
      my $xref = {};
      $count++;
      $xref->{ACCESSION} = $accession;
      $xref->{LABEL} = $label;
      $xref->{VERSION} = $version;
      $xref->{SPECIES_ID} = $species_id;
      $xref->{SEQUENCE_TYPE} = 
      $xref->{INFO_TYPE} = $info_type;
      $xref->{SOURCE_ID} = $xref_source_id;
      $xref->{DESCRIPTION} = $description;

      # Add sequence if there is some
      $get_sequence_sth->execute($xref_id);
      $get_sequence_sth->bind_columns(\$parsed_seq, \$type, \$status);
      while ($get_sequence_sth->fetch) {
	$xref->{SEQUENCE_TYPE} = $type;
        $xref->{STATUS} = $status;
	$xref->{SEQUENCE} = $parsed_seq;
      }

      # Look for synonyms
      $get_synonym_sth->execute($xref_id);
      $get_synonym_sth->bind_columns(\$synonym);
      while ($get_synonym_sth->fetch) {
        push (@{$xref->{SYNONYMS} }, $synonym);
      }

      # Look for direct xref
      $get_direct_sth->execute($xref_id);
      $get_direct_sth->bind_columns(\$stable_id, \$linkage_xref);
      while ($get_direct_sth->fetch) {
        my %direct;
	my $isoform;
	$direct{STABLE_ID} = $stable_id;
	$direct{ENSEMBL_TYPE} = 'Translation';
	$direct{LINKAGE_TYPE} = $linkage_xref;
	$direct{SOURCE_ID} = $xref_source_id;
        push (@{$xref->{DIRECT_XREFS}}, \%direct);
      }

      #Add any dependent xrefs
      $get_dependent_sth->execute($xref_id);
      $get_dependent_sth->bind_columns(\$dep_xref_id, \$dep_accession, \$dep_species_id, \$dep_source_id, \$linkage_source_id);
      my %dep;
      while ($get_dependent_sth->fetch) {
        if ($dep_species_id != $species_id) { next; }
	$dep{ACCESSION} = $dep_accession;
	$dep{LABEL} = $dep_accession;
	$dep{SOURCE_ID} = $dep_source_id;
	$dep{LINKAGE_SOURCE_ID} = $linkage_source_id;
	$get_synonym_sth->execute($dep_xref_id);
	$get_synonym_sth->bind_columns(\$synonym);
	while ($get_synonym_sth->fetch) {
          push (@{$dep{SYNONYMS} }, $synonym);
        }
	push @{$xref->{DEPENDENT_XREFS}}, \%dep;
      }

      push @xrefs, $xref;

      if ($count > 1000) {
        $self->upload_xref_object_graphs( \@xrefs, $dbi );
        $count = 0;
        undef @xrefs;
      }
    }
  }

  $self->upload_xref_object_graphs(\@xrefs, $dbi) if scalar(@xrefs) > 0;
  return 0; # successfull
}

1;
