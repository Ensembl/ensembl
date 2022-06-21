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

# Parse RefSeq data from central database to create species specific xrefs.

package XrefParser::RefSeqDatabaseParser;

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

  my $source_name = $self->get_source_name_for_source_id($source_id, $dbi);
  my @source_ids;

  if ($source_name =~ /RefSeq_dna/) {
    my $mrna_source_id = $self->get_source_id_for_source_name('RefSeq_mRNA','refseq', $dbi);
    push @source_ids, $mrna_source_id;
    my $pred_mrna_source_id = $self->get_source_id_for_source_name('RefSeq_mRNA_predicted','refseq', $dbi);
    push @source_ids, $pred_mrna_source_id;
    my $ncrna_source_id = $self->get_source_id_for_source_name('RefSeq_ncRNA', undef, $dbi);
    push @source_ids, $ncrna_source_id;
    my $pred_ncrna_source_id = $self->get_source_id_for_source_name('RefSeq_ncRNA_predicted', undef, $dbi);
    push @source_ids, $pred_ncrna_source_id;
  } elsif ($source_name =~ /RefSeq_peptide/) {
    my $peptide_source_id = $self->get_source_id_for_source_name('RefSeq_peptide', undef, $dbi);
    push @source_ids, $peptide_source_id;
    my $pred_peptide_source_id = $self->get_source_id_for_source_name('RefSeq_peptide_predicted', undef, $dbi);
    push @source_ids, $pred_peptide_source_id;
  }

  my $entrez_source_id = $self->get_source_id_for_source_name('EntrezGene', undef, $dbi);
  my $wiki_source_id = $self->get_source_id_for_source_name('WikiGene', undef, $dbi);

  # Retrieve existing NCBIGene xrefs
  my (%entrez)     = %{$self->get_acc_to_label("EntrezGene",$species_id, undef, $dbi)};

  # Get existing mrna, entrezgene and wikigene accession => xref_id
  my (%refseq_ids, %entrez_ids, %wiki_ids, $add_dependent_xref_sth);
  if ($source_name =~ /RefSeq_peptide/) {
    (%refseq_ids) = %{ $self->get_valid_codes("RefSeq_mRNA", $species_id, $dbi) };
    (%entrez_ids) = %{ $self->get_valid_codes("EntrezGene", $species_id, $dbi) };
    (%wiki_ids) = %{ $self->get_valid_codes("WikiGene", $species_id, $dbi) };
    $add_dependent_xref_sth = $dbi->prepare("INSERT INTO dependent_xref  (master_xref_id, dependent_xref_id, linkage_source_id) VALUES (?,?, $entrez_source_id)");
  }

  my $get_xref_sql = "SELECT xref_id, accession, version, label, description, info_type ".
  "FROM xref WHERE species_id = ? AND source_id = ?";
  my $get_xref_sth = $xref_source->prepare($get_xref_sql);
  my $get_dependent_sql = "SELECT x.xref_id, x.accession, x.version, x.label, x.description, x.source_id, x.species_id, dx.linkage_source_id FROM xref x, dependent_xref dx ".
  "WHERE dx.dependent_xref_id = x.xref_id and dx.master_xref_id = ?";
  my $get_dependent_sth = $xref_source->prepare($get_dependent_sql);
  my $get_sequence_sql = "SELECT sequence, sequence_type, status FROM primary_xref WHERE xref_id = ?";
  my $get_sequence_sth = $xref_source->prepare($get_sequence_sql);
  my $get_synonym_sql = "SELECT synonym FROM synonym WHERE xref_id = ?";
  my $get_synonym_sth = $xref_source->prepare($get_synonym_sql);
  my $get_pair_sql = "SELECT accession2 FROM pairs where accession1 = ?";
  my $get_pair_sth = $xref_source->prepare($get_pair_sql);
  my ($xref_id, $accession, $version, $label, $description, $info_type, $parsed_seq, $type, $status, $dep_xref_id, $dep_accession, $dep_version, $dep_label, $dep_description, $dep_source_id, $dep_species_id, $linkage_source_id, $synonym, $refseq_pair);

  my @xrefs;
  my $count = 0;
  foreach my $xref_source_id (@source_ids) {
    $get_xref_sth->execute($species_id, $xref_source_id);
    $get_xref_sth->bind_columns(\$xref_id, \$accession, \$version, \$label, \$description, \$info_type);
    while ($get_xref_sth->fetch()) {
      my $xref = {};
      $count++;
  
      $xref->{ACCESSION} = $accession;
      $xref->{LABEL} = $label;
      $xref->{VERSION} = $version;
      $xref->{SPECIES_ID} = $species_id;
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

      # Add pair information if there is some
      $get_pair_sth->execute($accession);
      $get_pair_sth->bind_columns(\$refseq_pair);
      while ($get_pair_sth->fetch) {
	$refseq_pair =~ s/\.[0-9]*//;
        $xref->{PAIR} = $refseq_pair;
      }

      # Look for synonyms
      $get_synonym_sth->execute($xref_id);
      $get_synonym_sth->bind_columns(\$synonym);
      while ($get_synonym_sth->fetch) {
        push (@{$xref->{SYNONYMS} }, $synonym);
      }

      # Add any dependent xrefs
      $get_dependent_sth->execute($xref_id);
      $get_dependent_sth->bind_columns(\$dep_xref_id, \$dep_accession, \$dep_version, \$dep_label, \$dep_description, \$dep_source_id, \$dep_species_id, \$linkage_source_id);
      while ($get_dependent_sth->fetch) {
        if ($dep_species_id != $species_id) { next; }
	if (defined $entrez{$dep_accession}) {
          my %dep;
	  $dep{ACCESSION} = $dep_accession;
	  $dep{LABEL} = $entrez{$dep_accession};
	  $dep{VERSION} = $dep_version;
	  $dep{DESCRIPTION} = $dep_description;
	  $dep{SOURCE_ID} = $entrez_source_id;
	  $dep{LINKAGE_SOURCE_ID} = $linkage_source_id;
	  push @{$xref->{DEPENDENT_XREFS}}, \%dep;

	  my %dep2;
	  $dep2{ACCESSION} = $dep_accession;
	  $dep2{LABEL} = $entrez{$dep_accession};
	  $dep2{VERSION} = $dep_version;
	  $dep2{DESCRIPTION} = $dep_description;
	  $dep2{SOURCE_ID} = $wiki_source_id;
	  $dep2{LINKAGE_SOURCE_ID} = $linkage_source_id;
	  push @{$xref->{DEPENDENT_XREFS}}, \%dep2;

	  # Add dependent xrefs for RefSeq mRNA as well where available
          # only after they are added in priorit 1
          if (defined $refseq_pair) {
            if ($refseq_ids{$refseq_pair}) {
              foreach my $refseq_id (@{ $refseq_ids{$refseq_pair} }) {
                foreach my $entrez_id (@{ $entrez_ids{$dep_accession} }) {
                  $add_dependent_xref_sth->execute($refseq_id, $entrez_id);
                }
                foreach my $wiki_id (@{ $wiki_ids{$dep_accession} }) {
                  $add_dependent_xref_sth->execute($refseq_id, $wiki_id);
                }
              }
            }
          }
        }
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
