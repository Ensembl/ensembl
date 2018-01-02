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

package XrefMapper::drosophila;
use strict;

use  XrefMapper::BasicMapper;
use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

use XrefMapper::BasicMapper qw(%stable_id_to_internal_id %object_xref_mappings %xref_to_source %xref_accessions %source_to_external_db);

my %genes_to_transcripts;
my %transcript_to_translation;
my %translation_to_transcript;
my %transcript_length;

sub gene_description_filter_regexps {

		return ();

}

sub set_methods{

  my $default_method = 'ExonerateGappedBest1';
  my %override_method_for_source = (
	   ExonerateGappedBest5 => ['RefSeq_mRNA','RefSeq_mRNA_predicted', 'RefSeq_ncRNA', 'RefSeq_ncRNA_predicted' ],
         );

  return $default_method, \%override_method_for_source;
}



# Special logic for drosophila display_xrefs:
#
# gene: flybase_name if present, else gadfly_gene_cgid
#
# transcript: flybase_name if present, else gadfly_transcript_cgid
sub xref_offset{
		my ($self, $val) = @_;

		if(defined($val)){
				$self->{'_xref_offset'} = $val;
		}
		return $self->{'_xref_offset'};
}

sub gene_description_sources {
  return (
          "FlyBaseName_gene",
          "FlyBaseCGID_gene",
         );
}

sub transcript_display_xref_sources {
  my $self     = shift;

  my @list = qw(FlyBaseName_transcript
                        FlyBaseCGID_transcript);

  my %ignore;

  return [\@list,\%ignore];

}

sub gene_display_xref_sources {
  my $self     = shift;

  my @list = qw(FlyBaseName_gene
                FlyBaseCGID_gene
                flybase_gene_id);

  my %ignore;


  return [\@list,\%ignore];

}


sub build_genes_to_transcripts {

  my ($self) = @_;

  my $sql = "SELECT gene_id, transcript_id, seq_region_start, seq_region_end FROM transcript";
  my $sth = $self->core->dbc->prepare($sql);
  $sth->execute();

  my ($gene_id, $transcript_id, $start, $end);
  $sth->bind_columns(\$gene_id, \$transcript_id, \$start, \$end);

  # Note %genes_to_transcripts is global
  while ($sth->fetch()) {
    push @{$genes_to_transcripts{$gene_id}}, $transcript_id;
    $transcript_length{$transcript_id} = $end- $start;
  }

}

sub load_translation_to_transcript{
  my ($self) = @_;

  my $sth = $self->core->dbc->prepare("SELECT translation_id, transcript_id FROM translation");
  $sth->execute();

  my ($translation_id, $transcript_id);
  $sth->bind_columns(\$translation_id, \$transcript_id);

  while ($sth->fetch()) {
    $translation_to_transcript{$translation_id} = $transcript_id;
    $transcript_to_translation{$transcript_id} = $translation_id if ($translation_id);
  }
}

# Want to use FlyBase transcript names, rather than deriving them from genes.
sub transcript_names_from_gene {
  return;
}

1;
