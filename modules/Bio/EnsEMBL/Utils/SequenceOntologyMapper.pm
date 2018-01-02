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

=pod


=head1 NAME

SequenceOntologyMapper - Translates EnsEMBL objects into Sequence Ontology terms

=head1 SYNOPSIS

use Bio::EnsEMBL::Utils::SequenceOntologyMapper

# get an Ensembl feature somehow in scalar $feature 
...
...

my $ontology_adaptor = $registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );
my $mapper = SequenceOntologyMapper->new($ontology_adaptor);

print $mapper->to_accession($feature), "\n";
print $mapper->to_name($feature), "\n";

=head1 DESCRIPTION

Basic mapper from Ensembl feature or related objects to Sequence Ontology 
(http://www.sequenceontology.org) terms.

The interface allows to map to SO accessions and names.

=cut

package Bio::EnsEMBL::Utils::SequenceOntologyMapper;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Cache;
use Bio::EnsEMBL::Utils::Exception;

my %gene_so_mapping = 
  (
# Protein coding gene biotype
   'protein_coding' 		=> 'SO:0001217',
   'IG_C_gene' 	  		=> 'SO:0001217',
   'IG_D_gene' 	  		=> 'SO:0001217',
   'IG_gene' 	  		=> 'SO:0001217',
   'IG_J_gene' 	  		=> 'SO:0001217',
   'IG_LV_gene'   		=> 'SO:0001217',
   'IG_M_gene'   		=> 'SO:0001217',
   'IG_V_gene' 	  		=> 'SO:0001217',
   'IG_Z_gene' 	  		=> 'SO:0001217',
   'mRNA' 	  		=> 'SO:0001217',
   'nontranslating_CDS' 	=> 'SO:0001217',
   'polymorphic'		=> 'SO:0001217',
   'polymorphic_pseudogene'	=> 'SO:0001217',
   'TR_C_gene' 	  		=> 'SO:0001217',
   'TR_D_gene' 	  		=> 'SO:0001217',
   'TR_gene' 	  		=> 'SO:0001217',
   'TR_J_gene' 	  		=> 'SO:0001217',
   'TR_V_gene' 	  		=> 'SO:0001217',
   
# Pseudogene biotype
   'IG_C_pseudogene' 			=> 'SO:0000336',
   'IG_D_pseudogene' 			=> 'SO:0000336',
   'IG_J_pseudogene' 			=> 'SO:0000336',
   'IG_pseudogene' 			=> 'SO:0000336',
   'IG_V_pseudogene' 			=> 'SO:0000336',
   'miRNA_pseudogene' 			=> 'SO:0000336',
   'misc_RNA_pseudogene' 		=> 'SO:0000336',
   'Mt_tRNA_pseudogene' 		=> 'SO:0000336',
   'ncbi_pseudogene' 			=> 'SO:0000336',
   'ncRNA_pseudogene' 			=> 'SO:0000336',
   'processed_pseudogene' 		=> 'SO:0000336',
   'pseudogene' 			=> 'SO:0000336',
   'rRNA_pseudogene' 			=> 'SO:0000336',
   'scRNA_pseudogene' 			=> 'SO:0000336',
   'snoRNA_pseudogene' 			=> 'SO:0000336',
   'snRNA_pseudogene' 			=> 'SO:0000336',
   'transcribed_processed_pseudogene' 	=> 'SO:0000336',
   'transcribed_unitary_pseudogene' 	=> 'SO:0000336',
   'transcribed_unprocessed_pseudogene' => 'SO:0000336',
   'translated_processed_pseudogene' 	=> 'SO:0000336',
   'translated_unprocessed_pseudogene' 	=> 'SO:0000336',
   'tRNA_pseudogene' 			=> 'SO:0000336',
   'TR_J_pseudogene' 			=> 'SO:0000336',
   'TR_pseudogene' 			=> 'SO:0000336',
   'TR_V_pseudogene' 			=> 'SO:0000336',
   'unitary_pseudogene' 		=> 'SO:0000336',
   'unprocessed_pseudogene' 		=> 'SO:0000336',

# ncRNA gene biotypes
   '3prime_overlapping_ncrna' 		=> 'SO:0001263',
   'ambiguous_orf' 			=> 'SO:0001263',
   'antisense' 				=> 'SO:0001263',
   'antisense_RNA' 			=> 'SO:0001263',
   'antitoxin' 				=> 'SO:0001263',
   'bidirectional_promoter_lncrna'	=> 'SO:0001263',
   'class_II_RNA' 			=> 'SO:0001263',
   'class_I_RNA' 			=> 'SO:0001263',
   'CRISPR' 				=> 'SO:0001263',
   'guide_RNA' 				=> 'SO:0001263',
   'known_ncrna' 			=> 'SO:0001263',
   'lincRNA' 				=> 'SO:0001263',
   'lncRNA' 				=> 'SO:0001263',
   'macro_lncRNA' 			=> 'SO:0001263',
   'miRNA' 				=> 'SO:0001263',
   'misc_RNA' 				=> 'SO:0001263',
   'Mt_rRNA' 				=> 'SO:0001263',
   'Mt_tRNA' 				=> 'SO:0001263',
   'ncRNA' 				=> 'SO:0001263',
   'ncrna_host' 			=> 'SO:0001263',
   'non_coding' 			=> 'SO:0001263',
   'piRNA' 				=> 'SO:0001263',
   'pre_miRNA' 				=> 'SO:0001263',
   'processed_transcript' 		=> 'SO:0001263',
   'retained_intron' 			=> 'SO:0001263',
   'ribozyme' 				=> 'SO:0001263',
   'RNase_MRP_RNA' 			=> 'SO:0001263',
   'RNase_P_RNA' 			=> 'SO:0001263',
   'rRNA' 				=> 'SO:0001263',
   'scaRNA' 				=> 'SO:0001263',
   'scRNA' 				=> 'SO:0001263',
   'sense_intronic' 			=> 'SO:0001263',
   'sense_overlapping' 			=> 'SO:0001263',
   'snlRNA' 				=> 'SO:0001263',
   'snoRNA' 				=> 'SO:0001263',
   'snRNA' 				=> 'SO:0001263',
   'sRNA' 				=> 'SO:0001263',
   'SRP_RNA' 				=> 'SO:0001263',
   'telomerase_RNA' 			=> 'SO:0001263',
   'tmRNA' 				=> 'SO:0001263',
   'tRNA' 				=> 'SO:0001263',
   'vaultRNA' 				=> 'SO:0001263',
   'Y_RNA' 				=> 'SO:0001263'
  );

my %transcript_so_mapping = 
  (

# mRNA biotypes
   'protein_coding'		=> 'SO:0000234',
   'mRNA'			=> 'SO:0000234',
   'nonsense_mediated_decay'	=> 'SO:0000234',
   'nontranslating_CDS'		=> 'SO:0000234',
   'non_stop_decay'		=> 'SO:0000234',
   'polymorphic_pseudogene'	=> 'SO:0000234',
   
# IG biotypes (SO:3000000 gene_segment)
   'IG_C_gene'			=> 'SO:0000478', # C_gene_segment
   'TR_C_gene'			=> 'SO:0000478', # C_gene_segment
   'IG_D_gene'			=> 'SO:0000458', # D_gene_segment
   'TR_D_gene'			=> 'SO:0000458', # D_gene_segment
   'IG_gene'			=> 'SO:3000000', # gene_segment
   'TR_gene'			=> 'SO:3000000', # gene_segment
   'IG_J_gene'			=> 'SO:0000470', # J_gene_segment
   'TR_J_gene'			=> 'SO:0000470', # J_gene_segment
   'IG_LV_gene'			=> 'SO:3000000', # gene_segment
   'IG_M_gene'			=> 'SO:3000000', # gene_segment
   'IG_V_gene'			=> 'SO:0000466', # V_gene_segment
   'TR_V_gene'			=> 'SO:0000466', # V_gene_segment
   'IG_Z_gene'			=> 'SO:3000000', # gene_segment

# Pseudogenic_transcript biotypes
   'pseudogene' 			=> 'SO:0000516',
   'disrupted_domain'			=> 'SO:0000516',
   'IG_C_pseudogene'			=> 'SO:0000516',
   'IG_D_pseudogene'			=> 'SO:0000516',
   'IG_J_pseudogene'			=> 'SO:0000516',
   'IG_pseudogene'			=> 'SO:0000516',
   'IG_V_pseudogene'			=> 'SO:0000516',
   'miRNA_pseudogene'			=> 'SO:0000516',
   'misc_RNA_pseudogene'		=> 'SO:0000516',
   'Mt_tRNA_pseudogene'			=> 'SO:0000516',
   'ncbi_pseudogene'			=> 'SO:0000516',
   'ncRNA_pseudogene'			=> 'SO:0000516',
   'processed_pseudogene'		=> 'SO:0000516',
   'rRNA_pseudogene'			=> 'SO:0000516',
   'scRNA_pseudogene'			=> 'SO:0000516',
   'snoRNA_pseudogene'			=> 'SO:0000516',
   'snRNA_pseudogene'			=> 'SO:0000516',
   'transcribed_processed_pseudogene'	=> 'SO:0000516',
   'transcribed_unitary_pseudogene'	=> 'SO:0000516',
   'transcribed_unprocessed_pseudogene'	=> 'SO:0000516',
   'translated_processed_pseudogene'	=> 'SO:0000516',
   'translated_unprocessed_pseudogene'	=> 'SO:0000516',
   'tRNA_pseudogene'			=> 'SO:0000516',
   'TR_J_pseudogene'			=> 'SO:0000516',
   'TR_pseudogene'			=> 'SO:0000516',
   'TR_V_pseudogene'			=> 'SO:0000516',
   'unitary_pseudogene'			=> 'SO:0000516',
   'unprocessed_pseudogene'		=> 'SO:0000516',

# ncRNA transcript biotypes
## Long non coding RNAs
   '3prime_overlapping_ncrna'		=> 'SO:0001877',
   'ambiguous_orf'			=> 'SO:0001877',
   'antisense'				=> 'SO:0001877',
   'antisense_RNA'			=> 'SO:0001877',
   'antitoxin'				=> 'SO:0001877',
   'bidirectional_promoter_lncrna'	=> 'SO:0001877',
   'lincRNA'				=> 'SO:0001877',
   'macro_lncRNA'			=> 'SO:0001877',
   'ncrna_host' 			=> 'SO:0001877',
   'non_coding'				=> 'SO:0001877',
   'processed_transcript'		=> 'SO:0001877',
   'retained_intron'			=> 'SO:0001877',
   'ribozyme'				=> 'SO:0001877',
   'sense_intronic'			=> 'SO:0001877',
   'sense_overlapping'			=> 'SO:0001877',
   
## Short non coding RNAs
   'class_II_RNA'			=> 'SO:0000989', # class_II_RNA
   'class_I_RNA'			=> 'SO:0000990', # class_I_RNA
   'guide_RNA'				=> 'SO:0000602', # guide_RNA
   'miRNA'				=> 'SO:0000276', # miRNA
   'known_ncRNA'			=> 'SO:0000655', # ncRNA
   'misc_RNA'				=> 'SO:0000655', # ncRNA
   'ncRNA'				=> 'SO:0000655', # ncRNA
   'piRNA'				=> 'SO:0001035', # piRNA
   'pre_miRNA'				=> 'SO:0001244', # pre_miRNA
   'RNase_MRP_RNA'			=> 'SO:0000385', # RNase_MRP_RNA
   'RNase_P_RNA'			=> 'SO:0000386', # RNase_P_RNA
   'rRNA' 				=> 'SO:0000252', # rRNA
   'Mt_rRNA' 				=> 'SO:0000252', # rRNA
   'scaRNA'				=> 'SO:0000013', # scRNA
   'scRNA'				=> 'SO:0000013', # scRNA
   'snoRNA' 				=> 'SO:0000275', # snoRNA
   'sRNA'				=> 'SO:0000274', # snRNA
   'snlRNA'				=> 'SO:0000274', # snRNA
   'snRNA'				=> 'SO:0000274', # snRNA
   'SRP_RNA'				=> 'SO:0000590', # SRP_RNA
   'telomerase_RNA'			=> 'SO:0000390', # telomerase_RNA
   'tmRNA'				=> 'SO:0000584', # tmRNA
   'tRNA'				=> 'SO:0000253', # tRNA
   'Mt_tRNA'				=> 'SO:0000253', # tRNA
   'vaultRNA'				=> 'SO:0002040', # vaultRNA_primary_transcript
   'vault_RNA'                          => 'SO:0002040', # vaultRNA_primary_transcript
   'Y_RNA'				=> 'SO:0000405', # Y_RNA
   
  );

my %utr_so_mapping =
  (
   'UTR'             => 'SO:0000203', # UTR
   'five_prime_utr'  => 'SO:0000204', # five_prime_UTR
   'three_prime_utr' => 'SO:0000205'  # three_prime_UTR
  );

my %region_so_mapping =
  (
   'chromosome'  => 'SO:0000340', # chromosome
   'supercontig' => 'SO:0000148', # supercontig
   'scaffold'    => 'SO:0000148', # supercontig
   'contig'      => 'SO:0000149'  # contig
  );

my %feature_so_mapping = 
  (
   'Bio::EnsEMBL::Feature' => 'SO:0000001', # region
   'Bio::EnsEMBL::Gene' => 'SO:0000704',    # gene
   'Bio::EnsEMBL::Transcript' => 'SO:0000673', # transcript
   'Bio::EnsEMBL::PredictionTranscript' => 'SO:0000673', # transcript
   'Bio::EnsEMBL::Exon' => 'SO:0000147',       # exon
   'Bio::EnsEMBL::PredictionExon' => 'SO:0000147',       # exon
   'Bio::EnsEMBL::UTR'  => 'SO:0000203',       # UTR
   'Bio::EnsEMBL::ExonTranscript' => 'SO:0000147', # Exon
   'Bio::EnsEMBL::CDS'   => 'SO:0000316',      # CDS
   'Bio::EnsEMBL::Slice' => 'SO:0000001',      # region
   'Bio::EnsEMBL::SimpleFeature' => 'SO:0001411', # biological_region
   'Bio::EnsEMBL::MiscFeature' => 'SO:0001411',	  # biological_region
   'Bio::EnsEMBL::RepeatFeature' => 'SO:0000657', # repeat region
   'Bio::EnsEMBL::Variation::VariationFeature' => 'SO:0001060', # sequence variant
   'Bio::EnsEMBL::Variation::StructuralVariationFeature' => 'SO:0001537', # structural variant
   'Bio::EnsEMBL::Compara::ConstrainedElement' => 'SO:0001009', #DNA_constraint_sequence ????
   'Bio::EnsEMBL::Funcgen::RegulatoryFeature' => 'SO:0005836', # regulatory_region
   'Bio::EnsEMBL::DnaDnaAlignFeature' => 'SO:0000347', # nucleotide_match
   'Bio::EnsEMBL::DnaPepAlignFeature' => 'SO:0000349', # protein_match
   'Bio::EnsEMBL::KaryotypeBand' => 'SO:0000341', # chromosome_band
  );


=head1 METHODS

=head2 new

    Constructor
    Arg [1]    : OntologyTermAdaptor from the EnsEMBL registry
    Returntype : Bio::EnsEMBL::SequenceOntologyMapper
    Exceptions : If no ontology term adaptor is passed as argument

=cut

sub new {
  my ($class, $oa) = @_;
  defined $oa or throw "No ontology term adaptor specified";

  my $self = 
    { 	
     ontology_adaptor => $oa,
     feat_to_acc => \%feature_so_mapping,
     gene_to_acc => \%gene_so_mapping,
     utr_to_acc => \%utr_so_mapping,
     region_to_acc => \%region_so_mapping,
     tran_to_acc => \%transcript_so_mapping
    };
 
  $self->{ontology_adaptor}->isa('Bio::EnsEMBL::DBSQL::OntologyTermAdaptor') or 
    throw "Argument is not an OntologyTermAdaptor object";

  tie my %cache, 'Bio::EnsEMBL::Utils::Cache', 100;
  $self->{cache} = \%cache;

  bless $self, $class;
  return $self;
}

=head2 to_accession

    Arg [0]    : Instance of Bio::EnsEMBL::Feature, subclass or
                 related Storable
    Description: translates a Feature type into an SO term accession
    Returntype : String; the SO accession
    Exceptions : if cannot map to SO term

=cut

sub to_accession {
  my $self = shift;
  my $feature = shift;

  my $so_accession;
  my $ref = ref($feature);
  
  my ($gene_to_acc, $tran_to_acc, $feat_to_acc, $utr_to_acc, $region_to_acc) = 
    ($self->{gene_to_acc}, $self->{tran_to_acc}, $self->{feat_to_acc}, $self->{utr_to_acc}, $self->{region_to_acc});
  
  if ($feature->isa('Bio::EnsEMBL::Gene') and 
      exists $gene_to_acc->{$feature->biotype}) {
    $so_accession = $gene_to_acc->{$feature->biotype};
  } elsif ($feature->isa('Bio::EnsEMBL::Transcript') and 
	   exists $tran_to_acc->{$feature->biotype}) {
    $so_accession = $tran_to_acc->{$feature->biotype};
  } elsif ($feature->isa('Bio::EnsEMBL::UTR') and
           exists $utr_to_acc->{$feature->type}) {
    $so_accession = $utr_to_acc->{$feature->type};
  } elsif ($feature->isa('Bio::EnsEMBL::Slice') and
           exists $region_to_acc->{$feature->coord_system_name}) {
    $so_accession = $region_to_acc->{$feature->coord_system_name};
  }

  if (not $so_accession and exists $feat_to_acc->{$ref}) {
    $so_accession = $feat_to_acc->{$ref};
  } else {
    $so_accession = $feature->SO_term()
      if $feature->can('SO_term');
  }

  throw sprintf "%s: mapping to sequence ontology accession not found", $ref
    unless $so_accession;
    
  return $so_accession;
}

=head2 to_name

    Arg [0]    : Instance of Bio::EnsEMBL::Feature, subclass or
                 related Storable
    Description: translates a Feature type into an SO term name
    Returntype : String; the SO term name
    Exceptions : if cannot map to an SO term

=cut

sub to_name {
  my $self = shift;
  my $feature = shift;

  my $so_name;
  my $so_accession = eval {
    $self->to_accession($feature);
  };
  if ($@) {
    $so_name = $feature->class_SO_term()
      if $feature->isa('Bio::EnsEMBL::Variation::BaseVariationFeature');
  } else {
    $so_name = $self->_fetch_SO_name_by_accession($so_accession);
  }

  throw sprintf "%s: mapping to sequence ontology name not found", ref($feature)
    unless $so_name;

  return $so_name;
}

=head2 gene_biotype_to_name

    Arg [0]    : Biotype string
    Description: translates a biotype into an SO term name
    Returntype : String; the SO term name
    Exceptions : if cannot map to an SO term

=cut

sub gene_biotype_to_name {
  my $self = shift;
  my $biotype = shift;

  if (exists $gene_so_mapping{$biotype}) {
    return $gene_so_mapping{$biotype};
  } else {
    throw "Biotype not found in gene SO mapping $biotype";
  }
}

=head2 transcript_biotype_to_name

    Arg [0]    : Biotype string
    Description: translates a biotype into an SO term name
    Returntype : String; the SO term name
    Exceptions : if cannot map to an SO term

=cut

sub transcript_biotype_to_name {
  my $self = shift;
  my $biotype = shift;

  if (exists $transcript_so_mapping{$biotype}) {
    return $transcript_so_mapping{$biotype};
  } else {
    throw "Biotype not found in transcript SO mapping $biotype";
  }
}


=head1 PRIVATE METHODS

=head2 _fetch_SO_name_by_accession

  Arg [0]    : String; Sequence Ontology accession
  Description: Returns the name linked to the given accession. These are
               internally cached for speed.
  Returntype : String; the name of the given accession
  Exceptions : None

=cut

sub _fetch_SO_name_by_accession {
  my ($self, $so_accession) = @_;
  my $so_name = $self->{cache}->{$so_accession};

  if(!$so_name) {
    my $so_term = $self->{'ontology_adaptor'}->fetch_by_accession($so_accession);
    $so_name = $so_term->name();
    $self->{cache}->{$so_accession} = $so_name;
  }

  return $so_name;
}

1;
