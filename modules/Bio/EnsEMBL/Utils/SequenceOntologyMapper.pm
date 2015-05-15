=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
   'protein_coding' 		=> 'SO:0001217', # protein_coding_gene
   'pseudogene' 			=> 'SO:0000336', # pseudogene
   'processed_transcript' 	=> 'SO:0001503', # processed_transcript
   'lincRNA' 				=> 'SO:0001641', # lincRNA_gene
   'polymorphic_pseudogene'=> 'SO:0000336',		 # pseudogene
   'Mt_tRNA' 				=> 'SO:0000088', # mt_gene
   'IG_D_gene' 			=> 'SO:0000510',	 # D_gene
   'snoRNA' 				=> 'SO:0001267', #snoRNA_gene
   'misc_RNA' 				=> 'SO:0000356', #RNA
   'miRNA' 				=> 'SO:0001265', #miRNA_gene
   'rRNA' 					=> 'SO:0001637', #rRNA_gene
   'snRNA'					=> 'SO:0001268', #snRNA_gene
   'snRNA_pseudogene'		=> 'SO:0000336', # pseudogene
   'tRNA_pseudogene'		=> 'SO:0000778', # pseudogenic_tRNA
   'rRNA_pseudogene'		=> 'SO:0000777', # pseudogenic_rRNA
   'TR_J_gene'				=> 'SO:0000470', # J_gene
   'TR_V_gene'				=> 'SO:0000466', # V_gene
   'TR_C_gene'				=> 'SO:0000478', # C_gene
   'ncRNA'					=> 'SO:0001263', # ncRNA_gene
   'tRNA'					=> 'SO:0001272', # tRNA_gene
   'retrotransposed'		=> 'SO:0000569', # retrotransposed
   ## heavily abbreviated
  );

my %transcript_so_mapping = 
  (
   'processed_transcript' 				=> 'SO:0001503', # processed_transcript
   'nonsense_mediated_decay' 			=> 'SO:0001621', # NMD_transcript_variant
   'retained_intron' 					=> 'SO:0000681', # aberrant_processed_transcript
   'transcribed_unprocessed_pseudogene'=> 'SO:0000516', # pseudogenic_transcript
   'processed_pseudogene' 				=> 'SO:0000043', # processed_pseudogene
   'unprocessed_pseudogene' 			=> 'SO:0000336', # pseudogene
   'unitary_pseudogene'				=> 'SO:0000336',        # pseudogene
   'pseudogene' 						=> 'SO:0000336', # pseudogene
   'transcribed_processed_pseudogene'	=> 'SO:0000043',                # processed_pseudogene
   'retrotransposed' 					=> 'SO:0000569', #retrotransposed
   'ncrna_host' 						=> 'SO:0000483', # nc_primary_transcript
   'polymorphic_pseudogene'			=> 'SO:0000336',        # pseudogene
   'lincRNA'							=> 'SO:0001463', # lincRNA
   'ncrna_host'						=> 'SO:0000483', # nc_primary_transcript
   '3prime_overlapping_ncrna'			=> 'SO:0000483',         # nc_primary_transcript
   'TR_V_gene'							=> 'SO:0000466', # V_gene_segment
   'TR_V_pseudogene'					=> 'SO:0000336', # pseudogene
   'TR_J_gene'							=> 'SO:0000470',
   'IG_C_gene'							=> 'SO:0000478',
   'IG_C_pseudogene'					=> 'SO:0000336', # pseudogene
   'TR_C_gene'							=> 'SO:0000478', # C_gene_segment
   'IG_J_pseudogene'					=> 'SO:0000336', # pseudogene
   'miRNA'								=> 'SO:0000276', #miRNA
   'miRNA_pseudogene'					=> 'SO:0000336', # pseudogene
   'disrupted_domain' 					=> 'SO:0000681', # aberrant_processed_transcript
   'rRNA' 								=> 'SO:0000252', #rRNA
   'rRNA_pseudogene'					=> 'SO:0000777', # pseudogenic_rRNA
   'scRNA_pseudogene'					=> 'SO:0000336', # pseudogene
   'snoRNA' 							=> 'SO:0000275', # snoRNA
   'snoRNA_pseudogene'					=> 'SO:0000336', # pseudogene
   'snRNA'								=> 'SO:0000274', # snRNA
   'snRNA_pseudogene'					=> 'SO:0000336',  # pseudogene
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
   'Bio::EnsEMBL::Exon' => 'SO:0000147',       # exon
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
