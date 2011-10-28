=pod

=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 NAME

BiotypeMapper - Translates EnsEMBL biotypes into Sequence Ontology terms and back

=head1 AUTHOR

Kieron Taylor, 2011 - ktaylor@ebi.ac.uk

=head1 SYNOPSIS

use Bio::EnsEMBL::Utils::BiotypeMapper

my $ontology_adaptor = $registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );
my $biotype_mapper = new BiotypeMapper($ontology_adaptor);

print $biotype_mapper->translate_feature_to_SO_term($feature);

=head1 DESCRIPTION

BiotypeMapper provides a series of nearest matches between EnsEMBL biotypes and
the Sequence Ontology (http://www.sequenceontology.org)

Mappings are imperfect due to the inexact correspondance of biotypes to 
several SO terms. The a best guess has been chosen in each case.

Reverse mappings from SO to biotype are vague, due to many-to-one relationships.
In this case a list of possible terms is given. 

=cut

package Bio::EnsEMBL::Utils::BiotypeMapper;

use strict;
use warnings;
use Carp;

my %gene_so_mapping = (
	'protein_coding' 		=> 'SO:0001217', # protein_coding_gene
	'pseudogene' 			=> 'SO:0000336', # pseudogene
	'processed_transcript' 	=> 'SO:0001503', # processed_transcript
	'lincRNA' 				=> 'SO:0001641', # lincRNA_gene
	'polymorphic_pseudogene'=> 'SO:0000336', # pseudogene
	'Mt_tRNA' 				=> 'SO:0000088', # mt_gene
	'IG_D_gene' 			=> 'SO:0000510', # D_gene
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

my %transcript_so_mapping = (
	'processed_transcript' 				=> 'SO:0001503', # processed_transcript
	'nonsense_mediated_decay' 			=> 'SO:0001621', # NMD_transcript_variant
	'retained_intron' 					=> 'SO:0000681', # aberrant_processed_transcript
	'transcribed_unprocessed_pseudogene'=> 'SO:0000516', # pseudogenic_transcript
	'processed_pseudogene' 				=> 'SO:0000043', # processed_pseudogene
	'unprocessed_pseudogene' 			=> 'SO:0000336', # pseudogene
	'unitary_pseudogene'				=> 'SO:0000336',
	'pseudogene' 						=> 'SO:0000336', # pseudogene
	'transcribed_processed_pseudogene'	=> 'SO:0000043', 
	'retrotransposed' 					=> 'SO:0000569', #retrotransposed
	'ncrna_host' 						=> 'SO:0000483',
	'polymorphic_pseudogene'			=> 'SO:0000336',
	'lincRNA'							=> 'SO:0001463',
	'ncrna_host'						=> 'SO:0000483',
	'3prime_overlapping_ncrna'			=> 'SO:0000483',
	'TR_V_gene'							=> 'SO:0000466',
	'TR_V_pseudogene'					=> 'SO:0000336',

	'TR_J_gene'							=> 'SO:0000470',
	'IG_C_gene'							=> 'SO:0000478',
	'IG_C_pseudogene'					=> 'SO:0000336',
	'TR_C_gene'							=> 'SO:0000478',
	'IG_J_pseudogene'					=> 'SO:0000336',
	'miRNA'								=> 'SO:0000276', #miRNA
	'miRNA_pseudogene'					=> 'SO:0000336',
	'disrupted_domain' 					=> 'SO:0000681', # aberrant_processed_transcript
	'rRNA' 								=> 'SO:0000252', #rRNA
	'rRNA_pseudogene'					=> 'SO:0000777', 
	'scRNA_pseudogene'					=> 'SO:0000336',
	'snoRNA' 							=> 'SO:0000275', # snoRNA
	'snoRNA_pseudogene'					=> 'SO:0000336',
	'snRNA'								=> 'SO:0000274', # snRNA
	'snRNA_pseudogene'					=> 'SO:0000336',

	);

my %feature_so_mapping = (
	'Bio::EnsEMBL::Gene' => 'SO:0000704', # gene
	'Bio::EnsEMBL::Transcript' => 'SO:0000673', # transcript
	'Bio::EnsEMBL::Slice' => 'SO:0000001', # region
	'Bio::EnsEMBL::Variation::VariationFeature' => 'SO:0001060', # sequence variant
	'Bio::EnsEMBL::Variation::StructuralVariationFeature' => 'SO:0001537', # structural variant
    'Bio::EnsEMBL::Compara::ConstrainedElement' => 'SO:0001009', #DNA_constraint_sequence ????
	'Bio::EnsEMBL::Funcgen::RegulatoryFeature' => 'SO:0001679', # transcription_regulatory_region
);

=head2 new

    Constructor
    Arg [1]    : OntologyAdaptor from the EnsEMBL registry
	Returntype : Bio::EnsEMBL::BiotypeMapper

=cut

sub new {
	my $class = shift;
	my $self = { 	
			ontology_adaptor => shift,
	};

	bless $self, $class;
    return $self;
}

=head2 translate_feature_to_SO_term

	Arg [0]    : Bio::EnsEMBL::Feature, subclass or related Storable
	Description: Translates a Feature type into an SO term. If the Feature is a
	             Gene or Transcript, then a further refinement of the type is made
				 via Biotype
	Returntype : String

=cut

sub translate_feature_to_SO_term {
	my $self = shift;
	my $feature = shift;
	my $so_accession;
	my $so_term;
	if (ref($feature) eq "Bio::EnsEMBL::Gene" and exists $gene_so_mapping{$feature->biotype}) {
		$so_accession = $gene_so_mapping{$feature->biotype};
	}
	elsif (ref($feature) eq "Bio::EnsEMBL::Transcription" and exists $transcript_so_mapping{$feature->biotype}) {
		$so_accession = $transcript_so_mapping{$feature->biotype};
	}
	else {
		$so_accession = $feature_so_mapping{ref($feature)};
	}
	if (defined($so_accession)) {
		$so_term = $self->{'ontology_adaptor'}->fetch_by_accession($so_accession);
	}
	else {
		carp "Ontology mapping not found for ".ref($feature)."\n";
		return "????????";
	}

	return $so_term->name;
}


=head2 translate_SO_to_biotype

	Arg [0]    : Sequence Ontology term, either in name or URI format
	Description: Returns the closest corresponding Ensembl biotypes to a given SO term
	Returntype : String containing a comma-separated list of biotypes
=cut

sub translate_SO_to_biotype {
	my $self = shift;
	my $translate_me = shift;

	my @so_names;
# look up text in ontology database
	if ($translate_me !~ /^SO:/) {
		my $so_terms = $self->{'ontology_adaptor'}->fetch_all_by_name($translate_me);
		@so_names = [];
		foreach my $term (@{$so_terms}) {
			push @so_names,$term->accession();
		}
	}
	else {
		push @so_names,$translate_me;
	}
# convert list of accessions into biotypes
	my @biotypes;
	foreach my $accession (@so_names) {
		foreach	my $key (keys %gene_so_mapping) {
			if ($gene_so_mapping{$key} eq $accession) {
				push @biotypes,$key;
			}
		}
		foreach my $key (keys %transcript_so_mapping) {
			if ($transcript_so_mapping{$key} eq $accession) {
				push @biotypes,$key;
			}
		}
		foreach my $key (keys %feature_so_mapping) {
			if ($feature_so_mapping{$key} eq $accession) {
				push @biotypes,$key;
			}
		}
	}

	return join (',',@biotypes);
}

1;
