=pod

=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
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

my $list_of_biotypes = $biotype_mapper->biotypes_belonging_to_group('protein-coding');

=head1 DESCRIPTION

BiotypeMapper provides a series of nearest matches between EnsEMBL biotypes and
the Sequence Ontology (http://www.sequenceontology.org). In addition, biotypes
are members of groupings, such as "short non-coding". This allows one to
conveniently select all the biotypes of a certain kind. 

SO Mappings are imperfect due to the inexact correspondance of biotypes to 
several SO terms. The a best guess has been chosen in each case.

Reverse mappings from SO to biotype are vague, due to many-to-one relationships.
In this case a list of possible terms is given. 

=cut

package Bio::EnsEMBL::Utils::BiotypeMapper;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Utils::Cache;

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
	'Bio::EnsEMBL::Exon' => 'SO:0000147',
	'Bio::EnsEMBL::Slice' => 'SO:0000001', # region
	'Bio::EnsEMBL::SimpleFeature' => 'SO:0001411', # biological_region
	'Bio::EnsEMBL::MiscFeature' => 'SO:0001411', # biological_region
	'Bio::EnsEMBL::RepeatFeature' => 'SO:0000657', # repeat region
	'Bio::EnsEMBL::Variation::VariationFeature' => 'SO:0001060', # sequence variant
	'Bio::EnsEMBL::Variation::StructuralVariationFeature' => 'SO:0001537', # structural variant
  'Bio::EnsEMBL::Compara::ConstrainedElement' => 'SO:0001009', #DNA_constraint_sequence ????
	'Bio::EnsEMBL::Funcgen::RegulatoryFeature' => 'SO:0005836', # regulatory_region
);

my %grouping_of_biotypes = (
    # Genebuilder/Havana categorisation
    'protein_coding' => [qw( protein_coding polymorphic_pseudogene   )],
    'pseudogene'     => [qw( pseudogene retrotransposed snRNA_pseudogene tRNA_pseudogene
                            miRNA_pseudogene Mt_tRNA_pseudogene rRNA_pseudogene
                            scRNA_pseudogene misc_RNA_pseudogene snoRNA_pseudogene
                            processed_pseudogene
                        )],
    'long_noncoding' => [qw( 3prime_overlapping_ncrna antisense lincRNA ncrna_host non_coding 
                            processed_transcript sense_intronic sense_overlapping
                        )],
    'short_noncoding'=> [qw( miRNA misc_RNA  Mt_tRNA 
                             rRNA snoRNA  snRNA 
                        )],
    # practical Ensembl core categories for fasta dumping
    'cdna'              => [qw( protein_coding polymorphic_pseudogene IG_V_gene TR_V_gene 
                                IG_J_gene TR_J_gene IG_D_gene IG_C_gene TR_C_gene pseudogene
                                retrotransposed IG_V_pseudogene TR_V_pseudogene 
                                IG_J_pseudogene IG_C_pseudogene processed_transcript
                                antisense ambiguous_orf transcribed_processed_pseudogene
                                disrupted_domain processed_pseudogene
                           )],
    'peptide_producing' => [qw( protein_coding polymorphic_pseudogene IG_V_gene TR_V_gene 
                                IG_J_gene TR_J_gene IG_D_gene IG_C_gene TR_C_gene IG_LV_gene
                                nonsense_mediated_decay
                           )],
    'ncrna'             => [qw( ncRNA miRNA miRNA_pseudogene misc_RNA misc_RNA_pseudogene Mt_tRNA 
                            Mt_tRNA_pseudogene Mt_rRNA rRNA rRNA_pseudogene scRNA_pseudogene 
                            snoRNA snoRNA_pseudogene snRNA snRNA_pseudogene tRNA_pseudogene
                            3prime_overlapping_ncrna lincRNA ncrna_host non_coding 
                             sense_intronic sense_overlapping tRNA
                            )],
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
	
	tie my %cache, 'Bio::EnsEMBL::Utils::Cache', 100;
  $self->{cache} = \%cache;

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
	my $so_name;
	my $ref = ref($feature);
	if ($feature->isa('Bio::EnsEMBL::Gene') && exists $gene_so_mapping{$feature->biotype}) {
		$so_accession = $gene_so_mapping{$feature->biotype};
	}
	elsif ($feature->isa('Bio::EnsEMBL::Transcript') && exists $transcript_so_mapping{$feature->biotype}) {
		$so_accession = $transcript_so_mapping{$feature->biotype};
	}
	elsif ($feature->isa('Bio::EnsEMBL::Variation::BaseVariationFeature')) {
	  $so_name = $feature->class_SO_term();
	}
	
	if (! $so_accession && ! $so_name && exists $feature_so_mapping{$ref}) {
	  $so_accession = $feature_so_mapping{$ref};
	}
	else {
		if($feature->can('SO_term')) {
		  $so_accession = $feature->SO_term();
		}
	}
	
	if ($so_accession) {
		$so_name = $self->fetch_SO_name_by_accession($so_accession);
	}
	
	throw ("Ontology mapping not found for ".ref($feature)) unless $so_name;
	
	return $so_name;
}

=head2 fetch_SO_name_by_accession

  Arg [0]    : Sequence Ontology accession
  Description: Returns the name linked to the given accession. These are
               internally cached for speed.
  Returntype : The name of the given accession.

=cut

sub fetch_SO_name_by_accession {
  my ($self, $so_accession) = @_;
  my $so_name = $self->{cache}->{$so_accession};
  if(!$so_name) {
    my $so_term = $self->{'ontology_adaptor'}->fetch_by_accession($so_accession);
    $so_name = $so_term->name();
    $self->{cache}->{$so_accession} = $so_name;
  }
  return $so_name;
}


=head2 translate_SO_to_biotype

	Arg [0]    : Sequence Ontology term, either in name or URI format
	Description: Returns the closest corresponding Ensembl biotypes to a given SO term
	Returntype : Listref Array of Strings containing possible biotypes
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

	return \@biotypes;
}

=head2 belongs_to_groups

    Arg [0]    : Biotype (string)
    Description: Returns the group names that include the given biotype
    Returntype : Listref of strings
=cut

sub belongs_to_groups {
    my $self = shift;
    my $member = shift;
    my @belongs_to;
    foreach my $group (keys %grouping_of_biotypes) {
        $group = lc($group);
        foreach my $biotype ( @{ $grouping_of_biotypes{$group} }) {
            if ($biotype eq $member) {push @belongs_to,$group;}
        }
    }
    return \@belongs_to;
}

=head2 group_members

    Arg [0]    : Biotype group name (string)
    Description: Returns a list of biotypes that belong in the group.
    Returntype : Listref of strings
=cut

sub group_members {
    my $self = shift;
    my $group = lc(shift);
    if (exists($grouping_of_biotypes{$group})) {
        my @biotypes = @{ $grouping_of_biotypes{$group} };
        return \@biotypes;
    }
    else {
        throw ("$group is not a valid group name for biotypes");
    }
}

=head2 member_of_group

    Arg [0]    : Biotype (string)
    Arg [1]    : Group to check (string) 
    Description: Returns true if a biotype is present in a group
    Returntype : Boolean
=cut

sub member_of_group {
    my $self = shift;
    my $biotype = shift;
    my $query_group = lc(shift);
    my @groups = @{ $self->belongs_to_groups($biotype) };
    while (my $group = lc(shift @groups)) {
        if ($group eq $query_group) {
            return 1;   
        }
    }
    return 0;
}

1;
