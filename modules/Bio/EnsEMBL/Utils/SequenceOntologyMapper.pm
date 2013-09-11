=pod

=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 NAME

SequenceOntologyMapper - Translates EnsEMBL objects into Sequence Ontology terms

=head1 SYNOPSIS

use Bio::EnsEMBL::Utils::SequenceOntologyMapper

# get an Ensembl feature somehow
...
...

my $ontology_adaptor = $registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );
my $mapper = SequenceOntologyMapper->new($ontology_adaptor);

print $mapper->translate($feature);

=head1 DESCRIPTION

Basic mapper from Ensembl feature objects to Sequence Ontology 
(http://www.sequenceontology.org) terms.

There's no reverse mapping as there doesn't seem to be any utility 
in it at the moment (would require to create empty feature objects).

=cut

package Bio::EnsEMBL::Utils::BiotypeMapper;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Cache;
use Bio::EnsEMBL::Utils::Exception;

=head1 METHODS

=head2 new

    Constructor
    Arg [1]    : OntologyTermAdaptor from the EnsEMBL registry
    Returntype : Bio::EnsEMBL::SequenceOntologyMapper
    Exceptions : If no ontology term adaptor is given

=cut

sub new {
  my $class = shift;
  my $self = 
    { 	
     ontology_adaptor => shift || throw "No ontology term adaptor specified",
     feat_to_acc => 
     {
      'Bio::EnsEMBL::Feature' => 'SO:0000001', # region
      'Bio::EnsEMBL::Gene' => 'SO:0000704', # gene
      'Bio::EnsEMBL::Transcript' => 'SO:0000673', # transcript
      'Bio::EnsEMBL::Exon' => 'SO:0000147', # exon
      'Bio::EnsEMBL::Slice' => 'SO:0000001', # region
      'Bio::EnsEMBL::SimpleFeature' => 'SO:0001411', # biological_region
      'Bio::EnsEMBL::MiscFeature' => 'SO:0001411', # biological_region
      'Bio::EnsEMBL::RepeatFeature' => 'SO:0000657', # repeat region
      'Bio::EnsEMBL::Variation::VariationFeature' => 'SO:0001060', # sequence variant
      'Bio::EnsEMBL::Variation::StructuralVariationFeature' => 'SO:0001537', # structural variant
      'Bio::EnsEMBL::Compara::ConstrainedElement' => 'SO:0001009', # DNA_constraint_sequence ????
      'Bio::EnsEMBL::Funcgen::RegulatoryFeature' => 'SO:0005836', # regulatory_region
     }
    };
 
  $self->{ontology_adaptor}->isa('Bio::EnsEMBL::DBSQL::OntologyTermAdaptor') or 
    throw "Argument is not an OntologyTermAdaptor object";

  tie my %cache, 'Bio::EnsEMBL::Utils::Cache', 100;
  $self->{cache} = \%cache;

  bless $self, $class;
  return $self;
}

=head2 translate

    Arg [0]    : Bio::EnsEMBL::Feature, subclass or related Storable
    Description: Translates a Feature type into an SO term
    Returntype : String; the SO term
    Exceptions : if argument is not an instance of Bio::EnsEMBL::Feature
                 or if cannot translate

=cut

sub translate {
    my $self = shift;
    my $feature = shift;

    my $so_accession;
    my $so_name;
    my $ref = ref($feature);

    $feature->isa('Bio::EnsEMBL::Feature') or 
      throw "Need an instance of Bio::EnsEMBL::Feature";
  
    my $mapping = $self->{feat_to_acc};
    if (exists $mapping->{$ref}) {
      $so_accession = $mapping->{$ref};
    } elsif ($feature->can('SO_term')) {
      $so_accession = $feature->SO_term();
    } else {
      $so_accession = $mapping->{'Bio::EnsEMBL::Feature'};
    }

    if ($so_accession) {
       $so_name = $self->_fetch_SO_name_by_accession($so_accession);
     } else {
       #
       # WARNING
       # 
       # There doesn't seem to be a class_SO_term method in class
       # BaseVariationFeature or its ancestors
       #
       $so_name = $feature->class_SO_term()
	 if $feature->isa('Bio::EnsEMBL::Variation::BaseVariationFeature');
     }
    
    throw sprintf "%s: sequence ontology mapping not found", $ref
      unless $so_name;
    
    return $so_name;
}


=head1 PRIVATE METHODS

=head2 _fetch_SO_name_by_accession

  Arg [0]    : Sequence Ontology accession
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
