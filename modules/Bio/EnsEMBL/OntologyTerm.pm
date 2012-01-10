=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::OntologyTerm

=head1 DESCRIPTION

An ontology term object, (most often) created by
Bio::EnsEMBL::DBSQL::GOTermAdaptor and used in querying for
transcripts, genes, and translations using the relevant adaptors and
methods.

=head1 METHODS

=cut

package Bio::EnsEMBL::OntologyTerm;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use base qw( Bio::EnsEMBL::Storable );

=head2 new

  Arg [-ACCESSION]  : String
                      The accession of the ontology term.

  Arg [-ONTOLOGY]   : String
                      The ontology that the term belongs to.

  Arg [-NAMESPACE]  : String
                      The namespace of the ontology term.

  Arg [-NAME]       : String
                      The name of the ontology term.

  Arg [-SUBSETS]     : (optional) Listref of strings
                      The subsets within the ontology to which this
                      term belongs.

  Arg [-DEFINITION] : (optional) String
                      The definition of the ontology term.

  Arg [-SYNONYMS]   : (optional) Listref of strings
                      The synonyms of this term.

  Arg               : Further arguments required for parent class
                      Bio::EnsEMBL::Storable.

  Description   : Creates an ontology term object.

  Example       :

    my $term = Bio::EnsEMBL::OntologyTerm->new(
      '-accession'  => 'GO:0021785',
      '-ontology'   => 'GO',
      '-namespace'  => 'biological_process',
      '-name'       => 'branchiomotor neuron axon guidance',
      '-definition' => 'The process in which a branchiomotor '
        . 'neuron growth cone is directed to a specific target site. '
        . 'Branchiomotor neurons are located in the hindbrain and '
        . 'innervate branchial arch-derived muscles that control jaw '
        . 'movements, facial expression, the larynx, and the pharynx.',
      '-synonyms' => [ 'BMN axon guidance',
                       'branchial motor axon guidance',
                       'special visceral motor neuron axon guidance' ]

        # ... other arguments required by Bio::EnsEMBL::Storable.
    );

  Return type   : Bio::EnsEMBL::OntologyTerm

=cut

sub new {
  my $proto = shift(@_);

  my $this = $proto->SUPER::new(@_);

  my ( $accession, $ontology, $namespace, $name, $definition, $subsets )
    = rearrange( [ 'ACCESSION',  'ONTOLOGY', 'NAMESPACE', 'NAME',
                   'DEFINITION', 'SUBSETS' ],
                 @_ );

  $this->{'accession'}  = $accession;
  $this->{'ontology'}   = $ontology;
  $this->{'namespace'}  = $namespace;
  $this->{'name'}       = $name;
  $this->{'definition'} = $definition;
  $this->{'subsets'}    = [ @{$subsets} ];

  $this->{'child_terms_fetched'}  = 0;
  $this->{'parent_terms_fetched'} = 0;

  return $this;
}

=head2 accession

  Arg           : None

  Description   : Returns the accession for the ontology term in question.

  Example       : my $accession = $term->accession();

  Return type   : String

=cut

sub accession {
  my ($this) = @_;
  return $this->{'accession'};
}

=head2 ontology

  Arg           : None

  Description   : Returns the ontology for the ontology term in question.

  Example       : my $ontology = $term->ontology();

  Return type   : String

=cut

sub ontology {
  my ($this) = @_;
  return $this->{'ontology'};
}

=head2 namespace

  Arg           : None

  Description   : Returns the namespace for the ontology term in question.

  Example       : my $acc = $term->namespace();

  Return type   : String

=cut

sub namespace {
  my ($this) = @_;
  return $this->{'namespace'};
}

=head2 name

  Arg           : None

  Description   : Returns the name for the ontology term in question.

  Example       : my $name = $term->name();

  Return type   : String

=cut

sub name {
  my ($this) = @_;
  return $this->{'name'};
}

=head2 definition

  Arg           : None

  Description   : Returns the definition for the ontology term in question.

  Example       : my $definition = $term->definition();

  Return type   : String

=cut

sub definition {
  my ($this) = @_;
  return $this->{'definition'};
}

=head2 synonyms

  Arg           : None

  Description   : Returns the list of synonyms defined for this term
                  (if any).

  Example       : my @synonyms = @{ $term->synonyms() };

  Return type   : Listref of strings

=cut

sub synonyms {
  my ($this) = @_;

  if ( !exists( $this->{'synonyms'} ) ) {
    $this->{'synonyms'} =
      $this->adaptor()->_fetch_synonyms_by_dbID( $this->dbID() );
  }

  return $this->{'synonyms'};
}

=head2 subsets

  Arg           : None

  Description   : Returns a list of subsets that this term is part
                  of.  The list might be empty.

  Example       : my @subsets = @{ $term->subsets() };

  Return type   : listref of strings

=cut

sub subsets {
  my ($this) = @_;
  return $this->{'subsets'};
}

=head2 children

  Arg           : (optional) List of strings
                  The type of relations to retrieve children for.

  Description   : Returns the children terms of this ontology term.

  Example       : my @children =
                    @{ $term->children( 'is_a', 'part_of' ) };

  Return type   : listref of Bio::EnsEMBL::OntologyTerm

=cut

sub children {
  my ( $this, @relations ) = @_;

  my @terms = @{ $this->adaptor()->fetch_all_by_parent_term($this) };

  if (@relations) {
    @terms = ();
    foreach my $relation (@relations) {
      if ( exists( $this->{'children'}{$relation} ) ) {
        push( @terms, @{ $this->{'children'}{$relation} } );
      }
    }
  }

  return \@terms;
}

=head2 descendants

  Arg           : None

  Description   : Returns the complete set of 'is_a' and 'part_of'
                  descendant terms of this ontology term, down to
                  and including any leaf terms.

  Example       : my @descendants = @{ $term->descendants() };

  Return type   : listref of Bio::EnsEMBL::OntologyTerm

=cut

sub descendants {
  my ($this) = @_;
  return $this->adaptor()->fetch_all_by_ancestor_term($this);
}

=head2 parents

  Arg           : (optional) List of strings
                  The type of relations to retrieve parents for.

  Description   : Returns the parent terms of this ontology term.

  Example       : my @parents =
                    @{ $term->parents( 'is_a', 'part_of' ) };

  Return type   : listref of Bio::EnsEMBL::OntologyTerm

=cut

sub parents {
  my ( $this, @relations ) = @_;

  my @terms = @{ $this->adaptor()->fetch_all_by_child_term($this) };

  if (@relations) {
    @terms = ();
    foreach my $relation (@relations) {
      if ( exists( $this->{'parents'}{$relation} ) ) {
        push( @terms, @{ $this->{'parents'}{$relation} } );
      }
    }
  }

  return \@terms;
}

=head2 ancestors

  Arg           : None

  Description   : Returns the complete set of 'is_a' and 'part_of'
                  ancestor terms of this ontology term, up to and
                  including the root term.

  Example       : my @ancestors = @{ $term->ancestors() };

  Return type   : listref of Bio::EnsEMBL::OntologyTerm

=cut

sub ancestors {
  my ($this) = @_;
  return $this->adaptor()->fetch_all_by_descendant_term($this);
}

1;
