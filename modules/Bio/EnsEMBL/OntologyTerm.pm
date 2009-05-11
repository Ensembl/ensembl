=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

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

  Arg [-NAMESPACE]  : String
                      The namespace of the ontology term.

  Arg [-NAME]       : String
                      The name of the ontology term.

  Arg [-DEFINITION] : (optional) String
                      The definition of the ontology term.

  Arg               : Further arguments required for parent class
                      Bio::EnsEMBL::Storable.

  Description   : Creates an ontology term object.

  Example       :

    my $term = Bio::EnsEMBL::OntologyTerm->new(
      '-accession'  => 'GO:0030326',
      '-namespace'  => 'biological_process',
      '-name'       => 'embryonic limb morphogenesis',
      '-definition' => 'The process, occurring in the embryo, '
        . 'by which the anatomical structures of the limb '
        . 'are generated and organized. '
        . 'Morphogenesis pertains to the creation of form. '
        . 'A limb is an appendage of an animal '
        . 'used for locomotion or grasping.'
        # ... other arguments required by Bio::EnsEMBL::Storable.
    );

  Return type   : Bio::EnsEMBL::OntologyTerm

=cut

sub new {
  my $proto = shift(@_);

  my $this = $proto->SUPER::new(@_);

  my ( $accession, $namespace, $name, $definition ) =
    rearrange( [ 'ACCESSION', 'NAMESPACE', 'NAME', 'DEFINITION' ], @_ );

  $this->{'accession'}  = $accession;
  $this->{'namespace'}  = $namespace;
  $this->{'name'}       = $name;
  $this->{'definition'} = $definition;

  $this->{'child_terms_fetched'}  = 0;
  $this->{'parent_terms_fetched'} = 0;

  return $this;
}

=head2 accession

  Arg           : (optional) String
                  The new accession of this ontology term.

  Description   : Accesses and/or assigns the accession for the
                  ontology term in question.

  Example       : my $accession = $term->accession();

  Return type   : String

=cut

sub accession {
  my ( $this, $value ) = @_;
  if ( defined($value) ) { $this->{'accession'} = $value }
  return $this->{'accession'};
}

=head2 namespace

  Arg           : (optional) String
                  The new namespace of this ontology term.

  Description   : Accesses and/or assigns the namespace for the
                  ontology term in question.

  Example       : my $acc = $term->namespace();

  Return type   : String

=cut

sub namespace {
  my ( $this, $value ) = @_;
  if ( defined($value) ) { $this->{'namespace'} = $value }
  return $this->{'namespace'};
}

=head2 name

  Arg           : (optional) String
                  The new name of this ontology term.

  Description   : Accesses and/or assigns the name for the ontology
                  term in question.

  Example       : my $name = $term->name();

  Return type   : String

=cut

sub name {
  my ( $this, $value ) = @_;
  if ( defined($value) ) { $this->{'name'} = $value }
  return $this->{'name'};
}

=head2 definition

  Arg           : (optional) String
                  The new definition of this ontology term.

  Description   : Accesses and/or assigns the definition for the
                  ontology term in question.

  Example       : my $definition = $term->definition();

  Return type   : String

=cut

sub defintion {
  my ( $this, $value ) = @_;
  if ( defined($value) ) { $this->{'defintion'} = $value }
  return $this->{'definition'};
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

  my @terms = @{ $this->adaptor()->fetch_by_parent_term($this) };

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

  Example       : my @desc = @{ $term->descendants() };

  Return type   : listref of Bio::EnsEMBL::OntologyTerm

=cut

sub descendants {
  my ($this) = @_;
  return $this->adaptor()->fetch_all_by_parent_term($this);
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

  my @terms = @{ $this->adaptor()->fetch_by_child_term($this) };

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

  Example       : my @desc = @{ $term->descendants() };

  Return type   : listref of Bio::EnsEMBL::OntologyTerm

=cut

sub ancestors {
  my ($this) = @_;
  return $this->adaptor()->fetch_all_by_child_term($this);
}

1;
