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

Bio::EnsEMBL::DBSQL::OntologyTermAdaptor

=head1 SYNOPSIS

  my $goa = $registry->get_adaptor( 'Multi', 'Ontology', 'GOTerm' );

  my $term             = $goa->fetch_by_accession('GO:0010885');

  my @child_terms      = @{ $goa->fetch_by_parent_term($term) };
  my @descendant_terms = @{ $goa->fetch_all_by_parent_term($term) };

  my @parent_terms   = @{ $goa->fetch_by_child_term($term) };
  my @ancestor_terms = @{ $goa->fetch_all_by_child_term($term) };

  my %ancestor_chart = %{ $goa->fetch_ancestor_chart($term) };

=head1 DESCRIPTION

An adaptor for fetching ontology terms, creates
Bio::EnsEMBL::OntologyTerm objects.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::OntologyTermAdaptor;

use strict;
use warnings;

use DBI qw( :sql_types );

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );

use Bio::EnsEMBL::OntologyTerm;

use base qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );

=head2 new

  Arg [1]       : Bio::EnsEMBL::DBSQL::DBAdaptor
                  Argument required for parent class
                  Bio::EnsEMBL::DBSQL::BaseAdaptor.

  Arg [2]       : String
                  The particular ontology that this ontology adaptor
                  deals with.

  Description   : Creates an ontology term adaptor.

  Example       :

    my $ot_adaptor =
      Bio::EnsEMBL::DBSQL::OntologyTermAdaptor->new( $dba, 'GO' );

  Return type   : Bio::EnsEMBL::DBSQL::OntologyTermAdaptor

=cut

sub new {
  my ( $proto, $dba, $ontology ) = @_;

  if ( !ref($dba) || !$dba->isa('Bio::EnsEMBL::DBSQL::DBAdaptor') ) {
    throw('First argument needs to be a '
        . 'Bio::EnsEMBL::DBSQL::DBAdaptor object' );
  }

  my $this = $proto->SUPER::new($dba);

  $this->{'ontology'} = $ontology;

  return $this;
}

=head2 ontology

  Arg           : None

  Description   : Returns the name of the ontology which this
                  adaptor is used to retrieve terms for.

  Example       :

    my $ontology = $ot_adaptor->ontology();

  Return type   : String

=cut

sub ontology {
  my ($this) = @_;
  return $this->{'ontology'};
}

=head2 fetch_by_accession

  Arg [1]       : String

  Description   : Fetches an ontology term given an accession.

  Example       :

    my $term = $ot_adaptor->fetch_by_accession('GO:0030326');

  Return type   : Bio::EnsEMBL::OntologyTerm

=cut

sub fetch_by_accession {
  my ( $this, $accession ) = @_;

  my $statement = q(
SELECT  term.term_id,
        term.name,
        term.definition,
        ontology.namespace
FROM    ontology,
        term
WHERE   ontology.name = ?
  AND   ontology.ontology_id = term.ontology_id
  AND   term.accession = ?);

  my $sth = $this->prepare($statement);
  $sth->bind_param( 1, $this->{'ontology'}, SQL_VARCHAR );
  $sth->bind_param( 2, $accession,          SQL_VARCHAR );

  $sth->execute();

  my ( $dbid, $name, $definition, $namespace );
  $sth->bind_columns( \( $dbid, $name, $definition, $namespace ) );

  $sth->fetch();
  my $term = Bio::EnsEMBL::OntologyTerm->new(
    '-dbid'       => $dbid,
    '-adaptor'    => $this,
    '-accession'  => $accession,
    '-namespace'  => $namespace,
    '-name'       => $name,
    '-definition' => $definition
  );
  $sth->finish();

  return $term;
} ## end sub fetch_by_accession

=head2 fetch_by_parent_term

  Arg [1]       : Bio::EnsEMBL::OntologyTerm
                  The term whose children terms should be fetched.

  Description   : Given a parent ontology term, returns a list of
                  its immediate children terms.

  Example       : my @children =
                    @{ $ot_adaptor->fetch_by_parent_term($term) };

  Return type   : listref of Bio::EnsEMBL::OntologyTerm

=cut

sub fetch_by_parent_term {
  my ( $this, $term ) = @_;

  if ( !ref($term) || !$term->isa('Bio::EnsEMBL::OntologyTerm') ) {
    throw('Argument needs to be a Bio::EnsEMBL::OntologyTerm object');
  }

  my @terms;

  if ( !$term->{'child_terms_fetched'} ) {
    my $statement = q(
SELECT  child_term.term_id,
        child_term.accession,
        child_term.name,
        child_term.definition,
        rt.name
FROM    term child_term,
        relation,
        relation_type rt
WHERE   relation.child_term_id = child_term.term_id
  AND   relation.parent_term_id = ?
  AND   relation.relation_type_id = rt.relation_type_id);

    my $sth = $this->prepare($statement);
    $sth->bind_param( 1, $term->dbID(), SQL_INTEGER );

    $sth->execute();

    my ( $dbid, $accession, $name, $definition, $relation );
    $sth->bind_columns(
      \( $dbid, $accession, $name, $definition, $relation ) );

    while ( $sth->fetch() ) {
      my $child_term = Bio::EnsEMBL::OntologyTerm->new(
        '-dbid'       => $dbid,
        '-adaptor'    => $this,
        '-accession'  => $accession,
        '-namespace'  => $term->{'namespace'},
        '-name'       => $name,
        '-definition' => $definition,
      );

      push( @terms,                              $child_term );
      push( @{ $term->{'children'}{$relation} }, $child_term );
    }

    $term->{'child_terms_fetched'} = 1;
  } else {
    foreach my $relation ( values( %{ $term->{'children'} } ) ) {
      push( @terms, @{$relation} );
    }
  }

  return \@terms;
} ## end sub fetch_by_parent_term

=head2 fetch_all_by_parent_term

  Arg [1]       : Bio::EnsEMBL::OntologyTerm
                  The term whose descendant terms should be fetched.

  Description   : Given a parent ontology term, returns a list of
                  all its descendant terms, down to and including
                  any leaf terms.  In GO, relations of the type
                  'is_a' and 'part_of' are followed, but not
                  'regulates' etc.

  Example       : my @descendants =
                    @{ $ot_adaptor->fetch_all_by_parent_term($term) };

  Return type   : listref of Bio::EnsEMBL::OntologyTerm

=cut

sub fetch_all_by_parent_term {
  my ( $this, $term ) = @_;

  if ( !ref($term) || !$term->isa('Bio::EnsEMBL::OntologyTerm') ) {
    throw('Argument needs to be a Bio::EnsEMBL::OntologyTerm object');
  }

  my $statement = q(
SELECT DISTINCT
        child_term.term_id,
        child_term.accession,
        child_term.name,
        child_term.definition
FROM    term child_term,
        closure
WHERE   closure.child_term_id = child_term.term_id
  AND   closure.parent_term_id = ?
  AND   closure.distance > 0
ORDER BY closure.distance, child_term.accession);

  my $sth = $this->prepare($statement);
  $sth->bind_param( 1, $term->dbID(), SQL_INTEGER );

  $sth->execute();

  my ( $dbid, $accession, $name, $definition );
  $sth->bind_columns( \( $dbid, $accession, $name, $definition ) );

  my @terms;
  while ( $sth->fetch() ) {
    push(
      @terms,
      Bio::EnsEMBL::OntologyTerm->new(
        '-dbid'       => $dbid,
        '-adaptor'    => $this,
        '-accession'  => $accession,
        '-namespace'  => $term->{'namespace'},
        '-name'       => $name,
        '-definition' => $definition,
      ) );
  }

  return \@terms;
} ## end sub fetch_all_by_parent_term

=head2 fetch_by_child_term

  Arg [1]       : Bio::EnsEMBL::OntologyTerm
                  The term whose parent terms should be fetched.

  Description   : Given a child ontology term, returns a list of
                  its immediate parent terms.

  Example       : my @parents =
                    @{ $ot_adaptor->fetch_by_child_term($term) };

  Return type   : listref of Bio::EnsEMBL::OntologyTerm

=cut

sub fetch_by_child_term {
  my ( $this, $term ) = @_;

  if ( !ref($term) || !$term->isa('Bio::EnsEMBL::OntologyTerm') ) {
    throw('Argument needs to be a Bio::EnsEMBL::OntologyTerm object');
  }

  my @terms;

  if ( !$term->{'parent_terms_fetched'} ) {
    my $statement = q(
SELECT  parent_term.term_id,
        parent_term.accession,
        parent_term.name,
        parent_term.definition,
        rt.name
FROM    term parent_term,
        relation,
        relation_type rt
WHERE   relation.child_term_id = ?
  AND   relation.parent_term_id = parent_term.term_id
  AND   relation.relation_type_id = rt.relation_type_id);

    my $sth = $this->prepare($statement);
    $sth->bind_param( 1, $term->dbID(), SQL_INTEGER );

    $sth->execute();

    my ( $dbid, $accession, $name, $definition, $relation );
    $sth->bind_columns(
      \( $dbid, $accession, $name, $definition, $relation ) );

    while ( $sth->fetch() ) {
      my $parent_term = Bio::EnsEMBL::OntologyTerm->new(
        '-dbid'       => $dbid,
        '-adaptor'    => $this,
        '-accession'  => $accession,
        '-namespace'  => $term->{'namespace'},
        '-name'       => $name,
        '-definition' => $definition,
      );

      push( @terms,                             $parent_term );
      push( @{ $term->{'parents'}{$relation} }, $parent_term );
    }

    $term->{'parent_terms_fetched'} = 1;
  } else {
    foreach my $relation ( values( %{ $term->{'parents'} } ) ) {
      push( @terms, @{$relation} );
    }
  }

  return \@terms;
} ## end sub fetch_by_child_term

=head2 fetch_all_by_child_term

  Arg [1]       : Bio::EnsEMBL::OntologyTerm
                  The term whose ancestor terms should be fetched.

  Description   : Given a child ontology term, returns a list of all
                  its ancestor terms, up to and including any root
                  term.  In GO, relations of the type 'is_a' and
                  'part_of' are followed, but not 'regulates' etc.

  Example       : my @ancestors =
                    @{ $ot_adaptor->fetch_all_by_parent_term($term) };

  Return type   : listref of Bio::EnsEMBL::OntologyTerm

=cut

sub fetch_all_by_child_term {
  my ( $this, $term ) = @_;

  if ( !ref($term) || !$term->isa('Bio::EnsEMBL::OntologyTerm') ) {
    throw('Argument needs to be a Bio::EnsEMBL::OntologyTerm object');
  }

  my $statement = q(
SELECT DISTINCT
        parent_term.term_id,
        parent_term.accession,
        parent_term.name,
        parent_term.definition
FROM    term parent_term,
        closure
WHERE   closure.child_term_id = ?
  AND   closure.parent_term_id = parent_term.term_id
  AND   closure.distance > 0
ORDER BY closure.distance, parent_term.accession);

  my $sth = $this->prepare($statement);
  $sth->bind_param( 1, $term->dbID(), SQL_INTEGER );

  $sth->execute();

  my ( $dbid, $accession, $name, $definition );
  $sth->bind_columns( \( $dbid, $accession, $name, $definition ) );

  my @terms;
  while ( $sth->fetch() ) {
    push(
      @terms,
      Bio::EnsEMBL::OntologyTerm->new(
        '-dbid'       => $dbid,
        '-adaptor'    => $this,
        '-accession'  => $accession,
        '-namespace'  => $term->{'namespace'},
        '-name'       => $name,
        '-definition' => $definition,
      ) );
  }

  return \@terms;
} ## end sub fetch_all_by_child_term


sub fetch_ancestor_chart {
  my ( $this, $term ) = @_;

  my $statement = q(
SELECT  subparent_term.accession,
        parent_term.accession,
        relation_type.name
FROM    closure,
        relation,
        relation_type,
        term subparent_term,
        term parent_term
WHERE   closure.child_term_id = ?
  AND   relation.parent_term_id = closure.parent_term_id
  AND   relation.child_term_id = closure.subparent_term_id
  AND   relation.relation_type_id = relation_type.relation_type_id
  AND   subparent_term.term_id = closure.subparent_term_id
  AND   parent_term.term_id = closure.parent_term_id
ORDER BY closure.distance);

  my $sth = $this->prepare($statement);
  $sth->bind_param( 1, $term->dbID(), SQL_INTEGER );

  $sth->execute();

  my ( $subparent_id, $parent_id, $relation );
  $sth->bind_columns( \( $subparent_id, $parent_id, $relation ) );

  my %chart;
  while ( $sth->fetch() ) {
    if ( !exists( $chart{$parent_id} ) ) {
      $chart{$parent_id} = {};
    }
    push( @{ $chart{$subparent_id}{$relation} }, $parent_id );
  }

  return \%chart;
} ## end sub fetch_ancestor_chart

#-----------------------------------------------------------------------
# Useful public methods that implement functionality not properly
# provided by the parent class Bio::EnsEMBL::DBSQL::BaseAdaptor.

sub fetch_by_dbID {
  my ( $this, $dbid ) = @_;

  my $statement = q(
SELECT  term.accession,
        term.name,
        term.definition,
        ontology.namespace
FROM    ontology,
        term
WHERE   ontology.ontology_id = term.ontology_id
  AND   term.term_id = ?);

  my $sth = $this->prepare($statement);
  $sth->bind_param( 1, $dbid, SQL_INTEGER );

  $sth->execute();

  my ( $accession, $name, $definition, $namespace );
  $sth->bind_columns( \( $accession, $name, $definition, $namespace ) );

  $sth->fetch();
  my $term = Bio::EnsEMBL::OntologyTerm->new(
    '-dbid'       => $dbid,
    '-adaptor'    => $this,
    '-accession'  => $accession,
    '-namespace'  => $namespace,
    '-name'       => $name,
    '-definition' => $definition
  );
  $sth->finish();

  return $term;
} ## end sub fetch_by_dbID

sub fetch_by_dbID_list {
  my ( $this, $dbids ) = @_;

  if ( !@{$dbids} ) { return [] }

  my $stmt = q(
SELECT  term.term_id,
        term.accession,
        term.name,
        term.definition,
        ontology.namespace
FROM    ontology,
        term
WHERE   ontology.ontology_id = term.ontology_id
  AND   term.term_id IN (%s));

  my $statement = sprintf(
    $stmt,
    join( ',',
      map { $this->dbc()->db_handle()->quote( $_, SQL_INTEGER ) }
        @{$dbids} ) );

  my $sth = $this->prepare($statement);

  $sth->execute();

  my ( $dbid, $accession, $name, $definition, $namespace );
  $sth->bind_columns(
    \( $dbid, $accession, $name, $definition, $namespace ) );

  my @terms;
  while ( $sth->fetch() ) {
    push(
      @terms,
      Bio::EnsEMBL::OntologyTerm->new(
        '-dbid'       => $dbid,
        '-adaptor'    => $this,
        '-accession'  => $accession,
        '-namespace'  => $namespace,
        '-name'       => $name,
        '-definition' => $definition
      ) );
  }

  return \@terms;
} ## end sub fetch_by_dbID_list

1;
