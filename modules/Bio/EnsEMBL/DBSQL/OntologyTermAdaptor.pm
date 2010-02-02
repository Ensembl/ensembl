
=head1 LICENSE

  Copyright (c) 1999-2010 The European Bioinformatics Institute and
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

  my $term = $goa->fetch_by_accession('GO:0010885');

  my @children    = @{ $goa->fetch_all_by_parent_term($term) };
  my @descendants = @{ $goa->fetch_all_by_ancestor_term($term) };

  my @parents   = @{ $goa->fetch_all_by_child_term($term) };
  my @ancestors = @{ $goa->fetch_all_by_descendant_term($term) };

  my %ancestor_chart = %{ $goa->_fetch_ancestor_chart($term) };

=head1 DESCRIPTION

An abstract adaptor class for fetching ontology
terms, creates Bio::EnsEMBL::OntologyTerm objects.
Specialized by Bio::EnsEMBL::DBSQL::GOTermAdaptor and
Bio::EnsEMBL::DBSQL::SOTermAdaptor

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::OntologyTermAdaptor;

use strict;
use warnings;

use DBI qw( :sql_types );

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

  Caller        : Bio::EnsEMBL::DBSQL::GOTermAdaptor
                  Bio::EnsEMBL::DBSQL::SOTermAdaptor

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
        term.subsets,
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

  my ( $dbid, $name, $definition, $subsets, $namespace );
  $sth->bind_columns(
    \( $dbid, $name, $definition, $subsets, $namespace ) );

  $sth->fetch();
  $subsets ||= '';
  my $term = Bio::EnsEMBL::OntologyTerm->new(
    '-dbid'       => $dbid,
    '-adaptor'    => $this,
    '-accession'  => $accession,
    '-namespace'  => $namespace,
    '-subsets'    => [ split( /,/, $subsets ) ],
    '-name'       => $name,
    '-definition' => $definition
  );
  $sth->finish();

  return $term;
} ## end sub fetch_by_accession

=head2 fetch_all_by_parent_term

  Arg [1]       : Bio::EnsEMBL::OntologyTerm
                  The term whose children terms should be fetched.

  Description   : Given a parent ontology term, returns a list of
                  its immediate children terms.

  Example       :

    my @children =
      @{ $ot_adaptor->fetch_all_by_parent_term($term) };

  Return type   : listref of Bio::EnsEMBL::OntologyTerm

=cut

sub fetch_all_by_parent_term {
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
        child_term.subsets,
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

    my ( $dbid, $accession, $name, $definition, $subsets, $relation );
    $sth->bind_columns(
      \( $dbid, $accession, $name, $definition, $subsets, $relation ) );

    while ( $sth->fetch() ) {
      $subsets ||= '';

      my $child_term = Bio::EnsEMBL::OntologyTerm->new(
        '-dbid'       => $dbid,
        '-adaptor'    => $this,
        '-accession'  => $accession,
        '-namespace'  => $term->{'namespace'},
        '-subsets'    => [ split( /,/, $subsets ) ],
        '-name'       => $name,
        '-definition' => $definition,
      );

      push( @terms,                              $child_term );
      push( @{ $term->{'children'}{$relation} }, $child_term );
    }

    $sth->finish();

    $term->{'child_terms_fetched'} = 1;
  } else {
    foreach my $relation ( values( %{ $term->{'children'} } ) ) {
      push( @terms, @{$relation} );
    }
  }

  return \@terms;
} ## end sub fetch_all_by_parent_term

=head2 fetch_all_by_ancestor_term

  Arg [1]       : Bio::EnsEMBL::OntologyTerm
                  The term whose descendant terms should be fetched.

  Description   : Given a parent ontology term, returns a list of
                  all its descendant terms, down to and including
                  any leaf terms.  Relations of the type 'is_a' and
                  'part_of' are followed.

  Example       :

    my @descendants =
      @{ $ot_adaptor->fetch_all_by_ancestor_term($term) };

  Return type   : listref of Bio::EnsEMBL::OntologyTerm

=cut

sub fetch_all_by_ancestor_term {
  my ( $this, $term ) = @_;

  if ( !ref($term) || !$term->isa('Bio::EnsEMBL::OntologyTerm') ) {
    throw('Argument needs to be a Bio::EnsEMBL::OntologyTerm object');
  }

  my $statement = q(
SELECT DISTINCT
        child_term.term_id,
        child_term.accession,
        child_term.name,
        child_term.definition,
        child_term.subsets
FROM    term child_term,
        closure
WHERE   closure.child_term_id = child_term.term_id
  AND   closure.parent_term_id = ?
  AND   closure.distance > 0
ORDER BY closure.distance, child_term.accession);

  my $sth = $this->prepare($statement);
  $sth->bind_param( 1, $term->dbID(), SQL_INTEGER );

  $sth->execute();

  my ( $dbid, $accession, $name, $definition, $subsets );
  $sth->bind_columns(
    \( $dbid, $accession, $name, $definition, $subsets ) );

  my @terms;

  while ( $sth->fetch() ) {
    $subsets ||= '';

    push(
      @terms,
      Bio::EnsEMBL::OntologyTerm->new(
        '-dbid'       => $dbid,
        '-adaptor'    => $this,
        '-accession'  => $accession,
        '-namespace'  => $term->{'namespace'},
        '-subsets'    => [ split( /,/, $subsets ) ],
        '-name'       => $name,
        '-definition' => $definition,
      ) );
  }

  $sth->finish();

  return \@terms;
} ## end sub fetch_all_by_ancestor_term

=head2 fetch_all_by_child_term

  Arg [1]       : Bio::EnsEMBL::OntologyTerm
                  The term whose parent terms should be fetched.

  Description   : Given a child ontology term, returns a list of
                  its immediate parent terms.

  Example       :

    my @parents = @{ $ot_adaptor->fetch_all_by_child_term($term) };

  Return type   : listref of Bio::EnsEMBL::OntologyTerm

=cut

sub fetch_all_by_child_term {
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
        parent_term.subsets,
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

    my ( $dbid, $accession, $name, $definition, $subsets, $relation );
    $sth->bind_columns(
      \( $dbid, $accession, $name, $definition, $subsets, $relation ) );

    while ( $sth->fetch() ) {
      $subsets ||= '';

      my $parent_term = Bio::EnsEMBL::OntologyTerm->new(
        '-dbid'       => $dbid,
        '-adaptor'    => $this,
        '-accession'  => $accession,
        '-namespace'  => $term->{'namespace'},
        '-subsets'    => [ split( /,/, $subsets ) ],
        '-name'       => $name,
        '-definition' => $definition,
      );

      push( @terms,                             $parent_term );
      push( @{ $term->{'parents'}{$relation} }, $parent_term );
    }

    $sth->finish();

    $term->{'parent_terms_fetched'} = 1;
  } else {
    foreach my $relation ( values( %{ $term->{'parents'} } ) ) {
      push( @terms, @{$relation} );
    }
  }

  return \@terms;
} ## end sub fetch_all_by_child_term

=head2 fetch_all_by_descendant_term

  Arg [1]       : Bio::EnsEMBL::OntologyTerm
                  The term whose ancestor terms should be fetched.

  Arg [2]       : (optional) String
                  The subset within the ontolgy to which the query
                  should be restricted.  The subset may be specified as
                  a SQL pattern, e.g., "%goslim%" (but "goslim%" might
                  not do what you expect), or as a specific subset name,
                  e.g., "goslim_goa".

  Arg [3]       : (optional) Boolean
                  If true (non-zero), only return the closest
                  term(s).  If this argument is true, and the
                  previous argument is left undefined, this method
                  will return the parent(s) of the given term.

  Description   : Given a child ontology term, returns a list of
                  all its ancestor terms, up to and including any
                  root term.  Relations of the type 'is_a' and
                  'part_of' are followed.  Optionally, only terms in
                  a given subset of the ontology may be returned,
                  and additionally one may ask to only get the
                  closest term(s) to the given child term.

  Example       :

    my @ancestors =
      @{ $ot_adaptor->fetch_all_by_descendant_term($term) };

  Return type   : listref of Bio::EnsEMBL::OntologyTerm

=cut

sub fetch_all_by_descendant_term {
  my ( $this, $term, $subset, $closest_only ) = @_;

  if ( !ref($term) || !$term->isa('Bio::EnsEMBL::OntologyTerm') ) {
    throw('Argument needs to be a Bio::EnsEMBL::OntologyTerm object');
  }

  $closest_only ||= 0;

  my $statement = q(
SELECT DISTINCT
        parent_term.term_id,
        parent_term.accession,
        parent_term.name,
        parent_term.definition,
        parent_term.subsets,
        closure.distance
FROM    term parent_term,
        closure
WHERE   closure.child_term_id = ?
  AND   closure.parent_term_id = parent_term.term_id
  AND   closure.distance > 0);

  if ( defined($subset) ) {
    if ( index( $subset, '%' ) != -1 ) {
      $statement .= q(
  AND   parent_term.subsets LIKE ?);
    } else {
      $statement .= q(
  AND   FIND_IN_SET(?, parent_term.subsets) > 0);
    }
  }

  $statement .= q(
ORDER BY closure.distance, parent_term.accession);

  my $sth = $this->prepare($statement);
  $sth->bind_param( 1, $term->dbID(), SQL_INTEGER );

  if ( defined($subset) ) {
    $sth->bind_param( 2, $subset, SQL_VARCHAR );
  }

  $sth->execute();

  my ( $dbid, $accession, $name, $definition, $subsets, $distance );
  $sth->bind_columns(
    \( $dbid, $accession, $name, $definition, $subsets, $distance ) );

  my @terms;
  my $min_distance;

  while ( $sth->fetch() ) {
    $subsets ||= '';
    $min_distance ||= $distance;

    if ( !$closest_only || $distance == $min_distance ) {
      push(
        @terms,
        Bio::EnsEMBL::OntologyTerm->new(
          '-dbid'       => $dbid,
          '-adaptor'    => $this,
          '-accession'  => $accession,
          '-namespace'  => $term->{'namespace'},
          '-subsets'    => [ split( /,/, $subsets ) ],
          '-name'       => $name,
          '-definition' => $definition,
        ) );
    } else {
      last;
    }
  }

  $sth->finish();

  return \@terms;
} ## end sub fetch_all_by_descendant_term

=head2 _fetch_ancestor_chart

  Arg [1]       : Bio::EnsEMBL::OntologyTerm
                  The term whose ancestor terms should be fetched.

  Description   : Given a child ontology term, returns a hash
                  structure containing its ancestor terms, up to and
                  including any root term.  Relations of the type
                  'is_a' and 'part_of' are included.

  Example       :

    my %chart = %{ $ot_adaptor->_fetch_ancestor_chart($term) };

  Return type   : A reference to a hash structure like this:

    {
      'GO:XXXXXXX' => {
        'term' =>           # ref to Bio::EnsEMBL::OntologyTerm object
        'is_a'    => [...], # listref of Bio::EnsEMBL::OntologyTerm
        'part_of' => [...], # listref of Bio::EnsEMBL::OntologyTerm
      },
      'GO:YYYYYYY' => {
        # Similarly for all ancestors,
        # and including the query term itself.
      }
    }

=cut

sub _fetch_ancestor_chart {
  my ( $this, $term ) = @_;

  if ( !ref($term) || !$term->isa('Bio::EnsEMBL::OntologyTerm') ) {
    throw('Argument needs to be a Bio::EnsEMBL::OntologyTerm object');
  }

  my $statement = q(
SELECT  subparent_term.term_id,
        parent_term.term_id,
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

  my %id_chart;
  my %acc_chart;

  while ( $sth->fetch() ) {
    if ( !exists( $id_chart{$parent_id} ) ) {
      $id_chart{$parent_id} = {};
    }
    push( @{ $id_chart{$subparent_id}{$relation} }, $parent_id );
  }

  $sth->finish();

  my @terms = @{ $this->fetch_all_by_dbID_list( [ keys(%id_chart) ] ) };

  foreach my $term (@terms) {
    $id_chart{ $term->dbID() }{'term'}       = $term;
    $acc_chart{ $term->accession() }{'term'} = $term;
  }

  foreach my $term (@terms) {
    my $accession = $term->accession();
    my $dbID      = $term->dbID();

    foreach my $relation ( keys( %{ $id_chart{$dbID} } ) ) {
      if ( $relation eq 'term' ) { next }

      foreach my $id ( @{ $id_chart{$dbID}{$relation} } ) {
        push(
          @{ $acc_chart{$accession}{$relation} },
          $id_chart{$id}{'term'} );
      }
    }
  }

  return \%acc_chart;
} ## end sub _fetch_ancestor_chart

#-----------------------------------------------------------------------
# Useful public methods that implement functionality not properly
# provided by the parent class Bio::EnsEMBL::DBSQL::BaseAdaptor.

sub fetch_by_dbID {
  my ( $this, $dbid ) = @_;

  my $statement = q(
SELECT  term.accession,
        term.name,
        term.definition,
        term.subsets,
        ontology.namespace
FROM    ontology,
        term
WHERE   ontology.ontology_id = term.ontology_id
  AND   term.term_id = ?);

  my $sth = $this->prepare($statement);
  $sth->bind_param( 1, $dbid, SQL_INTEGER );

  $sth->execute();

  my ( $accession, $name, $definition, $subsets, $namespace );
  $sth->bind_columns(
    \( $accession, $name, $definition, $subsets, $namespace ) );

  $sth->fetch();
  $subsets ||= '';
  my $term = Bio::EnsEMBL::OntologyTerm->new(
    '-dbid'       => $dbid,
    '-adaptor'    => $this,
    '-accession'  => $accession,
    '-namespace'  => $namespace,
    '-subsets'    => [ split( /,/, $subsets ) ],
    '-name'       => $name,
    '-definition' => $definition
  );
  $sth->finish();

  return $term;
} ## end sub fetch_by_dbID

sub fetch_all_by_dbID_list {
  my ( $this, $dbids ) = @_;

  if ( !@{$dbids} ) { return [] }

  my $stmt = q(
SELECT  term.term_id,
        term.accession,
        term.name,
        term.definition,
        term.subsets,
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

  my ( $dbid, $accession, $name, $definition, $subsets, $namespace );
  $sth->bind_columns(
    \( $dbid, $accession, $name, $definition, $subsets, $namespace ) );

  my @terms;

  while ( $sth->fetch() ) {
    $subsets ||= '';

    push(
      @terms,
      Bio::EnsEMBL::OntologyTerm->new(
        '-dbid'       => $dbid,
        '-adaptor'    => $this,
        '-accession'  => $accession,
        '-namespace'  => $namespace,
        '-subsets'    => [ split( /,/, $subsets ) ],
        '-name'       => $name,
        '-definition' => $definition
      ) );
  }

  $sth->finish();

  return \@terms;
} ## end sub fetch_all_by_dbID_list

sub fetch_all {
  my ($this) = @_;

  my $statement = q(
SELECT  term.term_id,
        term.accession,
        term.name,
        term.definition,
        term.subsets,
        ontology.namespace
FROM    ontology,
        term
WHERE   ontology.ontology_id = term.ontology_id
  AND   ontology.name = ?);

  my $sth = $this->prepare($statement);

  $sth->bind_param( 1, $this->{'ontology'}, SQL_VARCHAR );

  $sth->execute();

  my ( $dbid, $accession, $name, $definition, $subsets, $namespace );
  $sth->bind_columns(
    \( $dbid, $accession, $name, $definition, $subsets, $namespace ) );

  my @terms;

  while ( $sth->fetch() ) {
    $subsets ||= '';

    push(
      @terms,
      Bio::EnsEMBL::OntologyTerm->new(
        '-dbid'       => $dbid,
        '-adaptor'    => $this,
        '-accession'  => $accession,
        '-namespace'  => $namespace,
        '-subsets'    => [ split( /,/, $subsets ) ],
        '-name'       => $name,
        '-definition' => $definition
      ) );
  }

  $sth->finish();

  return \@terms;
} ## end sub fetch_all

1;
