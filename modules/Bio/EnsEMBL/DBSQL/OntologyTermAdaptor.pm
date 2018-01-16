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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::OntologyTermAdaptor

=head1 SYNOPSIS

  my $goa =
    $registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );

  my $term = $goa->fetch_by_accession('GO:0010885');

  my @children    = @{ $goa->fetch_all_by_parent_term($term) };
  my @descendants = @{ $goa->fetch_all_by_ancestor_term($term) };

  my @parents   = @{ $goa->fetch_all_by_child_term($term) };
  my @ancestors = @{ $goa->fetch_all_by_descendant_term($term) };

  my %ancestor_chart = %{ $goa->_fetch_ancestor_chart($term) };

=head1 DESCRIPTION

An abstract adaptor class for fetching ontology
terms, creates Bio::EnsEMBL::OntologyTerm objects.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::OntologyTermAdaptor;

use strict;
use warnings;

use DBI qw( :sql_types );

use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );

use Bio::EnsEMBL::OntologyTerm;

use base qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );

=head2 fetch_all_by_name

  Arg [1]       : String, name of term, or SQL pattern
  Arg [2]       : (optional) String, name of ontology
  Arg [3]       : (optional) Boolean, search through obsolete terms as well

  Description   : Fetches ontology term(s) given a name, a synonym, or a
                  SQL pattern like "%splice_site%"

  Example       :

    my ($term) =
      @{ $ot_adaptor->fetch_by_name( 'DNA_binding_site', 'SO' ) };

    # Will find terms in both SO and GO:
    my @terms = @{ $ot_adaptor->fetch_by_name('%splice_site%') };

  Return type   : listref of Bio::EnsEMBL::OntologyTerm

=cut

sub fetch_all_by_name {
  my ( $this, $pattern, $ontology, $include_obsolete ) = @_;

  my $statement = q(
SELECT DISTINCT
        term.term_id,
        term.accession,
        term.name,
        term.definition,
        term.subsets,
        term.is_root,
        term.is_obsolete,
        ontology.name,
        ontology.namespace,
        ontology.data_version
FROM    ontology
  JOIN  term USING (ontology_id)
  LEFT JOIN  synonym USING (term_id)
WHERE   ( term.name LIKE ? OR synonym.name LIKE ? ));

  if ( defined($ontology) ) {
    $statement .= " AND ontology.name = ?";
  }
  $statement .= " AND term.is_obsolete = 0" unless $include_obsolete;

  my $sth = $this->prepare($statement);
  $sth->bind_param( 1, $pattern, SQL_VARCHAR );
  $sth->bind_param( 2, $pattern, SQL_VARCHAR );

  if ( defined($ontology) ) {
    $sth->bind_param( 3, $ontology, SQL_VARCHAR );
  }

  $sth->execute();

  my ( $dbid, $accession, $name, $definition, $subsets, $is_root, $is_obsolete, $namespace, $ontology_version );
  $sth->bind_columns(
     \( $dbid, $accession, $name, $definition, $subsets, $is_root, $is_obsolete, $ontology, $namespace, $ontology_version ) );

  my @terms;

  while ( $sth->fetch() ) {
    $subsets ||= '';

    push @terms,
      Bio::EnsEMBL::OntologyTerm->new(
                               '-dbid'        => $dbid,
                               '-adaptor'     => $this,
                               '-accession'   => $accession,
                               '-is_root'     => $is_root,
                               '-is_obsolete' => $is_obsolete,
                               '-ontology'    => $ontology,
                               '-ontology_version' => $ontology_version,
                               '-namespace'   => $namespace,
                               '-subsets'     => [ split( /,/, $subsets ) ],
                               '-name'        => $name,
                               '-definition'  => $definition,
                               '-synonyms'    => $this->_fetch_synonyms_by_dbID($dbid)
    );

  }

  return \@terms;
} ## end sub fetch_all_by_name


=head2 fetch_by_accession

  Arg [1]       : String
  Arg [2]       : (optional) Boolean, search through obsolete terms as well

  Description   : Fetches an ontology term given an accession.

  Example       :

    my $term = $ot_adaptor->fetch_by_accession('GO:0030326');

  Return type   : Bio::EnsEMBL::OntologyTerm

=cut

sub fetch_by_accession {
  my ( $this, $accession, $include_obsolete ) = @_;

  my $statement = q(
SELECT  term.term_id,
        term.accession,
        term.name,
        term.definition,
        term.subsets,
        term.is_root,
        term.is_obsolete,
        ontology.name,
        ontology.namespace,
        ontology.data_version
FROM    ontology
  JOIN  term USING (ontology_id)
WHERE   term.accession = ?);

  $statement .= " AND term.is_obsolete = 0" unless $include_obsolete;

  my $sth = $this->prepare($statement);
  $sth->bind_param( 1, $accession, SQL_VARCHAR );

  $sth->execute();

  my ( $dbid, $name, $definition, $subsets, $is_root, $is_obsolete, $ontology, $namespace, $ontology_version );
  $sth->bind_columns(
      \( $dbid, $accession, $name, $definition, $subsets, $is_root, $is_obsolete, $ontology, $namespace, $ontology_version ) );

  $sth->fetch();
  $sth->finish();

  my $term;
  if (!$dbid) {
    $term = $this->fetch_by_alt_id($accession);
  } else {
    $subsets ||= '';

    $term =
      Bio::EnsEMBL::OntologyTerm->new(
                    '-dbid'       => $dbid,
                    '-adaptor'    => $this,
                    '-accession'  => $accession,
                    '-is_root'    => $is_root,
                    '-is_obsolete'=> $is_obsolete,
                    '-ontology'   => $ontology,
                    '-ontology_version' => $ontology_version,
                    '-namespace'  => $namespace,
                    '-subsets'    => [ split( /,/, $subsets ) ],
                    '-name'       => $name,
                    '-definition' => $definition,
                    '-synonyms' => $this->_fetch_synonyms_by_dbID($dbid)
      );
  }

  return $term;
} ## end sub fetch_by_accession


=head2 fetch_by_alt_id

  Arg [1]       : String

  Description   : Fetches an ontology term given an alt_id.

  Example       :

    my $term = $ot_adaptor->fetch_by_alt_id('GO:0019952');

  Return type   : Bio::EnsEMBL::OntologyTerm

=cut

sub fetch_by_alt_id {
  my ( $this, $accession ) = @_;

  my $statement = q(
SELECT  term.term_id,
        alt_id.accession,
        term.name,
        term.definition,
        term.subsets,
        term.is_root,
        term.is_obsolete,
        ontology.name,
        ontology.namespace,
        ontology.data_version
FROM    ontology
  JOIN  term USING (ontology_id)
  JOIN  alt_id USING (term_id)
WHERE   alt_id.accession = ?);

  my $sth = $this->prepare($statement);
  $sth->bind_param( 1, $accession, SQL_VARCHAR );

  $sth->execute();

  my ( $dbid, $name, $definition, $subsets, $is_root, $is_obsolete, $ontology, $namespace, $ontology_version );
  $sth->bind_columns(
      \( $dbid, $accession, $name, $definition, $subsets, $is_root, $is_obsolete, $ontology, $namespace, $ontology_version ) );

  $sth->fetch();

  # early exit in the event of bad $accession
  unless ($dbid) {return;}

  $subsets ||= '';
  my $term =
    Bio::EnsEMBL::OntologyTerm->new(
                    '-dbid'       => $dbid,
                    '-adaptor'    => $this,
                    '-accession'  => $accession,
                    '-ontology'   => $ontology,
                    '-ontology_version' => $ontology_version,
                    '-namespace'  => $namespace,
                    '-subsets'    => [ split( /,/, $subsets ) ],
                    '-name'       => $name,
                    '-definition' => $definition,
                    '-synonyms' => $this->_fetch_synonyms_by_dbID($dbid)
    );
  $sth->finish();

  return $term;
} ## end sub fetch_by_alt_id



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
  my ( $this, $term, $ontology, $include_obsolete ) = @_;

  assert_ref( $term, 'Bio::EnsEMBL::OntologyTerm' );

  my @terms;

  if ( !$term->{'child_terms_fetched'} ) {
    my $statement = q(
SELECT  child_term.term_id,
        child_term.accession,
        child_term.name,
        child_term.definition,
        child_term.subsets,
        child_term.is_root,
        child_term.is_obsolete,
        rt.name,
        ontology.data_version
FROM    term child_term
  JOIN  relation ON (relation.child_term_id = child_term.term_id)
  JOIN  relation_type rt USING (relation_type_id)
  JOIN  ontology ON (ontology.ontology_id = relation.ontology_id)
WHERE   relation.parent_term_id = ?
   AND  ontology.name = ?);

    $statement .= " AND child_term.is_obsolete = 0" unless $include_obsolete;

    my $sth = $this->prepare($statement);
    $sth->bind_param( 1, $term->dbID(), SQL_INTEGER );

    if (!defined $ontology) {
      $ontology = $term->{'ontology'};
    }
    $sth->bind_param( 2, $ontology, SQL_VARCHAR );


    $sth->execute();

    my ( $dbid, $accession, $name, $definition, $subsets, $is_root, $is_obsolete, $relation, $ontology_version );
    $sth->bind_columns(
      \( $dbid, $accession, $name, $definition, $subsets, $is_root, $is_obsolete, $relation, $ontology_version ) );

    while ( $sth->fetch() ) {
      $subsets ||= '';

      my $child_term =
        Bio::EnsEMBL::OntologyTerm->new(
                               '-dbid'        => $dbid,
                               '-adaptor'     => $this,
                               '-accession'   => $accession,
                               '-is_root'     => $is_root,
                               '-is_obsolete' => $is_obsolete,
                               '-ontology'    => $term->{'ontology'},
                               '-ontology_version' => $ontology_version,
                               '-namespace'   => $term->{'namespace'},
                               '-subsets'     => [ split( /,/, $subsets ) ],
                               '-name'        => $name,
                               '-definition'  => $definition, 
                               '-synonyms'    => $this->_fetch_synonyms_by_dbID($dbid)
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
  my ( $this, $term, $ontology ) = @_;

  assert_ref( $term, 'Bio::EnsEMBL::OntologyTerm' );

  my $statement = q(
SELECT DISTINCT
        child_term.term_id,
        child_term.accession,
        child_term.name,
        child_term.definition,
        child_term.subsets,
        child_term.is_root,
        child_term.is_obsolete,
        ontology.data_version,
        closure.distance
FROM    term child_term
  JOIN  closure ON (closure.child_term_id = child_term.term_id)
  JOIN  ontology ON (closure.ontology_id = ontology.ontology_id)
WHERE   closure.parent_term_id = ?
  AND   closure.distance > 0
  AND   closure.ontology_id = child_term.ontology_id
  AND   ontology.name = ?
ORDER BY closure.distance, child_term.accession);

  my $sth = $this->prepare($statement);
  $sth->bind_param( 1, $term->dbID(), SQL_INTEGER );
  if (!defined $ontology) {
    $ontology = $term->{'ontology'};
  }
  $sth->bind_param( 2, $ontology, SQL_VARCHAR );
  $sth->execute();

  my ( $dbid, $accession, $name, $definition, $subsets, $is_root, $is_obsolete, $ontology_version, $closure_distance );
  $sth->bind_columns(
                 \( $dbid, $accession, $name, $definition, $subsets, $is_root, $is_obsolete, $ontology_version, $closure_distance ) );

  my @terms;

  while ( $sth->fetch() ) {
    $subsets ||= '';

    push( @terms,
          Bio::EnsEMBL::OntologyTerm->new(
                               '-dbid'        => $dbid,
                               '-adaptor'     => $this,
                               '-accession'   => $accession,
                               '-is_root'     => $is_root,
                               '-is_obsolete' => $is_obsolete,
                               '-ontology'    => $term->{'ontology'},
                               '-ontology_version' => $ontology_version,
                               '-namespace'   => $term->{'namespace'},
                               '-subsets'     => [ split( /,/, $subsets ) ],
                               '-name'        => $name,
                               '-definition'  => $definition,
                               '-synonyms'    => $this->_fetch_synonyms_by_dbID($dbid)
    ) );
  }

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
  my ( $this, $term, $ontology ) = @_;

  assert_ref( $term, 'Bio::EnsEMBL::OntologyTerm' );

  my @terms;

  if ( !$term->{'parent_terms_fetched'} ) {
    my $statement = q(
SELECT  parent_term.term_id,
        parent_term.accession,
        parent_term.name,
        parent_term.definition,
        parent_term.subsets,
        parent_term.is_root,
        parent_term.is_obsolete,
        rt.name,
        ontology.data_version
FROM    term parent_term
  JOIN  relation ON (relation.parent_term_id = parent_term.term_id)
  JOIN  relation_type rt USING (relation_type_id)
  JOIN  ontology ON (ontology.ontology_id = relation.ontology_id)
WHERE   relation.child_term_id = ?
   AND  ontology.name = ?);

    my $sth = $this->prepare($statement);
    $sth->bind_param( 1, $term->dbID(), SQL_INTEGER );

    if (!defined $ontology) {
      $ontology = $term->{'ontology'};
    }
    $sth->bind_param( 2, $ontology, SQL_VARCHAR );

    $sth->execute();

    my ( $dbid, $accession, $name, $definition, $subsets, $is_root, $is_obsolete, $relation, $ontology_version );
    $sth->bind_columns(
      \( $dbid, $accession, $name, $definition, $subsets, $is_root, $is_obsolete, $relation, $ontology_version ) );

    while ( $sth->fetch() ) {
      $subsets ||= '';

      my $parent_term =
        Bio::EnsEMBL::OntologyTerm->new(
                               '-dbid'        => $dbid,
                               '-adaptor'     => $this,
                               '-accession'   => $accession,
                               '-is_root'     => $is_root,
                               '-is_obsolete' => $is_obsolete,
                               '-ontology'    => $term->{'ontology'},
                               '-ontology_version' => $ontology_version,
                               '-namespace'   => $term->{'namespace'},
                               '-subsets'     => [ split( /,/, $subsets ) ],
                               '-name'        => $name,
                               '-definition'  => $definition,
                               '-synonyms'    => $this->_fetch_synonyms_by_dbID($dbid)
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
} ## end sub fetch_all_by_child_term

=head2 fetch_all_by_descendant_term

  Arg [1]       : Bio::EnsEMBL::OntologyTerm
                  The term whose ancestor terms should be fetched.

  Arg [2]       : (optional) String
                  The subset within the ontolgy to which the query
                  should be restricted.  The subset may be specified as
                  a SQL pattern, e.g., "%goslim%" (but "goslim%" might
                  not do what you expect), or as a specific subset name,
                  e.g., "goslim_generic".

  Arg [3]       : (optional) Boolean
                  If true (non-zero), only return the closest
                  term(s).  If this argument is true, and the
                  previous argument is left undefined, this method
                  will return the parent(s) of the given term.
  
  Arg [4]       : (optional) Boolean
                  If true we will allow the retrieval of terms whose distance
                  to the current term is 0. If false then we will only return
                  those which are above the current term in the ontology

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
  my ( $this, $term, $subset, $closest_only, $allow_zero_distance, $ontology ) = @_;

  assert_ref( $term, 'Bio::EnsEMBL::OntologyTerm' );

  $closest_only ||= 0;

  my $statement = q(
SELECT DISTINCT
        parent_term.term_id,
        parent_term.accession,
        parent_term.name,
        parent_term.definition,
        parent_term.subsets,
        parent_term.is_root,
        parent_term.is_obsolete,
        closure.distance,
        ontology.data_version
FROM    term parent_term
  JOIN  closure ON (closure.parent_term_id = parent_term.term_id)
  JOIN  ontology ON (closure.ontology_id = ontology.ontology_id)
WHERE   closure.child_term_id = ?
  AND   closure.distance > ?
  AND   closure.ontology_id = parent_term.ontology_id
  AND   ontology.name = ?);

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
  my $query_distance = ($allow_zero_distance) ? -1 : 0;
  $sth->bind_param( 2, $query_distance, SQL_INTEGER );
  if (!defined $ontology) {
    $ontology = $term->{'ontology'};
  }
  $sth->bind_param( 3, $ontology, SQL_VARCHAR );

  if ( defined($subset) ) {
    $sth->bind_param( 4, $subset, SQL_VARCHAR );
  }

  $sth->execute();

  my ( $dbid, $accession, $name, $definition, $subsets, $is_root, $is_obsolete, $distance, $ontology_version );
  $sth->bind_columns(
      \( $dbid, $accession, $name, $definition, $subsets, $is_root, $is_obsolete, $distance, $ontology_version ) );

  my @terms;
  my $min_distance;

  while ( $sth->fetch() ) {
    $subsets ||= '';
    $min_distance ||= $distance;

    if ( !$closest_only || $distance == $min_distance ) {
      push( @terms,
            Bio::EnsEMBL::OntologyTerm->new(
                               '-dbid'        => $dbid,
                               '-adaptor'     => $this,
                               '-accession'   => $accession,
                               '-is_root'     => $is_root,
                               '-is_obsolete' => $is_obsolete,
                               '-ontology'    => $term->{'ontology'},
                               '-ontology_version' => $ontology_version,
                               '-namespace'   => $term->{'namespace'},
                               '-subsets'     => [ split( /,/, $subsets ) ],
                               '-name'        => $name,
                               '-definition'  => $definition,
                               '-synonyms'    => $this->_fetch_synonyms_by_dbID($dbid)
    ) );
    } else {
      $sth->finish();
      last;
    }
  }

  return \@terms;
} ## end sub fetch_all_by_descendant_term

sub _fetch_synonyms_by_dbID {
  my ( $this, $dbID ) = @_;

  my $statement = q(
SELECT  synonym.name
FROM    synonym
WHERE   synonym.term_id = ?);

  my $sth = $this->prepare($statement);
  $sth->bind_param( 1, $dbID, SQL_INTEGER );

  $sth->execute();

  my $synonym;
  $sth->bind_col( 1, \$synonym );

  my @synonyms;
  while ( $sth->fetch() ) {
    push( @synonyms, $synonym );
  }

  return \@synonyms;
}



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
  my ( $this, $term, $ontology ) = @_;

  assert_ref( $term, 'Bio::EnsEMBL::OntologyTerm' );

  my $statement = q(
SELECT  subparent_term.term_id,
        parent_term.term_id,
        relation_type.name
FROM    closure
  JOIN  relation
    ON (relation.parent_term_id = closure.parent_term_id
      AND relation.child_term_id = closure.subparent_term_id
      AND closure.ontology_id = relation.ontology_id)
  JOIN  relation_type USING (relation_type_id)
  JOIN  term subparent_term
    ON (subparent_term.term_id = closure.subparent_term_id)
  JOIN  term parent_term ON (parent_term.term_id = closure.parent_term_id)
  JOIN  ontology ON (ontology.ontology_id = closure.ontology_id)
WHERE   closure.child_term_id = ?
   AND  ontology.name = ?
ORDER BY closure.distance);

  my $sth = $this->prepare($statement);
  $sth->bind_param( 1, $term->dbID(), SQL_INTEGER );
  if (!defined $ontology) {
    $ontology = $term->{'ontology'};
  }
  $sth->bind_param( 2, $ontology, SQL_VARCHAR );

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
        push( @{ $acc_chart{$accession}{$relation} },
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
  my ( $this, $dbid, $include_obsolete ) = @_;

  my $statement = q(
SELECT  term.term_id,
        term.accession,
        term.name,
        term.definition,
        term.subsets,
        term.is_root,
        term.is_obsolete,
        ontology.name,
        ontology.namespace,
        ontology.data_version
FROM    ontology
  JOIN  term USING (ontology_id)
WHERE   term.term_id = ?);

  $statement .= " AND term.is_obsolete = 0" unless $include_obsolete;

  my $sth = $this->prepare($statement);
  $sth->bind_param( 1, $dbid, SQL_INTEGER );

  $sth->execute();

  my ( $accession, $name, $definition, $subsets, $is_root, $is_obsolete, $ontology,
       $namespace, $ontology_version );
  $sth->bind_columns(
      \( $dbid, $accession, $name, $definition, $subsets, $is_root, $is_obsolete, $ontology, $namespace, $ontology_version
      ) );

  $sth->fetch();
  
  unless ($accession) {return;}
  
  $subsets ||= '';

  my $term =
    Bio::EnsEMBL::OntologyTerm->new(
                    '-dbid'        => $dbid,
                    '-adaptor'     => $this,
                    '-accession'   => $accession,
                    '-is_root'     => $is_root,
                    '-is_obsolete' => $is_obsolete,
                    '-ontology'    => $ontology,
                    '-ontology_version' => $ontology_version,
                    '-namespace'   => $namespace,
                    '-subsets'     => [ split( /,/, $subsets ) ],
                    '-name'        => $name,
                    '-definition'  => $definition,
                    '-synonyms'    => $this->_fetch_synonyms_by_dbID($dbid)
    );
  $sth->finish();

  return $term;
} ## end sub fetch_by_dbID

sub fetch_all_by_dbID_list {
  my ( $this, $dbids, $include_obsolete ) = @_;

  if ( !@{$dbids} ) { return [] }

  my $stmt = q(
SELECT  term.term_id,
        term.accession,
        term.name,
        term.definition,
        term.subsets,
        term.is_root,
        term.is_obsolete,
        ontology.name,
        ontology.namespace,
        ontology.data_version
FROM    ontology
  JOIN  term USING (ontology_id)
WHERE   term.term_id IN (%s));

  my $statement = sprintf(
    $stmt,
    join(
      ',',
      map {
        $this->dbc()->db_handle()->quote( $_, SQL_INTEGER )
        } @{$dbids} ) );

  $statement .= " AND term.is_obsolete = 0" unless $include_obsolete;

  my $sth = $this->prepare($statement);

  $sth->execute();

  my ( $dbid, $accession, $name, $definition, $subsets, $is_root, $is_obsolete, $ontology,
       $namespace, $ontology_version );
  $sth->bind_columns( \( $dbid,    $accession, $name, $definition,
                         $subsets, $is_root, $is_obsolete, $ontology,  $namespace, $ontology_version ) );

  my @terms;

  while ( $sth->fetch() ) {
    $subsets ||= '';

    push( @terms,
          Bio::EnsEMBL::OntologyTerm->new(
                               '-dbid'        => $dbid,
                               '-adaptor'     => $this,
                               '-accession'   => $accession,
                               '-is_root'     => $is_root,
                               '-is_obsolete' => $is_obsolete,
                               '-ontology'    => $ontology,
                               '-ontology_version' => $ontology_version,
                               '-namespace'   => $namespace,
                               '-subsets'     => [ split( /,/, $subsets ) ],
                               '-name'        => $name,
                               '-definition'  => $definition,
                               '-synonyms'    => $this->_fetch_synonyms_by_dbID($dbid)
    ) );
  }

  return \@terms;
} ## end sub fetch_all_by_dbID_list


=head2 fetch_all_alt_ids

  Arg [1]       : String

  Description   : Fetches all alt_ids for a given ontology term

  Example       :

    my ($accessions) =
      @{ $ot_adaptor->fetch_all_alt_ids( 'GO:0000003' ) };

  Return type   : listref of accessions

=cut
sub fetch_all_alt_ids {
  my ($this, $accession) = @_;

  my $statement = q(
SELECT  alt_id.accession
FROM    term
  JOIN  alt_id USING (term_id)
 WHERE  term.accession = ?);

  my $sth = $this->prepare($statement);
  $sth->bind_param(1, $accession, SQL_VARCHAR);
  $sth->execute();

  my (@accessions, $alt_id);
  $sth->bind_columns( \$alt_id);

  while ( $sth->fetch() ) {
    push( @accessions, $alt_id);
  }

  return \@accessions;
}




=head2 fetch_all_roots

  Arg [1]       : (optional) String, name of ontology

  Description   : Fetches all roots for all ontologies
                  Optionally, can be restricted to a given ontology

  Example       :

    my ($terms) =
      @{ $ot_adaptor->fetch_all_roots( 'SO' ) };

    # Will find terms in EFO, SO and GO:
    my @terms = @{ $ot_adaptor->fetch_all_roots() };

  Return type   : listref of Bio::EnsEMBL::OntologyTerm

=cut

sub fetch_all_roots {
  my ($this, $ontology_name) = @_;

  my $statement = q(
SELECT  term.term_id,
        term.accession,
        term.name,
        term.definition,
        term.subsets,
        term.is_root,
        term.is_obsolete,
        ontology.name,
        ontology.namespace,
        ontology.data_version
FROM    ontology
  JOIN  term USING (ontology_id)
 WHERE  is_root = 1);

  if (defined $ontology_name) {
    $statement .= " AND ontology.name = ?";
  }

  my $sth = $this->prepare($statement);
  if (defined $ontology_name) {
    $sth->bind_param( 1, $ontology_name, SQL_VARCHAR );
  }
  $sth->execute();

  my ( $dbid, $accession, $name, $definition, $subsets, $is_root, $is_obsolete, $ontology,
       $namespace, $ontology_version );
  $sth->bind_columns( \( $dbid,    $accession, $name, $definition,
                         $subsets, $is_root, $is_obsolete, $ontology,  $namespace, $ontology_version ) );

  my @terms;

  while ( $sth->fetch() ) {
    $subsets ||= '';

    push( @terms,
          Bio::EnsEMBL::OntologyTerm->new(
                               '-dbid'        => $dbid,
                               '-adaptor'     => $this,
                               '-accession'   => $accession,
                               '-is_root'     => $is_root,
                               '-is_obsolete' => $is_obsolete,
                               '-ontology'    => $ontology,
                               '-ontology_version' => $ontology_version,
                               '-namespace'   => $namespace,
                               '-subsets'     => [ split( /,/, $subsets ) ],
                               '-name'        => $name,
                               '-definition'  => $definition,
                               '-synonyms'    => $this->_fetch_synonyms_by_dbID($dbid)

    ));
  }

  return \@terms;
}


sub fetch_all {
  my ($this, $include_obsolete) = @_;

  my $statement = q(
SELECT  term.term_id,
        term.accession,
        term.name,
        term.definition,
        term.subsets,
        term.is_root,
        term.is_obsolete,
        ontology.name,
        ontology.namespace,
        ontology.data_version
FROM    ontology
  JOIN  term USING (ontology_id));

  $statement .= " WHERE term.is_obsolete = 0" unless $include_obsolete;

  my $sth = $this->prepare($statement);
  $sth->execute();

  my ( $dbid, $accession, $name, $definition, $subsets, $is_root, $is_obsolete, $ontology,
       $namespace, $ontology_version );
  $sth->bind_columns( \( $dbid,    $accession, $name, $definition,
                         $subsets, $is_root, $is_obsolete, $ontology,  $namespace, $ontology_version ) );

  my @terms;

  while ( $sth->fetch() ) {
    $subsets ||= '';

    push( @terms,
          Bio::EnsEMBL::OntologyTerm->new(
                               '-dbid'        => $dbid,
                               '-adaptor'     => $this,
                               '-accession'   => $accession,
                               '-is_root'     => $is_root,
                               '-is_obsolete' => $is_obsolete,
                               '-ontology'    => $ontology,
                               '-ontology_version' => $ontology_version,
                               '-namespace'   => $namespace,
                               '-subsets'     => [ split( /,/, $subsets ) ],
                               '-name'        => $name,
                               '-definition'  => $definition,
                               '-synonyms'    => $this->_fetch_synonyms_by_dbID($dbid)
    ) );
  }

  return \@terms;
} ## end sub fetch_all

1;
