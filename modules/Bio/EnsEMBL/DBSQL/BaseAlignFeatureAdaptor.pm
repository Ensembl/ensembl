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

Bio::EnsEMBL::DBSQL::BaseAlignFeatureAdaptor - Abstract Base class for
AlignFeatureAdaptors

=head1 SYNOPSIS

Abstract class, should not be instantiated.  Implementation of abstract
methods must be performed by subclasses.

=head1 DESCRIPTION

This is a base adaptor for the align feature adaptors
DnaAlignFeatureAdaptor and ProteinAlignFeatureAdaptor.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::BaseAlignFeatureAdaptor;
use vars qw(@ISA @EXPORT);
use strict;

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);

@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});

=head2 fetch_all_by_Slice_and_hcoverage

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice from which to obtain align features.
  Arg [2]    : (optional) float $hcoverage 
               A lower bound for the hcoverage of feats to obtain.
  Arg [3]    : (optional) string $logic_name
               The logic name of the type of features to obtain.
  Example    : @feats = @{
                $adaptor->fetch_all_by_Slice_and_hcoverage( $slice,
                  50.0 ) };
  Description: Returns a listref of features created from the
               database which are on the Slice $slice and with a
               hcoverage greater than $hcoverage.  If logic name
               is defined, only features with an analysis of type
               $logic_name will be returned.
  Returntype : listref of Bio::EnsEMBL::BaseAlignFeatures
               in Slice coordinates
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_Slice_and_hcoverage {
  my ( $self, $slice, $hcoverage, $logic_name ) = @_;

  my $constraint;
  if ( defined($hcoverage) ) {
    $constraint = "hcoverage > $hcoverage";
  }

  return
    $self->fetch_all_by_Slice_constraint( $slice, $constraint,
                                          $logic_name );
}

=head2 fetch_all_by_Slice_and_external_db

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice from which to obtain align features.
  Arg [2]    : String $external_db_name
               Name of the external DB to which the align features
               should be restricted.
  Arg [3]    : (optional) string $logic_name
               The logic name of the type of features to obtain.
  Example    : @feats = @{
                  $adaptor->fetch_all_by_Slice_and_external_db( $slice,
                    'EMBL' ) };
  Description: Returns a listref of features created from the
               database which are on the Slice $slice and associated
               with external DB $external_db_name.  If logic name
               is defined, only features with an analysis of type
               $logic_name will be returned.
  Returntype : listref of Bio::EnsEMBL::BaseAlignFeatures
               in Slice coordinates
  Exceptions : thrown if $external_db_name is not defined or if
               the subclass does not return a table alias for the
               external_db table from _tables()
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_Slice_and_external_db {
  my ( $self, $slice, $external_db_name, $logic_name ) = @_;

  if ( !defined($external_db_name) ) {
    throw("Need name of external DB to restrict to");
  }

  my @join_tables = $self->_tables();

  my $edb_alias;
  foreach my $join_table (@join_tables) {
    my ( $table, $table_alias ) = @{$join_table};
    if ( $table eq 'external_db' ) {
      $edb_alias = $table_alias;
      last;
    }
  }

  if ( !defined($edb_alias) ) {
    throw("Can not find alias for external_db table");
  }

  my $constraint = sprintf( "%s.db_name = %s",
                            $edb_alias,
                            $self->dbc()->db_handle()
                              ->quote( $external_db_name, SQL_VARCHAR )
  );

  return
    $self->fetch_all_by_Slice_constraint( $slice, $constraint,
                                          $logic_name );
} ## end sub fetch_all_by_Slice_and_external_db

=head2 fetch_all_by_Slice_and_pid

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice from which to obtain align features.
  Arg [2]    : (optional) float $pid 
               A lower bound for the percentage identity of features
               to obtain.
  Arg [3]    : (optional) string $logic_name
               The logic name of the type of features to obtain.
  Example    : @feats =
                 @{ $adaptor->fetch_all_by_Slice_and_pid( $slice, 50.0 ) };
  Description: Returns a listref of features created from the
               database which are on the Slice $slice and with a
               percentage identity greater than $pid.  If logic name
               is defined, only features with an analysis of type
               $logic_name will be returned.
  Returntype : listref of Bio::EnsEMBL::BaseAlignFeatures
               in Slice coordinates
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Slice_and_pid {
  my ( $self, $slice, $pid, $logic_name ) = @_;

  #  #get the primary table alias
  #  my @tabs = $self->_tables;
  #  my $alias = $tabs[0]->[1];

  #  if(defined $pid) {
  #    $constraint = "${alias}.perc_ident > $pid";
  #  }

  my $constraint;
  if ( defined($pid) ) {
    $constraint = sprintf( "perc_ident > %s",
                           $self->dbc()->db_handle()
                             ->quote( $pid, SQL_FLOAT ) );
  }

  return
    $self->fetch_all_by_Slice_constraint( $slice, $constraint,
                                          $logic_name );
}


=head2 fetch_all_by_hit_name

  Arg [1]    : string $hit_name
               The hit_name of the features to obtain
  Arg [2]    : (optional) string $logic_name
               The analysis logic name of the type of features to
               obtain.
  Example    : @feats =
                 @{ $adaptor->fetch_all_by_hit_name( 'AK078491.1',
                   'vertrna' ); }
  Description: Returns a listref of features created from the
               database which correspond to the given hit_name.  If
               logic name is defined, only features with an analysis
               of type $logic_name will be returned.
  Returntype : listref of Bio::EnsEMBL::BaseAlignFeatures
  Exceptions : thrown if hit_name is not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_hit_name {
  my ( $self, $hit_name, $logic_name ) = @_;

  if ( !defined($hit_name) ) {
    throw("hit_name argument is required");
  }

  # Construct a constraint like 't1.hit_name = "123"'
  my @tabs = $self->_tables();
  my ( $name, $syn ) = @{ $tabs[0] };

  my $constraint = sprintf( "%s.hit_name = %s",
           $syn,
           $self->dbc()->db_handle()->quote( $hit_name, SQL_VARCHAR ) );

  if ( defined($logic_name) ) {
    # Add the $logic_name constraint
    $constraint =
      $self->_logic_name_to_constraint( $constraint, $logic_name );
  }

  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_hit_name_unversioned

  Arg [1]    : string $hit_name
               The beginning of the hit_name of the features to
               obtain, e.g. AA768786 would retrieve AA768786.1,
               AA768786.2 etc.
  Arg [2]    : (optional) string $logic_name
               The analysis logic name of the type of features to
               obtain.
  Example    : @feats =
                  @{ $adaptor->fetch_all_by_hit_name( $name,
                    $logic_name ) };
  Description: Returns a listref of features created from the
               database which start with the given hit_name.  If
               logic name is defined, only features with an analysis
               of type $logic_name will be returned.
  Returntype : listref of Bio::EnsEMBL::BaseAlignFeatures
  Exceptions : thrown if hit_name is not defined
  Caller     : general
  Status     : At risk

=cut

sub fetch_all_by_hit_name_unversioned {
  my ( $self, $hit_name, $logic_name ) = @_;

  if ( !defined($hit_name) ) {
    throw("hit_name argument is required");
  }
  $hit_name =~ s/_/\\_/;

  #construct a constraint like 't1.hit_name = "123"'
  my @tabs = $self->_tables;
  my ( $name, $syn ) = @{ $tabs[0] };

  my $constraint = sprintf( "%s.hit_name LIKE %s",
    $syn,
    $self->dbc()->db_handle()->quote( $hit_name . '.%', SQL_VARCHAR ) );

  if ( defined($logic_name) ) {
    # Add the $logic_name constraint
    $constraint =
      $self->_logic_name_to_constraint( $constraint, $logic_name );
  }

  return $self->generic_fetch($constraint);
}


##implemented by subclasses:
# store
# _tables
# _columns
# _obj_from_hashref



1;
