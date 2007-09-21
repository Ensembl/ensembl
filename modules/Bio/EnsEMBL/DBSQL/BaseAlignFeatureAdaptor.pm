#
# EnsEMBL module for Bio::EnsEMBL::BaseFeatureAdaptor
#
# Copyright (c) 2003 EnsEMBL
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::BaseAlignFeatureAdaptor - Abstract Base class for
AlignFeatureAdaptors

=head1 SYNOPSIS

Abstract class - should not be instantiated.  Implementation of
abstract methods must be performed by subclasses.

=head1 DESCRIPTION

This is a base adaptor for align adaptors.  Since DnaAlignFeatureAdaptor and
PepAlignFeatureAdaptor had almost the same functionality it made sense to 
streamline by creating this superclass.

=head1 CONTACT

Post questions/comments to the ensembl developer list: ensembl-dev@ebi.ac.uk

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
               a lower bound for the hcoverage of feats to obtain
  Arg [3]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Example    : @feats = $adaptor->fetch_all_by_Slice_and_hcoverage($slice, 50.0);
  Description: Returns a listref of features created from the database which
               are on the Slice $slice and with a hcoverage
               greater than $hcoverage.  If logic name is defined, only features
               with an analysis of type $logic_name will be returned.
  Returntype : listref of Bio::EnsEMBL::BaseAlignFeatures in Slice coordinates
  Exceptions : thrown if pid is not defined
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_Slice_and_hcoverage {
  my ($self,$slice,$hcoverage, $logic_name) = @_;
  my $constraint;


  if(defined $hcoverage){
    $constraint = "hcoverage > $hcoverage";
  }

  return $self->fetch_all_by_Slice_constraint($slice, $constraint,
                                              $logic_name);
}


=head2 fetch_all_by_Slice_and_pid

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice from which to obtain align features.
  Arg [2]    : (optional) float $pid 
               a lower bound for the percentage identity of feats to obtain
  Arg [3]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Example    : @feats = $adaptor->fetch_all_by_Slice_and_pid($slice, 50.0);
  Description: Returns a listref of features created from the database which
               are on the Slice $slice and with a percentage identity
               greater than $pid.  If logic name is defined, only features
               with an analysis of type $logic_name will be returned.
  Returntype : listref of Bio::EnsEMBL::BaseAlignFeatures in Slice coordinates
  Exceptions : thrown if pid is not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Slice_and_pid {
  my ($self,$slice,$pid, $logic_name) = @_;
  my $constraint;


#  #get the primary table alias
#  my @tabs = $self->_tables;
#  my $alias = $tabs[0]->[1];

#  if(defined $pid) {
#    $constraint = "${alias}.perc_ident > $pid";
#  }

  if(defined $pid){
    $constraint = "perc_ident > $pid";
  }

  return $self->fetch_all_by_Slice_constraint($slice, $constraint,
                                              $logic_name);
}


=head2 fetch_all_by_hit_name

  Arg [1]    : string $hit_name
               the hit_name of the features to obtain
  Arg [2]    : (optional) string $logic_name
               the analysis logic name of the type of features to obtain
  Example    : @feats = $adaptor->fetch_all_by_hit_name($name, $logic_name);
  Description: Returns a listref of features created from the database
               which correspond to the given hit_name. If logic name
               is defined, only features with an analysis of type
               $logic_name will be returned.
  Returntype : listref of Bio::EnsEMBL::BaseAlignFeatures
  Exceptions : thrown if hit_name is not defined
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_hit_name{
  my( $self, $hit_name, $logic_name ) = @_;
  throw("hit_name argument is required") if(! $hit_name);

  #construct a constraint like 't1.hit_name = "123"'
  my @tabs = $self->_tables;
  my ($name, $syn) = @{$tabs[0]};
  my $constraint = ( "${syn}.hit_name = '$hit_name'" );

  if( $logic_name ){
    # Add the $logic_name constraint
    $constraint = $self->_logic_name_to_constraint($constraint, $logic_name);
  }
  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_hit_name_unversioned

  Arg [1]    : string $hit_name
               the beginning of the hit_name of the features to obtain
               e.g. AA768786 would retrieve AA768786.1, AA768786.2 etc
  Arg [2]    : (optional) string $logic_name
               the analysis logic name of the type of features to obtain
  Example    : @feats = $adaptor->fetch_all_by_hit_name($name, $logic_name);
  Description: Returns a listref of features created from the database
               which start with the  given hit_name. If logic name
               is defined, only features with an analysis of type
               $logic_name will be returned.
  Returntype : listref of Bio::EnsEMBL::BaseAlignFeatures
  Exceptions : thrown if hit_name is not defined
  Caller     : general
  Status     : At risk

=cut

sub fetch_all_by_hit_name_unversioned {
  my( $self, $hit_name, $logic_name ) = @_;
  throw("hit_name argument is required") if(! $hit_name);

  #construct a constraint like 't1.hit_name = "123"'
  my @tabs = $self->_tables;
  my ($name, $syn) = @{$tabs[0]};
  my $constraint = ( "${syn}.hit_name LIKE '$hit_name.%'" );

  if( $logic_name ){
    # Add the $logic_name constraint
    $constraint = $self->_logic_name_to_constraint($constraint, $logic_name);
  }
  return $self->generic_fetch($constraint);
}



=head2 fetch_all_by_RawContig_and_pid

  Description: DEPRECATED use fetch_all_by_Slice_and_pid instead

=cut

sub fetch_all_by_RawContig_and_pid {
  my($self, $contig, $pid, $logic_name) = @_;

  my $constraint;

  #get the primary table alias
  my @tabs = $self->_tables;
  my $alias = $tabs[0]->[1];

  if(defined $pid) {
    $constraint = "${alias}.perc_ident > $pid";
  }

  return $self->fetch_all_by_RawContig_constraint($contig, 
						  $constraint, 
						  $logic_name);
}




##implemented by subclasses:
# store
# _tables
# _columns
# _obj_from_hashref



1;
