#
# BioPerl module for Bio::EnsEMBL::DBSQL::BaseAlignFeatureAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
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

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal 
methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::BaseAlignFeatureAdaptor;
use vars qw(@ISA);
use strict;

#Object preamble - inherits from Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);



=head2 fetch_by_contig_id_and_pid

  Arg [1]    : int $cid
               the unique identifier (dbID) for contig to obtain feats from
  Arg [2]    : float $pid 
               a lower bound for the percentage identifier of feats to obtain
  Arg [3]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Example    : @feats = $adaptor->fetch_by_contig_id_and_pid(1234, 50.0);
  Description: Returns a list of features created from the database which are 
               are on the contig defined by $cid and with a percentage id 
               greater than $pid.  If logic name is defined, only features
               with an analysis of type $logic_name will be returned. 
  Returntype : list of Bio::EnsEMBL::*AlignFeature in contig coordinates
  Exceptions : thrown if $pid is not defined
  Caller     : general

=cut

sub fetch_by_contig_id_and_pid {
  my($self, $cid, $pid, $logic_name) = @_;

  my $constraint;

  if(!defined $pid) {
    $self->throw("need a pid even if its 0\n");
  } else {
    $constraint = "perc_ident > $pid";
  }

  my @features = 
    $self->fetch_by_contig_id_constraint($cid, $constraint, $logic_name);
  
  return @features;
}


=head2 fetch_by_Slice_and_pid

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice from which to obtain align features.
  Arg [2]    : float $pid 
               a lower bound for the percentage identifier of feats to obtain
  Arg [3]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Example    : @feats = $adaptor->fetch_by_Slice_and_pid($slice, 50.0);
  Description: Returns a list of features created from the database which are 
               are on the Slice $slice and with a percentage id 
               greater than $pid.  If logic name is defined, only features
               with an analysis of type $logic_name will be returned. 
  Returntype : list of Bio::EnsEMBL::*AlignFeature in Slice coordinates
  Exceptions : thrown if pid is not defined
  Caller     : general

=cut

sub fetch_by_Slice_and_pid {
  my ($self,$slice,$pid, $logic_name) = @_;
  my $constraint;

  if(!defined $pid){
    $self->throw("need a pid even if its 0\n");
  }else{
    $constraint = "perc_ident > $pid";
  }

  return $self->fetch_by_Slice_constraint($slice, $constraint, $logic_name);
}  


=head2 fetch_by_assembly_location_and_pid

  Arg [1]    : int $start
               the start of the assembly region from which to obtain align
               features (in chromosomal coords).
  Arg [2]    : int $end
               the end of the assembly region from which to obtain align 
               features (in chromosomal coords).
  Arg [3]    : string $chr
               the name of the chromosome from which to obtain align features.
  Arg [4]    : string $type
               the type of assembly to be used
  Arg [4]    : float $pid
               a lower bound for the percentage identifier of feats to obtain
  Arg [5]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Example    : @f = $a->fetch_by_assembly_location_and_pid(1,10000,'9','NCBI30',50.0);
  Description: Returns a list of features created from the database which are 
               are in the assembly region defined by $start, $end, and $chr, 
               and with a percentage identity greater than $pid.  If 
               $logic_name is defined, only features with an analysis of type 
               $logic_name will be returned. 
  Returntype : list of Bio::EnsEMBL::*AlignFeature in chromosomal coordinates
  Exceptions : thrown if $pid is not defined
  Caller     : general

=cut

sub fetch_by_assembly_location_and_pid{
  my ($self,$start,$end,$chr,$type, $pid, $logic_name) = @_;
  my $constraint;

  if(!defined $pid){
    $self->throw("need a pid even if its 0\n");
  }else{
    $constraint = "perc_ident > $pid";
  }

  return $self->fetch_by_assembly_location_constraint($start, $end, $chr,$type,
						     $constraint, $logic_name);

}

##Abstract methods inherited from BaseFeatureAdaptor must still be
##implemented by subclasses:
# store
# _tablename
# _columns
# _obj_from_hashref


1;


