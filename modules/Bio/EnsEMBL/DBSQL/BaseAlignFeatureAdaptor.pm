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

Abstract class should not be instantiated.  Implementation of
abstract methods must be performed by subclasses.

=head1 DESCRIPTION

This is a base adaptor for align adaptors.  Since DnaAlignFeatureAdaptor and
PepAlignFeatureAdaptor had almost the same functionality it made sense to 
streamline by creating this superclass.

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::BaseAlignFeatureAdaptor;
use vars qw(@ISA);
use strict;

#Object preamble - inherits from Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);

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


