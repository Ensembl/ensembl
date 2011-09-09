=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

Bio::EnsEMBL::DnaPepAlignFeature - Ensembl specific dna-pep pairwise
alignment feature

=head1 SYNOPSIS

  See BaseAlignFeature

=cut 


package Bio::EnsEMBL::DnaPepAlignFeature;

use strict;

use Bio::EnsEMBL::BaseAlignFeature;
use Scalar::Util qw(weaken isweak);

use vars qw(@ISA);

@ISA = qw( Bio::EnsEMBL::BaseAlignFeature );


=head2 new_fast

  Arg [1]    : hashref $hashref
               A hashref which will be blessed into a PepDnaAlignFeature. 
  Example    : none
  Description: This allows for very fast object creation when a large number 
               of DnaPepAlignFeatures needs to be created.  This is a bit of 
               a hack but necessary when thousands of features need to be
               generated within a couple of seconds for web display. It is
               not recommended that this method be called unless you know what
               you are doing.  It requires knowledge of the internals of this
               class and its superclasses.  
  Returntype : Bio::EnsEMBL::DnaPepAlignFeature
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor
  Status     : Stable

=cut

sub new_fast {
  my ($class, $hashref) = @_;
  my $self = bless $hashref, $class;
  weaken($self->{adaptor})  if ( ! isweak($self->{adaptor}) );
  return $self;
}


=head2 _hit_unit

  Arg [1]    : none
  Example    : none
  Description: PRIVATE implementation of abstract superclass method.  Returns
               1 as the 'unit' used for the hit sequence. 
  Returntype : int
  Exceptions : none
  Caller     : Bio::EnsEMBL::BaseAlignFeature
  Status     : Stable


=cut

sub _hit_unit {
  return 1;
}


=head2 _query_unit

  Arg [1]    : none
  Example    : none
  Description: PRIVATE implementation of abstract superclass method.  Returns
               3 as the 'unit' used for the query sequence. 
  Returntype : int
  Exceptions : none
  Caller     : Bio::EnsEMBL::BaseAlignFeature
  Status     : Stable


=cut

sub _query_unit {
  return 3;
}




1;
