package Bio::EnsEMBL::PepDnaAlignFeature;

# EnsEMBL module for storing dna-dna pairwise alignments
#
# Cared for by Michele Clamp <michele@sanger.ac.uk>
#
# You may distribute this module under the same terms as perl itself
#

=head1 NAME

  Bio::EnsEMBL::PepDnaAlignFeature - Ensembl specific pep-dna pairwise alignment feature

=head1 SYNOPSIS

  See BaseAlignFeature

=cut 

use Bio::EnsEMBL::BaseAlignFeature;

use vars qw(@ISA);
use strict;

@ISA = qw( Bio::EnsEMBL::BaseAlignFeature );


=head2 new_fast

  Arg [1]    : hashref $hashref
               A hashref which will be blessed into a PepDnaAlignFeature. 
  Example    : none
  Description: This allows for very fast object creation when a large number 
               of PepDnaAlignFeatures needs to be created.  This is a bit of 
               a hack but necessary when thousands of features need to be
               generated within a couple of seconds for web display. It is
               not recommended that this method be called unless you know what
               you are doing.  It requires knowledge of the internals of this
               class and its superclasses.  
  Returntype : Bio::EnsEMBL::PepDnaAlignFeature
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor

=cut

sub new_fast {
  my ($class, $hashref) = @_;

  return bless $hashref, $class;
}


sub transform {
  my $self = shift;

  $self->throw( "PepDnaAlignFeatures cant be transformed as".
		" they are not on EnsEMBL coord system" );
}

sub _hit_unit {
  return 3;
}

sub _query_unit {
  return 1;
}

1;
