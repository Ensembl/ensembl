# Interface which ContextMapper will implement
# 
# Author: Arne Stabenau
#
# Copyright EMBL/EBI 2001
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL:ContextMapperI

=head1 SYNOPSIS

Just interface description

=head1 DESCRIPTION

A ContextMapper takes objects which implement Bio::EnsEMBL::RangeI and
translates them to another context. Most typical example is RawContig
to VirtualContig (Chromosome) context.

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk


=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::ContextMapperI;

=head2 change_context

 Title   : change_context
 Usage   : $contextMapper->change_context( $rangeObject, $newContext )
 Function: changes the rangeObject to the given new context. ContextMapper 
           are usually specialised in one context. 
 Example : 
 Returns : throws various exceptions. Changes happen in place.
  Args    : $rangeObject is an object with a Bio::EnsEMBL::Range
            $newContext is a context with name and type.

=cut

sub change_context {
  my ($self, $range, $context ) = @_;

  if( ! $range->isa( "Bio::EnsEMBL::RangeI" )) {
    die( "Your object is no RangeI" );
  }

  die( "Your context mapper does not have change_context" );
}

1;
