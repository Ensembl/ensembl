
#
# Ensembl module for Bio::EnsEMBL::SimpleFeature
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::SimpleFeature - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::SimpleFeature;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root

use Bio::EnsEMBL::SeqFeature;

@ISA = qw(Bio::EnsEMBL::SeqFeature);

# new() is inherieted from SeqFeature


=head2 display_label

  Title   : display_label
  Usage   : $obj->display_label($newval)
  Function: 
  Example : 
  Returns : value of display_label
  Args    : newvalue (optional)

=cut

sub display_label{
   my ($self, $value) = @_;

   if (defined $value) {
      $self->{'display_label'} = $value;
   }
   return $self->{'display_label'};
}


sub dbID{
  my ($self, $arg) = @_;

  if(defined $arg) {
    $self->{_database_id} = $arg;
  }
  return $self->{_database_id}; 
}


=head2 adaptor

  Arg [1]    : Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor $adaptor
  Example    : none
  Description: get/set for this objects Adaptor
  Returntype : Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor
  Exceptions : none
  Caller     : general, set from adaptor on store

=cut

sub adaptor {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'adaptor'} = $value;
   }
   return $self->{'adaptor'};

}


1;
