#
# EnsEMBL module for Bio::EnsEMBL::Clone
#
# Cared for by EnsEMBL (www.ensembl.org)
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Clone - DEPRECATED use Bio::EnsEMBL::Slice instead

=head1 DESCRIPTION

DEPRECATED - Create Slices on clone regions instead

=head1 CONTACT

Post questions to the EnsEMBL developer list: <ensembl-dev@ebi.ac.uk> 

=cut

package Bio::EnsEMBL::Clone;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::RawContig;

use Bio::EnsEMBL::Utils::Exception qw(throw deprecate);

@ISA = qw( Bio::EnsEMBL::Slice );


=head2 new

  Description: DEPRECATED.  Use Bio::EnsEMBL::Slice instead

=cut

sub new {
  my $class = shift;

  deprecate("Bio::EnsEMBL::Clone is a deprecated class\n" .
            "Use Bio::EnsEMBL::Slice instead");

  return $class->SUPER::new(@_);
}




=head2 get_all_Contigs

  Description: DEPRECATED. Use Slice::project instead

=cut


sub get_all_Contigs {
  my( $self ) = @_;

  deprecate('Use Slice::project instead');

  # Assume that we actually want to project to sequence level
  my $projection = $self->project('seqlevel');

  my @out;
  foreach my $segment (@$projection) {
    my $contig = $segment->[2];
    #bless slices into RawContigs for backwards compatibility
    bless $contig, "Bio::EnsEMBL::RawContig";
    push @out, $contig;
  }

  return \@out;
}



=head2 add_Contig

  Description: DEPRECATED.  There is currently no replacement for this method
               Possibly one should be added?

=cut

sub add_Contig {
  my ($self, $contig) = @_;

  ### XXX Should there be a replacement for this?
  deprecate('There is currently no replacement for this method');

  return 0;
}



=head2 delete_by_dbID

  Description: DEPRECATED. There is currently no replacement for this method
               Possibly there should be one?

=cut

sub delete_by_dbID {
  my ($self)=shift;

  ### XXX Should there be a replacement for this?
  deprecate('There is currently no replacement for this method');

  return 0;
}



=head2 get_RawContig_by_position

  Description: DEPRECATED. Use Slice::project instead

=cut

sub get_RawContig_by_position {
  my ($self, $pos) = @_;

  deprecate('Use Slice::project instead');

  throw("get_rawcontig_by_position error: Position must be > 0") if($pos < 1);

  my $projection = $self->project();
  foreach my $segment (@$projection) {
    my($start,$end,$contig) = @$segment;
    if($start <= $pos) {
      return bless($contig, 'Bio::EnsEMBL::RawContig');
    }
  }

  return undef;
}



=head2 htg_phase

  Arg [1]    : string $htg_phase
               0,1,2,3 representing how finished the clone is
  Example    : none
  Description: get/set for attribute htg_phase
               ( high throughput genome project phase ) 
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub htg_phase {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'htg_phase'} = $value;
    }
    return $obj->{'htg_phase'};
}



=head2 created

  Arg [1]    : string $created
  Example    : none
  Description: get/set for attribute created.
               Gives the unix time value of the created 
               datetime field, which indicates
               the first time this clone was put in ensembl
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub created {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'created'} = $value;
    }
    return $obj->{'created'};
}



=head2 modified

  Arg [1]    : string $modified
  Example    : none
  Description: get/set for attribute modified
               Gives the unix time value of the modified 
               datetime field, which indicates
               the last time this clone was modified in ensembl
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub modified {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'modified'} = $value;
    }
    return $obj->{'modified'};
}



=head2 version

  Arg [1]    : string $version
  Example    : none
  Description: get/set for attribute version
               this could contain an ensembl version for the clone.
               Usually we just use the EMBL one though. EnsEMBL version
               are currently not generated or maintained for clones.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub version{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'version'} = $value;
    }
    return $obj->{'version'};

}



=head2 embl_version

  Arg [1]    : string $embl_version
  Example    : none
  Description: get/set for attribute embl_version
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub embl_version {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'embl_version'} = $value;
    }
    return $obj->{'embl_version'};
}



=head2 embl_id

  Arg [1]    : string $embl_id
  Example    : none
  Description: get/set for attribute embl_id
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub embl_id {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'embl_id'} = $value;
    }
    return $obj->{'embl_id'};

}



=head2 id

  Args       : string $value (optional)
               The new name of this clone
  Example    : $name = $clone->id();
  Description: The name of hte clone.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub id {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_clone_id'} = $value;
    }
    return $obj->{'_clone_id'};

}



=head2 dbID

  Description: Deprecated. Use SliceAdaptor::get_seq_region_id instead

  Exceptions : none
  Caller     : general, set from adaptor on store

=cut

sub dbID{
  my $self = shift;
  deprecate('Use Bio::EnsEMBL::Slice instead of Bio::EnsEMBL::RawContig');

  return $self->adaptor->get_seq_region_id($self);
}





1;






