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

  my $projection = $self->project('seqlevel');
  foreach my $segment (@$projection) {
    my($start,$end,$contig) = @$segment;
    if($start <= $pos) {
      return bless($contig, 'Bio::EnsEMBL::RawContig');
    }
  }

  return undef;
}



=head2 htg_phase

  Description: DEPRECATED - use $slice->get_attribute('htg_phase') instead

=cut

sub htg_phase {
   my $self = shift;
   my ($htg_phase) = $self->get_attribute('htg_phase');
   return $htg_phase;
}



=head2 created

  Description: DEPRECATED - Created information no longer stored

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

  Description: DEPRECATE - Modified information no longer stored

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

  Description: DEPRECATED - use Slice::seq_region_name

=cut

sub version{ embl_version(@_) }



=head2 embl_version

  Description: DEPRECATED - use Slice::seq_region_name

=cut

sub embl_version {
  my $self = shift;
  my $acc_ver = $self->seq_region_name();

  #strip version off end of accession
  my $ver;
  (undef, $ver) = split(/\./, $acc_ver);
  return $ver;
}



=head2 embl_id

  description: DEPRECATED - use Slice::seq_region_name

=cut

sub embl_id {
  my $self = shift;
  my $acc = $self->seq_region_name();

  #strip off version
  ($acc) = split(/\./, $acc);
  return $acc;
}


#what is actually meant by clone->name is seq_region_name not name
sub name {
  my $self = shift;
  return $self->seq_region_name();
}


=head2 id

Description: DEPRECATED - use Slice::seq_region_name

=cut

sub id { embl_id(@_);}



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






