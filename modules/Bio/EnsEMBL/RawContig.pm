# EnsEMBL Contig object
#
# Copyright EMBL-EBI 2001
#
# cared for by:: Arne Stabenau
# Date : 04.12.2001
#

=head1 NAME

Bio::EnsEMBL::RawContig This class is DEPRECATED, use Bio::EnsEMBL::Slice
instead

=cut


package Bio::EnsEMBL::RawContig;

use vars qw( @ISA );
use strict;

use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Clone;
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);

@ISA = qw( Bio::EnsEMBL::Slice );

=head2 new

Description: DEPRECATED use Bio::EnsEMBL::Slice instead

=cut

sub new {
  my ( $class, @args ) = @_;

  deprecate("Bio::EnsEMBL::RawContig is a deprecated class.\n" .
            "Use Bio::EnsEMBL::Slice instead");

  return $class->SUPER::new(@args);
}




sub embl_offset {
  my $self = shift;

  deprecate('Use Bio::EnsEMBL::Slice::project instead');

  my @projection = @{$self->project('clone')};

  if(@projection != 1) {
    warning('This contig not associated with a single clone');
    return 1;
  }

  #return the start position into the clone
  return $projection[0]->[2]->start();
}




sub clone {
  my $self = shift;

  deprecate('Use Bio::EnsEMBL::Slice::project instead');

  my @projection = @{$self->project('clone')};

  if(@projection != 1) {
    warning('This contig not associated with a single clone');
    return undef;
  }

  my $clone = $projection[0]->[2];

  #get full clone instead of potentially partial clone
  $clone = $self->adaptor->fetch_by_region($clone->coord_system->name(),
                                           $clone->seq_region_name(),
                                           undef,
                                           undef,
                                           undef,
                                           $clone->coord_system->version());


  #rebless the slice as a Bio::EnsEMBL::Clone so old method calls still work
  return bless $clone, 'Bio::EnsEMBL::Clone';
}




=head2:  ctg2genomic

  Description: use Bio::EnsEMBL::Slice instead

=cut

sub ctg2genomic{
  # Map the internal ID onto the golden path
  my $self   = shift;
	my $start  = shift || 1;
	my $end    = shift || $self->length;
	my $strand = shift || 1;

  deprecate('Use Bio::EnsEMBL::Slice instead');

  my $db = $self->adaptor->db();
  my $aa = $db->get_AssemblyMapperAdaptor();
  my $csa = $db->get_CoordSystemAdaptor();

  my $top_level = $csa->fetch_top_level();
  my $ma = $aa->fetch_by_CoordSystems($top_level, $self->coord_system);

  return $ma->map_coordinates_to_assembly( $self->seq_region_name,
                                           $start,
                                           $end,
                                           $strand,
                                           $self->coord_system);
}



#
# The name that is wanted is actually the seq_region name
#
sub name {
  my $self = shift;

  return $self->seq_region_name();
}



1;
