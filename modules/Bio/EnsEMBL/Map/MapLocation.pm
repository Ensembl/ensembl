# EnsEMBL module for MapLocation
# Copyright EMBL-EBI/Sanger center 2003
#
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Map::MapLocation

=head1 SYNOPSIS


=head1 DESCRIPTION

Represents a location on a genetic map, yac map, radition hybrid map, etc.

=cut


package Bio::EnsEMBL::Map::MapLocation;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Root);


=head2 new

  Arg [1]    : (optional) string $name
  Arg [2]    : (optional) string $map_name
  Arg [3]    : (optional) Bio::EnsEMBL::Chromosome $chromosome
  Arg [4]    : (optional) string $position
  Arg [5]    : (optional) float $lod_score
  Example    : $map_location = Bio::EnsEMBL::Map::MapLocation('DS1234',
							      'genethon',
							       $chr,
							      '12.39',
							      50.12);
  Description: Creates a new MapLocation
  Returntype : Bio::EnsEMBL::Map::MapLocation
  Exceptions : none
  Caller     : general

=cut

sub new {
  my ($caller, $name, $map_name, $chromosome, $position, $lod_score) = @_;

  my $class = ref($caller) || $caller;

  return bless( {'map_name'   => $map_name,
                 'name'       => $name, 
		 'chromosome' => $chromosome, 
		 'position'   => $position, 
		 'lod_score'  => $lod_score}, $class );
}



=head2 map_name

  Arg [1]    : string $map_name
  Example    : $map_name = $map_location->map_name;
  Description: Retrieves the map_name
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub map_name {
  my $self = shift;
  
  if(@_) {
    $self->{'map_name'} = shift;
  }
  
  return $self->{'map_name'};
}



=head2 name

  Arg [1]    : (optional) string $name
  Example    : $name = $map_location->name;
  Description: A name associated with the marker at this position.  For
               example if this is a genethon map location the name will be
               the synonym of source genethon. 
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub name {
  my $self = shift;

  if(@_) {
    $self->{'name'} = shift;
  }

  return $self->{'name'};  
}


=head2 chromosome

  Arg [1]    : (optional) Bio::EnsEMBL::Chromosome
  Example    : $chr = $map_location->chromosome;
  Description: Getter/Setter for the chromosome of this map location
  Returntype : Bio::EnsEMBL::Chromosome
  Exceptions : thrown if an invalid arg is passed
  Caller     : general

=cut

sub chromosome {
  my $self = shift;

  if(@_) {
    my $chr = shift;
    if($chr && !(ref($chr) && $chr->isa('Bio::EnsEMBL::Chromosome'))) {
      $self->throw("arg must be a Bio::EnsEMBL::Chromosome not [$chr]");
    }

    $self->{'chromosome'} = $chr;
  }

  return $self->{'chromosome'};
}



=head2 position

  Arg [1]    : (optional) string $position
  Example    : $pos = $map_location->position;
  Description: Getter/Setter for the position of this map location
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub position {
  my $self = shift;

  if(@_) {
    $self->{'position'} = shift;
  }
  
  return $self->{'position'};
}



=head2 lod_score

  Arg [1]    : (optional) float $lod
  Example    : $lod = $map_location->lod_score;
  Description: Getter/Setter for lod score of this map location
  Returntype : float
  Exceptions : none
  Caller     : general

=cut

sub lod_score {
  my $self = shift;

  if(@_) {
    $self->{'lod_score'} = shift;
  }

  return $self->{'lod_score'};
}


1;
