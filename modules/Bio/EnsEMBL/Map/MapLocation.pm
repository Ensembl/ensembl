=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
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

Bio::EnsEMBL::Map::MapLocation

=head1 SYNOPSIS

=head1 DESCRIPTION

Represents a location on a genetic map, yac map, radition hybrid map,
etc.

=head1 METHODS

=cut


package Bio::EnsEMBL::Map::MapLocation;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Utils::Exception qw(deprecate);

=head2 new

  Arg [1]    : (optional) string $name
  Arg [2]    : (optional) string $map_name
  Arg [3]    : (optional) string $chromosome_name
  Arg [4]    : (optional) string $position
  Arg [5]    : (optional) float $lod_score
  Example    : $map_location = Bio::EnsEMBL::Map::MapLocation('DS1234',
                                                              'genethon',
                                                              'X',
                                                              '12.39',
                                                              50.12);
  Description: Creates a new MapLocation
  Returntype : Bio::EnsEMBL::Map::MapLocation
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub new {
  my ($caller, $name, $map_name, $chromosome_name, $position, $lod_score) = @_;

  my $class = ref($caller) || $caller;

  return bless( {'map_name'        => $map_name,
                 'name'            => $name,
                 'chromosome_name' => $chromosome_name,
                 'position'        => $position,
                 'lod_score'       => $lod_score}, $class );
}



=head2 map_name

  Arg [1]    : string $map_name
  Example    : $map_name = $map_location->map_name;
  Description: Getter/Setter for the map name
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub map_name {
  my $self = shift;
  $self->{'map_name'} = shift if(@_);
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
  Status     : stable

=cut

sub name {
  my $self = shift;
  $self->{'name'} = shift if(@_);
  return $self->{'name'};
}


=head2 chromosome_name

  Arg [1]    : (optional) string $chromosome_name
  Example    : $chr_name = $map_location->chromosome_name;
               $map_location->chromosome_name('X');
  Description: The name of the chromosome associated with this map location
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub chromosome_name{
  my $self = shift;
  $self->{'chromosome_name'} = shift if(@_);
  return $self->{'chromosome_name'};
}



=head2 position

  Arg [1]    : (optional) string $position
  Example    : $pos = $map_location->position;
  Description: Getter/Setter for the position of this map location
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub position {
  my $self = shift;
  $self->{'position'} = shift if(@_);
  return $self->{'position'};
}



=head2 lod_score

  Arg [1]    : (optional) float $lod
  Example    : $lod = $map_location->lod_score;
  Description: Getter/Setter for lod score of this map location
  Returntype : float
  Exceptions : none
  Caller     : general
  Status     : stable

=cut

sub lod_score {
  my $self = shift;
  $self->{'lod_score'} = shift if(@_);
  return $self->{'lod_score'};
}



=head2 chromosome

  Description: DEPRECATED use chromosome_name() instead

=cut

sub chromosome {
  my $self = shift;
  deprecate('use chromosome_name instead');

  if(@_) {
    my $chr = shift;
    if(ref($chr)) {
      $self->chromosome_name($chr->seq_region_name());
    } else {
      $self->chromosome_name($chr);
    }
  }

  #this object has no way to talk to db and thus no way to
  #get a chromosome object
  return $self->chromosome_name();
}




1;
