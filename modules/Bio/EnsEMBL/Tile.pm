# Ensembl module for Bio::EnsEMBL::Tile
#
# Cared for by Arne Stabenau <ensembl-dev@ebi.ac.uk>
#
# Copyright EnsEMBL 2002
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Tile - DEPRECATED - You should not have to use this class.
                     The Bio::EnsEMBL::Slice::project method should be able
                     to construct any tiling path that you need.


=head1 AUTHOR - Arne Stabenau

=head1 CONTACT

This modules is part of the Ensembl project http://www.ensembl.org

Questions can be posted to the ensembl-dev mailing list:
ensembl-dev@ebi.ac.uk

=cut


package Bio::EnsEMBL::Tile;
use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Utils::Exception qw(deprecate);

@ISA = qw(Bio::EnsEMBL::Root);



### Only new_fast() is used by the EnsEMBL API
### Otter uses new()
sub new {
    my $class = shift;

    deprecate("Bio::EnsEMBL::Tile is a deprecated class.\nIf you want a " .
              'tiling path it is best to use Bio::EnsEMBL::Slice::project');

    my $self = bless {}, $class;
    if (@_) {
        $self->throw("new does not take any arguments");
    }
    return $self;
}

=head2 new_fast

  Arg [1]    : Bio::EnsEMBL::Slice $assembled_seq
               The assembled sequence which assembly should be described in this Tile
  Arg [2]    : int $assembled_start
  Arg [3]    : int $assembled_end
               Where on the sequence the tile starts and ends
  Arg [4]    : Bio::EnsEMBL::RawContig $component_seq
               The sequence that is used to build up the assembled sequence.
  Arg [5]    : int $component_start
  Arg [6]    : int $component_end
               Which component piece makes up the assembled piece
  Arg [7]    : int $orientation (+1 or -1)
  Example    : none
  Description: Fast constructor for Tile object
  Returntype : Bio::EnsEMBL::Tile
  Exceptions : none
  Caller     : Slice::get_tiling_path()

=cut

sub new_fast {
  deprecate("Bio::EnsEMBL::Tile is a deprecated class.\nIf you want a " .
            'tiling path it is best to use Bio::EnsEMBL::Slice::project');

  # make it really fast
  return bless {
     'assembled_Seq'    => $_[1],
     'assembled_start'  => $_[2],
     'assembled_end'    => $_[3],
     'component_Seq'    => $_[4],
     'component_start'  => $_[5],
     'component_end'    => $_[6],
     'component_ori'    => $_[7],
# compatibility zone
     'start'            => $_[2],
     'end'              => $_[3],
     'strand'           => $_[7],
     'contig'           => $_[4]
  }, $_[0];
}


sub assembled_Seq {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{'assembled_Seq'} = $arg;
  }

  return $self->{'assembled_Seq'};
}

sub assembled_start {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{'assembled_start'} = $arg;
    $self->{'start'} = $arg;
  }

  return $self->{'assembled_start'};
}

sub assembled_end {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{'assembled_end'} = $arg;
    $self->{'end'} = $arg;
  }

  return $self->{'assembled_end'};
}

sub component_Seq {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{'component_Seq'} = $arg;
    $self->{'contig'} = $arg;
  }

  return $self->{'component_Seq'};
}

sub component_start {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{'component_start'} = $arg;
  }

  return $self->{'component_start'};
}

sub component_end {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{'component_end'} = $arg;
  }

  return $self->{'component_end'};
}

sub component_ori {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{'component_ori'} = $arg;
    $self->{'strand'} = $arg;
  }

  return $self->{'component_ori'};
}


1;
