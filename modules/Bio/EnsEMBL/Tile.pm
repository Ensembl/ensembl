# Ensembl module for Bio::EnsEMBL::Tile
#
# Cared for by Arne Stabenau <ensembl-dev@ebi.ac.uk>
#
# Copyright EnsEMBL 2002
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Tile - container objects for assembly style information

=head1 SYNOPSIS

   @tile_list = $slice->get_tiling_path()
       

=head1 DESCRIPTION
 
  Each tile contains information about an ungapped sequence mapping from one sequence containing object to another.

=head1 AUTHOR - Arne Stabenau

=head1 CONTACT

This modules is part of the Ensembl project http://www.ensembl.org

Questions can be posted to the ensembl-dev mailing list:
ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Tile;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root

use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Root);

sub new {
  my($class,@args) = @_;

  my $self = {};
  bless $self,$class;

  my ($assembled_seq,$assembled_start,$assembled_end,
      $component_seq, $component_start, $component_end,
      $component_ori ) =
    $self->_rearrange([qw(ASSEMBLED_SEQ
			  ASSEMBLED_START
			  ASSEMBLED_END
			  COMPONENT_SEQ
			  COMPONENT_START 
			  COMPONENT_END
			  COMPONENT_ORI
                         )],
		      @args);
  
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
  # make it really fast
  return bless {
     'assembled_Seq' => $_[1],
     'assembled_start' => $_[2],
     'assembled_end' => $_[3],
     'component_Seq' => $_[4],
     'component_start' => $_[5],
     'component_end' => $_[6],
     'component_ori' => $_[7],
# compatibility zone
     'start' => $_[2],
     'end' => $_[3],
     'strand' => $_[7],
     'contig' => $_[4]
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
