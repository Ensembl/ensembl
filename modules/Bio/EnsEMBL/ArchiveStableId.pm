# EnsEMBL module for ArchiveStableId
# Copyright EMBL-EBI/Sanger center 2003
#
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::ArchiveStableId

=head1 SYNOPSIS

ArchiveStableId objects are the main workunit for retrieving stable id archived information from
 EnsEMBL core database.


=head1 DESCRIPTION

 Attributes:
  type: Gene, Transcript, Translation, Exon, other, undef
  stable_id: eg. ENSG00000000001
  db_name: eg. homo_sapiens_core_12_31
  version: 1

 Methods:
  new:
  new_fast:
  get_all_direct_predecessors:
  get_all_direct_successors:
  get_all_predecessors:
  get_all_successors:

  get_components:
  

=cut



package Bio::EnsEMBL::ArchiveStableId;


use warnings;
use strict;
use Bio::EnsEMBL::Root;
use vars qw(@ISA);


@ISA = qw(Bio::EnsEMBL::Root);



sub new {
  my $class = shift;
  $class = ref( $class ) || $class;

  my $self = bless {}, $class;

  my ( $stable_id, $version, $db_name, $type, $adaptor ) = 
    $self->_rearrange( [ qw( STABLE_ID VERSION DB_NAME TYPE ADAPTOR ) ], @_ );  
  
  $self->{'stable_id'} = $stable_id;
  $self->{'version'} = $version;
  $self->{'db_name'} = $db_name;
  $self->{'type'} = $type;
  $self->{'adaptor'} = $adaptor;
}




sub new_fast {
  my $class = shift;
  
  $class = ref( $class ) || $class;

  my $self = bless {
		    'stable_id' => $_[0],
		    'version' => $_[1],
		    'db_name' => $_[2],
		    'type' => $_[3],
		    'adaptor' => $_[4]
		   }, $class;
  return $self;
}


sub get_all_direct_predecessors {
  my $self = shift;

  
}

sub get_all_direct_successors {
  my $self = shift;

}

sub get_all_predecessors {
  my $self = shift;

}

sub get_all_successor {
  my $self = shift;

}

sub get_components {
  my $self = shift;

}



# getter / setter attribute section

sub type {
  my $self = shift;
  if( @_ ) {
    $self->{'type'} = shift;
  }
  return $self->{'type'};
}

sub stable_id {
  my $self = shift;
  if( @_ ) {
    $self->{'stable_id'} = shift;
  }
  return $self->{'stable_id'};
}

sub db_name {
  my $self = shift;
  if( @_ ) {
    $self->{'db_name'} = shift;
  }
  return $self->{'db_name'};
}

sub adaptor {
  my $self = shift;
  if( @_ ) {
    $self->{'adaptor'} = shift;
  }
  return $self->{'adaptor'};
}


# lazy loading 
sub version {
  my $self = shift;
  if( @_ ) {
    $self->{'version'} = shift;
  } else {
    if( ! defined $self->{'version'} ) {
      # lazy loading
      print STDERR "Lazy loading of version not yet implemented\n";
    }
  }
}

