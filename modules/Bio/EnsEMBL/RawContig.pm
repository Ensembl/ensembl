# EnsEMBL RawContig object
#
# Copyright EMBL-EBI 2001
#
# cared for by:: Arne Stabenau
# Date : 04.12.2001
#


=head1 NAME

Bio::EnsEMBL::RawContig
  Contig object which represents part of an EMBL Clone.Mostly for database usage

=head1 SYNOPSIS

=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk

=head1 APPENDIX

=cut



package Bio::EnsEMBL::RawContig;

use vars qw( @ISA );
use strict;


use Bio::EnsEMBL::DBSQL::DBAdaptor;

@ISA = qw( Bio::Root::RootI );


# arguments
# 1: dbID database internal id
# 2: adaptor
# 3: clone this contig belongs to
# 4: sequence object, can be primarySeq
# 5: corder EMBL clone file position 
# 6: offset EMBL clone file position in bp from start
# name, length, international_name

sub new {
  my ( $class, @args ) = @_;

  my $self = {};
  bless $self, $class;
  
  

}


sub adaptor {
  my $self = shift;
  my $arg = shift;
  
  ( defined $arg ) &&
    ( $self->{_adaptor} = $arg );
  
  return $self->{_adaptor};
}

    


sub dbID {
  my $self = shift;
  my $arg = shift;
  
  ( defined $arg ) &&
    ( $self->{_dbID} = $arg );
  
  return $self->{_dbID};
}

    

sub name {
  my $self = shift;
  my $arg = shift;
  
  ( defined $arg ) &&
    ( $self->{_name} = $arg );
  
  return $self->{_name};
}

    

sub international_name {
  my $self = shift;
  my $arg = shift;
  
  ( defined $arg ) &&
    ( $self->{_international_name} = $arg );
  
  return $self->{_international_name};
}

    

sub offset {
  my $self = shift;
  my $arg = shift;
  
  ( defined $arg ) &&
    ( $self->{_offset} = $arg );
  
  return $self->{_offset};
}

    

sub corder {
  my $self = shift;
  my $arg = shift;
  
  ( defined $arg ) &&
    ( $self->{_corder} = $arg );
  
  return $self->{_corder};
}

    

sub clone {
  my $self = shift;
  my $arg = shift;
  
  ( defined $arg ) &&
    ( $self->{_clone} = $arg );
  
  return $self->{_clone};
}

    

sub length {
  my $self = shift;
  my $arg = shift;
  
  ( defined $arg ) &&
    ( $self->{_length} = $arg );
  
  return $self->{_length};
}

    


1;
