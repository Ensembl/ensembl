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
# 2: adaptor connecting to the database

# 3: name - a name you give this contig ...

# 4: sequence object, can be primarySeq
# 5: length of sequence (if you want to avoid a fetch on the seq
#    to find out ..

# 6: clone this contig belongs to
# 7: corder EMBL clone file position 
# 8: offset EMBL clone file position in bp from start

# 9: international_name

sub new {
  my ( $class, @args ) = @_;

  my $self = {};
  bless $self, $class;
  
  my ( $dbID, $adaptor, $name, $sequence, $length,
       $clone, $corder, $offset, $international_name ) = @args;

  unless (( defined $dbID && defined $adaptor ) ||
	  ( defined $clone && defined $sequence )) {
    return undef;
  }

  (defined $dbID) && $self->dbID( $dbID );
  (defined $adaptor) && $self->adaptor( $adaptor );
  (defined $clone) && $self->clone( $clone );
  (defined $sequence) && $self->sequence( $sequence );
  (defined $name) && $self->name( $name );
  (defined $length) && $self->length( $length );
  (defined $corder) && $self->corder( $corder );
  (defined $offset) && $self->offset( $offset );
  (defined $international_name) && $self->international_name
    ( $international_name );
}




# things to take over from DBSQL/RawContig ?
# fetch - done in DBSQL/RawContigAdaptor
# get_all_Genes - should be in GeneAdaptor
#    How to handle halve off genes?
# get_old_Genes - dead
# get_all_Exons - ExonAdaptor
#   What happens with halve on/off ?
# get_old_Exons - dead
# get_Genes_by_Type - see get_all_Genes
# has_Genes - Gene or ExonAdaptor ...
# primary_seq - seq should do it now
# db_primary_seq - Assume seq does it ?
# perl_primary_seq - and again
# get_all_SeqFeatures - Call the relevant Adaptors
# get_all_SimilarityFeatures_above_score
# get_all_SimilarityFeatures
# get_all_RepeatFeatures

# Are this used from Pipeline? Web uses StaticContig
# get_MarkerFeatures
# get_landmark_MarkerFeatures


# get_all_PredictionFeatures - use PredictionTranscriptAdaptor
# get_genscan_peptides - should work via above, so deprecated
# get_all_ExternalFeatures - which adaptors provide them??
# get_all_ExternalGenes - again, we need adaptors for them
# cloneid - dead, use clone->dbID or clone->name,embl_name etc
# old_chromosome - ?? I assume unused and dead
# seq_version - Maybe clone->embl_version ?
# embl_order - good name, corder to be renamed
# embl_offset- good name, offset to be renamed
# embl_accession - Is there one for our contigs ?
# id - should be name
# internal_id - dbID
# dna_id - hidden in seq
# seq_date - hidden in seq

# static_golden_path functionality 
# All things should lie on an Assembly object.
#  get info out by supplying the contig as parameter





############################
#
#  Attribute section       #
#                          #
############################


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
  
  if( defined $arg ) {
    $self->{_name} = $arg ;
  } else {
    if( ! defined $self->{_name} &&
      defined $self->adaptor() ) {
      $self->adaptor->fetch( $self );
    }
  }
  
  return $self->{_name};
}

sub international_name {
  my $self = shift;
  my $arg = shift;
  
  if( defined $arg ) {
    $self->{_international_name} = $arg ;
  } else {
    if( ! defined $self->{_international_name} &&
      defined $self->adaptor() ) {
      $self->adaptor->fetch( $self );
    }
  }
  
  return $self->{_international_name};
}

sub offset {
  my $self = shift;
  my $arg = shift;
  
  if( defined $arg ) {
    $self->{_offset} = $arg ;
  } else {
    if( ! defined $self->{_offset} &&
      defined $self->adaptor() ) {
      $self->adaptor->fetch( $self );
    }
  }
  
  return $self->{_offset};
}

sub corder {
  my $self = shift;
  my $arg = shift;
  
  if( defined $arg ) {
    $self->{_corder} = $arg ;
  } else {
    if( ! defined $self->{_corder} &&
      defined $self->adaptor() ) {
      $self->adaptor->fetch( $self );
    }
  }
  
  return $self->{_corder};
}

sub clone {
  my $self = shift;
  my $arg = shift;
  
  if( defined $arg ) {
    $self->{_clone} = $arg ;
  } else {
    if( ! defined $self->{_clone} &&
      defined $self->adaptor() ) {
      $self->adaptor->fetch( $self );
    }
  }
  
  return $self->{_clone};
}

sub length {
  my $self = shift;
  my $arg = shift;
  
  if( defined $arg ) {
    $self->{_length} = $arg ;
  } else {
    if( ! defined $self->{_length} &&
      defined $self->adaptor() ) {
      $self->adaptor->fetch( $self );
    }
  }
  
  return $self->{_length};
}

sub seq {
  my $self = shift;
  my $arg = shift;
  
  if( defined $arg ) {
    $self->{_seq} = $arg ;
  } else {
    if( ! defined $self->{_seq} &&
      defined $self->adaptor() ) {
      $self->adaptor->fetch( $self );
    }
  }
  
  return $self->{_seq};
}



1;
