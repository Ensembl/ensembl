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

  get_components:
  

=cut



package Bio::EnsEMBL::ArchiveStableId;


use warnings;
use strict;
use Bio::EnsEMBL::Root;
use vars qw(@ISA);


@ISA = qw(Bio::EnsEMBL::Root);



=head2 new

  Arg  1     : -stable_id $stable_id 
  Arg [ ]    : -version $version 
  Arg [ ]    : -db_name $db_name 
  Arg [ ]    : -adaptor $adaptor 
  Arg [ ]    : -type $type 
  "Gene", "Transcript", "Translation", "Exon"
  Example    : none
  Description: standard constructor with named arguments to create ArchiveStableId
  Returntype : Bio::EnsEMBL::ArchiveStableId
  Exceptions : none
  Caller     : Adaptor

=cut


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

  return $self;
}



=head2 new_fast

  Arg [1]    : string $stable_id
  Arg [2]    : int $version
  Arg [3]    : string $db_name
  Arg [4]    : string $type
  Arg [5]    : Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor $adaptor
  Example    : none
  Description: faster version of above constructor
  Returntype : Bio::EnsEMBL::ArchiveStableId
  Exceptions : none
  Caller     : general, Adaptor

=cut


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


=head2 get_all_predecessors

  Args       : none
  Example    : none
  Description: Retrieve a list of ArchiveStableIds that were mapped to this one. 
  Returntype : listref Bio::EnsEMBL::ArchiveStableId
  Exceptions : none
  Caller     : general

=cut


sub get_all_predecessors {
  my $self = shift;

  $self->adaptor->fetch_pre_by_arch_id( $self );
}

=head2 get_all_successors

  Args       : none
  Example    : none
  Description: Retrieve a list of ArchiveStableIds that this one was mapped to.
  Returntype : listref Bio::EnsEMBL::ArchiveStableId
  Exceptions : none
  Caller     : general

=cut

sub get_all_successors {
  my $self = shift;

  $self->adaptor->fetch_succ_by_arch_id( $self );
}



=head2 get_peptide

  Args       : none
  Example    : none
  Description: Retrieves the peptide string for this ArchiveStableId.
               Undef if this is not a Translation or cant be found in the database.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub get_peptide {
  my $self = shift;

  if( $self->type() eq "Translation" ) {
    return $self->adaptor->get_peptide( $self );
  } else { 
    return undef;
  }
}


=head2 get_all_by_transcript_archive_id

  Args       : none
  Example    : none
  Description: If this is a genes ArchiveStableId and found in the database, this 
    function gets the transcripts archiveStableIds from it. Returns undef otherwise.
  Returntype : listref Bio::EnsEMBL::ArchiveStableId
  Exceptions : empty if not a gene stable id or not in database
  Caller     : general

=cut


sub get_all_transcript_archive_ids {
  my $self = shift;

  if( $self->type() eq "Gene" ) {
    return $self->adaptor->fetch_all_by_gene_archive_id( $self );
  } else {
    return undef;
  }
}



=head2 get_translation_archive_id

  Args       : none
  Example    : none
  Description: Retrieves the Translation ArchiveStableId for this transcript stable id. 
    If not found or this is not a transcripts id return undef
  Returntype : Bio::EnsEMBL::ArchiveStableId
  Exceptions : undef if not in db or not a Transcript
  Caller     : general

=cut


sub get_translation_archive_id {
  my $self = shift;

  if( $self->type() eq "Transcript" ) {
    return [$self->adaptor->fetch_by_transcript_archive_id( $self )];
  } elsif( $self->type() eq "Gene" ) {
    my $transcripts = $self->adaptor->fetch_all_by_gene_archive_id( $self );
	my @peptides ;
	for (@$transcripts) {
		push @peptides , $self->adaptor->fetch_by_transcript_archive_id( $_ );
	}
	return \@peptides;
  } else {
    return undef;
  }
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
      if( defined $self->{'db_name'} && defined $self->{'adaptor'} ) {
	# lazy loading
	$self->{'adaptor'}->_lookup_version( $self );
      }
    }       
  }
  return $self->{'version'};
}

1;
