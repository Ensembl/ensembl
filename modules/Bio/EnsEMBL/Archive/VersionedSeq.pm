#
# EnsEMBL module for Bio::EnsEMBL::Archive::VersionedSeq
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Archive::VersionedSeq - Simple archive VersionedSeq object

=head1 SYNOPSIS

   my $seq = Bio::EnsEMBL::Archive::Seq->new(
					     -name => 'ENSE0006734',
					     -type => 'exon',
					     );

   my $versioned_seq = Bio::EnsEMBL::Archive::VersionedSeq->new(
				     -seq => $seq,
				     -version => 1,
				     -sequence => 'ATCGTAGAT',
				     -start_clone => 'AC000234',
                                     -seq_start => 23423,
				     -end_clone => 'AC123132',
				     -seq_end => 1243,
				     -release_version => 100,
				     );

=head1 DESCRIPTION

This object is a full representation of a certain version of a sequence in the archive. It respects teh PrimarySeqI interface as far as sequence handling goes, and it takes advantage of the adaptor schema for making the queries via mysql when it is attached to a database, rather than in memory.

This object also holds archive specific information such as the accession/coordinate pairs for start and end, the version number, and the release number in which it existed.

Finally it holds onto the simpler Archive::Seq object to point to the original ensembl identifier and creation date.

Multiple versioned objects will have the same underlying seq object.

=head1 CONTACT

Elia Stupka - elia@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Bio::EnsEMBL::Archive::VersionedSeq;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI
# (inherits for methods like throw and rearrange)

use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI Bio::PrimarySeqI);

sub new {
  my($class, @args) = @_;
  
  my $self = $class->SUPER::new(@args);
  
   my ($dbid,$archive_seq,$version,$start_clone,
       $seq_start,$end_clone,$seq_end,$sequence,
       $release,$adaptor) = 
       $self->_rearrange([qw(
	DBID
	ARCHIVE_SEQ
	VERSION
	START_CLONE
	SEQ_START
	END_CLONE
	SEQ_END
	SEQUENCE
	RELEASE_VERSION			           
	ADAPTOR
        )],@args);

  ($archive_seq && $archive_seq->isa('Bio::EnsEMBL::Archive::Seq')) || $self->throw("An Archive VersionedSeq object must have an attached Bio::EnsEMBL::Archive::Seq [got $archive_seq]");
  $version || $self->throw("An Archive VersionedSeq object must have a version");
  $start_clone || $self->throw("An Archive VersionedSeq object must have a start_clone");
  $seq_start || $self->throw("An Archive VersionedSeq object must have a seq_start");
  $end_clone || $self->throw("An Archive VersionedSeq object must have an end_clone");
  $seq_end || $self->throw("An Archive VersionedSeq object must have a seq_end");
  $release || $self->throw("An Archive VersionedSeq object must have a release number");

  #If the adaptor is attached, no need to set the sequence in memory
  if ($adaptor) {
      $dbid || $self->throw("Creating Archive VersionedSeq with adaptor but without a db_ID");
      $self->adaptor($adaptor);
      $self->db_ID($dbid);
  }
  else {
      $seq && $self->seq($sequence);
  }
  $self->archive_seq($archive_seq);
  $self->version($version);
  $self->start_clone($start_clone);
  $self->end_clone($end_clone);
  $self->seq_start($seq_start);
  $self->seq_end($seq_end);
  $self->release_version($release);

  return $self;
}

=head2 display_id

 Title   : display_id
 Usage   : $self->display_id($newval)
 Function: get/set method for the dna id
 Returns : value of display_id
 Args    : newvalue (optional) 

=cut

sub display_id {
    my ($self,$value) = @_;

    return $self->archive_seq->name;
} 

=head2 primary_id

 Title   : primary_id
 Usage   : $self->primary_id($newval);
 Function: get/set method for the primary id
 Returns : value of primary id
 Args    : None


=cut

sub primary_id {
    my ($self,$value) = @_;

    return $self->archive_seq->name;
} 

=head2 accession_number

 Title   : accession_number
 Usage   : $self->accession_number($newval)
 Function: get/set method for the dna id
 Returns : value of accession_number
 Args    : newvalue (optional) 

=cut

sub accession_number {
    my ($self) = @_;

    return $self->archive_seq->name;
}

=head2 seq

 Title   : seq
 Usage   : $string    = $self->seq()
 Function: Returns the sequence as a string of letters.
 Returns : A scalar
 Args    : none

=cut

sub seq {
   my $self = shift;

   if (! $self->{'seq'}) {
       if( @_ ) {
           my $value = shift;
           $self->{'seq'} = $value;
       }
       elsif ($self->adaptor) {
           my $value = $self->adaptor->seq($self->db_ID);
           $self->{'seq'} = $value;
       }
   }
   return $self->{'seq'};
}

=head2 subseq

 Title   : subseq
 Usage   : $substring = $self->subseq(10,40);
 Function: returns the subseq from start to end, where the first base
           is 1 and the number is inclusive, ie 1-2 are the first two
           bases of the sequence
           Start cannot be larger than end but can be equal
 Returns : a string
 Args    : start and end scalars

=cut

sub subseq{
   my ($self,$start,$end) = @_;

   if( $start > $end ){
       $self->throw("in subseq, start [$start] has to be greater than end [$end]");
   }
   
   if( $start <= 0 || $end > $self->length ) {
       $self->throw("You have to have start positive and length less than the total length of sequence [$start:$end] Total ".$self->length."");
   }

   if ($self->adaptor) {
       return $self->adaptor->subseq($self->db_ID,$start,$end);
   }
   else {
       # remove one from start, and then length is end-start
       $start--;
       return substr $self->{'seq'}, $start, ($end-$start);
   }
}


=head2 moltype

 Title   : moltype
 Usage   : if( $self->moltype eq 'dna' ) { /Do Something/ }
 Function: Returns the type of sequence 
 Returns : dna
 Args    : none


=cut

sub moltype{
    my $type = $self->archive->type;

    if ($type =~ /exon|transcript/) {
	return 'dna';
    }
    elsif ($type =~ /translation/) {
	return 'protein';
    }
    else {
	return undef;
    }
}

=head2 length

 Title   : length
 Usage   : $len = $seq->length()
 Function: Returns the length of the sequence
 Returns : scalar
 Args    : none


=cut

sub length {
    my $self = shift;
    
    if (! $self->{'length'}) {
        if ($self->adaptor) {
            my $value = $self->adaptor->length($self->db_ID);
            $self->{'length'} = $value;
        }
        else {
            my $value = CORE::length($self->{'seq'});
            $self->{'length'} = $value;
        }
    }
    return $self->{'length'};
}

=head2 desc

 Title   : desc
 Usage   : $self->desc($newval)
 Function: 
 Example : 
 Returns : value of desc
 Args    : newvalue (optional)


=cut

sub desc {
   my ($self,$value) = @_;
   if (! defined $self->{'desc'} ) {
       if( defined $value && $value ne '' ) {
	   $self->{'desc'} = $value;
       } 
       else {
	   $self->{'desc'} = $self->archive_seq->type;
       }
   }
   return $self->{'desc'};
}

=head2 archive_seq

 Title   : archive_seq
 Usage   : $obj->archive_seq($newval)
 Function: Getset for archive_seq value
 Returns : value of archive_seq
 Args    : newvalue (optional)


=cut

sub archive_seq{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'archive_seq'} = $value;
    }
    return $obj->{'archive_seq'};

}

=head2 version

 Title   : version
 Usage   : $obj->version($newval)
 Function: Getset for version value
 Returns : value of version
 Args    : newvalue (optional)


=cut

sub version{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'version'} = $value;
    }
    return $obj->{'version'};

}

=head2 start_clone

 Title   : start_clone
 Usage   : $obj->start_clone($newval)
 Function: Getset for start_clone value
 Returns : value of start_clone
 Args    : newvalue (optional)


=cut

sub start_clone{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'start_clone'} = $value;
    }
    return $obj->{'start_clone'};

}

=head2 start

 Title   : start
 Usage   : $obj->start($newval)
 Function: Getset for start value
 Returns : value of start
 Args    : newvalue (optional)


=cut

sub start{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'start'} = $value;
    }
    return $obj->{'start'};

}

=head2 end_clone

 Title   : end_clone
 Usage   : $obj->end_clone($newval)
 Function: Getset for end_clone value
 Returns : value of end_clone
 Args    : newvalue (optional)


=cut

sub end_clone{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'end_clone'} = $value;
    }
    return $obj->{'end_clone'};

}

=head2 end

 Title   : end
 Usage   : $obj->end($newval)
 Function: Getset for end value
 Returns : value of end
 Args    : newvalue (optional)


=cut

sub end{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'end'} = $value;
    }
    return $obj->{'end'};

}

=head2 release_version

 Title   : release_version
 Usage   : $obj->release_version($newval)
 Function: Getset for release_version value
 Returns : value of release_version
 Args    : newvalue (optional)


=cut

sub release_version{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'release_version'} = $value;
    }
    return $obj->{'release_version'};
}


=head2 adaptor

 Title   : adaptor
 Usage   : $obj->adaptor($newval)
 Function: Getset for adaptor value
 Returns : value of adaptor
 Args    : newvalue (optional)


=cut

sub adaptor{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'adaptor'} = $value;
    }
    return $obj->{'adaptor'};

}

=head2 db_ID

 Title   : db_ID
 Usage   : $obj->db_ID($newval)
 Function: Getset for db_ID value
 Returns : value of db_ID
 Args    : newvalue (optional)


=cut

sub db_ID{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'db_ID'} = $value;
    }
    return $obj->{'db_ID'};

}
