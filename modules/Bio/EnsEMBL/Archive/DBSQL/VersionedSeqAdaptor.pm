
#
# BioPerl module for VersionedSeqAdaptor
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

VersionedSeqAdaptor - DB Adaptor for ArchiveSeq objects

=head1 SYNOPSIS

    my $asad = Bio::EnsEMBL::Archive::DBSQL->new($db);
    my @aseqs = $asad->fetch_by_ensembl_id('ENSE00000023423');

=head1 DESCRIPTION

The VersionedSeqAdaptor contains all the SQL needed to fetch/store ArchiveSeq objects in the Archive Database

=head1 CONTACT

e-mail: elia@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...
package Bio::EnsEMBL::Archive::DBSQL::VersionedSeqAdaptor;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Archive::VersionedSeq;

#The method ->new is inherited from the BaseAdaptor
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 fetch_by_dbID

 Title   : fetch_by_dbID
 Usage   :
 Function:
 Example :
 Returns : Bio::EnsEMBL::Overlap::Overlap
 Args    :


=cut

sub fetch_by_dbID {
    my ($self,$id) = @_;

    $self->throw("I need a dbID") unless $id;
    my $statement = "SELECT seq_id,version,start_contig,start_coord,end_contig,end_coord,modified,release_number from versioned_seq where versioned_seq_id = $id";
    my $sth = $self->db->execute($statement);
    my($seq_id,$v,$start_c,$start,$end_c,$end,$mod,$rel) = $sth->fetchrow_array;
    my $sad = $self->db->get_SeqAdaptor;
    my $seq = $sad->fetch_by_dbID($seq_id);
    my $vseq = Bio::EnsEMBL::Archive::VersionedSeq->new(
					      -dbid => $id,
					      -archive_seq => $seq,
					      -version => $v,
					      -start_contig => $start_c,
					      -start => $start,
					      -end_contig => $end_c,
					      -end => $end,
					      -modified => $mod,
					      -release_number => $rel,
					      -adaptor => $self
					      );
    
    $statement = "SELECT relative_versioned_seq_id from versioned_seq_relatives where master_versioned_seq_id = $id";
    my $sth = $self->db->execute($statement);
    while (my ($rid) = $sth->fetchrow_array) {
	$vseq->add_relative($self->fetch_by_dbID($rid));
    }
    
    $statement = "SELECT new_versioned_seq_id from versioned_seq_history where old_versioned_seq_id = $id";
    my $sth = $self->db->execute($statement);
    while (my ($fid) = $sth->fetchrow_array) {
	$vseq->add_future_vseq($self->fetch_by_dbID($fid));
    }
    return $vseq;
}

=head2 store

 Title   : store
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub store {
   my ($self,$vseq) = @_;

   unless (defined $vseq || $vseq->isa('Bio::EnsEMBL::Archive::VersionedSeq')) {
       $self->throw("Cannot store $vseq, you have to provide a Bio::EnsEMBL::Archive::Vseq object to store!");
   }
   if ($vseq->db_ID) {
       $self->warn("Trying to store this versionedseq ".$vseq->name." twice");
       return $vseq->db_ID;
   }
   
   my $seqda = $self->db->get_SeqAdaptor;
   my $seq_id = $seqda->store($vseq->archive_seq);
   
   if (! $self->_exists($vseq)) {
       my $statement;
       if ($vseq->archive_seq->type ne 'gene') {
	   (! defined $vseq->seq) && $self->warn("Trying to store a versioned seq without sequence for a ".$vseq->archive_seq->type." sequence");
	   
	   $statement = "INSERT INTO versioned_seq(versioned_seq_id,seq_id,version,sequence,start_contig,start_coord,end_contig,end_coord,modified,release_number) values (NULL,$seq_id,".$vseq->version.",'".$vseq->seq."','".$vseq->start_contig."',".$vseq->start.",'".$vseq->end_contig."',".$vseq->end.",'".$vseq->modified."',".$vseq->release_number.")";
       }
       else {
	   $statement = "INSERT INTO versioned_seq(versioned_seq_id,seq_id,version,sequence,start_contig,start_coord,end_contig,end_coord,modified,release_number) values (NULL,$seq_id,".$vseq->version.",NULL,'".$vseq->start_contig."',".$vseq->start.",'".$vseq->end_contig."',".$vseq->end.",'".$vseq->modified."',".$vseq->release_number.")";
       }
       my $sth = $self->db->execute($statement);
       my $id = $sth->{'mysql_insertid'};
       $vseq->db_ID($id);
       
       foreach my $relative ($vseq->each_relative) {
	   my $rid = $self->store($relative);
	   $statement = "INSERT INTO versioned_seq_relatives(master_versioned_seq_id,relative_versioned_seq_id) values ($id,$rid)";
	   $self->db->execute($statement);
       }
       
       foreach my $future_vseq ($vseq->each_future_vseq) {
	   my $fid = $self->store($future_vseq);
	   $statement = "INSERT INTO versioned_seq_history(old_versioned_seq_id,new_versioned_seq_id) values ($id,$fid)";
	   $self->db->execute($statement);
       }
   } 
   $vseq->adaptor;
   return $vseq->db_ID;
}

=head2 _exists

 Title   : _exists
 Usage   :
 Function: Check if this Seq exists already in the database
 Example :
 Returns : 
 Args    :

=cut

sub _exists{
    my ($self,$vseq) = @_;

    my $statement = "select versioned_seq_id from versioned_seq where seq_id = '".$vseq->archive_seq->db_ID."' and version = ".$vseq->version;
    my $sth = $self->db->execute($statement);
    my($id) = $sth->fetchrow_array;
    if ($id) { 
	$vseq->db_ID($id);
	return 1;
    }
    else {
	return 0;
    }
}

#Methods below are to provide quick mysql version of the PrimarySeqI methods

=head2 seq

 Title   : seq
 Usage   : $string    = $obj->seq()
 Function: Returns the sequence as a string of letters.
 Returns : A scalar
 Args    : none

=cut

sub seq {
   my ($self,$id) = @_;
   
   $self->throw("I need an id") unless $id;

   my $sth=$self->prepare("SELECT sequence FROM versioned_seq WHERE versioned_seq_id = $id");
   $sth->execute();
   my($str) = $sth->fetchrow;
   return $str;
}

=head2 subseq

 Title   : subseq
 Usage   : $substring = $obj->subseq(10,40);
 Function: returns the subseq from start to end, where the first base
           is 1 and the number is inclusive, ie 1-2 are the first two
           bases of the sequence
           Start cannot be larger than end but can be equal
 Returns : a string
 Args    : start and end scalars

=cut

sub subseq{
    my ($self,$id,$start,$end) = @_;
    
    $self->throw("I need an id") unless $id;
    $self->throw("I need a start") unless $start;
    $self->throw("I need an end") unless $end;

    my $length= $end-$start+1;
    
    my $sth=$self->prepare("SELECT SUBSTRING(sequence,$start,$length) FROM versioned_seq WHERE versioned_seq_id = $id");
    $sth->execute(); 
    
    my($subseq) = $sth->fetchrow;
    
    return $subseq;
}

=head2 length

 Title   : length
 Usage   : $len = $seq->length()
 Function: Returns the length of the sequence
 Returns : scalar
 Args    : none


=cut

sub length {
    my ($self,$id) = @_;

    $self->throw("I need an id") unless $id;

    my $sth=$self->prepare("SELECT length(sequence) FROM versioned_seq WHERE versioned_seq_id = $id");
    $sth->execute(); 
    
    my ($length) = $sth->fetchrow;
    
    return $length;
}

1;








