
#
# BioPerl module for SeqAdaptor
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

SeqAdaptor - DB Adaptor for ArchiveSeq objects

=head1 SYNOPSIS

    my $asad = Bio::EnsEMBL::Archive::DBSQL->new($db);
    my @aseqs = $asad->fetch_by_ensembl_id('ENSE00000023423');

=head1 DESCRIPTION

The SeqAdaptor contains all the SQL needed to fetch/store Archive::Seq 
objects in the Archive Database

=head1 CONTACT

e-mail: elia@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...
package Bio::EnsEMBL::Archive::DBSQL::SeqAdaptor;

use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Archive::Seq;

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
    my $statement = "SELECT name,type,created from seq where seq_id = $id";
    my $sth = $self->db->execute($statement);
    my($name,$type,$created) = $sth->fetchrow_array;
    my $seq = Bio::EnsEMBL::Archive::Seq->new(
					      -dbid => $id,
					      -name => $name,
					      -type => $type,
					      -created => $created,
					      -adaptor => $self
					      );
    return $seq;
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
   my ($self,$seq) = @_;

   unless (defined $seq || $seq->isa('Bio::EnsEMBL::Archive::Seq')) {
       $self->throw("Cannot store $seq, you have to provide a Bio::EnsEMBL::Archive::Seq object to store!");
   }
   if ($seq->db_ID) {
       return $seq->db_ID;
   }
   
   if (! $self->_exists($seq)) {
       my $statement = "INSERT INTO seq(seq_id,name,type,created) values (NULL,'".$seq->name."','".$seq->type."','".$seq->created."')";
       
       my $sth = $self->db->execute($statement);
       my $id = $sth->{'mysql_insertid'};
       $seq->db_ID($id);
   } 
   $seq->adaptor;
   return $seq->db_ID;
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
    my ($self,$seq) = @_;

    my $statement = "select seq_id from seq where name = '".$seq->name."' and type = '".$seq->type."'";
    my $sth = $self->db->execute($statement);
    my($id) = $sth->fetchrow_array;
    if ($id) { 
	$seq->db_ID($id);
	return 1;
    }
    else {
	return 0;
    }
}
1;
