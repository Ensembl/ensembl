#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::RawContigAdaptor
#
# Cared for by Imre Vastrik <vastrik@ebi.ac.uk>
#
# Copyright Imre Vastrik
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::RawContigAdaptor - MySQL database adapter class for EnsEMBL Feature Objects

=head1 SYNOPSIS



=head1 DESCRIPTION



=head1 CONTACT



=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are 
usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::DBSQL::RawContigAdaptor;

use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object


use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use DBI;
use Bio::EnsEMBL::DBSQL::DummyStatement;
use Bio::EnsEMBL::DBSQL::DBPrimarySeq;

use Bio::EnsEMBL::Clone;
use Bio::EnsEMBL::RawContig;


@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

# new() is inherited from Bio::EnsEMBL::DBSQL::BaseAdaptor

sub get_internal_id_by_id
{
    my ($self, $id) = @_;
    my $sth = $self->db->prepare
    (
         "select contig_id from contig where name = '$id'"
    );
    my $res = $sth->execute;
    if(my $rowhash = $sth->fetchrow_hashref) {
	return $rowhash->{contig_id};
    } else {
	$self->warn("Could not find contig with id $id");
    }
}

sub get_id_by_contig_id
{
    my ($self, $contig_id) = @_;
    my $sth = $self->db->prepare
    (
         "select name from contig where contig_id = '$contig_id'"
    );
    my $res = $sth->execute;
    if(my $rowhash = $sth->fetchrow_hashref) {
	return $rowhash->{name};
    } else {
	$self->warn("Could not find contig with contig_id $contig_id");
    }
}

#  contig_id | name                | clone_id | length | offset | corder | dna_id | chromosome_id | international_name 

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  my $sth = $self->prepare( "SELECT contig_id, name, clone_id, length, 
                          offset, corder, dna_id, 
                          international_name
                   FROM contig
                   WHERE contig_id = $dbID" );
  $sth->execute();
  
  my ( $contig ) = $self->_contig_from_sth( $sth );

  return $contig;
}


sub fetch_all {
  my $self  = shift;
  my @res;

  my $sth = $self->prepare( "SELECT contig_id, dna_id
                             FROM contig " );
  $sth->execute();
  while( my $aref = $sth->fetchrow_arrayref() ) {
    my $dbPrimarySeq = Bio::EnsEMBL::DBSQL::DBPrimarySeq->new
      ( $aref->[1], $self->db() ); # ?
    
    my $contig = Bio::EnsEMBL::RawContig->new( $aref->[0], $self );
    $contig->seq( $dbPrimarySeq );
    push( @res, $contig );
  }
  return @res;
}



sub fetch_by_name {
  my $self = shift;
  my $name = shift;

  my $sth = $self->prepare( "SELECT contig_id, name, clone_id, length, 
                             offset, corder, dna_id, 
                             international_name
                             FROM contig
                             WHERE name = '$name'" );
  $sth->execute();
  
  my ( $contig ) = $self->_contig_from_sth( $sth );

  return $contig;
}


sub fetch_by_clone {
  my $self = shift;
  my $clone = shift;

  my $clone_id = $clone->dbID;
  
  my $sth = $self->prepare( "SELECT contig_id, name, clone_id, length, 
                             offset, corder, dna_id, chromosome_id, 
                             international_name
                             FROM contig
                             WHERE clone_id = $clone_id" );

  my @res = $self->_contig_from_sth( $sth );
  return \@res;
}


# mainly used from RawContig object. 
# Argument is a ready contig which needs its attributes filled

sub fetch {
  my $self = shift;
  my $contig = shift;
  
  my $dbID = $contig->dbID();

  my $sth = $self->prepare( "SELECT contig_id, name, clone_id, length, 
                             offset, corder, dna_id, international_name
                             FROM contig
                             WHERE contig_id = $dbID" );
  $sth->execute();
  
  my $aref = $sth->fetchrow_arrayref();
  if( defined $aref ) {
    $self->_fill_contig_from_arrayref( $contig, $aref );
  } else {
    $self->throw( "Couldnt fetch contig, unexpected .." );
  }
}


sub _contig_from_sth {
  my $self = shift;
  my $sth = shift;

  if( !defined $sth ) {
      $self->throw("Bad internal error - no statement!");
  }

  my @res = ();

  $sth->execute();
  while( my $aref = $sth->fetchrow_arrayref() ) {
    
    my $contig = Bio::EnsEMBL::RawContig->new( $aref->[0], $self );
    $self->_fill_contig_from_arrayref( $contig, $aref );

    push( @res, $contig );
  }

  return @res;
}



sub _fill_contig_from_arrayref {
  my $self   = shift;
  my $contig = shift;
  my $aref   = shift;

  if( !defined $aref ) {
      $self->throw("Bad internal error - no array ref");
  }

  my ( $contig_id, $name, $clone_id, $length, $offset, $corder, $dna_id,
       $international_name ) = @$aref;

    
  my $dbPrimarySeq = Bio::EnsEMBL::DBSQL::DBPrimarySeq->new
    ( $dna_id, $self->db() ); # ?

  # my $contig = Bio::EnsEMBL::RawContig->new($contig_id,$self);


  (defined $dbPrimarySeq) && $contig->sequence( $dbPrimarySeq );
  (defined $length) && $contig->length( $length );
  (defined $name) && $contig->name( $name );
  (defined $offset) && $contig->offset( $offset );
  (defined $corder) && $contig->corder( $corder );
  (defined $international_name) && $contig->international_name( $international_name );

  # maybe these should be lazy fetched...
  $contig->_clone_id($clone_id);

  #my $clone = $self->db->get_CloneAdaptor->fetch_by_dbID($clone_id);
  #$contig->clone( $clone );
  

  return $contig;
}


sub store{
  my($self, $contig, $clone_id) = @_;

  $self->throw("$contig is not a Bio::EnsEMBL::DB::ContigI - cannot insert contig for clone $clone_id")
    unless $contig->isa('Bio::EnsEMBL::DB::ContigI');   
  my $dna = $contig->primary_seq  || $self->throw("No sequence in contig object");
  $dna->id                        || $self->throw("No contig id entered.");
  $clone_id                          || $self->throw("No clone_id entered.");


  $self->_insertSequence($dna->seq, $contig->seq_date);
  my $international_name = $contig->international_name;
  if(!$international_name){
    $international_name = 'NULL';
  } 
  my $sql = "insert into contig(name,
                                dna_id,
                                length,
                                clone_id,
                                offset,
                                corder,
                                international_name)
              values('".$contig->id."', 
                    LAST_INSERT_ID(), 
		    ".$contig->primary_seq->length." ,
                    ".$clone_id." ,
                    ".$contig->embl_offset." ,
                    ".$contig->order." , 
                    '".$international_name."')";
  
  my $sth = $self->prepare($sql);
  my $rv = $sth->execute();
  $self->throw("Failed to insert contig ".$contig->id."\n") unless $rv;
       
    
    $sth = $self->prepare("select last_insert_id()");
    $sth->execute;
    my ($id) = $sth->fetchrow
        or $self->throw("Failed to get last insert id");
    #can no longer do this as get_all_SeqFeatures no longer exists
    #if a contig is written to the database
    # this is a nasty hack. We should have a cleaner way to do this.
    #my @features = $contig->get_all_SeqFeatures;
    #print(STDERR "Contig $contigid - $id\n"); 
    # write sequence features. We write all of them together as it
    # is more efficient
    #$self->get_Feature_Obj->write($contig, @features);
    
    return 1;
}


sub _insertSequence{

   my ($self, $sequence, $date) = @_;
    
    $sequence =~ tr/atgcn/ATGCN/;
    

  
    
    my $statement = $self->prepare("
        insert into dna(sequence,created) 
        values(?, FROM_UNIXTIME(?))
        "); 
        
    my $rv = $statement->execute($sequence, $date); 
    
    $self->throw("Failed to insert dna $sequence") unless $rv;   


}



1;
