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

Bio::EnsEMBL::DBSQL::RawContigAdaptor - MySQL database adapter class 
for EnsEMBL RawContig Objects

=head1 SYNOPSIS

$contig_adaptor = $database_adaptor->get_RawContigAdaptor();
$contig = $contig_adaptor->fetch_by_dbID(1234); 

=head1 DESCRIPTION

Allows for the retrieval and storage of RawContig objects from the database.

=head1 CONTACT

Post questions to the EnsEMBL developer list <ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::DBSQL::RawContigAdaptor;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use DBI;


use Bio::EnsEMBL::Clone;
use Bio::EnsEMBL::RawContig;

use Bio::EnsEMBL::Utils::Cache; # CPAN LRU Cache module
my $RAW_CONTIG_CACHE_SIZE = 20;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 new

  Arg [1]    : list of arguments @args
  Example    : $contig_adaptor = new Bio::EnsEMBL::RawContigAdaptor($db);
  Description: Creates a new RawContigAdaptor
  Returntype : Bio::EnsEMBL::RawContig
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::DBConnection

=cut

sub new {
  my($class, @args) = @_;

  #call super class constructor
  my $self = $class->SUPER::new(@args);

  #Initialize caching data structures
  tie(%{$self->{_raw_contig_cache}}, 
      'Bio::EnsEMBL::Utils::Cache', 
      $RAW_CONTIG_CACHE_SIZE);

  return $self;
}



=head2 fetch_by_dbID

  Arg [1]    : int $dbID
               the database identifier for the contig to retrieve
  Example    : $contig = $raw_contig_adaptor->fetch_by_dbID(1234);
  Description: Retreives a RawContig object from the database.  RawContigs
               retrieved by this method are in fact lazy loaded and not
               populated at all.  Once data is actually requested from a 
               RawContig object the attributes are pulled from the database.
               This is far faster for many requests which do not require 
               populated Contigs.  However, if a list of fully filled contigs
               is going to be required, then it may be much faster to use
               the fetch_filled_by_dbIDs method. 
  Returntype : Bio::EnsEMBL::RawContig 
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  # retrieve the raw contig from the cache
  my $raw_contig = $self->{_raw_contig_cache}->{$dbID}; 
     
  unless($raw_contig) {
    # the cache did not contain a raw contig, create a new empty 'lazy loaded'
    # contig under the assumption that the contig data will not be 
    # needed.  If the data is needed then it can be loaded from the DB by
    # other db adaptor functions

    $raw_contig = new Bio::EnsEMBL::RawContig($dbID, $self);

    #store the new raw contig in the cache
    $self->{_raw_contig_cache}->{$dbID} = $raw_contig;
  }

  return $raw_contig;
}


=head2 fetch_all

  Arg [1]    : none
  Example    : @contigs = @{$raw_contig_adaptor->fetch_all()};
  Description: Retrieves all of the contig objects in the database.  This will
               be a very slow request on a full-sized database
  Returntype : listref of Bio::EnsEMBL::RawContig
  Exceptions : none
  Caller     : general

=cut

sub fetch_all {
  my $self  = shift;
  my @res;

  my $sth = $self->prepare( "SELECT contig_id, name, clone_id, length,
                             embl_offset, dna_id
                             FROM contig " );
  $sth->execute();

  return $self->_contig_from_sth( $sth );
}

=head2 fetch_by_name

  Arg [1]    : string $name
               the name of the contig to retrieve
  Example    : $contig = $rca->fetch_by_name('AC004501.1.1.39507');
  Description: Retrieves a contig from the database using the contig name
  Returntype : Bio::EnsEMBL::RawContig
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_name {
  my $self = shift;
  my $name = shift;

  my $sth = $self->prepare( "SELECT contig_id, name, clone_id, length, 
                             embl_offset, dna_id
                             FROM contig
                             WHERE name = '$name'" );
  $sth->execute();
  
  my ( $contig ) = @{$self->_contig_from_sth( $sth )};

  return $contig;
}


=head2 fetch_filled_by_dbIDs

  Arg [1]    : list of ints @contig_ids
               The list of contigs to retrieve
  Example    : @contigs = $raw_contig_adaptor->fetch_filled_by_dbIDs(120, 121);
  Description: Returns a hashref of RawContigs retrieved from the database.  
               The contigs are fully filled, and their attached clone is 
               retrieved fully filled as well. 
  Returntype : hasref of RawContigs, with dbIDs as keys
  Exceptions : none
  Caller     : Bio::EnsEMBL::Slice::get_tiling_path

=cut

sub fetch_filled_by_dbIDs {
  my $self = shift;
  my @contig_ids = @_;
  my %result = ();

  if(scalar(@contig_ids) == 0) {
    return ();
  }

  my $sth = $self->prepare( "SELECT co.contig_id, co.name, co.clone_id, 
                                    co.length, co.embl_offset, co.dna_id,
                                    cl.embl_acc, cl.embl_version,
                                    cl.name, cl.version, cl.htg_phase,
                                    UNIX_TIMESTAMP(cl.created), UNIX_TIMESTAMP(cl.modified)
                             FROM contig co, clone cl
                             WHERE co.contig_id IN ( " .
			           join( ", ", @contig_ids ) . ")  
                             AND co.clone_id = cl.clone_id"  );



  $sth->execute();
  
  while( my $aref = $sth->fetchrow_arrayref() ) {
    my($contig_id, $contig_name, $contig_clone_id, $contig_length, 
       $contig_embl_offset, $contig_dna_id, $clone_embl_acc, 
       $clone_embl_version, $clone_name, $clone_version, $clone_htg_phase, 
       $clone_created, $clone_modified) = @$aref;

    my $contig = Bio::EnsEMBL::RawContig->new( $contig_id, $self );
    $self->_fill_contig_from_arrayref( $contig,$aref );
    my $clone = Bio::EnsEMBL::Clone->new
      (
       $self->db->get_CloneAdaptor(),
       $contig_clone_id, $clone_name,
       $clone_embl_acc, $clone_version, ,
       $clone_embl_version, $clone_htg_phase,
       $clone_created, $clone_modified
      );
    $contig->clone( $clone );

    $result{ $contig->dbID() } = $contig;

    $self->{_raw_contig_cache}->{$contig->dbID()} = $contig;
  }

  return \%result;
}


=head2 fetch_all_by_Clone

  Arg [1]    : Bio::EnsEMBL::Clone $clone
               The clone object that the contigs are desired from 
  Example    : $contigs = raw_contig_adaptor->fetch_all_by_Clone($clone);
  Description: Returns a list reference of contig objects on a particular
               clone.
  Returntype : list reference of Bio::EnsEMBL::RawContig objects
  Exceptions : none
  Caller     : general

=cut

sub fetch_all_by_Clone {
  my $self = shift;
  my $clone = shift;

  my $clone_id = $clone->dbID;
  
  my $sth = $self->prepare( "SELECT contig_id, name, clone_id, length, 
                             embl_offset, dna_id
                             FROM contig
                             WHERE clone_id = $clone_id" );

  return $self->_contig_from_sth( $sth );
}





=head2 fetch_attributes

  Arg [1]    : Bio::EnsEMBL::RawContig $contig 
  Example    : $raw_contig_adaptor->fetch_attributes($contig);
  Description: Fills in the attributes of a lazy-loaded RawContig object.
               RawContig.  
  Returntype : none
  Exceptions : thrown if the attributes for a contig could not be filled
  Caller     : Bio::EnsEMBL::RawContig

=cut

sub fetch_attributes {
  my $self = shift;
  my $contig = shift;
  
  my $dbID = $contig->dbID();

  unless($dbID && $contig->adaptor) {
    #do nothing, we do not know where to get the attributes from
    return;
  }

  my $sth = $self->prepare( "SELECT contig_id, name, clone_id, length, 
                             embl_offset, dna_id
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



=head2 store

  Arg [1]    : Bio::EnsEMBL::RawContig $contig
               The contig to store in the database.
  Arg [2]    : Bio::EnsEMBL::Clone $clone 
               The clone the contig is on
  Example    : my $contig_id = $raw_contig_adaptor->store($contig, $clone);
  Description: Stores a contig and its associated DNA sequence in the database,
               returns the database id of the new database record.  Attaches
               this adaptor to the the contig object and sets its dbID on 
               success.
  Returntype : int
  Exceptions : thrown if $contig arg is not defined, $clone arg is not a clone
               or if the database insertion fails. 
  Caller     : Bio::EnsEMBL::Clone::store

=cut

sub store {
  my ($self, $contig, $clone) = @_;

  unless($contig && $contig->isa('Bio::EnsEMBL::RawContig')) {
    $self->throw('contig arg is not a Bio::EnsEMBL::RawContig');
  }
  unless($clone && $clone->isa('Bio::EnsEMBL::Clone')) {
    $self->throw('clone arg is not a Bio::EnsEMBL::Clone');
  } 

  my $date     = $clone->created();
  my $dna      = $contig->seq();
  my $dna_id   = $self->db->get_SequenceAdaptor->store($dna, $date);
  
  my $sth = $self->prepare("INSERT INTO contig (name, 
                                                clone_id, 
                                                length, 
                                                embl_offset, 
                                                dna_id) 
                            VALUES ( ?, ?, ?, ?, ? )");
  my $rv = $sth->execute($contig->name(), 
                         $clone->dbID(),
                         $contig->length(),
                         $contig->embl_offset(),
                         $dna_id);

  unless($rv) {
    $self->throw("Failed to insert contig " . $contig->name . "into database");
  }
  
  $sth->finish;

  $sth = $self->prepare("SELECT last_insert_id()");
  $sth->execute();
  my ($contig_id) = $sth->fetchrow_array();

  $sth->finish;

  $contig->dbID($contig_id);
  $contig->adaptor($self);

  return $contig_id;
}


=head2 _contig_from_sth

  Arg [1]    : DBI Statementhandle $sth
  Example    : @contigs = $self->_contig_from_sth($sth);
  Description: PRIVATE creates a contig from an executed  DBI statement handle.
  Returntype : listref of Bio::EnsEMBL::RawContig
  Exceptions : thrown if the statement handle is not defined
  Caller     : internal

=cut

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

  return \@res;
}



=head2 _fill_contig_from_arrayref

  Arg [1]    : DBI array reference $aref
  Arg [2]    : Bio::EnsEMBL::RawContig $contig
  Example    : $contig = $self->_fill_contig_from_arrayref($arr_ref); 
  Description: PRIVATE Populates an empty RawContig object from an 
               array ref to an SQL query
  Returntype : Bio::EnsEMBL::RawContig
  Exceptions : thrown if 
  Caller     : internal

=cut

sub _fill_contig_from_arrayref {
  my $self   = shift;
  my $contig = shift;
  my $aref   = shift;

  if( !defined $aref ) {
      $self->throw("Bad internal error - no array ref");
  }

  my ( $contig_id, $name, $clone_id, $length, $offset, $dna_id) = @$aref;

    
  (defined $length) && $contig->length( $length );
  (defined $name) && $contig->name( $name );
  (defined $offset) && $contig->embl_offset( $offset );

  # maybe these should be lazy fetched...
  $contig->_clone_id($clone_id);

  #my $clone = $self->db->get_CloneAdaptor->fetch_by_dbID($clone_id);
  #$contig->clone( $clone );
  

  return $contig;
}


=head2 deleteObj

  Arg [1]    : none
  Example    : none
  Description: Cleans up object references contained within this object so
               that proper garbage collection can occur.
  Returntype : none
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::DBConnection::deleteObj

=cut

sub deleteObj {
  my $self = shift;

  #print STDERR "\t\tRawContigAdaptor::deleteObj\n";

  $self->SUPER::deleteObj();

  #flush internal cache
  %{$self->{'_raw_contig_cache'}} = ();
}



=head2 remove

  Arg [1]    : Bio::EnsEMBL::RawContig $contig
  Example    : $rawcontig_adaptor->remove($contig)
  Description: This removes a contig (itself) plus any attached features and dna.
               The method loops over available 5 available feature adaptors and
               removes features on them individually.
  Returntype : none
  Exceptions : Throw if unable to get one of the feature adaptors.
               Throw if unable to identify the dna to remove.
               Throw if unable to remove record from dna table.
               Throw if unable to remove contig.
  Caller     : Bio::EnsEMBL::DBSQL::CloneAdaptor

=cut

sub remove {
  my ($self, $contig) = @_;

  # The list of feature adaptors to be looped over
  my @adaptor_list = ($self->db()->get_SimpleFeatureAdaptor(),
		      $self->db()->get_RepeatFeatureAdaptor(), 
		      $self->db()->get_PredictionTranscriptAdaptor(),
		      $self->db()->get_ProteinAlignFeatureAdaptor(),
		      $self->db()->get_DnaAlignFeatureAdaptor());
  
  foreach my $adaptor ( @adaptor_list ) {
    $adaptor->remove_by_RawContig($contig);
  }

  # Delete DNA as long as we aren't using a remote DNA database.
  if ($self->db ne $self->db->dnadb) {
    $self->warn("Using a remote dna database - not deleting dna\n");
  } else {

    # Get the dna_id for the delete
    my $sth = $self->prepare("SELECT dna_id FROM contig where contig_id = ?");
    $sth->execute($contig->dbID);
    $self->throw("Failed to find any dna for dna_id '$contig->dbID'")
      unless $sth->rows;

    # Do the delete
    my $dna_id = $sth->fetchrow_array;
    $sth = $self->prepare("DELETE FROM dna WHERE dna_id = $dna_id");
    $sth->execute;
    my $flag_found=0;
    if($sth->rows){
      $flag_found=1;
    }
    # now try dnac - eval since table might not be present..
    $sth = $self->prepare("DELETE FROM dnac WHERE dna_id = $dna_id");
    eval{
      $sth->execute;
    };
    if($@){
    }else{
      if($sth->rows){
	$flag_found=1;
      }
    }
    print "got here\n";
    if($flag_found==0){
      $self->throw("Failed to delete dna for dna_id '$dna_id'")
    }
  }

  # Remove the contig
  my $sth = $self->prepare("DELETE FROM contig WHERE contig_id = ?");
  $sth->execute($contig->dbID);
  $self->throw("Failed to delete contigs for contig_id '$contig->dbID'")
    unless $sth->rows;

}



sub fetch_all_names{
  my ($self) = @_;

  my $sql = "select name from contig";
  my $sth = $self->prepare($sql);
  $sth->execute;
  my @names;
  while(my ($name) = $sth->fetchrow){
    push(@names, $name);
  }
  
  return \@names;
}

1;
