
#
# BioPerl module for DBArchive::Obj
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBArchive::Obj - Object representing the EnsEMBL Archive DB

=head1 SYNOPSIS

    $db = new Bio::EnsEMBL::DBArchive::Obj( -user => 'root', -db => 'pog' , -host => 'caldy' , -driver => 'mysql' );

    $clone = $db->write_seq('3452');

    $contig = $db->get_seq('3452');

    $gene  = $db->get_seq_by_clone('X45667');

    

=head1 DESCRIPTION

This object represents an archive database that is implemented somehow
(you shouldn\'t care much as long as you can get the object). The
archive database holds a slice of data for older versions of proteins,
genes, and exons. It comprises three methods for writing and
retrieving sequences from the database. The purpose of this object is
to allow versioning in EnsEMBL, holding only the most recent of an
entry in the main DBSQL database, and storing here only the relevant
information of older versions.


=head1 CONTACT

Elia Stupka - EBI (elia@ebi.ac.uk)

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBArchive::Obj;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::RootI;
use DBI;
use Bio::Seq;

@ISA = qw(Bio::Root::RootI);


sub new {
  my($class,@args) = @_;

  my $self = {};
  bless $self,$class;



  #print "Got",join(',',@args),"\n";
  my ($db,$host,$driver,$user,$password,$debug,$readonly) = 
      $self->_rearrange([qw(DBNAME
			    HOST
			    DRIVER
			    USER
			    PASS
			    DEBUG
			    READONLY
			    )],@args);
  #print "Got $db as db and $user as user\n";

  $db || $self->throw("Database object must have a database name");
  $user || $self->throw("Database object must have a user");

  $self->_readonly($readonly);
  if($self->_readonly){
      $self->warn("Archive Database accessed in READONLY mode");
  }

  if( $debug ) {
      $self->_debug($debug);
  } else {
      $self->_debug(0);
  }
  
  if( ! $driver ) {
      $driver = 'mysql';
  }
  if( ! $host ) {
      $host = 'localhost';
  }
  my $dsn = "DBI:$driver:database=$db;host=$host";

  if( $debug && $debug > 10 ) {
      $self->_db_handle("dummy dbh handle in debug mode $debug");
  } else {
      
      my $dbh = DBI->connect("$dsn","$user","$password",{RaiseError => 1});
      $dbh || $self->throw("Could not connect to database $db user $user using [$dsn] as a locator");
      
      if( $self->_debug > 3 ) {
	  $self->warn("Using connection $dbh");
      }
     
      $self->_db_handle($dbh);
  }

  return $self; # success - we hope!
}

=head2 get_new_id_from_old_id

 Title   : get_new_id_from_old_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_new_id_from_old_id{
   my ($self,$type,$old_id) = @_;

   if( !defined $old_id ) {
       $self->throw("Must give type and old_id");
   }

   my $sth = $self->prepare("select old_id,new_id from deleted_id where old_id = '$old_id' and id_type = '$type'");
   my $res = $sth->execute();
   my ($old,$new) = $sth->fetchrow_array();
   
   if( !defined $new || $new eq "" ) {
       return "__DELETED__";
   } else {
       return $new;
   }

}


=head2 get_seq

 Title   : get_seq
 Usage   : $db->get_seq (id, version)
 Function: Gets a sequence object out of the Archive database
 Example : $db->get_seq (ENSP0000012,1.2)
 Returns : $seq object
 Args    : id, version


=cut

sub get_seq{
    my ($self,$seqid,$seqversion) = @_;
    
    $seqid || $self->throw("Attempting to get a sequence with no id");
    $seqversion || $self->throw("Attempting to get a sequence without a version number");
    
    #For now $seqtype is not passed to this method, because each type uses a different
    #pre-tag in the id, this might change later...
    #$seqtype || $self->throw("Attempting to get a sequence without a sequence type");
        
    # get the sequence object
    my $sth = $self->prepare("select id,version,sequence from sequence where (id = '$seqid' && version = '$seqversion')");
    my $res = $sth->execute();
    my @out = $self->_create_seq_obj($sth);
    return $out[0];
}

=head2 get_seq_by_id

 Title   : get_seq_by_id
 Usage   : $db->get_seq (id)
 Function: Gets a sequence object for each version for a given id out of the Archive database
 Example : $db->get_seq_by_id (ENST00000007087)
 Returns : array of $seq objects
 Args    : id


=cut

sub get_seq_by_id{
    my ($self,$seqid) = @_;
    
    $seqid || $self->throw("Attempting to get a sequence with no id");
    
    # get the sequence object
    my $sth = $self->prepare("select id,version,sequence from sequence where id = '$seqid'");
    my $res = $sth->execute();
    my @out = $self->_create_seq_obj($sth);
    return @out;
}

=head2 get_seq_by_clone_version

 Title   : get_seq_by_clone_version
 Usage   : $db->get_seq (clone_id, clone_version, seq_type)
 Function: Gets all the sequence objects for a given clone_id, clone_version and sequence 
           type out of the Archive database
 Example : $db->get_seq_by_clone ('AL021546','1','exon')
 Returns : array of $seq objects
 Args    : clone_id


=cut

sub get_seq_by_clone_version{
    my ($self,$clone_id, $clone_version, $seq_type) = @_;
    my $where_clause;

    $clone_id || $self->throw("Attempting to get a sequence with no clone id");
    $clone_version || $self->throw("Attempting to get a sequence with no clone version");
    $seq_type || $self->throw("Attempting to get a sequence with no sequence type");

    if ($clone_version eq 'all') {
	$where_clause = "where (clone_id = '$clone_id' && seq_type='$seq_type')";
    }
    else {
	$where_clause = "where (clone_id = '$clone_id' && clone_version = '$clone_version' && seq_type='$seq_type')";
    }

    # get the sequence objects
    my $sth = $self->prepare("select id,version,sequence from sequence $where_clause");
    my $res = $sth->execute();
    my @out = $self->_create_seq_obj($sth);
    return @out;
}

=head2 get_seq_by_gene_version

 Title   : get_seq_by_gene_version
 Usage   : $db->get_seq (gene_id, version)
 Function: If version is specified, gets all the sequence objects for a given gene_id 
           and gene_version out of the Archive database. If version is equal to 'all',
           then gets all the sequence objects for a given gene_id
 Example : $db->get_seq_by_gene ('AL021546','1'), or $db->get_seq_by_gene ('AL021546','all')
 Returns : array of $seq objects
 Args    : gene_id, version


=cut

sub get_seq_by_gene_version{
    my ($self,$gene_id, $gene_version, $seq_type) = @_;
    my $where_clause;
    
    $gene_id || $self->throw("Attempting to get a sequence with no gene id");
    $gene_version || $self->throw("Attempting to get a sequence with no gene version");
    $seq_type || $self->throw("Attempting to get a sequence with no sequence type"); 
    
    if ($gene_version eq 'all') {
	$where_clause = "where (gene_id = '$gene_id' && seq_type = '$seq_type')";
    }
    
    else {
	$where_clause = "where (gene_id = '$gene_id' && gene_version = '$gene_version' && seq_type = '$seq_type')";
    }

    # get the sequence object
    my $sth = $self->prepare("select id,version,sequence from sequence $where_clause");
    my $res = $sth->execute();
    my @out = $self->_create_seq_obj($sth);
    return @out;
}

=head2 write_deleted_id

 Title   : write_deleted_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub write_deleted_id{
    my ($self,$type,$old_id,$old_v,$new_id) = @_;
    
    $old_v || $self->throw("Must give type[got $type], old_id [got $old_id] and version [got $old_v] at least");
    $new_id ||= "NULL";
    my $sth = $self->prepare("insert into deleted_id (old_id,old_version,new_id,id_type) values('$old_id',$old_v,'$new_id','$type')");
    my $res = $sth->execute();
}


=head2 write_dead_geneid

 Title   : write_dead_geneid
 Usage   : $arcdb->write_dead_geneid($geneid)
 Function: write the id of a dead gene to the db 
 Example : $arcdb->write_dead_geneid('ENSG0000003456')
 Returns : nothing
 Args    : gene id string


=cut

sub write_dead_geneid{
    my ($self,$geneid) = @_;
    my $query="insert into dead_genes (id) values ('$geneid')";
    $self->_execute($query);
}

=head2 write_seq

 Title   : write_seq
 Usage   : $db->write_seq (seq,gene_id,gene_version,type,clone_id,clone_version)
 Function: Writes an entry in the archive database
 Example : $db->get_seq_by_id (ENSP0000012)
 Returns : array of $seq objects
 Args    : seq object, version, type, gene_id, gene_version, clone_id, clone_version
           Note that the id of the seq object contains the id of the
           db entry.


=cut

sub write_seq{
    my ($self, $seq, $version, $type, $gene_id,$gene_version,$cid,$cv) = @_;
   
    $seq || $self->throw("Attempting to write a sequence without a sequence object!");
    $type || $self->throw("Attempting to write a sequence without a sequence type!");
    $version || $self->throw("Attempting to write a sequence without a sequence version number!");
    $seq ||  $self->throw("Attempting to write a sequence without a sequence!");
    $gene_id || $self->throw("Attempting to write a sequence without a gene id!");
    $gene_version || $self->throw("Attempting to write a sequence without a gene version number!");
    $cid || $self->throw("Attempting to write a sequence without a clone id!");
    $cv || $self->throw("Attempting to write a sequence without a clone version number!");
    
    my $query="insert into sequence (id,version,seq_type,gene_id,gene_version,sequence,clone_id,clone_version) values ('".$seq->id()."','$version','$type','$gene_id','$gene_version','".$seq->seq."','".$cid."','".$cv."')";

    $self->_execute($query);
}

=head2 delete_seq

 Title   : delete_seq
 Usage   : $db->delete_seq (id, version)
 Function: Deletes a sequence entry from the Archive database
 Example : $db->delete_seq (ENSP0000012,1.2)
 Returns : 
 Args    : id, version


=cut

sub delete_seq{
    my ($self,$seqid,$seqversion) = @_;
    
    $seqid || $self->throw("Attempting to delete a sequence with no id");
    $seqversion || $self->throw("Attempting to delete a sequence without a version number");

    if ($self->_debug < 10) { 
	$self->throw ("Attempting to delete a sequence not in debug 10 mode!");
    }

    # delete the sequence entry
    my $query="delete from sequence where (id = '$seqid' && version = '$seqversion')";

    $self->_execute($query);
}

=head2 prepare

 Title   : prepare
 Usage   : $sth = $dbobj->prepare("select seq_start,seq_end from feature where analysis = \" \" ");
 Function: prepares a SQL statement on the DBI handle

           If the debug level is greater than 10, provides information into the
           DummyStatement object
 Example :
 Returns : A DBI statement handle object
 Args    : a SQL string


=cut

sub prepare{
   my ($self,$string) = @_;

   if( ! $string ) {
       $self->throw("Attempting to prepare an empty SQL query!");
   }

   if( $self->_debug > 10 ) {
       print STDERR "Prepared statement $string\n";
       my $st = Bio::EnsEMBL::DBSQL::DummyStatement->new();
       $st->_fileh(\*STDERR);
       $st->_statement($string);
       return $st;
   }

   # should we try to verify the string?

   return $self->_db_handle->prepare($string);
}

=head2 _debug

 Title   : _debug
 Usage   : $obj->_debug($newval)
 Function: 
 Example : 
 Returns : value of _debug
 Args    : newvalue (optional)


=cut

sub _debug{
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'_debug'} = $value;
    }
    return $self->{'_debug'};
    
}


=head2 _readonly

 Title   : _readonly
 Usage   : $obj->_readonly($newval)
 Function: 
 Example : 
 Returns : value of _readonly
 Args    : newvalue (optional)


=cut

sub _readonly{
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'_readonly'} = $value;
    }
    return $self->{'_readonly'};
}


=head2 _db_handle

 Title   : _db_handle
 Usage   : $obj->_db_handle($newval)
 Function: 
 Example : 
 Returns : value of _db_handle
 Args    : newvalue (optional)


=cut

sub _db_handle{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_db_handle'} = $value;
    }
    return $self->{'_db_handle'};

}

=head2 _create_seq_obj

 Title   : _create_seq_obj
 Usage   : $obj->_create_seq_obj ($sth)
 Function: 
 Example : 
 Returns : seq object
 Args    : $sth


=cut

sub _create_seq_obj{
    my ($self,$sth) = @_;

    my $seq = Bio::Seq->new;
    my @out;
    
    while( my $rowhash = $sth->fetchrow_hashref) {
	my $type;
	my $id = $rowhash->{'id'};
	$id .= ".";
	$id .= $rowhash->{'version'};
	if ($rowhash->{'seq_type'} eq 'protein') {
	    $type = 'amino';
	}
	else {
	    $type = 'dna';
	}
	$seq = Bio::Seq->new(
			     -seq=>$rowhash->{'sequence'},
			     -id=>$id,
			     -desc=>'Sequence from the EnsEMBL Archive database',
			     -type=>$type,
			     );
	push @out, $seq;
    }
    
    #Sort array of sequence objects by id
    @out = sort { my $aa = $a->id; $aa =~ s/^[^.]*.//g; my $bb = $b->id; $bb =~ s/^[^.]*.//g; return $aa <=> $bb } @out;
    return @out;
}


=head2

 Title   : get_new_stable_ids
 Usage   : @ids = $obj->get_new_stable_ids('exon',6);
 Function: 
 Example : 
 Returns : new id, ENS-style, as an array
 Args    : id type, one of 'exon','transcript','gene','translation'


=cut

my %stubhash = ( 'exon' => 'ENSE', 'transcript' => 'ENST','gene' => 'ENSG','translation' => 'ENSP' );

sub get_new_stable_ids {
    my ($self,$table,$number) = @_;

    if( !defined $number ) {
	$self->throw("Must call as table,number");
    }

    if( !defined $stubhash{$table} ) {
	$self->throw("Does not have $table as a valid stable id table");
    }
    
    my $stub= $stubhash{$table};
    $table .= "_stable";
    my @out;

    my $query="lock table $table write";

    $self->_execute($query);

    # wrap critical region in an eval so we can catch errors and release table

    eval {

	my $query = "select max(external_id) as id from $table where 1";
	
	my $sth   = $self->prepare($query);
	my $res   = $sth->execute;
	my $row   = $sth->fetchrow_hashref;
	my $id    = $row->{id};
	
	if (!defined($id) || $id eq "") {
	    $id = $stub . "00000000000";
	}
	
	if ($id =~ /\D+(\d+)$/) {
	    
	    my $newid  = $1;
	    my $i;
	    
	    foreach $i ( 1..$number ) {

		$newid++;
		
		
		if (length($newid) > 11) {
		    if ($newid =~ /^0/) {
			$newid =~ s/^0//;
		    } else {
			$self->throw("Can't truncate number string to generate new id [$newid]");
		    }
		}
		my $c = $stub . $newid;
		my $query = "insert into $table (internal_id,external_id,created) values (NULL,'$c',NOW())";

		$self->_execute($query);
			
		push(@out,$c);
	    }
	    
	    
	} else {
	    $self->throw("[$id] does not look like an object id (e.g. ENST00000019784)");
	}
    };

    my $error = undef;

    if( $@ ) {
	$error = $@;
    }

    $query = "unlock tables";
    $self->_execute($query);

    if( defined $error ) {
	$self->throw("Problem in making IDs. Unlocked tables. \n\n Error $@");
    }

    return @out;
    
}
=head2 _execute

 Title   : _execute
 Usage   :
 Function: Internal SQL prepare and execute function
           which does nothing if in readonly mode
 Example :
 Returns : 
 Args    :


=cut

sub _execute{
   my ($self,$query) = @_;

   if($self->_readonly){
       #print "READONLY: $query\n";
   }else{
       my $usth   = $self->prepare($query);
       $usth->execute;
   }
}



=head2 DESTROY

 Title   : DESTROY
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub DESTROY{
   my ($obj) = @_;

   if( $obj->{'_db_handle'} ) {
       $obj->{'_db_handle'}->disconnect;
       $obj->{'_db_handle'} = undef;
   }
}


