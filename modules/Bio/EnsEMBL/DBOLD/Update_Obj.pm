#
# EnsEMBL module for Bio::EnsEMBL::DBOLD::Update_Obj
#
# Cared for by Elia Stupka <elia@ebi.ac.uk>
#
# Copyright Elia Stupka
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBOLD::Update_Obj - MySQL database adapter class for EnsEMBL update system

=head1 SYNOPSIS

  use Bio::EnsEMBL::DBOLD::Obj;
  use Bio::EnsEMBL::DBOLD::Update_Obj;

  $db = new Bio::EnsEMBL::DBOLD::Obj( -user => 'root', -db => 'pog' , -host => 'caldy' , -driver => 'mysql' );
  my $update_obj=Bio::EnsEMBL::Update_Obj->new($obj);

  # Get the last update time - offset
  $update_obj->get_last_update_offset();

=head1 DESCRIPTION

This is one of the objects contained in Bio:EnsEMBL::DBOLD::Obj, dealing with
the update system, such identifying last update, getting updated objects, ghosts, etc.

The Obj object represents a database that is implemented somehow (you shouldn\'t
care much as long as you can get the object). 

=head1 CONTACT

Elia Stupka: elia@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are 
usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::DBOLD::Update_Obj;

use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::EnsEMBL::DBOLD::Obj;
use Bio::EnsEMBL::DB::Update_ObjI;
use Bio::EnsEMBL::Ghost;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;

use DBI;

use Bio::EnsEMBL::DBOLD::DummyStatement;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::Update_ObjI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub new {
  my($class,$db_obj) = @_;

  my $self = {};
  bless $self,$class;


  
  $db_obj || $self->throw("Database Gene object must be passed a db obj!");
  $self->_db_obj($db_obj);

  return $self; # success - we hope!
}

=head2 donor_locator
    
 Title   : get_donor_locator
 Usage   : $obj->get_donor_locator; 
 Function: Reads the meta table of the database to get the donor_database_locator
 Example : get_donor_locator
 Returns : locator string
 Args    : none


=cut

sub get_donor_locator {
    my ($self) = @_;
    
    my $sth     = $self->_db_obj->prepare("select donor_database_locator from meta");
    my $res     = $sth->execute();
    my $rowhash = $sth->fetchrow_hashref;
    my $donor   = $rowhash->{'donor_database_locator'};

    ($donor eq "") && $self->throw ("No value stored for database locator in meta table!");

    return $rowhash->{'donor_database_locator'};
}

=head2 get_last_update_offset

 Title   : get_last_update_offset
 Usage   : $obj->get_last_update_offset; 
 Function: Reads the meta table of the database to get the last_update time - offset time
 Example : get_last_update_offset
 Returns : UNIX TIME of last update - offset time
 Args    : none

=cut

sub get_last_update_offset{
    my ($self) = @_;
    
    #Get the last update time
    my $sth     = $self->_db_obj->prepare("select UNIX_TIMESTAMP(max(time_started)) from db_update where status = 'COMPLETE'");
    my $res     = $sth ->execute();
    my $rowhash = $sth->fetchrow_hashref;
    my $last    = $rowhash->{'UNIX_TIMESTAMP(max(time_started))'};
    if ($last eq "") { 
	$self->warn ("No value stored for last_update in db_update table! Setting it to 1!");
        $last=1801;
    }

    #Now get the offset time from the meta table, which is in time format
    #$sth     = $self->_db_obj->prepare("select UNIX_TIMESTAMP(offset_time) from meta");
    #$res     = $sth->execute();
    #$rowhash = $sth->fetchrow_hashref();
    #my $offset  = $rowhash->{'UNIX_TIMESTAMP(offset_time)'};

    #Hardcoded for the moment, because '00:30:00' does not covert to 1800 in MySQL
    my $offset = 1800;

    my $last_offset = $last - $offset;

    return $last_offset;
}

=head2 get_now_offset

 Title   : get_now_offset
 Usage   : $obj->get_now_minus_offset; 
 Function: Gets the current time from the point of view of the database, substracts the
           offset time found in the meta table and gives back unix time of now-offset
 Example : get_now_offset
 Returns : UNIX TIME of now - offset_time
 Args    : none


=cut

sub get_now_offset{
    my ($self) = @_;

    #Get the offset time from the meta table, which is in time format
    my $sth     = $self->_db_obj->prepare("select offset_time from meta");
    my $res     = $sth->execute();
    my $rowhash = $sth->fetchrow_hashref();
    my $offset = $rowhash->{'offset_time'};

    #Perform the subtraction in mysql
    $sth         = $self->_db_obj->prepare("select UNIX_TIMESTAMP(DATE_SUB(now(), INTERVAL '$offset' HOUR_SECOND))");
    $res         = $sth->execute();
    $rowhash     = $sth->fetchrow_hashref();

    my $now_offset = $rowhash->{"UNIX_TIMESTAMP(DATE_SUB(now(), INTERVAL '$offset' HOUR_SECOND))"};
     
    return $now_offset;
}

=head2 replace_last_update
    
 Title   : replace_last_update(@$now_offset)
 Usage   : $obj->replace_last_update($now_offset)
 Function: Replaces the time in the last update field of the meta table with the now_offset time of the recipient
 Example : 
 Returns : nothing
 Args    : 

=cut

sub replace_last_update {
    my ($self, $now_offset) = @_;
    
    $now_offset || $self->throw("Trying to replace last update without a now-offset time\n");
    
    my $last_offset= $self->get_last_update_offset;
    
    my $sth = $self->_db_obj->prepare("select FROM_UNIXTIME(".$now_offset.")");
    $sth->execute();
    my $rowhash = $sth->fetchrow_arrayref();
    $now_offset = $rowhash->[0];
    
    $sth = $self->_db_obj->prepare("select FROM_UNIXTIME(".$last_offset.")");
    $sth->execute();
    $rowhash = $sth->fetchrow_arrayref();
    $last_offset = $rowhash->[0];
    
    $sth = $self->_db_obj->prepare("update db_update set time_finished=now(),status='COMPLETE' where status='STARTED'");
    $sth->execute;
}

=head2 current_update
    
 Title   : current_update
 Usage   : $obj->current_update
 Function: Checks whether the database is in the middle of an update
 Example : 
 Returns : 0,1
 Args    : 

=cut

sub current_update {
    my ($self) = @_;

    my $sth = $self->_db_obj->prepare("select id from db_update where status = 'STARTED'");
    my $res = $sth ->execute;

    my $id = 0;
    my $rowhash = $sth->fetchrow_hashref;
       $id      = $rowhash->{'id'};

    return $id;
}

=head2 start_update
    
 Title   : start_update
 Usage   : my $id = $obj->start_update
 Function: Enters a new updating process in the db_update table
 Example : 
 Returns : int
 Args    : 

=cut

sub start_update {
    my ($self,$start,$end) = @_;

    my $cid = $self->current_update;
    
    $self->throw("No start time defined") unless defined($start);
    $self->throw("No end time defined")   unless defined($end);

    $self->throw("Update already in progresss id [$cid]") unless ! $cid;
    
    
    my $sth = $self->_db_obj->prepare("insert into db_update" . 
			     "(id,time_started,status,modified_start,modified_end) " . 
			     " values(NULL,now(),'STARTED',FROM_UNIXTIME($start),FROM_UNIXTIME($end))");
    $sth->execute;
    
    $sth = $self->_db_obj->prepare("select last_insert_id()");
    $sth->execute;
    my $rowhash = $sth->fetchrow_hashref;
    my $id      = $rowhash->{'last_insert_id()'};
    
    return $id;
}

=head2 finish_update
    
 Title   : finish_update
 Usage   : my $id = $obj->finish_update
 Function: Completes the current update process
 Example : 
 Returns : nothing
 Args    : None

=cut

sub finish_update {
    my ($self) = @_;

    my $cid = $self->current_update;

    $self->throw("No current updating process. Can't finish") unless $cid;

    my $sth = $self->_db_obj->prepare("select id,time_started,modified_end,modified_start from db_update where id = $cid");
    my $res = $sth->execute;

    my $rowhash = $sth->fetchrow_hashref;

    my $id              = $rowhash->{'id'};
    my $time_started    = $rowhash->{'time_started'};
    my $modified_end    = $rowhash->{'modified_end'};
    my $modified_start  = $rowhash->{'modified_start'};

    # Should get stats on each table here as well.

    $sth = $self->_db_obj->prepare("replace into db_update(id,time_started,time_finished,modified_start,modified_end,status)".
			     "values($id,\'$time_started\',now(),\'$modified_start\',\'$modified_end\','COMPLETE'");
    $sth->execute;

}

=head2 get_updated_Clone_id
    
 Title   : get_updated_Clone_id
 Usage   : $obj->get_updated_Clone_id ($recipient_last_update, $recipient_now)
 Function: Gets all the objects that have been updated (i.e.change in 
 Example : $obj->get_updated_Objects (973036800,973090800)
 Returns : database objects (clones and genes)
 Args    : $recipient_last_update, $recipient_now

=cut

sub get_updated_Clone_id {
    my ($self, $last_offset, $now_offset) = @_;
    
    $last_offset || $self->throw("Attempting to get updated objects without the recipient db last update time");
    $now_offset  || $self->throw("Attempting to get updated objects without the recipient db current time");

    ($last_offset>$now_offset) && $self->throw("Last update more recent than now-offset time, serious trouble");

    my $sth      = $self->_db_obj->prepare("select FROM_UNIXTIME(".$last_offset.")");
    my $res      = $sth->execute();
    my $rowhash  = $sth->fetchrow_arrayref();
    $last_offset = $rowhash->[0];               print STDERR "Last= $last_offset\n";

    $sth        = $self->_db_obj->prepare("select FROM_UNIXTIME(".$now_offset.")");
    $res        = $sth->execute();
    $rowhash    = $sth->fetchrow_arrayref();
    $now_offset = $rowhash->[0];                print STDERR "now: $now_offset\n";

    $sth = $self->_db_obj->prepare("select id from clone where stored > '".$last_offset." - 00:30:00' and stored <= '".$now_offset."'");
    $res = $sth->execute;

    my @clones;

    while( my $rowhash = $sth->fetchrow_hashref) {
	push(@clones,$rowhash->{'id'});
    }

    return @clones;
}

=head2 get_updated_Objects
    
 Title   : get_updated_Objects
 Usage   : $obj->get_updated_Objects ($recipient_last_update, $recipient_now)
 Function: Gets all the objects that have been updated (i.e.change in 
	   version number) between the current time - offset time given by
           the recipient database and the last update time stored in its meta table 
 Example : $obj->get_updated_Objects (973036800,973090800)
 Returns : database objects (clones and genes)
 Args    : $recipient_last_update, $recipient_now

=cut

sub get_updated_Objects{
    my ($self, $last_offset, $now_offset) = @_;
    
    $last_offset || $self->throw("Attempting to get updated objects without the recipient db last update time");
    $now_offset  || $self->throw("Attempting to get updated objects without the recipient db current time");

    ($last_offset>$now_offset) && $self->throw("Last update more recent than now-offset time, serious trouble");

    my $sth      = $self->_db_obj->prepare("select FROM_UNIXTIME(".$last_offset.")");
    my $res      = $sth->execute();
    my $rowhash  = $sth->fetchrow_arrayref();
    $last_offset = $rowhash->[0];       print STDERR "Last: $last_offset\n";

    $sth        = $self->_db_obj->prepare("select FROM_UNIXTIME(".$now_offset.")");
    $res        = $sth->execute();
    $rowhash    = $sth->fetchrow_arrayref();
    $now_offset = $rowhash->[0];        print STDERR "now: $now_offset\n";

    $sth = $self->_db_obj->prepare("select id from clone where stored > '".$last_offset." - 00:30:00' and stored <= '".$now_offset."'");
    $res = $sth->execute;

    my @out;
    my @clones;

    while( my $rowhash = $sth->fetchrow_hashref) {
	push(@clones,$rowhash->{'id'});
    }
    
    foreach my $cloneid (@clones) {
	my $clone = new Bio::EnsEMBL::DBOLD::Clone( -id    => $cloneid,
						    -dbobj => $self->_db_obj );
   
	$clone->fetch();
	push @out, $clone;
    }	
    
    $sth = $self->_db_obj->prepare("select id from gene where stored > '".$last_offset."' and stored <= '".$now_offset."'");
    $sth->execute;

    my @genes;

    while( $rowhash = $sth->fetchrow_hashref) {
	push(@genes,$rowhash->{'id'});
    }
    
    #Get all gene objects for the ids contained in @clones, and push them in @out
    foreach my $geneid (@genes) {
	my $gene_obj=Bio::EnsEMBL::DBOLD::Gene_Obj->new($self->_db_obj);
	my $gene = $gene_obj->get($geneid);
	push @out, $gene;
    }	
    return @out;
}

=head2 get_updated_Ghosts
    
 Title   : get_updated_Ghosts
 Usage   : $obj->get_updated_Ghosts ($recipient_last_update, $recipient_now_offset)
 Function: Gets all the ghosts for objects that have been deleted (i.e.permanently from 
	   the donor db) between the current time - offset time given by
           the recipient database and the last update time stored in its meta table 
 Example : $obj->get_updated_Ghosts (973036800,973090800)
 Returns : ghost objects
 Args    : $recipient_last_update, $recipient_now_offset

=cut

sub get_updated_Ghosts{
    my ($self, $last_offset, $now_offset) = @_;
    my @out;
    
    $last_offset || $self->throw("Attempting to get updated objects without the recipient db last update time");
    $now_offset  || $self->throw("Attempting to get updated objects without the recipient db current time");

    ($last_offset>$now_offset) && $self->throw("Last update more recent than now-offset time, serious trouble");
    
    my $sth      = $self->_db_obj->prepare("select FROM_UNIXTIME(".$last_offset.")");
    my $res      = $sth->execute();
    my $rowhash  = $sth->fetchrow_arrayref();
    $last_offset = $rowhash->[0];
    
    $sth        = $self->_db_obj->prepare("select FROM_UNIXTIME(".$now_offset.")");
    $res        = $sth->execute();
    $rowhash    = $sth->fetchrow_arrayref();
    $now_offset = $rowhash->[0];
    
    #Get all ghosts that have deleted times between last and now-offset
    $sth = $self->_db_obj->prepare("select id,version,obj_type,deleted,stored from ghost where stored > '".$last_offset."' and stored <= '".$now_offset."'");
    $res = $sth->execute();

    while(my $rowhash = $sth->fetchrow_hashref()) {
	my $ghost = Bio::EnsEMBL::Ghost->new();

	$ghost->id      ($rowhash->{'id'});
	$ghost->version ($rowhash->{'version'});
	$ghost->obj_type($rowhash->{'obj_type'});
	$ghost->deleted ($rowhash->{'deleted'});
	$ghost->_stored ($rowhash->{'stored'});

	push @out, $ghost;
    }

    return @out;
}

=head2 get_Ghost
    
 Title   : get_Ghost
 Usage   : $obj->get_Ghost ($ghost_id,$ghost_version,$ghost_obj_type)
 Function: Gets a ghost by id, version,obj_type  
 Example : $obj->get_Ghost ('test','1','transcript')
 Returns : ghost objects
 Args    : ghost id, version and object type

=cut

sub get_Ghost{
    my ($self, $g_id, $g_obj_type) = @_;
    my @out;
    
    $g_id                        || $self->throw("Attempting to get a ghost object without an id");
    $g_obj_type                  || $self->throw("Attempting to get a ghost object without an object type");

    my $sth     = $self->_db_obj->prepare("select id,version,obj_type,deleted,stored from ghost where id='" . $g_id . "' and obj_type = '" . $g_obj_type . "'");
    my $res     = $sth->execute();
    my $rv      = $sth->rows     || $self->throw("Ghost not found in database!");
    my $rowhash = $sth->fetchrow_hashref();

    my $ghost   = Bio::EnsEMBL::Ghost->new();

    $ghost->id      ($rowhash->{'id'});
    $ghost->version ($rowhash->{'version'});
    $ghost->obj_type($rowhash->{'obj_type'});
    $ghost->deleted ($rowhash->{'deleted'});
    $ghost->_stored ($rowhash->{'stored'});

    return $ghost;
}

=head2 write_Ghost
    
 Title   : write_Ghost
 Usage   : $obj->write_Ghost ($ghost)
 Function: Writes a ghost to the database  
 Example : $obj->write_Ghost ($ghost)
 Returns : 
 Args    : ghost object

=cut

sub write_Ghost{
    my ($self, $ghost) = @_;
    
    $ghost                             || $self->throw("Attempting to write a ghost without a ghost object");
    $ghost->isa("Bio::EnsEMBL::Ghost") || $self->throw("$ghost is not an EnsEMBL ghost - not dumping!");
    
    my $sth = $self->_db_obj->prepare("insert into ghost (id, version, obj_type,deleted,stored) values('".
			     $ghost->id          . "','" . 
			     $ghost->version     . "','" .
			     $ghost->obj_type    . "','" .
			     $ghost->deleted     . "',now())");

    $sth->execute();
    return 1;
}

=head2 archive_Gene
    
 Title   : archive_Gene
 Usage   : $obj->archive_gene($gene,$arcdb)
 Function: Deletes a gene and all its transcripts and exons, 
           and archives partial info in the archive db passed on.
 Example : 
 Returns : nothing
 Args    : $gene, $arcdb (archive database object)


=cut

sub archive_Gene {
   my ($self,$gene,$arc_db) = @_;

   my $sth;
   my $res;

   foreach my $transcript ($gene->each_Transcript) {
       
       my $seq = $transcript->dna_seq;
          $seq->id($transcript->id);
       
       #Temporary, since versions not stored yet...
       !$transcript->version && $transcript->version(1);
       !$gene      ->version && $gene      ->version(1);
       
       # Store transcript and protein translation in the archive database
       $arc_db->write_seq($seq, $transcript->version, 'transcript', $gene->id, $gene->version);
       
       $seq = $transcript->translate;

       $arc_db->write_seq($seq, $transcript->version, 'protein', $gene->id, $gene->version);

       #Delete transcript rows
       $sth = $self->_db_obj->prepare("delete from transcript where id = '".$transcript->id."'");
       $res = $sth->execute;
       
       foreach my $exon ($gene->each_unique_Exon) {
	   #Get out info needed to write into archive db
	   $seq = $exon->primary_seq;
	   $seq->id($exon->id);
	   
	   #Temporary, since versions not stored yet...
	   !$exon->version && $exon->version(1);
	   
	   #Write into archive db
	   $arc_db->write_seq($seq, $exon->version, 'exon', $gene->id, $gene->version);
       }

       foreach my $exon ($transcript->each_Exon) {
	   
	   #Delete exon_transcript and exon rows
	   $sth = $self->_db_obj->prepare("delete from exon_transcript where transcript = '" . $transcript->id . "'");
	   $res = $sth->execute;

	   $sth = $self->_db_obj->prepare("delete from exon where id = '" . $exon->id . "'");
	   $res = $sth->execute;
       }
   }
   
   # delete gene rows
   $sth = $self->_db_obj->prepare("delete from gene where id = '".$gene->id."'");
   $res = $sth->execute;
}   

=head2 _db_obj

 Title   : _db_obj
 Usage   : $obj->_db_obj($newval)
 Function: 
 Example : 
 Returns : value of _db_obj
 Args    : newvalue (optional)


=cut

sub _db_obj{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_db_obj'} = $value;
    }
    return $self->{'_db_obj'};

}
