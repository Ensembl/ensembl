#
# EnsEMBL module for Bio::EnsEMBL::AceDB::Update_Obj
#
# Cared for by Simon Kay <simon@sanger.ac.uk>
#
# Copyright Simon Kay
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::AceDB::Update_Obj - AceDB database adapter class for EnsEMBL update system

=head1 SYNOPSIS

  use Bio::EnsEMBL::AceDB::Obj;
  use Bio::EnsEMBL::AceDB::Update_Obj;

  $db = new Bio::EnsEMBL::AceDB::Obj( -host => 'wormsr1', -port => 100100);
  my $update_obj = Bio::EnsEMBL::AceDB::Update_Obj->new($obj);

  # Get the last update time - offset
  $update_obj->get_last_update_offset();

=head1 DESCRIPTION

This is one of the objects contained in the Bio:EnsEMBL::AceDB package, dealing with
the update system, such identifying last update, getting updated objects, ghosts, etc.

The Obj object represents a database that is implemented somehow (you shouldn\'t
care much as long as you can get the object). 

=head1 CONTACT

Simon Kay: sjk@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are 
usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::AceDB::Update_Obj;

use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object and Bio::EnsEMBL::DB::Update_ObjI

use Bio::Root::Object;
use Bio::EnsEMBL::AceDB::Obj;
use Bio::EnsEMBL::DB::Update_ObjI;
use Bio::EnsEMBL::Ghost;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;

use DBI;

use Bio::EnsEMBL::DBSQL::DummyStatement;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::Update_ObjI);

# new() is inherited from Bio::Root::Object
# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,$db_obj) = @_;

  my $make = $self->SUPER::_initialize;
  
  $db_obj || $self->throw("Database Update_Obj must be passed a db obj!");
  $self->_db_obj($db_obj);

  return $make; 
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
    ($last eq "") && $self->warn ("No value stored for last_update in db_update table!
Setting it to zero!");
    $last=0;

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
 Function: Gets the current time from the point of view of the database which for an AceDB,
            is assumed to be the value of the time function minus the database
           offset time which is just 30 minutes (1800 sec).
 Example : get_now_offset
 Returns : UNIX TIME of now - offset_time
 Args    : none


=cut

sub get_now_offset{
    my ($self) = @_;

    return time - 1800;
}

=head2 replace_last_update
    
 Title   : replace_last_update(@$now_offset)
 Usage   : $obj->replace_last_update($now_offset)
 Function: Replaces the time in the last update field of the meta table with the now_offset 
            time of the recipient. This field does not exist for an AceDB so this method does nothing.
 Example : 
 Returns : nothing
 Args    : 

=cut

sub replace_last_update {
    my ($self, $now_offset) = @_;
}

=head2 current_update
    
 Title   : current_update
 Usage   : $obj->current_update
 Function: Checks whether the database is in the middle of an update. This can't be
            found for an AceDB so 0 is just returned
 Example : 
 Returns : 0
 Args    : 

=cut

sub current_update {
    my ($self) = @_;

    return 0;
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
 Usage   : $obj->get_updated_Clone_id()
 Function: Gets all the objects ids that have been updated. This can't be determined in an AceDB
            so all the object ids are returned
 Example : $obj->get_updated_Objects()
 Returns : database objects (clones and genes)
 Args    : 

=cut

sub get_updated_Clone_id {
   my ($self) = @_;
 
   return $self->_db_obj->get_all_Clone_id;
}


=head2 get_updated_Ghosts
    
 Title   : get_updated_Ghosts
 Usage   : $obj->get_updated_Ghosts ($recipient_last_update, $recipient_now_offset)
 Function: Ghosts are not stored in AceDB so undef is returned. 
 Example : $obj->get_updated_Ghosts (973036800,973090800)
 Returns : undef
 Args    : $recipient_last_update, $recipient_now_offset

=cut

sub get_updated_Ghosts{
    return undef;
}

=head2 get_Ghost
    
 Title   : get_Ghost
 Usage   : $obj->get_Ghost ($ghost_id,$ghost_version,$ghost_obj_type)
 Function: Ghosts are not stored in AceDB so undef is returned. 
 Example : $obj->get_Ghost ('test','1','transcript')
 Returns : undef
 Args    : ghost id, version and object type

=cut

sub get_Ghost{
    return undef;
}

=head2 write_Ghost
    
 Title   : write_Ghost
 Usage   : $obj->write_Ghost ($ghost)
 Function: Writes a ghost to the database. This is not implemented with AceDB so an error is thrown.  
 Example : $obj->write_Ghost ($ghost)
 Returns : 
 Args    : ghost object

=cut

sub write_Ghost{ 
    my ($self) = @_;   
    $self->throw("Writing ghosts is not implemented with AceDB");
}

=head2 archive_Gene
    
 Title   : archive_Gene
 Usage   : $obj->archive_gene($gene,$arcdb)
 Function: Deletes a gene and all its transcripts and exons and archives partial info.
            This is not implemented with AceDB so an error is thrown.
 Example : 
 Returns : nothing
 Args    : $gene, $arcdb (archive database object)


=cut

sub archive_Gene {
    my ($self) = @_;
   $self->throw("Archiving is not implemented with AceDB");
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
