#
# bioperl module for database updates
#
# Cared for by Michele Clamp  <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::UpdateManager

=head1 SYNOPSIS

    my $fromlocator =  "Bio::EnsEMBL::TimDB::Obj";   # Timdb is different of course
    my $tolocator   =  "Bio::EnsEMBL::DBSQL::Obj/host=obi-wan:port=410000;dbname=pogtest;user=root;pass=''";
    my $arclocator  =  "Bio::EnsEMBL::DBArchive::Obj/host=obi-wan:port=410000;dbname=arctest;user=root;pass=''";

    my $manager     =   new Bio::EnsEMBL::Analysis::UpdateManager(-fromlocator => $fromlocator,
								  -tolocator   => $tolocator,
								  -arclocator  => $arclocator,
								  -fromtime    => 90000000,
								  -totime      => 90001000,
								  );

       $manager->chunksize(20);            # Transfer objects in chunks of 20
       $manager->nowrite(1);               # Don't write to recipient database (testing)
       $manager->verbose(1);               # Extra debugging info.

       $manager->update;                   # Start the updating process


=head1 DESCRIPTION

Manages the updating of one database to another.  Fires off the updates in chunks
by forking child processes to keep from running out of memory.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::UpdateManager;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object;

# Needed for TimDB
$SIG{INT}=sub { my $sig=shift;die "exited after SIG$sig";};

use Bio::EnsEMBL::TimDB::Obj;
use Bio::EnsEMBL::DBLoader;
use Bio::Root::Object;

# Inherits from the base bioperl object
@ISA = qw(Bio::Root::Object);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  
  my $make = $self->SUPER::_initialize;

  my ($fromlocator,$tolocator,$arclocator,$fromtime,$totime) = $self->_rearrange([qw(FROMLOCATOR
										     TOLOCATOR
										     ARCLOCATOR
										     FROMTIME
										     TOTIME
										     )],@args);

  $fromlocator || $self->throw("No donor database locator specified\n");
  $tolocator   || $self->throw("No recipient database locator specified\n");
  $arclocator  || $self->throw("No archive database locator specified\n");
  
  
  $self->fromlocator($fromlocator);
  $self->tolocator  ($tolocator);
  $self->arclocator ($arclocator);

  $self->fromtime   ($fromtime);
  $self->totime     ($totime);

  return $make; # success - we hope!
}


=head2 getchunk

  Title   : getchunk
  Usage   : $self->getchunk(@clonearray);
  Function: Gets a chunk of clones from an array
  Returns : int
  Args    : int

=cut

sub getchunk {
    my ($self,$current,@clones) = @_;

    my @chunk = splice(@clones,$current,$self->chunksize);

    return @chunk;
}

=head2 fromlocator

  Title   : fromlocator
  Usage   : $self->fromlocator($str);
  Function: Get/Set method for the donor database locator string
  Returns : string
  Args    : string

=cut


sub fromlocator {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_fromlocator} = $arg;
    }
    return $self->{_fromlocator};
}



=head2 tolocator

  Title   : tolocator
  Usage   : $self->tolocator($str);
  Function: Get/Set method for the recipient database locator string
  Returns : string
  Args    : string

=cut

sub tolocator {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_tolocator} = $arg;
    }
    
    return $self->{_tolocator};
}

=head2 arclocator

  Title   : arclocator
  Usage   : $self->arclocator($str);
  Function: Get/Set method for the archive database locator string
  Returns : string
  Args    : string

=cut

sub arclocator {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_arclocator} = $arg;
    }
    
    return $self->{_arclocator};
}


=head2 connect

  Title   : connect
  Usage   : my $db = $self->connect($locator);
  Function: Connects to a database using a locator string
  Returns : Bio::EnsEMBL::DB::ObjI
  Args    : String

=cut

sub connect {
    my ($self,$locator,@args) = @_;
    
    my $db;
    if ($locator eq "Bio::EnsEMBL::TimDB::Obj") {
	$db = "$locator"->new(\@args);
    } else {
	$db = new Bio::EnsEMBL::DBLoader($locator);
    }
    return $db;
}


=head2 chunksize

  Title   : chunksize
  Usage   : my $chunksize = $self->chunksize
  Function: Get/Set method for the number of clones to update in one chunk
  Returns : int
  Args    : int

=cut

sub chunksize {
    my ($self,$arg) = @_;

    my $default_chunk = 20;

    if (defined($arg)) {
	$self->{_chunksize} = $arg;
    } 
    return $self->{_chunksize} || 20;
}


=head2 retries

  Title   : retries
  Usage   : my $num = $self->retries
  Function: Get/Set method for the number of retries to do
  Returns : int
  Args    : int

=cut

sub retries {
    my ($self,$arg) = @_;
    
    my $default_tries = 1;

    if (defined($arg)) {
	$self->{_retries} = $arg;
    } elsif (!defined($self->{_retries})) {
	$self->{_retries} = $default_tries;
    }

    return $self->{_retries};
}

=head2 nowrite

  Title   : nowrite
  Usage   : $self->nowrite(1)
  Function: No writing is done to the recipient database
  Returns : Nothing
  Args    : 0,1

=cut

sub nowrite {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_nowrite} = $arg;
    }

    return $self->{_nowrite};
}

=head2 verbose

  Title   : verbose
    Usage   : $self->verbose(1);
  Function: Turns on/off verbose message printing
  Returns : Nothing
  Args    : 0,1

=cut

sub verbose {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_verbose} = $arg;
    }

    return $self->{_verbose};
}


=head2 fromtime 

  Title   : fromtime
  Usage   : my $time = $self->fromtime
  Function: Get/Set method for earliest modified time to update
  Returns : unix time string
  Args    : unix time string

=cut

sub fromtime {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_fromtime} = $arg;
    } 
    return $self->{_fromtime} || 1;
}

=head2 totime 

  Title   : totime
  Usage   : my $time = $self->totime
  Function: Get/Set method for latest modified time to update
  Returns : unix time string
  Args    : unix time string

=cut

sub totime {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_totime} = $arg;
    } 

    return $self->{_totime} || time;
}

=head2 get_updated_objects

  Title   : get_updated_objects
  Usage   : my @clones = $self->get_updated_objects
  Function: Fetches updated objects within the time slice in the object
  Returns : @string
  Args    : None

=cut

sub get_updated_objects {
    my ($self,$fromdb) = @_;

    $self->throw("Can't connect to donor database") unless $fromdb;

    my @clones = $fromdb->get_updated_Clone_id($self->fromtime,$self->totime);
    
    return @clones;
}


=head2 check_update_status

  Title   : check_update_status
  Usage   : $self->check_update_status
  Function: 
  Returns : int
  Args    : None

=cut

sub check_update_status {
    my ($self) = @_;

    my $tdb = $self->connect($self->tolocator);

    if ($tdb->current_update) {
	$self->throw("Update already running in recipient database. Can't start update");
    }

    my @clones;
    push(@clones,'AC000072');

    my $fdb = $self->connect($self->fromlocator,@clones);

    if ($self->fromlocator ne "Bio::EnsEMBL::TimDB::Obj") {

	if ($fdb->current_update) {
	    $self->throw("Update running in donor database.  Can't transfer objects");
	}
    } 

    
    return ($fdb,$tdb);
}
	    
=head2 update

  Title   : update
  Usage   : $self->update
  Function: Updates from the donor to the recipient database in chunks
  Returns : Nothing
  Args    : None

=cut
  
sub update {
    my ($self) = @_;


    my ($fdb,$tdb) = $self->check_update_status;
    my $id         = $tdb ->start_update unless $self->nowrite;
    my @clone_id   = $self->get_updated_objects($fdb);

    my $num_clones = scalar(@clone_id);
    my $current    = 0;

    while ($current < $num_clones) {
	print(STDERR "Starting new fork chunk $current/$num_clones\n");

	if (my $pid = fork) {
	    # Parent here
	    # Wait for the child to finish
	    my $status = waitpid $pid,0;
	    
	    print(STDERR "Exit status for child [$current] = $?\n");

	    $current += $self->chunksize;

	} elsif (defined($pid)) {
           $SIG{INT}=sub { my $sig=shift;die "exited after SIG$sig";};

	    # child here

	    my @clones = $self->getchunk($current,@clone_id);           print(STDERR  "In child. Transferring @clones\n");
	    my $fromdb = $self->connect ($self->fromlocator,@clones);   print(STDERR  "Connected to donor database\n");
	    my $todb   = $self->connect ($self->tolocator); 	        print(STDERR  "Connected to recipient database\n");
	    my $arcdb  = $self->connect ($self->arclocator); 	        print(STDERR  "Connected to archive database\n");

#	    eval {
		$self->transfer_chunk($fromdb,$todb,$arcdb,@clones);
#	    };

	    # must exit child. Big trouble otherwise
#	    exit($@);
	    exit(0);
	} else {
	    $self->throw("Couldn't fork a new process");
	}
    }

    if (!$self->nowrite) {
	my $todb = $self->connect($self->tolocator);  print(STDERR  "Connected to recipient database\n");
  	   $todb->replace_last_update($self->totime);
    }
}

=head2 transfer_chunk

  Title   : transfer_chunk
  Usage   : $self->transfer_chunk($fromdb,$todb,$arcdb,@clones);
  Function: Updates from the donor to the recipient database in chunks
  Returns : Nothing
  Args    : None.

=cut


sub transfer_chunk {
    my ($self,$fromdb,$todb,$arcdb,@clones) = @_;

    foreach my $id (@clones) {
	
	my $object=$fromdb->get_Clone($id);
	
	# Check if it is a clone object
	if ($object->isa("Bio::EnsEMBL::DB::CloneI")) {
	    $self->write_clone($todb,$arcdb,$object);
	}
	
	# These won't happen - updated clones only are returned from TimDB
	
	# Check if it is a gene
	elsif ($object->isa("Bio::EnsEMBL::Gene")) {
	    $self->write_gene($todb,$arcdb,$object);
	}
	
	# Check if it is an exon
	elsif ($object->isa("Bio::EnsEMBL::Exon")) {
	    $self->write_exon($todb,$object);
	}
    }
    
    
}

=head2 write_exon

  Title   : write_exon
  Usage   : $self->write_exon($db,$object);
  Function: Writes exon into the recipient database
  Returns : Nothing
  Args    : Bio::EnsEMBL::DB::ObjI,Bio::EnsEMBL::Exon

=cut


sub write_exon {
    my ($self,$db,$object) = @_;

    $self->verbose && print STDERR "Got exon with id ".$object->id."\n";

    my $rec_exon;

    #Check if it is already present in recipient
    eval {
	$rec_exon = $db->get_Exon($object->id);
    };


    if ( $rec_exon ) {

	$self->verbose && print STDERR "New Exon, writing it in the database\n";
	$self->nowrite || $db->write_Exon($object);

    } else {

	if ($object->version > $rec_exon->version) {

	    $self->verbose && print STDERR "Exon with new version, updating the database\n";
	    
	    unless ($self->nowrite) {
		$db->delete_Exon($object->id);
		$db->write_Exon ($object);
	    }

	} elsif ($rec_exon->version > $object->version) {
	    print STDERR "Something is seriously wrong, found a gene in the recipient database with version number higher than that of the donor database!!!\n";
	}  else {
	    $self->verbose && print STDERR "Exons with the same version, databases kept unchanged\n";
	}
    }
}


=head2 write_clone

  Title   : write_clone
  Usage   : $self->write_clone($db,$object);
  Function: Writes clone into the recipient database
  Returns : Nothing
  Args    : Bio::EnsEMBL::DB::ObjI,Bio::EnsEMBL::DB::CloneI

=cut


sub write_clone {
    my ($self,$db,$arcdb,$object) = @_;

    my $rec_clone;

    $self->verbose && print STDERR "Got clone with id ".$object->id."\n";
    
    eval {
	$rec_clone = $db->get_Clone($object->id);
    };

    if (! defined($rec_clone) ) { 
	$self->verbose &&  print STDERR "New Clone, writing it in the database\n";
	$db  ->write_Clone($object)  unless $self->nowrite;

	foreach my $gene ($object->get_all_Genes('evidence')) {
	    $self->verbose &&  print STDERR "Getting all genes via clone get_all_Genes method\n";
	    $self->write_gene($db,$arcdb,$gene,'1');
	}

    } else {
	print("Object 1 [$object] [$rec_clone]\n");
	if ($object->version > $rec_clone->version) {
	    $self->verbose && print STDERR "Clone with new version, updating the database, and deleting the old version\n";

	    unless ($self->nowrite) {
		$db->delete_Clone($rec_clone->id);
		$db->write_Clone ($object);
	    }

	    foreach my $gene ($object->get_all_Genes('evidence')) {

		if ($self->verbose) {
		    print STDERR "Getting all genes via clone get_all_Genes method\n";
		    print STDERR "Got gene ".$gene->id."\n";
		}

		$self->write_gene($db,$arcdb,$gene,'1');
	    }

	} elsif ($rec_clone->version > $object->version) {
	    print STDERR "ERROR: Something is seriously wrong, found a clone in the recipient database with version number higher than that of the donor database!!!\n";
	    exit;
	}  else {
	    $self->verbose && print STDERR "Clone versions equal, not modifying database\n";
	}
    }
}


=head2 write_gene

  Title   : write_gene
  Usage   : $self->write_gene($db,$object);
  Function: Writes gene into the recipient database
  Returns : Nothing
  Args    : Bio::EnsEMBL::DB::ObjI,Bio::EnsEMBL::Gene

=cut


sub write_gene {
    my ($self,$db,$arc_db,$don_gene,$clone_level) = @_;

    my $rec_gene;

    $self->verbose && print STDERR "Got gene with id ".$don_gene->id.", and version ".$don_gene->version."\n";    

    eval {
	$rec_gene = $db->get_Gene($don_gene->id,'evidence');
    };
    
    
    if ( ! defined($rec_gene) ) {
	$self->verbose && print STDERR "New Gene, writing it in the database\n";
	$self->nowrite || $db->write_Gene($don_gene);
    } else {
	if ($don_gene->version > $rec_gene->version) {
	    
	    $self->verbose && print STDERR "Gene with new version, updating the database, and archiving old version\n";

	    unless ($self->nowrite) {
		$db->archive_Gene($rec_gene,$arc_db);
		$db->write_Gene  ($don_gene);
	    }
	    
	}  elsif ($rec_gene->version > $don_gene->version) {
	    print STDERR "Something is seriously wrong, found a gene in the recipient database with version number higher than that of the donor database!!!\n";
	} else {
	    if ($clone_level) {
		$self->verbose && print STDERR "Genes with the same version, deleting recipient gene and writing one from donor without archiving\n";  

		unless ($self->nowrite) {
		    $db->delete_Gene($rec_gene->id);
		    $db->write_Gene ($don_gene);
		}
	    } else {
		$self->verbose && print STDERR "Genes with the same version, nothing needs to be done\n"; 
	    }
	}
    }
}   


