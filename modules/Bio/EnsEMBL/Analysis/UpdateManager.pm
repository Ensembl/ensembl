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
    my $manager     =   new Bio::EnsEMBL::Analysis::UpdateManager(-fromlocator => $fromlocator,
								  -tolocator   => $tolocator,
								  -fromtime    => 90000000,
								  -totime      => 90001000,
								  );
       $manager->chunksize(20);
       $manager->update;


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

$| = 1;

# Needed for TimDB
$SIG{INT}=sub {my $sig=shift;die "exited after SIG$sig";};


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

  my ($fromlocator,$tolocator,$fromtime,$totime) = $self->_rearrange([qw(FROMLOCATOR
									 TOLOCATOR
									 FROMTIME
									 TOTIME
									 )],@args);

  $fromlocator || $self->throw("No donor database locator specified\n");
  $tolocator   || $self->throw("No recipient database locator specified\n");
  
  
  $self->fromlocator($fromlocator);
  $self->tolocator  ($tolocator);
  $self->fromtime   ($fromtime);
  $self->totime     ($totime);

  return $make; # success - we hope!
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


=head2 disconnect

  Title   : disconnect
  Usage   : $self->disconnect($db)
  Function: Disconnects from a database
  Returns : nothing
  Args    : Bio::EnsEMBL::DB::ObjI

=cut

sub disconnect {
    my ($self,$db) = @_;

    # This not good really
    $db->DESTROY;

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
    } elsif (!defined($self->{_chunksize})) {
	$self->{_chunksize} = $default_chunk;
    }

    return $self->{_chunksize};
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
    if (!defined($self->{_fromtime})) {
	$self->{_fromtime} = 1;
    }
    return $self->{_fromtime};
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

    if (!defined($self->{_totime})) {
	$self->{_totime} = time;
    }

    return $self->{_totime};
}

=head2 get_updated_objects

  Title   : get_updated_objects
  Usage   : my @clones = $self->get_updated_objects
  Function: Fetches updated objects within the time slice in the object
  Returns : @string
  Args    : None

=cut

sub get_updated_objects {
    my ($self) = @_;

    my $fromdb = $self->connect($self->fromlocator);
    $self->throw("Can't connect to donor database") unless $fromdb;

    my @clones = $fromdb->get_updated_Clone_id($self->fromtime,$self->totime);
    
    $self->disconnect($fromdb);
    return @clones;
}

=head2 getchunk

  Title   : getchunk
  Usage   : $self->getchunk($current,@clones);
  Function: Gets the next slice of clone ids from the whole array
  Returns : Nothing
  Args    : None

=cut

sub getchunk {
    my ($self,$current,@clones) = @_;

    my @chunk = splice(@clones,0,$self->chunksize);

    return @chunk;
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

    my @clone_id   = $self->get_updated_objects;
    my $num_clones = scalar(@clone_id);
    my $current    = 0;
    
    while ($current < $num_clones) {
	print(STDERR "Starting new fork chunk $current\n");

	if (my $pid = fork) {
	    # Parent here
	    # Wait for the child to finish
	    my $status = waitpid $pid,0;
	    
	    if ($status) {
		$current += $self->chunksize;
	    } else {
		print(STDERR "Exit status for child [$current] = $status\n");
	    }
	} elsif (defined($pid)) {
	    # child here

	    my @clones = $self->getchunk($current,@clone_id);           print(STDERR "In child. Transferring @clones\n");
	    my $fromdb = $self->connect ($self->fromlocator,@clones);   print(STDERR "Connected to donor database\n");
	    my $todb   = $self->connect ($self->tolocator); 	        print(STDERR "Connected to recipient database\n");

	    $self->transfer_chunk($fromdb,$todb,@clones);
	    
	    # must exit child. Big trouble otherwise
	    exit(0);
	}
    }
}

=head2 transfer_chunk

  Title   : transfer_chunk
  Usage   : $self->transfer_chunk($fromdb,$todb,@clones);
  Function: Updates from the donor to the recipient database in chunks
  Returns : Nothing
  Args    : None

=cut


sub transfer_chunk {
    my ($self,$fromdb,$todb,@clones) = @_;

    eval {
	foreach my $id (@clones) {
	    
	    my $object=$fromdb->get_Clone($id);
	    my $type;
	    
	    # Check if it is a clone object
	    if ($object->isa("Bio::EnsEMBL::DB::CloneI")) {
		$type = "Clone";
		print STDERR "Got clone with id ".$object->id."\n";
	    }
	    
	    # Check if it is a gene
	    elsif ($object->isa("Bio::EnsEMBL::Gene")) {
		print STDERR "Gene level: got gene with id ".$object->id.
		    ", and version ".$object->version."\n";
	    }
	    
	    # Check if it is an exon
	    elsif ($object->isa("Bio::EnsEMBL::Exon")) {
		print STDERR "Got exon with id ".$object->id."\n";
	    }
	}
    };
    
    if ($@) {
	warn("ERROR: clone(s) not updated @clones\n");
    }
    
}


