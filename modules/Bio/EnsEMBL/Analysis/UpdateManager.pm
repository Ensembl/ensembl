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
        my $args = "";
	$db = "$locator"->new(-nogene  => $self->nogene,
                              -freeze  => $self->freeze,
                              \@args);
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

=head2 usefile

  Title   : usefile
  Usage   : $self->usefile(clone.list);
  Function: Uses a list of clones, rather than get_update_Objects
  Returns : Nothing
  Args    : filename

=cut

sub usefile {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_usefile} = $arg;
    }

    return $self->{_usefile};
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

=head2 nogene 

  Title   : nogene
  Usage   : 
  Function: Flag used when connecting to TimDB which
            says whether to get genes or not 
  Returns : 0,1 
  Args    : 0,1

=cut

sub nogene {
    my ($self,$arg) = @_;

    if (defined($arg)) {
        $self->{_nogene} = $arg;
    }

    return $self->{_nogene} || 0;
}

=head2 freeze

  Title   : freeze 
  Usage   :
  Function: Flag used when connecting to TimDB which
            says which sequence freeze version to use
  Returns : int 
  Args    : int 

=cut

sub freeze {
    my ($self,$arg) = @_;

    if (defined($arg)) {
        $self->{_freeze} = $arg;
    }

    return $self->{_freeze} || 0;
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
    my @clones;

    $self->throw("Can't connect to donor database") unless $fromdb;
    
    if ($self->usefile) {
	print STDERR "Using ".$self->usefile." as a list of clones to update!\n";
	my $file = $self->usefile();
	open(IN,"<$file");
	while(<IN>){
	    if(/^(\S+)/){
		push(@clones,$1);
		print STDERR "Clone to update is $1\n";
	    }
	}
    } else {
	eval {
	    @clones = $fromdb->get_Update_Obj->get_updated_Clone_id($self->fromtime,$self->totime);
	};
	if ($@) {
	    print "Could not call get_updated_Clone_id from the donor database:\n$@!";
	}
    }
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
    
    my $t_update_obj=$tdb->get_Update_Obj();

    if ($t_update_obj->current_update) {
	$self->warn("Update already running in recipient database. Can't start update");
    }

    my @clones;
    #push(@clones,'AC000072');

    my $fdb = $self->connect($self->fromlocator,@clones);

    if ($self->fromlocator ne "Bio::EnsEMBL::TimDB::Obj") {

	my $f_update_obj = $fdb->get_Update_Obj;
	if ($f_update_obj->current_update) {
	    $self->warn("Update running in donor database, watch out!");
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

    my $tuobj= $tdb->get_Update_Obj;
    my $fuobj= $fdb->get_Update_Obj;

    my $id         = $tuobj ->start_update($self->fromtime,$self->totime) unless $self->nowrite;
    my @clone_id   = $self->get_updated_objects($fdb);
    
    # test
    # @clone_id      = ('AL132766');

    my $num_clones = scalar(@clone_id);
    my $current    = 0;

    $self->verbose && print STDERR "$num_clones to check for updates\n";

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

	   my @clones = $self->getchunk($current,@clone_id);           
           print(STDERR  "In child [$current]. Transferring @clones\n");
	   $self->verbose() && print STDERR "Connecting to donor database...\n"; 
	   my $fromdb = $self->connect ($self->fromlocator,@clones);   
           print(STDERR  "Connected to donor database\n");
	    $self->verbose() && print STDERR "Connecting to recipient database...\n"; 
	   my $todb   = $self->connect ($self->tolocator); 	        
           print(STDERR  "Connected to recipient database\n");

	   my $arcdb;

	   if ($self->arclocator eq "none") {
	       $arcdb = "none";
	   } else {
	       print STDERR "Connecting to archive database...\n" if $self->verbose;
	       $arcdb  = $self->connect ($self->arclocator); 	        
	       print STDERR  "Connected to archive database\n";
	   }
	   
	   eval {
	       $self->transfer_chunk($fromdb,$todb,$arcdb,@clones);
	   };

	   # must exit child. Big trouble otherwise
	   if ($@) {
	        warn($@);
                exit(1);
	    }
            exit(0);
           
       } else {
	   $self->throw("Couldn't fork a new process");
       }
    }
        
    
    if (!$self->nowrite) {
	my $todb = $self->connect($self->tolocator);  print(STDERR  "Replacing last update time with current time\n");
	my $tuobj = $todb->get_Update_Obj;
	$tuobj->replace_last_update($self->totime);
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
        my $object;
	eval {
	    $object=$fromdb->get_Clone($id);
	};
	
	if ($@) {
	    warn("Could not fetch clone: $@ Skipping clone\n");
	    next;
	}
	
        eval {
            # Check if it is a clone object
            if ($object->isa("Bio::EnsEMBL::DB::CloneI")) {
                $self->write_clone($todb,$arcdb,$object);
            }
        };
        if ($@) {
            warn("ERROR: problems in updating clone $id: $@ Deleting it from recipient database to preserve data integrity\n");
                
                # Check if it is a clone object
                if ($object->isa("Bio::EnsEMBL::DB::CloneI")) {
                    my $clone = new Bio::EnsEMBL::DBSQL::Clone( -id => $id, -dbobj => $todb);
                    $clone->delete();
                }
            
                # Check if it is a gene
                elsif ($object->isa("Bio::EnsEMBL::Gene")) {
                    $todb->gene_Obj->delete($id);
                }
            
                # Check if it is an exon
                elsif ($object->isa("Bio::EnsEMBL::Exon")) {
                   $self->gene_Obj->delete_Exon($id);
                }
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
    my %old_genes;
    my $rec_clone;

    $self->verbose && print STDERR "Got clone with id " . $object->id . "\n";
   
    # Trys calling $clone->fetch to see if the clone is present in recipient DB
    eval {        
        my $clone = new Bio::EnsEMBL::DBSQL::Clone( -id    => $object->id,
						-dbobj => $db );                                                
        $rec_clone = $clone->fetch();        
    };

    # If the clone isn't in the recipient DB ask the DB object to write it there
    if (! defined($rec_clone) ) { 
        $self->verbose &&  print STDERR "New Clone, writing it in the database\n";
        unless ($self->nowrite) {
            $db  ->write_Clone($object);
        }
       
        unless ($self->nogene) { 
        foreach my $gene ($object->get_all_Genes('evidence')) {
            $self->verbose && print STDERR "New Gene, writing it in the database\n";
	    $self->nowrite || $db->gene_Obj->write($gene);
        }
        }
    } 
    
    else {
        print("Object 1 [$object] [$rec_clone]\n");
        $self->throw("Can't write clone as version is not known") unless $object->version;
         
        if ($object->version > $rec_clone->version) {
            $self->verbose && print STDERR "Clone with new version, updating the database, and deleting the old version\n";
            if ($self->nogene == 0) {
                print STDERR "Getting all genes from donor clone\n";
                my @new_genes=$object->get_all_Genes('evidence');
                print STDERR "Getting all genes from recipient clone\n";    
                foreach my $old_gene ($rec_clone->get_all_Genes('evidence')) {
                    $old_genes{$old_gene->id} = $old_gene;
                }
                unless ($self->nowrite) {
                    $rec_clone->delete();
                    $db->write_Clone ($object);
                }
                foreach my $gene (@new_genes) {
                    $self->write_gene_hash($db,$arcdb,$gene,%old_genes,'1');
                }
            }
        }       
        elsif ($rec_clone->version > $object->version) {
            print STDERR "ERROR: Something is seriously wrong, found a clone in the recipient database with version number higher than that of the donor database!!!\n";
            exit;
        }  
        
        else {
            $self->verbose && print STDERR "Clone versions equal, not modifying database\n";
        }
    }
}


=head2 write_gene_hash

  Title   : write_gene_hash
  Usage   : $self->write_gene_hash($db,$object);
  Function: Writes gene into the recipient database
  Returns : Nothing
  Args    : Bio::EnsEMBL::DB::ObjI,Bio::EnsEMBL::Gene

=cut


sub write_gene_hash {
    my ($self, $db, $arc_db, $don_gene, %old_genes, $clone_level) = @_;

    my $rec_gene;
    
    $self->verbose && print STDERR "Got gene with id " . $don_gene->id . 
        ", and version " . ($don_gene->version ? $don_gene->version : "undefined") . "\n";    
    
    if(!$old_genes{$don_gene->id} ) {
	$self->verbose && print STDERR "New Gene, writing it in the database\n";
	$self->nowrite || $db->gene_Obj->write($don_gene);
    } else {
	if ($don_gene->version > $old_genes{$don_gene->id}->version) {
	    
	    unless ($self->nowrite) {
		if ($arc_db ne "none") {
		    $self->verbose && print STDERR "Gene with new version, updating the database, and archiving old version\n";
		    $db->archive_Gene($old_genes{$don_gene->id},$arc_db);
		}
		else {
		    $self->verbose && print STDERR "Gene with new version, updating the database, and deleting old version\n";
		    $db->delete_Gene($old_genes{$don_gene->id}->id);
		}
		$db->gene_Obj->write($don_gene);
	    }
	    
	}  elsif ($old_genes{$don_gene->id}->version > $don_gene->version) {
	    print STDERR "Something is seriously wrong, found a gene in the recipient database with version number higher than that of the donor database!!!\n";
	} else {
	    if ($clone_level) {
		$self->verbose && print STDERR "Genes with the same version, deleting recipient gene and writing one from donor without archiving\n";  

		unless ($self->nowrite) {
		    my $gene_obj=Bio::EnsEMBL::DBSQL::Gene_Obj->new($db);
		    $gene_obj->delete($old_genes{$don_gene->id}->id);
		    $gene_obj->write($don_gene);
		}
	    } else {
		$self->verbose && print STDERR "Genes with the same version, nothing needs to be done\n"; 
	    }
	}
    }
} 

