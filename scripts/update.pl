#!usr/local/bin/perl
=head1 NAME

Update

=head1 SYNOPSIS
 
  update.pl

=head1 DESCRIPTION

This script updates a recipient database by checking its donor database

=head1 OPTIONS

    -dbtype    Database tpye (only used for TimDB)

    -host      host name for database (gets put as host= in locator)

    -port      For RDBs, what port to connect to (port= in locator)

    -dbname    For RDBs, what name to connect to (dbname= in locator)

    -dbuser    For RDBs, what username to connect as (dbuser= in locator)

    -dbpass    For RDBs, what password to use (dbpass= in locator)

    -module    Module name to load (Defaults to Bio::EnsEMBL::DBSQL::Obj)

    -arcpass   password for the archive database

    -help      Displays script documentation with PERLDOC

=cut

use strict;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Ghost;
use Bio::SeqIO;
use Getopt::Long;
use vars qw(@ISA);

@ISA = qw(Bio::Root::Object);

my $dbtype = 'rdb';
my $host   = 'localhost';
my $port   = '410000';
my $dbname = 'eliatest2';
my $dbuser = 'root';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::Obj';
my $help;
my $arcpass = undef;

&GetOptions( 
	     'dbtype:s'   => \$dbtype,
	     'host:s'     => \$host,
	     'port:n'     => \$port,
	     'dbname:s'   => \$dbname,
	     'dbuser:s'   => \$dbuser,
	     'dbpass:s'   => \$dbpass,
	     'module:s'   => \$module,
	     'arcpass:s'  => \$arcpass,
	     'h|help'     => \$help
	     );


if ($help) {
    exec('perldoc', $0);
}

my $arcname = 'archive';
my $arctype = 'rdb';
my $archost = 'localhost';
my $arcport = '410000';
my $arcuser = 'root';
my $arcmodule = 'Bio::EnsEMBL::DBArchive::Obj';
$arcpass || die "You forgot to give the password for the archive database! Sorry, no access!\n";

my $arclocator = "$arcmodule/host=$archost;port=$arcport;dbname=$arcname;user=$arcuser;pass=$arcpass";
my $arc_db = Bio::EnsEMBL::DBLoader->new($arclocator);
    
print STDERR "\nConnecting to recipient...\n";
my $recipient_locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $rec_db =  Bio::EnsEMBL::DBLoader->new($recipient_locator);

print STDERR "\nConnecting to donor...\n";
my $donor_locator = $rec_db->get_donor_locator;
my $don_db =  Bio::EnsEMBL::DBLoader->new($donor_locator);

my $last = $rec_db->get_last_update;
print STDERR "\nLast update: $last\n";

my $now = $rec_db->get_now_offset;
print STDERR "Time now at recipient - offset update time: $now\n";

if ($last > $now) {
    print STDERR "Time of last update more recent than now-offset, exiting!\n";
    exit;
}
print "\nTransferring objects updated between last update and now-offset from donor to recipient...\n";

my @object_array = $don_db->get_updated_Objects($last, $now);

#Get updated and new objects (clones and genes)
foreach my $object (@object_array) {
    my $type;
    #Check if it is a gene
    if ($object->isa("Bio::EnsEMBL::Gene")) {
	print STDERR "Got a gene with id ".$object->id."\n";
	$type = "Gene";
	my $rec_gene;
	#Check if it is already present in recipient
	eval {
	    $rec_gene = $rec_db->get_Gene($object->id);
	};
	
	#If not present in recipient, write it in
	if ( $@ ) {
	    print STDERR "New Gene, writing it in the database\n";
	    $rec_db->write_Gene($object);
	}
	#If present in recipient, check donor and recipient version
	else {
	    #If donor gene version greater than recipient gene version, update
	    if ($object->version > $rec_gene->version) {
		print "Gene with new version, updating the database, and archiving old version\n";
		$rec_db->archive_Gene($rec_gene,$arc_db);
		$rec_db->write_Gene($object);
	    }
	    
	    #If donor gene version is less than the recipient gene version, error 
	    #NOTE: Better catching to be implemented
	    elsif ($rec_gene->version > $object->version) {
		print STDERR "Something is seriously wrong, found a gene in the recipient database with version number higher than that of the donor database!!!\n";
	    }
	    #If versions equal, nothing needs to be done
	    else {
		print STDERR "Genes with the same version, databases kept unchanged\n";
	    }
	}
    }
    elsif ($object->isa("Bio::EnsEMBL::DB::CloneI")) {
	$type = "Clone";
	print "Got a Clone with id ".$object->id."\n";
	
	my $rec_clone;
	#Check if it is already present in recipient
	eval {
	    $rec_clone = $rec_db->get_Clone($object->id);
	};
	
        #If not present in recipient, write it in
	if ( $@ ) { 
	    print STDERR "New Clone, writing it in the database\n";
	    $rec_db->write_Clone($object);
	}
	#If present in recipient, check donor and recipient version
	else {
	    #If donor clone version greater than recipient clone version, replace
	    #NOTE: We permanently delete clones, without archiving them!
	    if ($object->version > $rec_clone->version) {
		print "Clone with new version, updating the database, and deleting the old version\n";
		$rec_db->delete_Clone($rec_clone->id);
		$rec_db->write_Clone($object);
	    }
	    
	    #If donor Clone version is less than the recipient Clone version, error 
	    #NOTE: Better catching to be implemented
	    elsif ($rec_clone->version > $object->version) {
		print STDERR "ERROR: Something is seriously wrong, found a clone in the recipient database with version number higher than that of the donor database!!!\n";
	    }
	  	  
	    #If versions equal, nothing needs to be done
	    else {
		print STDERR "Clone versions equal, database left unchanged\n";
	    }
	}
    }

    else {
	print STDERR "ERROR: Got an object of unkown type, exiting\n";
    }
}

my @object_array = $don_db->get_Ghosts_by_deleted($last, $now);
foreach my $ghost (@object_array) {
    print STDERR "Got ghost ".$ghost->id."\n";
    
    #Check if the ghost is already present in the recipient database
    eval {
	my $try=$rec_db->get_Ghost($ghost->id,$ghost->version,$ghost->obj_type);
    };
   
    #If not present in recipient, archive objects, and write the ghost
    if ( $@ ) {
	if ($ghost->obj_type eq 'clone') {
	    $rec_db->delete_Clone($ghost->id);
	}
	elsif ($ghost->obj_type eq 'gene') {
	    my $gene = $rec_db->get_Gene($ghost->id);
	    $rec_db->archive_Gene($gene,$arc_db);
	}
	#Finally write ghost in recipient
	$rec_db->write_Ghost($ghost);
    }

    #Nothing needs to be done if the ghost is already present
}
$rec_db->replace_last_update($now);


