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

    -verbose   Gets all the print STDERR for testing purposes

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
my $verbose;

&GetOptions( 
	     'dbtype:s'   => \$dbtype,
	     'host:s'     => \$host,
	     'port:n'     => \$port,
	     'dbname:s'   => \$dbname,
	     'dbuser:s'   => \$dbuser,
	     'dbpass:s'   => \$dbpass,
	     'module:s'   => \$module,
	     'arcpass:s'  => \$arcpass,
	     'h|help'     => \$help,
	     'v|verbose'  => \$verbose
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

$verbose && print "\nConnecting to archive database...\n";
my $arclocator = "$arcmodule/host=$archost;port=$arcport;dbname=$arcname;user=$arcuser;pass=$arcpass";
my $arc_db = Bio::EnsEMBL::DBLoader->new($arclocator);
    
$verbose && print "\nConnecting to recipient database...\n";
my $recipient_locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $rec_db =  Bio::EnsEMBL::DBLoader->new($recipient_locator);

$verbose && print "\nConnecting to donor database...\n";
my $donor_locator = $rec_db->get_donor_locator;
my $don_db =  Bio::EnsEMBL::DBLoader->new($donor_locator);

my $last_offset = $rec_db->get_last_update_offset;
$verbose && print "\nLast update-offset: $last_offset\n";

my $now_offset = $rec_db->get_now_offset;
$verbose && print "Time now-offset at recipient: $now_offset\n";

if ($last_offset > $now_offset) {
    print "Time of last_offset update more recent than now-offset, exiting!\n";
    exit;
}

$verbose && print "\nTransferring new ghosts from donor to recipient...\n";

my @object_array = $don_db->get_updated_Ghosts($last_offset, $now_offset);

#Sort with genes before exons (hopefully!)
@object_array = sort { if ( $a->obj_type eq 'gene' ) { return 1; } else { return $a->obj_type cmp $b->obj_type } } @object_array;

foreach my $ghost (@object_array) {
    $verbose && print "Got ghost with id ".$ghost->id."\n";
    
    #Check if the ghost is already present in the recipient database
    eval {
	$rec_db->get_Ghost($ghost->id,$ghost->obj_type);
    };
    
    #If not present in recipient, archive objects, and write the ghost
    if ( $@ ) {
	#If ghost of clone, delete clone
	if ($ghost->obj_type eq 'clone') {
	    $verbose && print "New ghost for deleted clone ".$ghost->id.", deleting clone from database\n";
	    $rec_db->delete_Clone($ghost->id);
	}
	#If ghost of gene, archive gene
	elsif ($ghost->obj_type eq 'gene') {
	    $verbose && print "New ghost for deleted gene ".$ghost->id.", archiving gene from recipient to archive\n";
	    my $gene = $rec_db->get_Gene($ghost->id);
	    $rec_db->archive_Gene($gene,$arc_db);
	}
	#If ghost of exon, delete exon
	elsif ($ghost->obj_type eq 'exon') {
            $verbose && print "New ghost for deleted exon ".$ghost->id.", deleting exon from database\n";
	    $rec_db->delete_Exon($ghost->id);
	}
	#Finally write ghost in recipient
	$verbose && print "Writing ghost in recipient database\n";
	$rec_db->write_Ghost($ghost);
    }

    #Nothing needs to be done if the ghost is already present
}

$verbose && print "\nTransferring updated and new objects from donor to recipient...\n";
#Get updated and new objects (clones and genes)
@object_array = [];
@object_array = $don_db->get_updated_Objects($last_offset, $now_offset);
#Should sort with clones first! Not implemented yet!

foreach my $object (@object_array) {
    my $type;
    
    #Check if it is a clone object
    if ($object->isa("Bio::EnsEMBL::DB::CloneI")) {
	$type = "Clone";
	$verbose && print "Got clone with id ".$object->id."\n";
	
	my $rec_clone;
	#Check if it is already present in recipient
	eval {
	    $rec_clone = $rec_db->get_Clone($object->id);
	};
	
        #If not present in recipient, write the clone in, and get all its genes out
	if ( $@ ) { 
	    $verbose &&  print "New Clone, writing it in the database\n";
	    $rec_db->write_Clone($object);
	    
	    #Get all the genes from this clone and check them, write them, archive them accordingly
	    foreach my $gene ($object->get_all_Genes()) {
		$verbose &&  print "Getting all genes via clone method\n";
		$verbose &&  print "Gene ".$gene->id." has version ".$gene->version."\n";
		&_place_gene($gene,'1');
	    }
	}
	
	#If present in recipient, check donor and recipient version
	else {
	    #If donor clone version greater than recipient clone version, replace
	    #NOTE: We permanently delete clones, without archiving them!
	    if ($object->version > $rec_clone->version) {
		$verbose && print "Clone with new version, updating the database, and deleting the old version\n";
		$rec_db->delete_Clone($rec_clone->id);
		$rec_db->write_Clone($object);
		#Get all the genes from this clone
		foreach my $gene ($object->get_all_Genes()) {
		    &_place_gene($gene,'1');
		}
	    }
	    
	    #If donor Clone version is less than the recipient Clone version, error 
	    #NOTE: Better catching to be implemented
	    elsif ($rec_clone->version > $object->version) {
		print "ERROR: Something is seriously wrong, found a clone in the recipient database with version number higher than that of the donor database!!!\n";
		exit;
	    }
	  	  
	    #If versions equal, nothing needs to be done
	    else {
		$verbose && print "Clone versions equal, not modifying database\n";
	    }
	}
    }

    #Check if it is a gene
    elsif ($object->isa("Bio::EnsEMBL::Gene")) {
	$verbose && print "In get_updated objects, gene level: got gene with id ".$object->id.", and version ".$object->version."\n";
	$type = "Gene";
	&_place_gene($object);
    }
    
    #Check if it is an exon
    elsif ($object->isa("Bio::EnsEMBL::Exon")) {
	$verbose && print "Got exon with id ".$object->id."\n";
	$type = "Exon";
	my $rec_exon;
	#Check if it is already present in recipient
	eval {
	    $rec_exon = $rec_db->get_Exon($object->id);
	};
	#If not present in recipient, write it in
	if ( $@ ) {
	    $verbose && print "New Exon, writing it in the database\n";
	    $rec_db->write_Exon($object);
	}
	#If exon present in recipient, check donor and recipient version
    else {
	#If donor exon version greater than recipient exon version, update
	if ($object->version > $rec_exon->version) {
	    $verbose && print "Exon with new version, updating the database\n";
	    $rec_db->delete_Exon($object->id);
	    $rec_db->write_Exon($object);
	}
	
	#If donor gene version is less than the recipient gene version, error 
	#NOTE: Better catching to be implemented
	elsif ($rec_exon->version > $object->version) {
	    print "Something is seriously wrong, found a gene in the recipient database with version number higher than that of the donor database!!!\n";
	}
	#If versions equal, nothing needs to be done
	else {
	    $verbose && print "Genes with the same version, databases kept unchanged\n";
	}
    }
	
    }
    else {
	print "ERROR: Got an object of unkown type with id ".$object->id."\n";	
    }
}


$rec_db->replace_last_update($now_offset);

sub _place_gene {
    my ($don_gene,$clone_level) = @_;
    my $rec_gene;

    #Check if the gene is present in the recipient
    eval {
	$rec_gene = $rec_db->get_Gene($don_gene->id);
    };
    
    #If gene not present in recipient, write it in
    if ( $@ ) {
	$verbose && print "New Gene, writing it in the database\n";
	$verbose &&  print "In _place_gene: Gene ".$don_gene->id." has version ".$don_gene->version."\n";
	$rec_db->write_Gene($don_gene);
    }
    #If gene present in recipient, check donor and recipient version
    else {
	#If donor gene version greater than recipient gene version, update
	if ($don_gene->version > $rec_gene->version) {
	    $verbose && print "Gene with new version, updating the database, and archiving old version\n";
	    $rec_db->archive_Gene($rec_gene,$arc_db);
	    $rec_db->write_Gene($don_gene);
	}
	
	#If donor gene version is less than the recipient gene version, error 
	#NOTE: Better catching to be implemented
	elsif ($rec_gene->version > $don_gene->version) {
	    print "Something is seriously wrong, found a gene in the recipient database with version number higher than that of the donor database!!!\n";
	}
	#If versions equal, nothing needs to be done
	else {
	    if ($clone_level) {
		$verbose && print "Genes with the same version, deleting recipient gene and writing one from donor without archiving\n";  
		$rec_db->delete_Gene($rec_gene->id);
		$rec_db->write_gene($don_gene);
	    }
	}
    }
}   
    
