#!/usr/local/bin/perl
=head1 NAME

Update

=head1 SYNOPSIS
 
  update.pl

=head1 DESCRIPTION

This script updates a recipient database by checking its donor database

=head1 OPTIONS

    -host      host name for database (gets put as host= in locator)

    -port      For RDBs, what port to connect to (port= in locator)

    -dbname    For RDBs, what name to connect to (dbname= in locator)

    -dbuser    For RDBs, what username to connect as (dbuser= in locator)

    -dbpass    For RDBs, what password to use (dbpass= in locator)

    -help      Displays script documentation with PERLDOC
    
    -nowrite   Runs entire script without writing in recipient

    -verbose   Gets all the print STDERR for testing purposes

=cut

use strict;
use Bio::EnsEMBL::TimDB::Obj;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Ghost;
use Bio::SeqIO;
use Getopt::Long;
use vars qw(@ISA);

@ISA = qw(Bio::Root::Object);

my $dbtype = 'rdb';
my $host   = 'localhost';
my $port   = '410000';
my $dbname = 'ensembl';
my $dbuser = 'ensro';
my $module = "Bio::EnsEMBL::DBSQL::Obj";
my $dbpass = undef;
my $help;
my $nowrite;
my $verbose;

&GetOptions( 
	     'dbtype:s'  => \$dbtype,
	     'host:s'    => \$host,
	     'port:n'    => \$port,
	     'dbname:s'  => \$dbname,
	     'dbuser:s'  => \$dbuser,
	     'dbpass:s'  => \$dbpass,
	     'module=s'  => \$module,
	     'h|help'    => \$help,
	     'nowrite'   => \$nowrite,
	     'v|verbose' => \$verbose
	     );


# signal handler
$SIG{INT}=sub {my $sig=shift;die "exited after SIG$sig";};

if ($help) {
    exec('perldoc', $0);
}

my $last_offset = 1;         $verbose && print "\nLast update-offset: $last_offset\n";
my $now_offset  = 949329046; $verbose && print "Time now-offset at recipient: $now_offset\n";

if ($last_offset > $now_offset) {
    print "Time of last_offset update more recent than now-offset, exiting!\n";
    exit;
}

$verbose && print "\nConnecting to recipient database...\n";
    
my $recipient_locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $rec_db =  Bio::EnsEMBL::DBLoader->new($recipient_locator);

$verbose && print "\nConnecting to donor database...\n\n";
# dummy list, to stop it loading everything
# object only used to get this list, else it will be missing the genes
my @clone_id;
{
    my @clones;
    my $don_db = Bio::EnsEMBL::TimDB::Obj->new(\@clones);

    # Get updated and new objects (clones and genes)
    @clone_id = $don_db->get_updated_Clone_id($last_offset, $now_offset);

    $verbose && print "\n".scalar(@clone_id)." clones found for updating\n";
}

$verbose && print "\nTransferring updated and new objects from donor to recipient...\n";

my $slice_size=2;
my @slice_array;

while(@slice_array=splice(@clone_id,0,$slice_size)){

    print "\nActive clones: ".join(',',@slice_array)."\n\n";

    {
	my @clone_array=@slice_array;
	my $don_db=Bio::EnsEMBL::TimDB::Obj->new(\@clone_array);

	# loop over clones
	foreach my $id (@slice_array) {
	    my $object=$don_db->get_Clone($id);
	    my $type;
	    
	    # Check if it is a clone object
	    if ($object->isa("Bio::EnsEMBL::DB::CloneI")) {
		$type = "Clone";
		$verbose && print "Got clone with id ".$object->id."\n";
	    }
	    
	    # Check if it is a gene
	    elsif ($object->isa("Bio::EnsEMBL::Gene")) {
		$verbose && print "Gene level: got gene with id ".$object->id.
		    ", and version ".$object->version."\n";
	    }
	    
	    # Check if it is an exon
	    elsif ($object->isa("Bio::EnsEMBL::Exon")) {
		$verbose && print "Got exon with id ".$object->id."\n";
	    }
	}
	$don_db->DESTROY;
    }
}


