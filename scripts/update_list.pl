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


if ($help) {
    exec('perldoc', $0);
}

my $don_db;

$verbose && print "\nConnecting to recipient database...\n";
    
my $recipient_locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $rec_db =  Bio::EnsEMBL::DBLoader->new($recipient_locator);

$verbose && print "\nConnecting to donor database...\n";
$don_db = Bio::EnsEMBL::TimDB::Obj->new(); 

my $last_offset = 1;         $verbose && print "\nLast update-offset: $last_offset\n";
my $now_offset  = 949329046; $verbose && print "Time now-offset at recipient: $now_offset\n";

if ($last_offset > $now_offset) {
    print "Time of last_offset update more recent than now-offset, exiting!\n";
    exit;
}

$verbose && print "\nTransferring updated and new objects from donor to recipient...\n";

#Get updated and new objects (clones and genes)
my  @object_array = $don_db->get_updated_Objects($last_offset, $now_offset);

foreach my $object (@object_array) {
    my $type;
    
    #Check if it is a clone object
    if ($object->isa("Bio::EnsEMBL::DB::CloneI")) {
	$type = "Clone";
	$verbose && print "Got clone with id ".$object->id."\n";
    }

    #Check if it is a gene
    elsif ($object->isa("Bio::EnsEMBL::Gene")) {
	$verbose && print "Gene level: got gene with id ".$object->id.", and version ".$object->version."\n";
    }
    
    #Check if it is an exon
    elsif ($object->isa("Bio::EnsEMBL::Exon")) {
	$verbose && print "Got exon with id ".$object->id."\n";
    }
}


