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

use Bio::EnsEMBL::Analysis::UpdateManager;

use strict;
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
my $slice;

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
	     'slice:s'   => \$slice,
	     'v|verbose' => \$verbose
	     );



if ($help) {
    exec('perldoc', $0);
}

my $last_offset = 949000000;  
my $now_offset  = time; 

$| = 1;

if ($last_offset > $now_offset) {
    print "Time of last_offset update more recent than now-offset, exiting!\n";
    exit;
}

my $from_locator     = "Bio::EnsEMBL::TimDB::Obj";
my $to_locator       = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";

my $update_manager   = new Bio::EnsEMBL::Analysis::UpdateManager(-fromlocator => $from_locator,
								 -tolocator   => $to_locator,
								 -fromtime    => $last_offset,
								 -totime      => $now_offset,
								 );

$update_manager->chunksize(20);
$update_manager->update;




