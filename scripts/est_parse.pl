#!/usr/local/bin/perl

=head1 NAME - EST dump

    Parses files in EST format (Gene Transfer Format)

=head1 SYNOPSIS - 

    ESTparse -dbname ensembl -parsefile genes.EST

=head1 DESCRIPTION

    This script parses EST files and writes the genes extracted to a database.
    The database is specified using the usual EnsEMBL options, described below.

    The actual parsing happens in the Bio::EnsEMBL::Utils::EST_handler module,
    which also handles the dumping of EST files.

    If the print option is specified, then the genes are not written to db, 
    but printed to STDOUT (mainly for testing)

=head1 OPTIONS

    -dbtype    database type (only needed for TimDB)

    -host      host name for the database (gets put as host= in locator) 

    -port      for RDBs, what port to connect to for the "to" database (port= in locator)

    -dbname    for RDBs, what database name to connect to (dbname= in locator)

    -dbuser    for RDBs, what username to connect as to the database (dbuser= in locator)

    -dbpass    for RDBs, what password to use to connect to to the database (dbpass= in locator)

    -module    module name to load to (Defaults to Bio::EnsEMBL::DBSQL::Obj)

    -parsefile name of the EST file to parse

    -print     prints gene structures to STDOUT

    -help      displays this documentation with PERLDOC

=cut

use Bio::EnsEMBL::Utils::EST_parser;
use Bio::EnsEMBL::DBLoader;
use strict;
use Getopt::Long;

#Database options
my $dbtype = 'rdb';
my $host   = 'localhost';
my $port   = '410000';
my $dbname = 'matloob_freeze17';
my $dbuser = 'root';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::Obj';

#Other options
my $parsefile;
my $print;
my $help;

&GetOptions( 
	     'dbtype:s'   => \$dbtype,
	     'host:s'     => \$host,
	     'port:n'     => \$port,
	     'dbname:s'   => \$dbname, 
	     'dbuser:s'   => \$dbuser,
	     'dbpass:s'   => \$dbpass,
	     'module:s'   => \$module,
       	     'parsefile:s'=> \$parsefile,
	     'print'      => \$print,
	     'h|help'     => \$help
	     );

my $ESTh=Bio::EnsEMBL::Utils::EST_parser->new();

open (PARSE,"$parsefile") || die("Could not open $parsefile for EST reading$!");
    
my @features=$ESTh->parse_file(\*PARSE);

if ($print) {
    $ESTh->print_ests;
}

 #DB writing option not yet implemented
 #Mapping of coordinates still needs to be done

else {
     my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
     my $db =  Bio::EnsEMBL::DBLoader->new($locator);
     my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($db);
     my $int_id=$features[0]->seqname;
     print STDERR "Got $int_id\n";
     my $contig=$db->get_Contig_by_international_id($int_id);
     $feature_obj->write($contig,@features);
}

