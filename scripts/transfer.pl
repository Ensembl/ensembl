#!/usr/local/bin/perl

=head1 NAME - transfer

 Transfers clones from one database to another

=head1 SYNOPSIS - 

    transfer dJ271M21
   
    clone2embl -fdbtype ace -tdbtype rdb dJ271M21 

    clone2embl -fdbtype rdb -host mysql.server.somewhere dJ271M21

=head1 OPTIONS

    -fdbtype   "from" database type (only needed for TimDB)
    -tdbtype   "to" database type (only needed for TimDB)

    -fhost     host name for "from"database (gets put as host= in locator)
    -thost     host name for "to" database (gets put as host= in locator) 

    -fport     for RDBs, what port to connect to for the "from" database (port= in locator)
    -tport     for RDBs, what port to connect to for the "to" database (port= in locator)

    -fdbname   for RDBs, what "from" database name to connect to (dbname= in locator)
    -tdbname   for RDBs, what "to" database name to connect to (dbname= in locator)

    -fdbuser   for RDBs, what username to connect as to the "from" database (dbuser= in locator)
    -tdbuser   for RDBs, what username to connect as to the "to" database (dbuser= in locator)

    -fdbpass   for RDBs, what password to use to connect to to the "from" database (dbpass= in locator)
    -tdbpass   for RDBs, what password to use to connect to to the "to" database (dbpass= in locator)

    -fmodule   module name to load from (Defaults to Bio::EnsEMBL::DBOLD::Obj)
    -tmodule   module name to load to (Defaults to Bio::EnsEMBL::DBSQL::Obj)

    -getall    all clones from the database [not applicable to timdb]

    -update    Used to update a clone, and store the old version in the archive database

    -arcpass   password for the archive database

    -usefile   read in on stdin a list of clones, one clone per line

    -start     start point in list of clones (useful with -getall)

    -end       end point in list of clones (useful with -getall)

    -help      Displays this documentation with PERLDOC

=head1 EXAMPLE CLONES

    dJ271M21/AL031983   T  single contig, mainly forward strand genes, but one reverse

    dJ718J7                single contig, single reverse strand gene with partial transcripts

    C361A3/AL031717     T  unfinished sanger, 3 contigs, 2 ordered by gene predictions

    C367G8/Z97634       T  unfinished sanger, 2 contigs, not ordered by gene predictions

    AP000228               External finished clone

    AC005776               External finished clone

    AC010144               External unfinished clone

=cut


use strict;

use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::TimDB::Obj;
use Bio::AnnSeqIO;
use Getopt::Long;

my $fdbtype = 'rdb';
my $fhost   = 'localhost';
my $fport   = '410000';
my $fdbname = 'ensdev';
my $fdbuser = 'ensembl';
my $fdbpass = undef;

my $tdbtype = 'rdb';
my $thost   = 'localhost';
my $tport   = '410000';
my $tdbname = 'ensdev';
my $tdbuser = 'ensembl';
my $tdbpass = undef;

my $usefile = 0;
my $use_embl = 0;
my $cstart = 0;
my $cend;
my $getall = 0;
my $update = 0;
my $arcpass = undef;
my $help;
my $fmodule = 'Bio::EnsEMBL::DBSQL::Obj';
my $tmodule = 'Bio::EnsEMBL::DBOLD::Obj';

&GetOptions( 
	     'fdbtype:s' => \$fdbtype,
	     'fhost:s'   => \$fhost,
	     'fport:n'   => \$fport,
	     'fdbname:s' => \$fdbname, 
	     'fdbuser:s' => \$fdbuser,
	     'fdbpass:s' => \$fdbpass,
	     'fmodule:s' => \$fmodule,

	     'tdbtype:s' => \$tdbtype,
	     'thost:s'   => \$thost,
	     'tport:n'   => \$tport,
	     'tdbname:s' => \$tdbname,
	     'tdbuser:s' => \$tdbuser,
	     'tdbpass:s' => \$tdbpass,
	     'tmodule:s' => \$tmodule,
	     
	     'embl'      => \$use_embl,
	     'getall'    => \$getall,
	     'update'    => \$update,
	     'arcpass:s' => \$arcpass,
	     'usefile'   => \$usefile,
	     'start:i'   => \$cstart,
	     'end:i'     => \$cend,
	     'h|help'    => \$help
	     );

my $from_db;
my $to_db;
my @clone;

  

if ($help){
    exec('perldoc', $0);
}

if( $usefile == 1 ) {
    while( <> ) {
	my ($en) = split;
	push(@clone,$en);
    }
} else {
    @clone = @ARGV;
}

if ( $tdbtype =~ 'timdb' ) {
    $to_db = Bio::EnsEMBL::TimDB::Obj->new(\@clone,0,0,1);
} else {
    my $locator = "$tmodule/host=$thost;port=$tport;dbname=$tdbname;user=$tdbuser;pass=$tdbpass";
    $to_db =  Bio::EnsEMBL::DBLoader->new($locator);
}

if ( $fdbtype =~ 'timdb' ) {
    $from_db = Bio::EnsEMBL::TimDB::Obj->new(\@clone,0,0,1);
} else {
    my $locator = "$fmodule/host=$fhost;port=$fport;dbname=$fdbname;user=$fdbuser;pass=$fdbpass";
    print STDERR "LOCATOR IS: $locator\n";
    $from_db = Bio::EnsEMBL::DBLoader->new($locator); 
}

if ( $getall == 1 ) {
    @clone = $from_db->get_all_Clone_id();
    print STDERR scalar(@clone)." clones found in DB\n";
}

if( defined $cend ) {
    print STDERR "splicing $cstart to $cend\n";
    my @temp = splice(@clone,$cstart,($cend-$cstart));
    @clone = @temp;
}

#This section is for the update mode, i.e. when a new version of a clone 
#needs to be stored, some of the data from the old version gets transferred to 
#the archive db, and the clone is deleted from the database, and the new one is written in

if ($update) {
    print STDERR "Update mode: storing old version in archivedb, deleting it from $fdbname, and storing new version in $fdbname!\n";
    
    my $arcname = 'archive';
    my $arctype = 'rdb';
    my $archost = 'localhost';
    my $arcport = '410000';
    my $arcuser = 'root';
    my $arcmodule = 'Bio::EnsEMBL::DBArchive::Obj';
    $arcpass || die "You forgot to give the password for the archive database! Sorry, no access!\n";
    
    my $arclocator = "$arcmodule/host=$archost;port=$arcport;dbname=$arcname;user=$arcuser;pass=$arcpass";
    my $arc_db = Bio::EnsEMBL::DBLoader->new($arclocator);
    
    
    
    foreach my $clone (@clone) {
	
	#First we need to get the information we need to store from the old version of the clone
	#Note: we are getting the clone from to_db, where the new version will be written
	
	print STDERR "Loading $clone\n"; 
		
	#Then we need to store the information for all transcripts, exons and proteins in the archive db
	foreach my $gene ($clone->get_all_Genes) {
	    #Delete genes,transcipts,translations, and exons from $to_db and store partial info in archive db
	    $to_db=archive_Gene($gene,$clone,$arc_db);
	}
	
	#Delete clone from $to_db
	my $clone = $to_db->delete_Clone($clone,$arc_db);
    }

    #Finally we proceed as normal, transferring the clone from from_db to to_db (out of the update loop)
}

foreach my $clone_id ( @clone ) {
    print STDERR "Loading $clone_id\n";
    eval {
	my $clone = $from_db->get_Clone($clone_id);
	
	$to_db->write_Clone($clone);
	
	foreach my $gene ( $clone->get_all_Genes() ) {
	    $to_db->write_Gene($gene);
	}
    };
    if ( $@ ) {
	
	print STDERR "Could not transfer clone $clone_id\n$@\n";
    }
}

close(ERROR);















