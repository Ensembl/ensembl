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
my $fhost   = 'sol28';
my $fport   = '410000';
my $fdbname = 'ensdev';
my $fdbuser = 'ensembl';
my $fdbpass = undef;

my $tdbtype = 'rdb';
my $thost   = 'sol28';
my $tport   = '410000';
my $tdbname = 'ensdev';
my $tdbuser = 'ensembl';
my $tdbpass = undef;

my $usefile = 0;
my $use_embl = 0;
my $cstart = 0;
my $cend;
my $getall = 0;
my $help;
my $fmodule = 'Bio::EnsEMBL::DBOLD::Obj';
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
	     
	     'embl'     => \$use_embl,
	     'getall'    => \$getall,
	     'usefile'   => \$usefile,
	     'start:i'   => \$cstart,
	     'end:i'     => \$cend,
	     'h|help'    => \$help
	     );

my $from_db;
my $to_db;
my @clone;


open(ERROR,">transfer.error\n");

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
	
	print ERROR "Clone $clone_id\n$@\n";
    }
}

close(ERROR);






