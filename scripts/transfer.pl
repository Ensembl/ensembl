#!/usr/local/bin/perl

use strict;

#use Bio::EnsEMBL::AceDB::Obj;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::TimDB::Obj;
use Bio::AnnSeqIO;

use Getopt::Long;

my $fdbtype = 'timdb';
my $fhost   = 'croc';
my $fport   = '410000';
my $fdbname = 'ensdev';

my $tdbtype = 'rdb';
my $thost   = 'croc';
my $tport   = '410000';
my $tdbname = 'ensdev';

my $usefile = 0;
my $use_embl = 0;
my $tdbuser = 'ensembl';

&GetOptions( 
	     'fembl'     => \$use_embl, 
	     'fdbtype:s' => \$fdbtype,
	     'fhost:s'   => \$fhost,
	     'fport:n'   => \$fport,
	     'tdbtype:s' => \$tdbtype,
	     'fdbname:s' => \$fdbname,
	     'thost:s'   => \$thost,
	     'tport:n'   => \$tport,
	     'tdbname:s' => \$tdbname,
	     'tdbuser:s' => \$tdbuser,
	     'usefile'   => \$usefile,
	     );

my $from_db;
my $to_db;

my @clone;

if( $usefile == 1 ) {
    while( <> ) {
	my ($en) = split;
	push(@clone,$en);
    }
} else {
    @clone = @ARGV;
}

open(ERROR,">transfer.error\n");


if( $tdbtype =~ 'ace' ) {
    $to_db = Bio::EnsEMBL::AceDB::Obj->new( -user => $tdbuser, -host => $thost, -port => $tport);
} elsif ( $tdbtype =~ 'rdb' ) {
    $to_db = Bio::EnsEMBL::DBSQL::Obj->new( -user => $tdbuser, -db => $tdbname , -host => $thost );
} else {
    die("$tdbtype is not a good type (should be ace, rdb)");
}

if( $fdbtype =~ 'ace' ) {
    $from_db = Bio::EnsEMBL::AceDB::Obj->new( -host => $fhost, -port => $fport);
} elsif ( $fdbtype =~ 'rdb' ) {
    $from_db = Bio::EnsEMBL::DBSQL::Obj->new( -user => 'root', -db => $fdbname , -host => $fhost );
} elsif ( $fdbtype =~ 'timdb' ) {
    $from_db = Bio::EnsEMBL::TimDB::Obj->new(\@clone,0,0,1);
} else {
    die("$fdbtype is not a good type (should be ace, rdb or timdb)");
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

