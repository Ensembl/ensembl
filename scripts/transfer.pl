#!/usr/local/bin/perl

use strict;

use Bio::EnsEMBL::AceDB::Obj;
use Bio::EnsEMBL::DB::Obj;
use Bio::EnsEMBL::TimDB::Obj;
use Bio::AnnSeqIO;

use Getopt::Long;

my $fdbtype = 'timdb';
my $fhost   = 'croc';
my $fport   = '410000';

my $tdbtype = 'rdb';
my $thost   = 'croc';
my $tport   = '410000';


&GetOptions( 
	     'fdbtype:s' => \$fdbtype,
	     'fhost:s'   => \$fhost,
	     'fport:n'   => \$fport,
	     'tdbtype:s' => \$tdbtype,
	     'thost:s'   => \$thost,
	     'tport:n'   => \$tport,
	     );

my $from_db;
my $to_db;


if( $fdbtype =~ 'ace' ) {
    $from_db = Bio::EnsEMBL::AceDB::Obj->new( -host => $fhost, -port => $fport);
} elsif ( $fdbtype =~ 'rdb' ) {
    $from_db = Bio::EnsEMBL::DB::Obj->new( -user => 'root', -db => 'ensdev' , -host => $fhost );
} elsif ( $fdbtype =~ 'timdb' ) {
    $from_db = Bio::EnsEMBL::TimDB::Obj->new(1);
} else {
    die("$fdbtype is not a good type (should be ace, rdb or timdb)");
}

if( $tdbtype =~ 'ace' ) {
    $to_db = Bio::EnsEMBL::AceDB::Obj->new( -host => $thost, -port => $tport);
} elsif ( $tdbtype =~ 'rdb' ) {
    $to_db = Bio::EnsEMBL::DB::Obj->new( -user => 'root', -db => 'ensdev' , -host => $thost );
} else {
    die("$tdbtype is not a good type (should be ace, rdb or timdb)");
}

foreach my $clone_id ( @ARGV ) {
    my $clone = $from_db->get_Clone($clone_id);
    
    $to_db->write_Clone($clone);
    
    foreach my $gene ( $clone->get_all_Genes() ) {
	$to_db->write_Gene($gene);
    }
}


