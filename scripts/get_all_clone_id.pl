#!/usr/local/bin/perl

use strict;

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::TimDB::Obj;
use Getopt::Long;
use Time::Local;

my $dbtype = 'rdb';
my $host   = 'croc';
my $port   = '0';
my $dbname = 'ensdev';
my $update_date = undef;
my $tmap = 0;

&GetOptions( 'dbtype:s'  => \$dbtype,
	     'host:s'    => \$host,
	     'port:n'    => \$port,
	     'dbname:s'  => \$dbname,
	     'update:s'  => \$update_date,
	     'tmap'       => \$tmap,
	     );


my @empty;
my $db;

if( $dbtype =~ 'ace' ) {
    $db = Bio::EnsEMBL::AceDB::Obj->new( -host => $host, -port => $port);
} elsif ( $dbtype =~ 'rdb' ) {
    $db = Bio::EnsEMBL::DBSQL::Obj->new( -user => 'ensembl', -db => $dbname , -host => $host );
} elsif ( $dbtype =~ 'timdb' ) {
    # clones required are passed to speed things up - cuts down on parsing of flat files
    $db = Bio::EnsEMBL::TimDB::Obj->new(\@empty,0);
} else {
    die("$dbtype is not a good type (should be ace, rdb or timdb)");
}

my @clones;

if( $update_date ) {
    my ($year,$month,$day) = split(/:/,$update_date);

    my $time = timelocal(0,0,0,$day,$month,$year);
    print STDERR "Using time $time\n";

    @clones = $db->get_updated_Clone_id($time);
} else {
    @clones = $db->get_all_Clone_id();
}

if( $tmap ) {
    my @temp;

    foreach my $clone ( @clones ) {
	my ($id) = $db->get_id_acc($clone);
	push(@temp,$id);
    }

    @clones = @temp;
}

foreach my $clone ( @clones ) {
    print "$clone\n";
}

