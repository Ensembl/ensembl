#!/usr/local/bin/perl

use strict;

use Bio::EnsEMBL::DBLoader;
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
my $dbuser = 'ensro';
my $dbpass = undef;
my $port = 410000;
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
    my $locator = "Bio::EnsEMBL::DBSQL::Obj/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
    $db = Bio::EnsEMBL::DBLoader->new($locator);

} elsif ( $dbtype =~ 'timdb' ) {
    # clones required are passed to speed things up - cuts down on parsing of flat files
    $db = Bio::EnsEMBL::TimDB::Obj->new(-freeze => 1,-nogene => 1,\@empty,0);
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

