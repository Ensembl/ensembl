#!/usr/local/bin/perl

use strict;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBArchive::Obj;
use Bio::EnsEMBL::DBSQL::StaticGoldenPathAdaptor;
use Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor;
use Bio::EnsEMBL::Pipeline::GeneComp;

my $fpc=shift(@ARGV);
my $logfile=shift(@ARGV);
my $mapfile=shift(@ARGV);
open (LOG,">$logfile");
open (MAP,">$mapfile");
open (FILE,"</work2/elia/contigmap");
my %map;
while (<FILE>) {
    chomp;
    $_ =~ /(\S+)\s(\S+)/;
    $map{$1}=$2;
}
close (FILE);
my $cross=Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor->new(-dbname=>'cross110',
							-host=>'ecs1c',
							-user=>'ensadmin');
my $db=$cross->new_dbobj();
my $olddb=$cross->old_dbobj();


$db->static_golden_path_type('UCSC');
print STDERR "New db= $db\n";
print STDERR "Old db= $olddb\n";

my $st=$db->get_StaticGoldenPathAdaptor;

print STDERR "Building Virtual Contig for ctg $fpc...\n";
print LOG "Building Virtual Contig for ctg $fpc...\n";
my $vc=$st->fetch_VirtualContig_by_fpc_name($fpc);

my $arcdb = Bio::EnsEMBL::DBArchive::Obj->new(-dbname=>'archive',
					      -host=>'ecs1c',
					      -user=>'ensadmin',
					      -readonly => 0);

my $gc =  Bio::EnsEMBL::Pipeline::GeneComp->new(-vc => $vc,
						-archive => $arcdb,
						-log => \*LOG,
						-map => \*MAP,
						-hashref => \%map
						);

$gc->map();
    

