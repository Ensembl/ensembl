#!/usr/local/bin/perl

use strict;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBArchive::Obj;
use Bio::EnsEMBL::DBSQL::StaticGoldenPathAdaptor;
use Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor;
use Bio::EnsEMBL::Pipeline::GeneComp;

my $chr=shift(@ARGV);
my $logfile=shift(@ARGV);
my $mapfile=shift(@ARGV);
open (LOG,">$logfile");
open (MAP,">$mapfile");
my $cross=Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor->new(-dbname=>'cross_oct07',-host=>'ecs1a',-user=>'ensadmin');
my $db=$cross->new_dbobj();
my $olddb=$cross->old_dbobj();

$db->static_golden_path_type('UCSC');
print STDERR "New db= $db\n";
print STDERR "Old db= $olddb\n";

my $st=$db->get_StaticGoldenPathAdaptor;

print STDERR "Building Virtual Contig for chromosome $chr...\n";
print LOG "Building Virtual Contig for $chr...\n";
my $vc=$st->fetch_VirtualContig_by_chr_name($chr);
#print $vc->primary_seq->seq;

my $arcdb = Bio::EnsEMBL::DBArchive::Obj->new(-dbname=>'archive',-host=>'ecs1c',-user=>'ensadmin');
#No final db writing, only logging
#my $finaldb = Bio::EnsEMBL::DBSQL::Obj->new(-dbname=>'freeze05_final',-host=>'ecs1c',-user=>'ensadmin');

my $gc =  Bio::EnsEMBL::Pipeline::GeneComp->new(-vc => $vc,
						-archive => $arcdb,
						-log => \*LOG,
						-map => \*MAP);

$gc->map();
    

