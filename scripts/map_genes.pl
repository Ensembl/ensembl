#!/usr/local/bin/perl
use strict;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBArchive::Obj;
use Bio::EnsEMBL::DBSQL::StaticGoldenPathAdaptor;
use Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor;
use Bio::EnsEMBL::Pipeline::GeneComp;

my $fpc=shift(@ARGV);
my $logfile=shift(@ARGV);
open (LOG,">>$logfile");
my $cross=Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor->new(-dbname=>'cross_sept05',-host=>'ecs1c',-user=>'ensadmin');
my $db=$cross->new_dbobj();
my $olddb=$cross->old_dbobj();

$db->static_golden_path_type('UCSC');
#print STDERR "db= $db\n";


my $st=$db->get_StaticGoldenPathAdaptor;

#print STDERR "st= $st\n";
print STDERR "Building Virtual Contig for $fpc...\n";
print LOG "Building Virtual Contig for $fpc...\n";
my $vc=$st->fetch_VirtualContig_by_fpc_name($fpc);
my $arcdb = Bio::EnsEMBL::DBArchive::Obj->new(-dbname=>'archive',-host=>'ecs1c',-user=>'ensadmin');
my $finaldb = Bio::EnsEMBL::DBSQL::Obj->new(-dbname=>'freeze05_final',-host=>'ecs1c',-user=>'ensadmin');
my $gc =  Bio::EnsEMBL::Pipeline::GeneComp->new(-vc => $vc,
						-archive => $arcdb,
						-finaldb => $finaldb,
						-log => \*LOG);

$gc->map();
    

