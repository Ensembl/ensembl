#!/usr/local/bin/perl
use strict;
use Bio::EnsEMBL::DBOLD::Obj;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBSQL::StaticGoldenPathAdaptor;
use Bio::EnsEMBL::ExternalData::SNPSQL::DBAdapter;

my $chr=shift(@ARGV);
my $ensembldb=Bio::EnsEMBL::DBSQL::Obj->new(-dbname=>'ensembl080',-host=>'ensrv5',-user=>'ensadmin');
my $snpdb = Bio::EnsEMBL::ExternalData::SNPSQL::DBAdapter->new( -dbname=>'snp080', 
						       -user=>'ensro',  
						       -host=>'ensrv5');
$ensembldb->static_golden_path_type('UCSC');
$ensembldb->add_ExternalFeatureFactory($snpdb);
my $old_st=$ensembldb->get_StaticGoldenPathAdaptor;
print STDERR "Building Virtual Contig for $chr from old db...\n";
my $old_vc=$old_st->fetch_VirtualContig_by_chr_start_end($chr,20000000,20100000);

print STDERR "Getting features\n";
my @old_features=$old_vc->get_all_ExternalFeatures();

my $size=scalar(@old_features);
print "$chr: got $size features\n";
foreach my $feature (@old_features) {
    print $feature->id,"\t",$feature->upStreamSeq,"\t",
    $feature->dnStreamSeq,"\t",$feature->alleles,"\t", 
    $feature->clone_name,"\t",$feature->start_in_clone_coord,"\t",
    $feature->start,"\t",$_,"\n";
}




