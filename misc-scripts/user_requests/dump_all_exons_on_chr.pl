#!/usr/local/bin/perl
use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::StaticGoldenPathAdaptor;
use Bio::SeqIO;

my $chr=shift(@ARGV);

my $ensembldb=Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname=>'ensembl100',-host=>'kaka.sanger.ac.uk',-user=>'anonymous');
$ensembldb->static_golden_path_type('UCSC');

my $st=$ensembldb->get_StaticGoldenPathAdaptor;

print STDERR "Building Virtual Contig for $chr from db...\n";
my $vc=$st->fetch_VirtualContig_by_chr_name($chr);

print STDERR "Getting exons\n";
my @genes=$vc->get_all_Genes();
my $seqio = Bio::SeqIO->newFh('-format' => 'Fasta');
foreach my $gene (@genes) {
    foreach my $exon ($gene->each_unique_Exon) {
	my $seq = $exon->seq;
	my $desc = "$chr start: ".$exon->start." end:".$exon->end; 
	$seq->display_id($exon->id);
	$seq->desc($desc);
	print $seqio $seq;
    }
}

