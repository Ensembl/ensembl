#!/usr/local/bin/perl -w

use strict;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;

my $seqio = new Bio::SeqIO('-format' => 'fasta',
			  '-fh' => \*STDOUT);
my $host = 'localhost';
my $port   = '3360';
my $dbname = 'ens100';
my $dbuser = 'root';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::DBAdaptor';
my $upstreamlen = 2000;
 
my $locator = "$module/host=$host;port=;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);

my @gene_ids = $db->gene_Obj->get_all_Gene_id();
$db->static_golden_path_type('UCSC');
my $sa = $db->get_StaticGoldenPathAdaptor();
foreach my $geneid ( @gene_ids )  {
 my ($chr,$base) = $sa->get_Gene_chr_bp($geneid);
 my $vc = $sa->fetch_VirtualContig_by_chr_start_end($chr, $base - $upstreamlen, $base);
 $vc->desc($geneid);
 $vc->display_id($geneid);
 $seqio->write_seq($vc); 
}
