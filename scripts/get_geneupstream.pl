#!/usr/local/bin/perl -w

use strict;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Getopt::Long;

my $host = 'localhost';
my $port   = '3360';
my $dbname = 'ens100';
my $dbuser = 'root';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::DBAdaptor';
my $upstreamlen = 2000;

&GetOptions( 
	     'host:s'          => \$host,
	     'port:n'          => \$port,
	     'db|dbname:s'     => \$dbname,
	     'user|dbuser:s'   => \$dbuser,
	     'p|dbpass:s'      => \$dbpass,
	     'm|module:s'      => \$module,
	     'u|upstream:i'    => \$upstreamlen,
	     );

my $locator = "$module/host=$host;port=;dbname=$dbname;user=$dbuser;pass=$dbpass";

my $seqio = new Bio::SeqIO('-format' => 'fasta',
			   '-fh' => \*STDOUT);

my $db =  Bio::EnsEMBL::DBLoader->new($locator);

my @gene_ids = $db->gene_Obj->get_all_Gene_id();
$db->static_golden_path_type('UCSC');
my $sa = $db->get_StaticGoldenPathAdaptor();
foreach my $geneid ( @gene_ids )  {
    my $gene = $db->gene_Obj->get($geneid);
    # just want first one
    # going to assume here that the first exon is 
    # on the same strand as the rest of the gene
    # which I know is not always true in
    # the Ensembl world because of misassembled sequence
    # in GP
    my ($exon) = $gene->each_unique_Exon();
    my ($chr,$base) = $sa->get_Gene_chr_bp($geneid);
    my $vc;
    my $strand = $exon->strand;
    if(  $strand == -1 ) {
	# gene is on the rev strand, go forwards 
	# (in the global chromosome sense)
	# '$upstream' bases
	$vc = $sa->fetch_VirtualContig_by_chr_start_end($chr, $base, 
							$base + $upstreamlen);
    } else {
	# gene is on fwd strand, go backwards (in global chromsome sense)
	# '$upstream' bases
	$vc = $sa->fetch_VirtualContig_by_chr_start_end($chr, 
							$base - $upstreamlen, 
							$base);
    }    
    $vc->id($geneid);
    $vc->desc("Strand:$strand Chrom:$chr basepair:$base" );
    $seqio->write_seq($vc); 
}
