use strict;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBSQL::StaticGoldenPathAdaptor;
use Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor;
use Bio::EnsEMBL::Pipeline::GeneComp;

my $logfile=shift(@ARGV);
open (LOG,">$logfile");
my $db=Bio::EnsEMBL::DBSQL::Obj->new(-dbname=>'july_dna',-host=>'ecs1c',-user=>'ensadmin');
my $cross=Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor->new(-dbname=>'crossmatch',-host=>'ecs1c',-user=>'ensadmin');
my $olddb=Bio::EnsEMBL::DBSQL::Obj->new(-dbname=>'june_dna',-host=>'ecs1c',-user=>'ensadmin');
$cross->_old_dbobj($olddb);

$db->_crossdb($cross);

$db->static_golden_path_type('UCSC');
print STDERR "db= $db\n";


my $st=$db->get_StaticGoldenPathAdaptor;

print STDERR "st= $st\n";

print STDERR "Building Virtual Contig...\n";
my $vc=$st->fetch_VirtualContig_by_chr_start_end('chr1',1,100000);

my (%temp_old,$mapped,$new,$untransf) = Bio::EnsEMBL::Pipeline::GeneComp::map_temp_Exons_to_real_Exons($vc,\*LOG);
my @keys = keys (%temp_old);
my $size = scalar @keys;

print STDERR "Got $size exons mapped in GeneComp\n";

my ($deadgeneid,$deadtranscriptid) = Bio::EnsEMBL::Pipeline::GeneComp::map_temp_Genes_to_real_Genes($vc,\*LOG,%temp_old);
my @d_genes=@$deadgeneid;
my @d_trans=@$deadtranscriptid;
print STDERR "Dead genes:\n";
foreach my $dg (@d_genes) {
    print STDERR "$dg\n";
}
print STDERR "Dead transcripts:\n";
foreach my $dt (@d_trans) {
    print "$dt\n";
}
			   
    

