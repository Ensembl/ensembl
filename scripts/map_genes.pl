use strict;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBArchive::Obj;
use Bio::EnsEMBL::DBSQL::StaticGoldenPathAdaptor;
use Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor;
use Bio::EnsEMBL::Pipeline::GeneComp;

my $chr=shift(@ARGV);
my $logfile=shift(@ARGV);
open (LOG,">$logfile");
my $cross=Bio::EnsEMBL::DBSQL::CrossMatchDBAdaptor->new(-dbname=>'crossmatch',-host=>'ecs1c',-user=>'ensadmin');
my $db=$cross->new_dbobj();
my $olddb=$cross->old_dbobj();

$db->static_golden_path_type('UCSC');
print STDERR "db= $db\n";


my $st=$db->get_StaticGoldenPathAdaptor;

print STDERR "st= $st\n";

print STDERR "Building Virtual Contig...\n";
my $vc=$st->fetch_VirtualContig_by_chr_name($chr);
my $arcdb = Bio::EnsEMBL::DBArchive::Obj->new(-dbname=>'test_archive',-host=>'ecs1c',-user=>'ensadmin');
my $finaldb = Bio::EnsEMBL::DBSQL::Obj->new(-dbname=>'test_final',-host=>'ecs1c',-user=>'ensadmin');
my $gc =  Bio::EnsEMBL::Pipeline::GeneComp->new(-vc => $vc,
						-archive => $arcdb,
						-finaldb => $finaldb,
						-log => \*LOG);

$gc->map();
#my (%temp_old,$mapped,$new,$untransf) = Bio::EnsEMBL::Pipeline::GeneComp::map_temp_Exons_to_real_Exons($vc,\*LOG);
#my @keys = keys (%temp_old);
#my $size = scalar @keys;

#print STDERR "Got $size exons mapped in GeneComp\n";

#my ($deadgeneid,$deadtranscriptid) = Bio::EnsEMBL::Pipeline::GeneComp::map_temp_Genes_to_real_Genes($vc,\*LOG,%temp_old);
#my @d_genes=@$deadgeneid;
#my @d_trans=@$deadtranscriptid;
#print STDERR "Dead genes:\n";
#foreach my $dg (@d_genes) {
#    print STDERR "$dg\n";
#}
#print STDERR "Dead transcripts:\n";
#foreach my $dt (@d_trans) {
#    print "$dt\n";
#}
			   
    

