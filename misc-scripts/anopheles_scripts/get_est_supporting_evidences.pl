use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Utils::Eprof('eprof_start','eprof_end','eprof_dump');

my $host      = 'kaka.sanger.ac.uk';
my $dbuser    = 'anonymous';
my $dbname    = 'anopheles_gambiae_estgene_11_2';
my $dbpass    = '';
my $path      = 'MOZ2';


print STDERR "Connecting to $host, $dbname\n";


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $host,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					   );

my $slice_adapt = $db->get_SliceAdaptor();


my $query1 = "select name from contig";

my $sth1 = $db->prepare($query1);
$sth1->execute();

while (my $contig_name = $sth1->fetchrow_array) {
    #print STDERR "NEW_CONTIG\n";

    my $slice = $slice_adapt->fetch_by_contig_name($contig_name);
    
    my @genes = @{$slice->get_all_Genes()};
    
    my $features = $slice->get_all_ProteinAlignFeatures();

    foreach my $gene(@genes) {
	
	my @transcripts = @{$gene->get_all_Transcripts()};

	foreach my $trans(@transcripts) {
	    
	    my $exons = $trans->get_all_Exons();
	    
	  FEATURE:
	    for my $feature ( @$features ) {
		print STDERR "FEAT: ".$feature->analysis->dbID."\n";
	      EXON:
		for my $exon ( @$exons ) {
		    if( $feature->start() > $exon->start &&
			$feature->end() < $exon->end ) {
			print $exon->dbID."\tprotein_align_feature\t".$feature->dbID."\n"
			}
		}
	    }
	}
    }
}
