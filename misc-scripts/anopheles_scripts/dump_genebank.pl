 use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::ProteinAdaptor;

use Bio::EnsEMBL::Utils::Eprof('eprof_start','eprof_end','eprof_dump');

my $host      = 'kaka.sanger.ac.uk';
my $dbuser    = 'anonymous';
my $dbname    = 'anopheles_gambiae_core_10_2';
my $dbpass    = '';
my $path      = 'MOZ2';

print STDERR "Connecting to $host, $dbname\n";


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $host,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					   );

my $clone_adapt = $db->get_CloneAdaptor();
my $gene_adapt = $db->get_GeneAdaptor();
my $slice_adapt = $db->get_SliceAdaptor();

open (OUT,">test.out");

my $query1 = "select clone_id,name from clone where name = 'AAAB01008961'";
my $sth1 = $db->prepare($query1);
$sth1->execute();

while(my ($id,$clone_name) = $sth1->fetchrow_array) {

    my $slice = $slice_adapt->fetch_by_clone_accession($clone_name);
    
    # my $clone_seq = $slice->seq;
    my $chr_name = $slice->chr_name;

    # print STDERR ">gnl|WGS:AAAB|$clone_name [organism=Anopheles gambiae str. PEST] [tech=wgs] [chromosome=$chr_name]\n$clone_seq\n";

    print OUT ">Feature gnl|WGS:AAAB|$clone_name\n";
    foreach my $gene(@{$slice->get_all_Genes}) {

	my $gene_id = $gene->dbID;
	my $start = $gene->start;
	my $end = $gene->end;
	my $name = $gene->stable_id;

	my $ebi;
	my $cel_gene;
	my $cel_transcript;
	my $cel_translation;

	if ($gene->strand == 1) {
	    print OUT "$start\t$end\tgene\n";
	}
	else {
	    print OUT "$end\t$start\n";
	}
	print OUT "\t\t\tlocus_tag\t$name\n";
	
	foreach my $tr(@{$gene->get_all_Transcripts}) {
	    my $count_ex = 0;
	    my $count_lation = 0;
	    my $tr_name = $tr->stable_id;
	    my $tr_dbID = $tr->dbID;
	    foreach my $ex(@{$tr->get_all_Exons}) {
		$count_ex++;

		my $ex_start = $ex->start;
		my $ex_end = $ex->end;
		if ($count_ex == 1) {
		    print OUT "$ex_start\t$ex_end\tmRNA\n";
		}
		else {
		    print OUT "$ex_start\t$ex_end\n";
		}
	    }

	    
	    print OUT "\t\t\tproduct\t$tr_name\n";
	    
	    my $translation = $tr->translation;
	    my $translation_name = $translation->stable_id;

	    	my $tr_dbID = $translation->dbID;

		my $query3 = "select x.display_label from xref x, external_db e, object_xref o where o.ensembl_id = $tr_dbID and o.xref_id = x.xref_id and x.external_db_id = e.external_db_id and (e.db_name = 'Celera_Pep' or e.db_name = 'Anopheles_paper')";
		my $sth3 = $db->prepare($query3);
		$sth3->execute();
		
		my $protein_id = $sth3->fetchrow;
	    

	    foreach my $lation (@{$tr->get_all_translateable_Exons}) {
		$count_lation++;
		
		my $tr_start = $lation->start;
		my $tr_end = $lation->end;
	
		if ($count_lation == 1) {
		    print OUT "$tr_start\t$tr_end\tCDS\n";
		}
		else {
		    print OUT "$tr_start\t$tr_end\n";
		}
	    }
	    print OUT "\t\t\tproduct\t$translation_name\n";
	    if ($protein_id) {
		print OUT "\t\t\tprotein_id\t$protein_id\n";
	    }
	    print OUT "\t\t\tevidence\tnot_experimental\n";
	}
    }
}

