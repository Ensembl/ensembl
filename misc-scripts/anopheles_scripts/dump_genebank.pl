use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::ProteinAdaptor;

use Bio::EnsEMBL::Utils::Eprof('eprof_start','eprof_end','eprof_dump');

my $host      = 'localhost';
my $dbuser    = 'manu';
my $dbname    = 'anopheles_gambiae_core_10_2';
my $dbpass    = '';
my $path      = 'MOZ2';

my %utr;
my %utr_tmp;

my %map;

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

open (MAP,"/Users/emmanuelmongin/work/ncbi_dump/AAAB01.output.p2g") || die;
open (OUT,">/Users/emmanuelmongin/work/ncbi_dump/test.out") || die;

while(<MAP>) {
    chomp;
    my @array = split;
    $map{$array[0]} = $array[1];
}

close(MAP);

#Get all of the scaffolds
my $query1 = "select clone_id,name from clone where name = 'AAAB01008961'";
my $sth1 = $db->prepare($query1);
$sth1->execute();

while(my ($id,$clone_name) = $sth1->fetchrow_array) {

    my $slice = $slice_adapt->fetch_by_clone_accession($clone_name);
    
    my $chr_name = $slice->chr_name;

    # print STDERR ">gnl|WGS:AAAB|$clone_name [organism=Anopheles gambiae str. PEST] [tech=wgs] [chromosome=$chr_name]\n$clone_seq\n";

    print OUT ">Feature gnl|WGS:AAAB|$clone_name\n";
    
    my @genes = @{$slice->get_all_Genes};

    foreach my $gene(@genes) {
	
	my $gene_id = $gene->dbID;
	
	my @transcripts = @{$gene->get_all_Transcripts};

	my $tr_start = $transcripts[0]->start;
	
	foreach my $tr(@transcripts) {
	    
	    $tr_start = $tr->start;
	}

	my $new_gene = &checks($gene);

	my @new_transcripts = @{$new_gene->get_all_Transcripts};

#print the gene coordinates

	my $new_gene_dbID = $new_gene->dbID;
	
	my $start = $gene->start;
	my $end = $gene->end;
	my $name = $gene->stable_id;

	my $ebi;
	my $cel_gene;
	my $cel_transcript;
	my $cel_translation;

	my $cdna_start;
	my $cdna_end;
	my $coding_start;
	my $coding_end;

#	print OUT "UTR COUNT ". $utr{$new_gene_dbID}->{'count'}."\n";
#	print OUT "UP ".$utr{$new_gene_dbID}->{'up'}."\n";
#	print OUT "DOWN ".$utr{$new_gene_dbID}->{'down'}."\n";

	if ($utr{$new_gene_dbID}->{'count'} == 2) {
	    if ($gene->strand == 1) {
		print OUT "$start\t$end\tgene\n";
	    }
	    else {
		print OUT "$end\t$start\tgene\n";
	    }
	}
	
	elsif (($utr{$new_gene_dbID}->{'count'} == 1)&&($utr{$new_gene_dbID}->{'up'})) {
	    if ($gene->strand == 1) {
		print OUT "$start\t>$end\tgene\n";
	    }
	    else {
		print OUT "$end\t>$start\tgene\n";
	    }
	}
	
	elsif (($utr{$new_gene_dbID}->{'count'} == 1)&&($utr{$new_gene_dbID}->{'up'})) {
	    if ($gene->strand == 1) {
		print OUT "<$start\t$end\tgene\n";
	    }
	    else {
		print OUT "<$end\t$start\tgene\n";
	    }
	}
	
	elsif ($utr{$new_gene_dbID}->{'count'} == 0) {
	    if ($gene->strand == 1) {
		print OUT "<$start\t>$end\tgene\n";
	    }
	    else {
		print OUT "<$end\t>$start\tgene\n";
	    }
	}
	print OUT "\t\t\tlocus_tag\t$name\n";

	foreach my $new_tr(@new_transcripts) {
	    &print_transcript_coordinates($new_tr);
	    &print_translation_coordinates($new_tr,$db);
	}
    }
}

sub checks {
    my ($gene) = @_;
    
    my $gene_dbID = $gene->dbID;
    
    my @transcripts = @{$gene->get_all_Transcripts};

    my @new_transcripts;
    
    foreach my $tr(@transcripts) { 
	my $tr_dbID = $tr->dbID;
	my $c_start = $tr->cdna_coding_start;
	my $c_end = $tr->cdna_coding_end;
	my $spl_sq = $tr->spliced_seq;
	
	my $start = $tr->start;
	my $end = $tr->end;
	
	print STDERR "START: $start\tEND: $end\n";
	
	my $pep = $tr->translate->seq;

	#Check if the transcipt has a 5' UTR
	if (($c_start > 1)&&($c_end == length($spl_sq))) {
	    #Transcript has a 5' UTR
	  
	    #Check for the stop codon methionine
	    if ($pep !~ /^M/) {
		#No methionine, in that case the cdna start should be equal to the coding start
	

		$tr = &chop_5($tr);

		push(@new_transcripts,$tr);
		
	    }
	    else {
		#The transcript has a UTR and a methionine, we keep the transcript as it is
		$utr_tmp{$gene_dbID}->{'up'} = 1;
		$utr_tmp{$gene_dbID}->{'count'}++;
		
		push(@new_transcripts,$tr);
	    }
	}
	
	#Check if the transcript has a 3' UTR
	
	if (($c_start == 1)&&($c_end < length($spl_sq))) {
	    #The transcript has a 3' UTR

	    if ($pep !~ /\*$/) {
		#THere is no stop codon but a UTR, chop the UTR

		$tr = &chop_3($tr);
		push(@new_transcripts,$tr);
	    }
	    else {
		#There is a stop codon with a UTR, leave the transcript as it is
		$utr_tmp{$gene_dbID}->{'down'} = 1;
		$utr_tmp{$gene_dbID}->{'count'}++;
			
		push(@new_transcripts,$tr);
	    }
	}

	if (($c_end == length($spl_sq))&&($c_start == 1)) {
	    #Nothig to do in that case, no UTRs
	    push(@new_transcripts,$tr);
	}
	if  (($c_start > 1)&&($c_end < length($spl_sq))) {
	    
	    if (($pep !~ /\*$/)&&($pep !~ /^M/)) {
		$tr = &chop_5($tr);
		$tr = &chop_3($tr);
		push(@new_transcripts,$tr);
	    }
	    elsif(($pep !~ /\*$/)&&($pep =~ /^M/)) {
		$tr = &chop_3($tr);
		push(@new_transcripts,$tr);
	    }
	    elsif (($pep =~ /\*$/)&&($pep !~ /^M/)) {
		$tr = &chop_5($tr);
		push(@new_transcripts,$tr);
	    }
	    elsif (($pep =~ /\*$/)&&($pep =~ /^M/)) {
		
		$utr_tmp{$gene_dbID}->{'up'} = 1;
		$utr_tmp{$gene_dbID}->{'down'} = 1;
		$utr_tmp{$gene_dbID}->{'count'}++;
		$utr_tmp{$gene_dbID}->{'count'}++;
		push(@new_transcripts,$tr);
	    }
	    else {
		die;
	    }
		
	}
	
	if ($utr_tmp{$gene_dbID}->{'count'} > $utr{$gene_dbID}->{'count'}) {
	    $utr{$gene_dbID} = $utr_tmp{$gene_dbID};
	    $utr_tmp{$gene_dbID} = undef;
	    
	}
    }

    $gene->{'_transcript_array'} =[];

    foreach my $new_tr(@new_transcripts) {
	
	$gene->add_Transcript($new_tr);
    }
    
    return $gene;

}

sub chop_5 {
    #chop 5' UTR, return the choped transcript
    my ($tr) = @_;
    
    my @translateable = @{$tr->get_all_translateable_Exons};
    my $nb1 = scalar(@translateable);
    
    my @new_exons;
    
    my $new_first_exon = $translateable[0];
    my $first_exon_id = $new_first_exon->dbID;

    my $tr_start = $new_first_exon->start();
    
    my $tr_end;
    my @all_exons =  @{$tr->get_all_Exons};
    
    my $seen;
    foreach my $e(@all_exons) {
	my $ex_dbid = $e->dbID;
	
	if ($ex_dbid == $first_exon_id) {
	    $seen = 1;
	}
	if ($seen) {
	    push(@new_exons,$e);
	}
    }
    
    $tr->flush_Exons();
		
    foreach my $ne (@new_exons) {
	$tr->add_Exon($ne);
    }
    
    my $nb_ex = scalar(@new_exons);
    
    my $translation = $tr->translation;

    $translation->start_Exon($new_exons[0]);
    $translation->end_Exon($new_exons[$nb_ex-1]);

    return $tr;
}

sub chop_3 {
    #chop 3' UTR, return the choped transcript
    my ($tr) = @_;
    
    my @translateable = @{$tr->get_all_translateable_Exons};
    
    my $nb1 = scalar(@translateable) - 1;
    
    my @new_exons;
    
    my $new_last_exon = $translateable[$nb1];
    my $last_exon_id = $new_last_exon->dbID;
    
    my $tr_end;
    my @all_exons =  @{$tr->get_all_Exons};
    
    my $seen;
    foreach my $e(@all_exons) {
	my $ex_dbid = $e->dbID;
	if($seen) {
	    next;
	}
	elsif ($ex_dbid == $last_exon_id) {
	    push(@new_exons,$e);
	    $seen = 1;
	}
	elsif ($ex_dbid != $last_exon_id) {
	    push(@new_exons,$e);
	}
	
    }

    my $nb_ex = scalar(@new_exons);
    
    my $translation = $tr->translation;
    
    $translation->start_Exon($new_exons[0]);
    $translation->end_Exon($new_exons[$nb_ex-1]);

    return $tr;
}

sub print_transcript_coordinates {
    my ($tr) = @_;
    
    my $coding_start = $tr->cdna_coding_start;
    my $coding_end = $tr->cdna_coding_end;
    
    my $cdna_length = length($tr->seq->seq);
    
    my $count_ex = 0;
    my $count_lation = 0;
    my $tr_name = $tr->stable_id;
    my $tr_dbID = $tr->dbID;
    
    my @exons = @{$tr->get_all_Exons};
    my $nb_ex = scalar(@exons);
    
    my $strand = $exons[$nb_ex-1]->strand;
    
    foreach my $ex(@exons) {
	$count_ex++;
	
	my $ex_start = $ex->start;
	my $ex_end = $ex->end;
	
#one exon case
	if ($nb_ex == 1) {
	    if (($coding_start > 1)&& ($coding_end != $cdna_length)) {
		if ($strand == 1) {
		    print OUT "$ex_start\t$ex_end\tmRNA\n";
		}
		else {
		    print OUT "$ex_end\t$ex_start\tmRNA\n";
		}
	    }
	    
	    elsif (($coding_start == 1)&& ($coding_end != $cdna_length)) {
		if ($strand == 1) {
		    print OUT "<$ex_start\t$ex_end\tmRNA\n";
		}
		else {
		    print OUT "<$ex_end\t$ex_start\tmRNA\n";
		}
	    }
	    
	    elsif (($coding_start > 1)&& ($coding_end == $cdna_length)) {
		if ($strand == 1) {
		    print OUT "$ex_start\t>$ex_end\tmRNA\n";
		}
		else {
		    print OUT "$ex_end\t>$ex_start\tmRNA\n";
		}
	    }
	    
	    elsif (($coding_start == 1)&& ($coding_end == $cdna_length)) {
		if ($strand == 1) {
		    print OUT "<$ex_start\t>$ex_end\tmRNA\n";
		}
		else {
		    print OUT "<$ex_end\t>$ex_start\tmRNA\n";
		}
	    }
	}
	
#Multiple exons case	
	elsif ($nb_ex > 1) {
	    
#First exon    
	    if ($count_ex == 1) {
		#5' UTR
		if ($coding_start > 1) { 
		    if ($strand == 1) {
			print OUT "$ex_start\t$ex_end\tmRNA\n";
		    }
		    else {
			print OUT "$ex_end\t$ex_start\tmRNA\n";
		    }
		}
#no 5' UTR partial sequence
		elsif  ($coding_start == 1) {
		    if ($strand == 1) {
			print OUT "<$ex_start\t$ex_end\tmRNA\n";
		    }
		    else {
			print OUT "<$ex_end\t$ex_start\tmRNA\n";
		    }
		}
	    }
#End first exon
	    
#Last exon
	    elsif ($count_ex == $nb_ex) {
		if ($coding_end == $cdna_length)  {
		    if ($strand == 1) {
			print OUT "$ex_start\t>$ex_end\n";
		    }
		    else {
			print OUT "$ex_end\t>$ex_start\n";
		    }
		}
		elsif($coding_end != $cdna_length)  {
		    if ($strand == 1) {
			print OUT "$ex_start\t>$ex_end\n";
		    }
		    else {
			print OUT "$ex_end\t>$ex_start\n";
		    }
		}
	    }
#End last exon
	    
#other exons
	    else {
		if ($strand == 1) {
		    print OUT "$ex_start\t$ex_end\n";
		}
		else {
		    print OUT "$ex_end\t$ex_start\n";
		}
	    }
	}
    }
    print OUT "\t\t\tproduct\t$tr_name\n";
}

sub print_translation_coordinates {
    my ($tr,$db) = @_;
    
    my @exons = @{$tr->get_all_Exons};
    my $nb_ex = scalar(@exons);
    
    my $strand = $exons[$nb_ex-1]->strand;

    my $count_lation;
    my $translation = $tr->translation;
    my $translation_name = $translation->stable_id;
    
    my $seq = $tr->translate->seq;
    
    #    print STDERR "SEQ: $seq\n";
    my ($first) = $seq =~ /(^\S)/;
    my ($last) = $seq =~ /(\S$)/;
#	    print STDERR "FIRST: $first\tLAST: $last\n";
    my $tr_dbID = $translation->dbID;
    
    my ($cel_id,$old_cel_id,$ebi_id,$old_ebi_id,$symbol) = &fetch_2update($tr_dbID,$db);
    
    
    my @translateable = @{$tr->get_all_translateable_Exons};
    my $nb1 = scalar(@translateable);
    
    foreach my $lation (@translateable) {
	$count_lation++;
	
	my $tr_start = $lation->start;
	my $tr_end = $lation->end;
	
#only one exon
	if ($nb1 == 1) {
#both methionine and stop codon are included
	    if (($first =~ /M/)&&($last =~ /\*/)) {
#takes in account the strand			
		if ($strand == 1) {
		    print OUT "$tr_start\t$tr_end\tCDS\n";
		}
		else {
		    print OUT "$tr_end\t$tr_start\tCDS\n";
		}
	    }
#Only methionine included		    
	    elsif (($first =~ /M/) && ($last !~ /\*/)) {
		if ($strand == 1) {
		    print OUT "$tr_start\t>$tr_end\tCDS\n";	
		}
		else {
		    print OUT "$tr_end\t>$tr_start\tCDS\n";
		}
	    }
#Only stop codon included
	    elsif (($first !~ /M/) && ($last =~ /\*/)) {
		if ($strand == 1) {
		    print OUT "<$tr_start\t$tr_end\tCDS\n";	
		}
		else {
		    print OUT "<$tr_end\t$tr_start\tCDS\n";
		}
	    }
	    
#None feature is included
	    elsif (($first !~ /M/) && ($last !~ /\*/)) {
		if ($strand == 1) {
		    print OUT "<$tr_start\t>$tr_end\tCDS\n";
		}
		else {
		    print OUT "<$tr_end\t>$tr_start\tCDS\n";
		}
	    }
	    else {
		die;
	    }
	}
#End one exon case
	
	elsif ($nb1 > 1) {
#Multiple exons case
#First exon
	    if ($count_lation == 1) {
#Has a methionine
		if ($first =~ /M/) { 
		    if ($strand == 1) {
			print OUT "$tr_start\t$tr_end\tCDS\n";
		    }
		    else {
			print OUT "$tr_end\t$tr_start\tCDS\n";
		    }
		}
#Partial feature		    
		else {
		    if ($strand == 1) {
			print OUT "<$tr_start\t$tr_end\tCDS\n";
		    }
		    else {
			print OUT "<$tr_end\t$tr_start\tCDS\n";
		    }
		}
	    }
#end first exon case
	    
#last exon case
	    elsif ($count_lation == $nb1) {
#last exon partial feature (no stop)
		if ($first !~ /\*/) {
		    if ($strand == 1) {
			print OUT "$tr_start\t>$tr_end\n";
		    }
		    else {
			print OUT "$tr_end\t>$tr_start\n";
		    }
		}
#last exon has a stop
		else {
#takes in account the strand
		    if ($strand == 1) {
			print OUT "$tr_start\t$tr_end\n";
		    }
		    else {
			print OUT "$tr_end\t$tr_start\n";
		    }
		}
	    }
#end of last exon case
	    
#other exons
	    else {
		if ($strand == 1) {
		    print OUT "$tr_start\t$tr_end\n";
		}
		else {
		    print OUT "$tr_end\t$tr_start\n";
		}
	    }
	}
    }
    print OUT "\t\t\tproduct\t$translation_name\n";
    
    if ($cel_id) {
	print OUT "\t\t\tprotein_id\t$map{$cel_id}\n";
	print OUT "\t\t\tprotein_id\t$cel_id\n";
	print OUT "\t\t\tprotein_id\t$old_cel_id\n";
    }
    if ($ebi_id) {
	print OUT "\t\t\tprotein_id\t$map{$ebi_id}\n";
	print OUT "\t\t\tprotein_id\t$ebi_id\n";
	print OUT "\t\t\tprotein_id\t$old_ebi_id\n";
    }
    if ($symbol) {
	print OUT "\t\t\tprotein_id\t$symbol\n";
    }
    my @interpro = &get_protein_annotation($tr_dbID,$db);

    foreach my $ipr(@interpro) {
	print OUT "\t\t\tinterpro\t$ipr\n";
    }
    print OUT "\t\t\tevidence\tnot_experimental\n";
    
}


sub fetch_2update {
    my ($tr_dbID,$db) = @_;
    
    my $query = "select x.display_label from xref x, external_db e, object_xref o where o.ensembl_id = $tr_dbID and o.xref_id = x.xref_id and x.external_db_id = e.external_db_id and e.db_name = 'Celera_Pep'";
    my $sth = $db->prepare($query);
    $sth->execute();
    
    my $cel_id = $sth->fetchrow;
    my $old_cel_id = "gnl|WGS:AAAB|".$cel_id;

     my $query1 = "select x.display_label from xref x, external_db e, object_xref o where o.ensembl_id = $tr_dbID and o.xref_id = x.xref_id and x.external_db_id = e.external_db_id and e.db_name = 'anopheles_paper'";
    
    my $sth1 = $db->prepare($query1);
    $sth1->execute();
    
    my $ebi_id = $sth1->fetchrow;
    
    if ($ebi_id) {
	my ($nid) = $ebi_id =~ /(\d+)$/;
	$ebi_id = "ebiP".int($nid);
    }

#    print STDERR "EBI: $ebi_id\n";

    my $old_ebi_id = "gnl|WGS:AAAB|".$ebi_id;
    
     my $query2 = "select x.display_label from xref x, external_db e, object_xref o where o.ensembl_id = $tr_dbID and o.xref_id = x.xref_id and x.external_db_id = e.external_db_id and e.db_name = 'anophele_symbol'";
    
    my $sth2 = $db->prepare($query2);
    $sth2->execute();
    
    my $symbol = $sth2->fetchrow;

    return ($cel_id,$old_cel_id,$ebi_id,$old_ebi_id,$symbol);
    
}

sub get_protein_annotation {
    my ($tr_dbID,$db) = @_;
    my @interpro;
    
    my $query = "select distinct(i.interpro_ac) from interpro i, protein_feature pf where pf.hit_id = i.id and translation_id = $tr_dbID";
    
    my $sth = $db->prepare($query);
    $sth->execute;
    while (my $ipr = $sth->fetchrow) {
	push (@interpro,$ipr);
    }
    return @interpro;
}
