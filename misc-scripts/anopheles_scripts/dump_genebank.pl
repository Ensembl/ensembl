use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::ProteinAdaptor;
use Bio::EnsEMBL::Exon;


my $host      = 'ecs2d';
my $dbuser    = 'ensro';
my $dbname    = 'anopheles_gambiae_core_17_2a';
my $dbpass    = '';

my %scafmap;
my %utr;
my %utr_tmp;
my %ebimap;

my %map;


GetOptions(
	   'dbname:s'    => \$dbname,
	   'dbhost:s'    => \$host,
	   'dbuser:s'    => \$dbuser,
	   'dbpass:s'    => \$dbpass,
	   );

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


open (MAP,"/acari/work1/mongin/anopheles_17_2a_ncbisub/input/AAAB01.output.p2g") || die;

open (SCAFMAP,"/acari/work1/mongin/anopheles_17_2a_ncbisub/input/accessions") || die;
open (EBIMAP,"/acari/work1/mongin/anopheles_17_2a_ncbisub/input/ebi_id_mapping.txt") || die;

while(<EBIMAP>) {
    chomp;
    my ($tr,$ac) = split;
    $ebimap{$tr} = $ac;
}

while(<SCAFMAP>) {
    chomp;
    my ($new,$a,$old) = split;
    $scafmap{$new} = $old;
}

while(<MAP>) {
    chomp;
    my @array = split;
    $map{$array[0]} = $array[1];
}

close(MAP);


#Get all of the distinct scaffolds
my $query1 = "select distinct(c.name), a.superctg_ori from clone c, assembly a where a.superctg_name = c.name";

my $sth1 = $db->prepare($query1);
$sth1->execute();

while (my ($clone_name,$ori) = $sth1->fetchrow_array) {

    open (OUT,">/acari/work1/mongin/anopheles_17_2a_ncbisub/output/$clone_name.tbl") || die;
    open (SEQ,">/acari/work1/mongin/anopheles_17_2a_ncbisub/output/$clone_name.fsa") || die;

    my $slice = $slice_adapt->fetch_by_clone_accession($clone_name);
    
    if ($ori == -1) {
	$slice = $slice->invert();
    }
    
    print STDERR "CLONE: $clone_name\n";

    my $chr_name = $slice->chr_name;
    
    my $clone_seq = $slice->seq;

    my $length_clone = length($clone_seq);

    my $old_clone_name = $scafmap{$clone_name};

    if ($chr_name !~ /UNKN/) {
	print SEQ ">gnl|WGS:AAAB|$old_clone_name|gb|$clone_name [organism=Anopheles gambiae str. PEST] [tech=wgs] [chromosome=$chr_name]\n$clone_seq\n";
    }
    else {
	print STDERR "HERE\n";
	print SEQ ">gnl|WGS:AAAB|$old_clone_name|gb|$clone_name [organism=Anopheles gambiae str. PEST] [tech=wgs]\n$clone_seq\n";
    }
    print OUT ">Feature gnl|WGS:AAAB|$old_clone_name|gb|$clone_name\n";
    
    my @genes = @{$slice->get_all_Genes};
    
   

    foreach my $gene(@genes) {

	 print STDERR "ID: ".$gene->stable_id."\tSTRAND: ".$gene->strand."\tSTART: ".$gene->start."\tEND: ".$gene->end."\n";

	if (($gene->start < 0) || ($gene->start > $length_clone) || ($gene->end < 0) || ($gene->end > $length_clone)) {
	    next;
	}
	 
	 else {
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
	    
#Check for the gene the presence of UTRs\
    	    
	    if (($utr{$new_gene_dbID}->{'up'}==1)&&($utr{$new_gene_dbID}->{'down'}==1)) {
		if ($gene->strand == 1) {
		    print OUT "$start\t$end\tgene\n";
		}
		else {
		    print OUT "$end\t$start\tgene\n";
		}
	    }
	    
	    elsif (($utr{$new_gene_dbID}->{'up'} == 1)&&($utr{$new_gene_dbID}->{'down'}!=1)) {
		if ($gene->strand == 1) {
		    print OUT "$start\t>$end\tgene\n";
		}
		else {
		    print OUT "$end\t>$start\tgene\n";
		}
	    }
	    
	    elsif (($utr{$new_gene_dbID}->{'up'} != 1)&&($utr{$new_gene_dbID}->{'down'}==1)) {
		if ($gene->strand == 1) {
		    print OUT "<$start\t$end\tgene\n";
		}
		else {
		    print OUT "<$end\t$start\tgene\n";
		}
	    }
	    
	    elsif (($utr{$new_gene_dbID}->{'up'} != 1)&&($utr{$new_gene_dbID}->{'down'}!=1)) {
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

    close(OUT);
    close(SEQ);
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
	my $spl_seq;
	#eval {
	$spl_seq = $tr->spliced_seq;
	 #   };
	#if ($@){
	#    next;
	#}

	my $tl_start = $tr->translation->start;
	my $tl_end = $tr->translation->end;

	my $pep = $tr->translate->seq;
	
	my $new_cdna_start;
	my $new_cdna_end;
	    
	if ($pep =~ /^M/) {
	    $tl_start = $c_start;
	    $new_cdna_start = 1;
	} else {
	    $tl_start = 1;
	    $new_cdna_start = $c_start;
	}

	$tl_end = $c_end - $c_start + $tl_start;
	if ($pep =~ /\*$/) {
	    $new_cdna_end = length( $spl_seq );
	} else {
	    $new_cdna_end = $c_end;
	}

	my @exon_coords = $tr->cdna2genomic($new_cdna_start,$new_cdna_end);

	my @new_exons;

	for my $exon_coord ( @exon_coords ) {
	    if ($exon_coord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
		
#		print STDERR $exon_coord->start."\t".$exon_coord->end."\t".$exon_coord->strand."\t".$exon_coord->id."\n";
		
		my $new_exon = new Bio::EnsEMBL::Exon($exon_coord->start,$exon_coord->end,$exon_coord->strand);
		$new_exon->contig($exon_coord->id);
		$new_exon->phase(0);
		$new_exon->end_phase(0);
		
		push (@new_exons,$new_exon);
		#print STDERR "NEW EX: $new_exon\n";
	    }
	    
	}

	  
	my $tl_start_exon;
	my $tl_end_exon;
	my $seen;
	
	#print STDERR "TRANSLATION: $tl_start\t$tl_end\n";

	foreach my $exon  (@new_exons ) {
	    if( ($tl_start > $exon->length) && (! defined $seen) ) {
	#	print STDERR "Translation_start: $tl_start\n";
		$tl_start -= $exon->length;
	    } elsif (! defined $seen) {
		$seen = 1;
		$tl_start_exon = $exon;
	    }
	    if(  $tl_end > $exon->length) {
		$tl_end -= $exon->length;
	    } else {
		$tl_end_exon = $exon;
		last;
	    }
	}

	$tr->flush_Exons();
		
	foreach my $ne (@new_exons) {
	    $tr->add_Exon($ne);
	    
	}
	
	$tr->{'cdna_coding_start'} =undef;
	$tr->{'cdna_coding_end'} =undef;

	#$tr->cdna_coding_start($new_cdna_start);
	#$tr->cdna_coding_end($new_cdna_end);
	
#	print STDERR "CDNA start: ".$tr->cdna_coding_end."\t".$tl_start."\n";

	#print STDERR "EX: ".$tl_end_exon."\n";

	my $translation = $tr->translation;
	$translation->start_Exon($tl_start_exon);
	$translation->end_Exon($tl_end_exon);
	$translation->start($tl_start);
	$translation->end($tl_end);

	
	$tr->{'translation'} = [];
	$tr->translation($translation);
	#print STDERR "TR: .".$tr->translation->start."\n";

	push(@new_transcripts,$tr);
	
    }


    $gene->{'_transcript_array'} =[];
	
    foreach my $new_tr(@new_transcripts) {
	my $cdna_length = length($new_tr->seq->seq);
	my $coding_start = $new_tr->cdna_coding_start;
	my $coding_end = $new_tr->cdna_coding_end;
	my $gene_dbid = $gene->dbID;
	my $tl_start = $new_tr->translation->start;
#	print STDERR "Coding start: $coding_start\n";

	#print STDERR "SEQ: ".$new_tr->translate->seq."\n";

	#my $coding_length = $coding_end - $coding_start + 1;
	
	#print STDERR "Coding start: $coding_start\tcoding end: $coding_end\tLength: $cdna_length\n";

	if ($tl_start > 1) {
	    $utr{$gene_dbid}->{'up'} = 1;
	}
	
	 if ($coding_end != $cdna_length) {
	     $utr{$gene_dbid}->{'down'} = 1;
	 }
	    
	$gene->add_Transcript($new_tr);
    }
    
    return $gene;

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
	#print OUT "$coding_start\t$coding_end\t$cdna_length\n";
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
			print OUT "$ex_start\t$ex_end\n";
		    }
		    else {
			print OUT "$ex_end\t$ex_start\n";
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
    
    #   print STDERR "SEQ: $seq\n";
    my ($first) = $seq =~ /(^\S)/;
    my ($last) = $seq =~ /(\S$)/;
    #print OUT "FIRST: $first\tLAST: $last\n";
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
		if ($last !~ /\*/) {
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
	print OUT "\t\t\tprotein_id\tgnl|WGS:AAAB|$cel_id|gb|$map{$cel_id}\n";
	print OUT "\t\t\tprot_desc\tgnl|WGS:AAAB|$cel_id|gb|$map{$cel_id}\n";
    }
    if ($ebi_id) {
	print OUT "\t\t\tprotein_id\tgnl|WGS:AAAB|$ebi_id|gb|$map{$ebi_id}\n";
	print OUT "\t\t\tprot_desc\tgnl|WGS:AAAB|$ebi_id|gb|$map{$ebi_id}\n";
    }
			   if ((!defined $cel_id)&&(!defined $ebi_id)) {
			       print OUT "\t\t\tprotein_id\tgnl|WGS:AAAB|$translation_name\n";
			   }
    if ($symbol) {
	print OUT "\t\t\tgene\t$symbol\n";
    }
    my @interpro = &get_protein_annotation($tr_dbID,$db);

    foreach my $ipr(@interpro) {
	print OUT "\t\t\tdb_xref\tInterpro:$ipr\n";
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

#    my $query1 = "select ts.stable_id from transcript_stable_id ts, xref x, external_db e, object_xref o, transcript t where o.ensembl_id = t.translation_id and t.transcript_id = ts.transcript_id and o.xref_id = x.xref_id and x.external_db_id = e.external_db_id and e.db_name = 'anopheles_paper' and o.ensembl_id = $tr_dbID";

    my $query1 = "select tr.stable_id from transcript_stable_id tr, transcript t, translation_stable_id ts where ts.translation_id = '$tr_dbID' and ts.translation_id = t.translation_id and t.transcript_id = tr.transcript_id";

    my $sth1 = $db->prepare($query1);
    $sth1->execute();

    my $transcript_id = $sth1->fetchrow;
    
    

    my $new_ebi_id = $ebimap{$transcript_id};;
    
    
    #if($ebi_id) {
#	my ($nid) = $ebi_id =~ /(\d+)$/;
#	$new_ebi_id = "ebiP".int($nid);
#    }
    

#    print STDERR "EBI: $ebi_id\n";

    my $old_ebi_id = "gnl|WGS:AAAB|".$new_ebi_id;
    
     my $query2 = "select x.display_label from xref x, external_db e, object_xref o where o.ensembl_id = $tr_dbID and o.xref_id = x.xref_id and x.external_db_id = e.external_db_id and e.db_name = 'anophele_symbol'";
    
    my $sth2 = $db->prepare($query2);
    $sth2->execute();
    
    my $symbol = $sth2->fetchrow;

    return ($cel_id,$old_cel_id,$new_ebi_id,$old_ebi_id,$symbol);
    
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






















