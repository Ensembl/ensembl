
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

my %transcript_map;
my %utr;


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

open (OUT,">/Users/emmanuelmongin/work/ncbi_dump/test.out");

#Get all of the scaffolds
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

#	print STDERR "ID: $gene_id\n";
	&checks($gene_id);
	
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

	if ($gene->strand == 1) {
	    if (($utr{$gene_id}->{'up'})&&($utr{$gene_id}->{'down'})) {
		print OUT "$start\t$end\tgene\n";
	    }
	    elsif ((!defined $utr{$gene_id}->{'up'})&&($utr{$gene_id}->{'down'})) {
		print OUT "<$start\t$end\tgene\n";
	    }
	    elsif (($utr{$gene_id}->{'up'})&&(! defined$utr{$gene_id}->{'down'})) {
		print OUT "$start\t>$end\tgene\n";
	    }
	    elsif ((! defined $utr{$gene_id}->{'up'})&&(! defined$utr{$gene_id}->{'down'})) {
		print OUT "<$start\t>$end\tgene\n";
	    }
	}
	else {
	   if (($utr{$gene_id}->{'up'})&&($utr{$gene_id}->{'down'})) {
		print OUT "$end\t$start\tgene\n";
	    }
	    elsif ((!defined $utr{$gene_id}->{'up'})&&($utr{$gene_id}->{'down'})) {
		print OUT "<$end\t$start\tgene\n";
	    }
	    elsif (($utr{$gene_id}->{'up'})&&(! defined$utr{$gene_id}->{'down'})) {
		print OUT "$end\t>$start\tgene\n";
	    }
	    elsif ((! defined $utr{$gene_id}->{'up'})&&(! defined$utr{$gene_id}->{'down'})) {
		print OUT "<$end\t>$start\tgene\n";
	    }  
	}
	print OUT "\t\t\tlocus_tag\t$name\n";
	
	my @transcripts = @{$gene->get_all_Transcripts};

	foreach my $tr(@transcripts) {
	    
	    my $coding_start = $tr->cdna_coding_start;
	    my $coding_end = $tr->cdna_coding_end;
	    
	    my $cdna_length = length($tr->seq->seq);
	    
#	    print OUT "START: $coding_start\n";
#	    print OUT "END: $coding_end\t$cdna_length\n";

	    my $count_ex = 0;
	    my $count_lation = 0;
	    my $tr_name = $tr->stable_id;
	    my $tr_dbID = $tr->dbID;
	    
	    my @exons = @{$tr->get_all_Exons};
	    my $nb_ex = scalar(@exons);
	    
	    my $strand = $exons[$nb_ex-1]->strand;
	    print OUT "STRAND: $strand\n";
	    
	    foreach my $ex(@exons) {
		$count_ex++;
		
		my $ex_start = $ex->start;
		my $ex_end = $ex->end;

#one exon case
		if ($nb_ex == 1) {
#UTRs on both sides
		    if (($coding_start > 1)&& ($coding_end != $cdna_length)) {
			if ($strand == 1) {
			    if ($transcript_map{$tr_dbID} eq 'no') {
				print OUT "$ex_start\t$ex_end\tmRNA\n";
			    }
			    else {
				print OUT "$ex_start\t$ex_end+3\tmRNA\n";
			    }
			}
			else {
			    if ($transcript_map{$tr_dbID} eq 'no') {
				print OUT "$ex_end\t$ex_start\tmRNA\n";
			    }
			    else {
				print OUT "$ex_end\t$ex_start-3\tmRNA\n";
			    }
			}
		    }
#no 5' UTR
		    elsif (($coding_start > 1)&& ($coding_end == $cdna_length)) {
			if ($strand == 1) {
			     if ($transcript_map{$tr_dbID} eq 'no') {
				 print OUT "$ex_start\t>$ex_end\tmRNA\n";
			     }
			     else {
				 print OUT "$ex_start\t>$ex_end+3\tmRNA\n";
			     }
			}
			else {
			    if ($transcript_map{$tr_dbID} eq 'no') {
				print OUT "$ex_end\t>$ex_start\tmRNA\n";
			    }
			    else {
				print OUT "$ex_end\t>$ex_start-3\tmRNA\n";
			    }
			}
		    }
#no 3' UTR		    
		    elsif (($coding_start == 1)&& ($coding_end != $cdna_length)) {
			if ($strand == 1) {
			    if ($transcript_map{$tr_dbID} eq 'no') {
				print OUT "<$ex_start\t$ex_end\tmRNA\n";
			    }
			    else {
				print OUT "<$ex_start\t$ex_end+3\tmRNA\n";
			    }
			}
			else {
			    if ($transcript_map{$tr_dbID} eq 'no') { 
				print OUT "<$ex_end\t$ex_start\tmRNA\n";
			    }
			    else {
				print OUT "<$ex_end\t$ex_start-3\tmRNA\n";
			    }
			}
			
		    }
		    else {
			if ($strand == 1) {
			    if ($transcript_map{$tr_dbID} eq 'no') {
				print OUT "<$ex_start\t>$ex_end\tmRNA\n";
			    }
			    else {
				print OUT "<$ex_start\t>$ex_end+3\tmRNA\n";
			    }
			}
			else {
			    if ($transcript_map{$tr_dbID} eq 'no') {
				print OUT "<$ex_end\t>$ex_start\tmRNA\n";
			    }
			    else {
				print OUT "<$ex_end\t>$ex_start-3\tmRNA\n";
			    }
 
			}
		    }
		}
#end of 1 exon case
	
#multiple exon case
#first exon
		elsif ($count_ex == 1) {
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

#last exon case		
		elsif ($count_ex == $nb_ex) {
#no 3' UTR-> partial sequence		    
		if ($coding_end == $cdna_length)  {
		    if ($strand == 1) {
			if ($transcript_map{$tr_dbID} eq 'no') {
			    print OUT "$ex_start\t>$ex_end\n";
			}
			else {
			    print OUT "$ex_start\t>$ex_end+3\n";
			}
		    }
		    else {
			if ($transcript_map{$tr_dbID} eq 'no') {
			    print OUT "$ex_end\t>$ex_start\n";
			}
			else {
			    print OUT "$ex_end\t>$ex_start-3\n";
			}
		    }
		}
#UTR 3' end
		else {
		    if ($strand == 1) {
			if ($transcript_map{$tr_dbID} eq 'no') {
			    print OUT "$ex_start\t$ex_end\n";
			}
			else {
			    print OUT "$ex_start\t$ex_end+3\n";
			}
		    }
		    else {
			if ($transcript_map{$tr_dbID} eq 'no') {
			    print OUT "$ex_end\t$ex_start\n";
			}
			else {
			    print OUT "$ex_end\t$ex_start-3\n";
			}
		    }
		}
	    }
#end last exon case

		else {
		    if ($strand == 1) {
			print OUT "$ex_start\t$ex_end\n";
		    }
		    else {
			print OUT "$ex_end\t$ex_start\n";
		    }
		}


	    }

	    
	    print OUT "\t\t\tproduct\t$tr_name\n";
	    
	    my $translation = $tr->translation;
	    my $translation_name = $translation->stable_id;

	    my $seq = $tr->translate->seq;
	    
	#    print STDERR "SEQ: $seq\n";
	    my ($first) = $seq =~ /(^\S)/;
	    my ($last) = $seq =~ /(\S$)/;
#	    print STDERR "FIRST: $first\tLAST: $last\n";
	    my $tr_dbID = $translation->dbID;

	    my $query3 = "select x.display_label from xref x, external_db e, object_xref o where o.ensembl_id = $tr_dbID and o.xref_id = x.xref_id and x.external_db_id = e.external_db_id and (e.db_name = 'Celera_Pep' or e.db_name = 'Anopheles_paper')";
		my $sth3 = $db->prepare($query3);
		$sth3->execute();
		
		my $protein_id = $sth3->fetchrow;
	    
	    my @translateable = @{$tr->get_all_translateable_Exons};
	    my $nb1 = scalar(@translateable);

	    foreach my $lation (@translateable) {
		$count_lation++;
		
		my $tr_start = $lation->start;
		my $tr_end = $lation->end;
	
#only one exon
		if ($nb1 == 1) {
#both methionine and stop codon are included
		    if (($first =~ /M/) && ($transcript_map{$tr_dbID} eq 'after')) {
#takes in account the strand			
			if ($strand == 1) {
			    print OUT "$tr_start\t$tr_end+3\tCDS\n";
			}
			else {
			     print OUT "$tr_end-3\t$tr_start\tCDS\n";
			}
		    }
#Only methionine included		    
		    elsif (($first =~ /M/) && ($transcript_map{$tr_dbID} eq 'no')) {
			if ($strand == 1) {
			    print OUT "$tr_start\t>$tr_end\tCDS\n";	
			}
			else {
			    print OUT "$tr_end\t>$tr_start\tCDS\n";
			}
		    }
#None feature is included
		    else {
			if ($strand == 1) {
			    print OUT "<$tr_start\t>$tr_end\tCDS\n";
			}
			else {
			    print OUT "<$tr_end\t>$tr_start\tCDS\n";
			}
		    }
		    
		}
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
		    if ($transcript_map{$tr_dbID} eq 'no') {
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
			    print OUT "$tr_start\t$tr_end+3\n";
			}
			else {
			    print OUT "$tr_end\t$tr_start-3\n";
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
	    print OUT "\t\t\tproduct\t$translation_name\n";
	    if ($protein_id) {
		print OUT "\t\t\tprotein_id\t$protein_id\n";
	    }
	    print OUT "\t\t\tevidence\tnot_experimental\n";
	}
    }
}


sub checks {
    my ($gene_id) = @_;
        
    my $gene = $gene_adapt->fetch_by_dbID($gene_id,1);
    
    foreach my $tr(@{$gene->get_all_Transcripts}) { 
	
	my $tr_dbID = $tr->dbID;
	my $spl_sq = $tr->spliced_seq;
	my $c_start = $tr->cdna_coding_start;
	my $c_end = $tr->cdna_coding_end;
	
	my $start = $tr->start;
	my $end = $tr->end;

	my $pep = $tr->translate->seq;
	
#check if at least a transcript of a gene has 5' utr
	if ($c_start > 1) {
	 $utr{$gene_id}->{'up'} = 1;  
	}
	#print STDERR $tr->stable_id()," ";
	if ($pep =~ /\*$/) {
	    $transcript_map{$tr_dbID} = "included";
	    #print STDERR "STOP CODON INCLUDED\n";
	}
	
	else {
	    if ($c_end < length($spl_sq)) {
		$utr{$gene_id}->{'down'} = 1;
		#print STDERR "UTR - ";
		my $codon = substr($spl_sq,$c_end,3);
		if ($codon =~ /TAG|TGA|TAA/i) {
		    $transcript_map{$tr_dbID} = "after";
		    #print STDERR "STOP AFTER TRANSCRIPT\n";
		}
		else {
		    $transcript_map{$tr_dbID} = "no";
		    #print STDERR "NOSTOP AFTER TRANSCRIPT\n";
		}
	    }
	    else {
		#print STDERR "NO UTR - ";
		my $end_exon = $tr->end_Exon;
		my $strand = $end_exon->strand;
		my $last_codon;
		if ($strand == 1) {
		    $last_codon = $end_exon->contig->subseq($end_exon->end+1,$end_exon->end+3,1);
		}
		
		else {
		    $last_codon = $end_exon->contig->subseq($end_exon->start-3,$end_exon->start-1,-1);
		}
		
		if ($last_codon =~ /TAG|TGA|TAA/i) {
		    #print STDERR "STOP AFTER TRANSCRIPT\n";
		     $transcript_map{$tr_dbID} = "after";
		}
		else {
		    #print STDERR "NOSTOP AFTER TRANSCRIPT\n";
		     $transcript_map{$tr_dbID} = "no";
		}
		
	    }
	    
	}
    }
}
