#!/usr/local/ensembl/bin/perl

=head1 NAME

  dump_translations.pl

=head1 SYNOPSIS
 
  dump_translations.pl

=head1 DESCRIPTION

dump_translations.pl dumps out the translations of all the genes in a database specified in GeneConf.pm
It\'s a stripped down version of gene2flat.

=head1 OPTIONS

=cut

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;

my $host = 'ecs2d';
my $dbuser = 'ensro';
my $dbname = 'homo_sapiens_core_10_30';
my $dbpass;
my $path;
my $port;


my $type = 'ensembl';

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $host,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					    '-port'   => $port
					    );


print STDERR "Connecting to $host, $dbname\n";

#my $query = "select t.transcript_id,t.translation_id from transcript t, gene g where t.gene_id = g.gene_id and g.type = '$type'"; 
#my $sth = $db->prepare($query);
#$sth->execute();

my $gene_adaptor = $db->get_GeneAdaptor;

my $gene_id_list = $gene_adaptor->list_geneIds;



while(my $gene_id = shift(@$gene_id_list)) {
    
    #my ($transc_dbid, $transl_dbid) = @array;
    
    my $gene = $gene_adaptor->fetch_by_dbID($gene_id,1);
    
    foreach my $tr(@{$gene->get_all_Transcripts}) { 
	
	my $spl_sq = $tr->spliced_seq;
	my $c_start = $tr->cdna_coding_start;
	my $c_end = $tr->cdna_coding_end;
	
	my $pep = $tr->translate->seq;
	
	print STDERR $tr->stable_id()," ";
	if ($pep =~ /\*$/) {
	    print STDERR "STOP CODON INCLUDED\n";
	}
	
	else {
	    if ($c_end < length($spl_sq)) {
		print STDERR "UTR - ";
		my $codon = substr($spl_sq,$c_end,3);
		if ($codon =~ /TAG|TGA|TAA/i) {
		    print STDERR "STOP AFTER TRANSCRIPT\n";
		}
		else {
		    print STDERR "NOSTOP AFTER TRANSCRIPT\n";
		}
	    }
	    else {
		print STDERR "NO UTR - ";
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
		    print STDERR "STOP AFTER TRANSCRIPT\n";
		}
		else {
		    print STDERR "NOSTOP AFTER TRANSCRIPT\n";
		}
		
	    }
	    
	}
    }

}




