#!/usr/local/bin/perl

=head1 NAME - gtf dump

    Parses files in GTF format (Gene Transfer Format)

=head1 SYNOPSIS - 

    gtfparse -dbname ensembl -parsefile genes.gtf

=head1 DESCRIPTION

    This script parses GTF files and writes the genes extracted to a database.
    The database is specified using the usual EnsEMBL options, described below.

    The actual parsing happens in the Bio::EnsEMBL::Utils::GTF_handler module,
    which also handles the dumping of GTF files.

    If the print option is specified, then the genes are not written to db, 
    but printed to STDOUT (mainly for testing)

=head1 OPTIONS

    -dbtype    database type (only needed for TimDB)

    -host      host name for the database (gets put as host= in locator) 

    -port      for RDBs, what port to connect to for the "to" database (port= in locator)

    -dbname    for RDBs, what database name to connect to (dbname= in locator)

    -dbuser    for RDBs, what username to connect as to the database (dbuser= in locator)

    -dbpass    for RDBs, what password to use to connect to to the database (dbpass= in locator)

    -module    module name to load to (Defaults to Bio::EnsEMBL::DBSQL::Obj)

    -parsefile name of the GTF file to parse

    -print     prints gene structures to STDOUT

    -help      displays this documentation with PERLDOC

=cut

use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::Utils::GTF_handler;
use Bio::EnsEMBL::DBLoader;
use strict;
use Getopt::Long;
use Bio::EnsEMBL::GeneComparison::GeneComparisonStats;

#Database options
my $dbtype = 'rdb';
my $host   = 'ensrv4.sanger.ac.uk';
my $port   = '410000';
my $dbname = 'ensembl_freeze17_michele';
my $dbuser = 'ensro';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::Obj';

#Other options
my $parse;
my $parse2;
my $print;
my $help;
my $check;
my $compare;
my $longest;

&GetOptions( 
	     'dbtype:s'   => \$dbtype,
	     'host:s'     => \$host,
	     'port:n'     => \$port,
	     'dbname:s'   => \$dbname, 
	     'dbuser:s'   => \$dbuser,
	     'dbpass:s'   => \$dbpass,
	     'module:s'   => \$module,
       	     'parse:s'    => \$parse,
	     'parse2:s'   => \$parse2,
	     'print'      => \$print,
	     'check'      => \$check,
	     'compare'    => \$compare,
	     'longest'    => \$longest,
	     'h|help'     => \$help
	     );

my $gtfh=Bio::EnsEMBL::Utils::GTF_handler->new();
open (PARSE,"$parse") || die("Could not open $parse for gtf reading$!");
my @gtf_genes=$gtfh->parse_file(\*PARSE);
my $g_n=scalar @gtf_genes;
print STDERR "Got $g_n genes from file $parse\n";
if ($print) {
     foreach my $gene (@gtf_genes) {
	print STDOUT "Gene ".$gene->id."\n";
	foreach my $trans ($gene->each_Transcript) {
	    print STDOUT "  Transcript ".$trans->id."\n";
	    foreach my $exon ($gene->each_unique_Exon) {
		print STDOUT "   exon ".$exon->id."\n";
	    }
	    print STDOUT "   translation start ".$trans->translation->start."\n";
	    print STDOUT "   translation end ".$trans->translation->end."\n";
	}
    }
 }

#DB writing option not yet implemented
#Mapping of coordinates still needs to be done
elsif ($check) {
    my $inputstream = Bio::SeqIO->new(-file => "ctg12382.fa",-format => 'Fasta');
    my $seq = $inputstream->next_seq();

    foreach my $gene (@gtf_genes) {
	foreach my $trans ($gene->each_Transcript) {
	    
	    print STDERR "Translation start is ".$trans->translation->start." in exon ".$trans->translation->start_exon_id."\n";
	    print STDERR "Translation end is ".$trans->translation->end." in exon ".$trans->translation->end_exon_id."\n";
		
	    foreach my $exon ($trans->each_Exon) {
		my $start=$exon->start;
		my $end=$exon->end;
		$exon->attach_seq($seq);
		my $eseq=$exon->seq;
		
		my $end=$exon->end+2;
		my $m_end=$end-1;
		my $start=$exon->start-2;
		my $m_start=$start+1;
		my $subseq=$seq->subseq($start,$m_start);
		my $subseq2=$seq->subseq($m_end,$end);
		print STDERR "The start of exon ".$exon->id." is ".$exon->start."\n";
		print STDERR "The sequence from $start to $m_start is $subseq\n";
		print STDERR "The end of exon ".$exon->id." is ".$exon->end."\n";
		print STDERR "The sequence from $m_end to $end is $subseq2\n";
	    }
	    #foreach my $trans ($gene->each_Transcript) {
	    #my $tseq=$trans->translate();
	    #print STDERR "Sequence direct from fa file:\n";
	    #print STDERR $tseq->seq."\n";
	    #}
	}
    }
    
}
elsif ($compare) {
    $parse =~/ctg(\d+).gtf/;
    my $parse2 = "../neomorphic_known/ctg$1.mrg.gtf";
    my $gtfh=Bio::EnsEMBL::Utils::GTF_handler->new();
    open (PARSE2,"$parse2") || die("Could not open $parse2 for ensembl gtf reading$!");
    
    my @standard_genes=$gtfh->parse_file(\*PARSE2);
    my $n=scalar @standard_genes;
    print STDERR "Got $n genes from Neomoprhic file $parse2\n";
    #foreach my $gene (@standard_genes) {
	#foreach my $exon ($gene->each_unique_Exon) {
	    #print STDERR "Neomorphic gene has exon: ".$exon->id."\n";
	#}
    #}
    my $stats=Bio::EnsEMBL::GeneComparison::GeneComparisonStats->new(-standard=>\@standard_genes, -predictor=>\@gtf_genes);
 my %genes=$stats->get_OverlapMap;
    my $knownoverlap=0;
    foreach my $gene_id (keys(%genes)){
	my $n=scalar(@{$genes{$gene_id}});
	
	if ($n) {
	    print STDERR "$gene_id overlaps $n gene(s): ";
	    $knownoverlap++;
	    foreach my $predictor (@{$genes{$gene_id}}) {
		print STDERR "$predictor ";
	    }
	    print STDERR "\n";
	}
    }
    print STDERR "In total, $knownoverlap known genes overlap ensembl genes\n";
    my @missed_genes=$stats->getMissedGenes();
    my $missed= scalar(@missed_genes);
    print STDERR "$missed known genes were missed completely: ";
    foreach my $id (@missed_genes) {
	print STDERR "$id ";
    }
    print STDERR "\n";
    my @matched_ensembl= $stats->getExactlyMatchedPredictorGenes();
    my $n_mens= scalar @matched_ensembl;
    print STDERR "$n_mens ensembl predictions were exactly matched: ";
    foreach my $id (@matched_ensembl) {
	print STDERR "Matched ensembl gene $id ";
    }
    print STDERR "\n";
    my @matched_known=$stats->getExactlyMatchedStandardGenes();
    my $n_mknown=scalar(@matched_known);
    print STDERR "$n_mknown known genes were exactly matched: ";
    foreach my $id (@matched_known) {
	print STDERR "$id ";
    }
    print STDERR "\n";
    my %overlaps = $stats->getOverlapStats;
    foreach my $overlap (keys %overlaps) {
	print STDERR "Gene overlap: ". $overlap. " Proportion of bases ".
	    $overlaps{$overlap}. "\n";                    
	
    }
    #$stats->make_merges;
    
    my $string=$stats->getGeneComparisonStats;
    print STDERR "Gene comparison stats:\n";
    print STDERR $string;
   
}

elsif ($longest) {
    my (%ens_long,%neo_long);
    my %length;
    open (LENGTH_FILE,"ens.length") || die("Could not open ens.length for ensembl pep length reading$!");
    while (<LENGTH_FILE>) {
	if (/(\w+\.\w+)\t(\d+)/) {
	    my $length=$2;
	    my $pep=$1;
	    $pep =~ s/TMPP\_/SEPT20T\./g;
	    $length{$pep}=$length;
	}
    }
    open (LENGTH_FILE,"neo.length") || die("Could not open neo.length for ensembl pep length reading$!");
    while (<LENGTH_FILE>) {
	if (/(.+)\t(\d+)/) {
	    $length{$1}=$2;
	}
    }
    foreach my $gene (@gtf_genes) {
	my $longest=0;
	my $start;
	my $end;
	my $fpc;
	my $trans_id;
	my $gene_id=$gene->id;
	print STDERR "Analysing gene ".$gene->id."\n";
	foreach my $trans ($gene->each_Transcript) {
	    
	    my $pep_length=$length{$trans->id};
	    print STDERR "           transcript ".$trans->id." length: $pep_length\n";
	    if ($pep_length > $longest) {
		$longest=$pep_length;
		$trans_id=$trans->id;
	    }
	    $start=$trans->first_exon->start;
	    $end=$trans->last_exon->end;
	    $fpc=$trans->first_exon->contig_id;
	}
	print STDERR "Longest transcript is $trans_id, and is $longest long\n";

	if ($trans_id =~ /TMPP/) {
	    my $type='ENS';
	    $ens_long{$trans_id}=[$gene_id,$type,$fpc,$start,$end];
	}
	elsif ($trans_id =~ /\w+\.ctg/) {
	    my $type='NEO';
	    $neo_long{$trans_id}=[$gene_id,$type,$fpc,$start,$end];
	}
    }
    
    open (PEP_FILE,"ens_sept25.pep") || die("Could not open ens_sept25.pep for ensembl pep length reading$!");
    print STDERR "Reading ensembl pep file\n";
    my $in = Bio::SeqIO->new(-fh   => \*PEP_FILE, -format=> 'Fasta');
    my $out = Bio::SeqIO->new(-fh => \*STDERR, -format => 'Fasta');
    while ( my $seq = $in->next_seq() ) {
	my $seqid=$seq->id;
	$seqid =~ s/TMPP\_/SEPT20T\./g;
	foreach my $trans_id (keys (%ens_long)) {
	    print STDERR "Trans id $trans_id - Seq id $seqid\n";
	    if ($seqid eq $trans_id) {
		my @string = @{$neo_long{$trans_id}};
		my $id=$string[0];
		$seq->display_id($id);
		my $desc=$string[1].": $trans_id FPC:".$string[2]." FPC_start: ".$string[3]." FPC_end: ".$string[4];
		$seq->desc($desc);
		$out->write_seq($seq);
	    }
	}
    }

    open (PEP_FILE,"neo_pred.pep") || die("Could not open neo_pred.pep for ensembl pep length reading$!");
    print STDERR "Reading neomorphic pep file\n";
    my $in = Bio::SeqIO->new(-fh   => \*PEP_FILE, -format=> 'Fasta');
    my $out = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'Fasta');
    while ( my $seq = $in->next_seq() ) {
	foreach my $trans_id (keys (%neo_long)) {
	    #print STDERR "Seq id: ".$seq->id." trans_id: $trans_id\n";
	    if ($seq->id =~ /$trans_id/) {
		my @string = @{$neo_long{$trans_id}};
		my $id=$string[0];
		$seq->display_id($id);
		my $desc=$string[1].": $trans_id FPC:".$string[2]." FPC_start: ".$string[3]." FPC_end: ".$string[4];
		$seq->desc($desc);
		$out->write_seq($seq);
	    }
	}
    }
}
    
else {
    my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
        
    my $db =  Bio::EnsEMBL::DBLoader->new($locator);
    my $gene_obj=Bio::EnsEMBL::DBSQL::Gene_Obj->new($db);
    foreach my $gene (@gtf_genes) {
	print STDERR "Gene id: ".$gene->id."\n";
	my @exons=$gene->each_unique_Exon;
	my $fpc=$exons[0]->contig_id;
	print STDERR "Got seqname $fpc\n";
	$db->static_golden_path_type('UCSC');
	my $sgp_adaptor = $db->get_StaticGoldenPathAdaptor();
	my $vc = $sgp_adaptor->fetch_VirtualContig_by_fpc_name($fpc);
	foreach my $exon ($gene->each_unique_Exon) {
	    $exon->contig_id($vc->id);
	    $exon->attach_seq($vc->primary_seq);
	    print STDERR "Got exon seq: ".$exon->seq."\n";
	}
	my $newgene = $vc->convert_Gene_to_raw_contig($gene);
	
	print STDERR "Dumping gene ".$gene->id."\n";
	foreach my $trans ($gene->each_Transcript) {
	    my $out = Bio::SeqIO->new(-fh => \*STDOUT, -format => 'Fasta');
	    my $seq= $trans->dna_seq;
	    my $desc="FPC: $fpc ";
	    $seq->desc($desc);
	    $out->write_seq($seq);
	    }
    }
    #$gene_obj->write($newgene);
}

