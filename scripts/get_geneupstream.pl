#!/usr/local/bin/perl -w

# $Id$

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

# this is more efficient SQL query, but you won't see results until entire
# set is returned
#foreach my $gene ( $db->gene_Obj->get_array_supporting('none',@gene_ids ) )
foreach my $geneid (@gene_ids )  
{
    my $gene = $db->gene_Obj->get($geneid);
#    my $geneid = $gene->id;
    foreach my $trans ( $gene->each_Transcript ) {
	my @exon = $trans->each_Exon;
	my $fe = $exon[0];
	my ($chr,$gene_start,$cdna_start) = find_trans_start($trans);	

	# hmm, we skip ones that are not correct
	if( ! $chr || ! $gene_start ) {
	    print STDERR "skipping $geneid because location is Chr:$chr Start:$gene_start\n"; 
	    next;
	}
	# Just want first exon.  Going to assume here that the first exon
	# is on the same strand as the rest of the gene which I know
	# is not always true in the Ensembl world because of
	# misassembled sequence in GP
	
	my $vc;
	my $strand = $fe->strand;
	if(  $strand == -1 ) {
	    # gene is on the rev strand, go forwards 
	    # (in the global chromosome sense)
	    # '$upstream' bases
	    $vc = $sa->fetch_VirtualContig_by_chr_start_end
		($chr, $gene_start, 
		 $gene_start + $upstreamlen);

	    # lets return all the sequence on the same strand in the
	    # same direction , so need to reverse complement;

	    $vc->revcom;

# this will fetch the 1st 10 bases of the gene

#	    $vc = $sa->fetch_VirtualContig_by_chr_start_end($chr, 
#	    $gene_start- 10,
#	    $gene_start);
#	    $vc = $vc->revcom;

	} else {
	    # gene is on fwd strand, go backwards (in global chromsome sense)
	    # '$upstream' bases

# this will fetch the $upstream bases from the gene
	    $vc = $sa->fetch_VirtualContig_by_chr_start_end
		($chr, 
		 $gene_start - $upstreamlen, 
		 $gene_start);

	    # this will fetch the 1st 10 bases of the gene
#	    $vc = $sa->fetch_VirtualContig_by_chr_start_end($chr, 
#							    $gene_start, 
#							    $gene_start + 10);

	}    
	$vc->id($geneid);
	$vc->desc("Strand:$strand Chrom:$chr basepair:$gene_start UpstreamBases:$upstreamlen" );	
	$seqio->write_seq($vc); 
    }
}

# taken from gene2flat.pl

sub  find_trans_start {
    my ($trans) = @_;

    my $start_pos;
    my $trans_pos;

    my $contig; 
    foreach my $exon ($trans->each_Exon) {
	if ($exon->id eq $trans->translation->start_exon_id) {
	    $contig = $exon->contig_id;
	    if ($exon->strand == 1) {
		$start_pos = $exon->start;
		$trans_pos = $exon->start + $trans->translation->start - 1;
	    } else {
		$start_pos = $exon->end;
		$trans_pos = $exon->end - $trans->translation->start + 1;
	    }
	}
    }
    if (!defined($start_pos)) {
	print STDERR "Couldn't find start exon for " . $trans->id . "\n";
	die;
    }
    my $query = "select chr_name,chr_start,chr_end,raw_start,raw_end,raw_ori from static_golden_path,contig where raw_id = internal_id and contig.id = \'" . $contig . "\'";
    my $sth  = $db->prepare($query);
    my $res  = $sth->execute;

    my $row = $sth->fetchrow_hashref;

    my $chr = $row->{chr_name};
    my $chr_start = $row->{chr_start};
    my $chr_end   = $row->{chr_end};
    my $raw_start = $row->{raw_start};
    my $raw_end   = $row->{raw_end}; 
    my $raw_ori   = $row->{raw_ori};

    my $gene_start; 
    my $cdna_start;

    if ($raw_ori == 1) {
	$gene_start = $chr_start + ($trans_pos - $raw_start);
	$cdna_start = $chr_start + ($start_pos - $raw_start);
    } else {
	$cdna_start = $chr_end   - ($start_pos - $raw_start);
	$gene_start = $chr_end   - ($trans_pos - $raw_start);
    } 

    return ($chr,$gene_start,$cdna_start);
}
