#!/usr/local/bin/perl
use Bio::SeqIO;
use Bio::PrimarySeq;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

$| = 1;

my $length = shift @ARGV;

my $db=Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname=>'homo_sapiens_core_110',-host=>'ensrv3',-user=>'ensadmin');
my $sth = $db->prepare("select start_exon,id from translation");
$sth->execute;
my %et;
my @start_exons;
while(my ($se,$id) = $sth->fetchrow_array) {
    $et{$se}=$id;
}

foreach my $se (keys %et) {
    print STDERR "Dumping sequence upstream of ".$et{$se}."\n";
    my $sth = $db->prepare("select e.strand,if(e.strand = 1,if(e.seq_start > $length,substring(sequence,e.seq_start-$length,$length),substring(sequence,$length-e.seq_start,e.seq_start)),if(c.length-e.seq_end > $length,substring(sequence,e.seq_end,$length),substring(sequence,e.seq_end,c.length))), if(e.strand = 1,if(e.seq_start > $length,length(substring(sequence,e.seq_start-$length,$length)),length(substring(sequence,$length-e.seq_start,e.seq_start))),if(c.length-e.seq_end > $length,length(substring(sequence,e.seq_end,$length)),length(substring(sequence,e.seq_end,c.length)))),c.internal_id from exon e, contig c, dna d where e.contig = c.internal_id and c.dna = d.id and e.id = '".$se."'");
    $sth->execute;
    TRANS: while (my ($strand,$string, $cl,$contig) = $sth->fetchrow_array) {
	my ($sth2,$sth3,$seqstr,$contig2);
	my @gaps = (50002,100002,200002);
	if ($cl < $length) {
	    my $diff = $length-$cl;
	    my $string2;
	    if ($strand == 1) {
		$sth2 = $db->prepare("select sgp2.raw_id from static_golden_path sgp1, static_golden_path sgp2 where sgp1.raw_id = $contig and sgp2.chr_end < sgp1.chr_start-102 and sgp2.chr_end < sgp1.chr_start and sgp1.chr_name = sgp2.chr_name");
		$sth2->execute;
		($contig2) = $sth2->fetchrow_array;
		my $c=0;
		while (!$contig2) {
		    if ($gaps[$c]) {
			$sth2 = $db->prepare("select sgp2.raw_id from static_golden_path sgp1, static_golden_path sgp2 where sgp1.raw_id = $contig and sgp2.chr_end < sgp1.chr_start-$gaps[$c] and sgp2.chr_end < sgp1.chr_start and sgp1.chr_name = sgp2.chr_name");
			$sth2->execute;
			($contig2) = $sth2->fetchrow_array;
			$c++;
		    }
		    else {
			last;
		    }
		}
		if ($contig2) { 
		    $sth3 =  $db->prepare("select substring(sequence,c.length-$diff,$diff), length(substring(sequence,c.length-$diff,$diff)),c.internal_id from contig c, dna d where c.dna = d.id and c.internal_id = $contig2");
		    $sth3->execute;
		    ($string2) = $sth3->fetchrow_array;
		}
	    }
	    else {
		$sth2 = $db->prepare("select sgp2.raw_id from static_golden_path sgp1, static_golden_path sgp2 where sgp1.raw_id = $contig and sgp2.chr_start < sgp1.chr_end+102 and sgp2.chr_start > sgp1.chr_end and sgp1.chr_name = sgp2.chr_name");
		$sth2->execute;
		($contig2) = $sth2->fetchrow_array;
		my $c=0;
		while (!$contig2) {
		    if ($gaps[$c]) {
			$sth2 = $db->prepare("select sgp2.raw_id from static_golden_path sgp1, static_golden_path sgp2 where sgp1.raw_id = $contig and sgp2.chr_start < sgp1.chr_end+$gaps[$c] and sgp2.chr_start > sgp1.chr_end and sgp1.chr_name = sgp2.chr_name");
			$sth2->execute;
			($contig2) = $sth2->fetchrow_array;
			$c++;
		    }
		    else {
			last;
		    }
		}
		if ($contig2) {
		    $sth3 =  $db->prepare("select substring(sequence,1,$diff), length(substring(sequence,1,$diff)),c.internal_id from contig c, dna d where c.dna = d.id and c.internal_id = $contig2");
		    $sth3->execute;
		    ($string2) = $sth3->fetchrow_array;
		}
	    }
	    $seqstr = $string.$string2;
	}
	else {
	    $seqstr = $string;
	}
	my $seq;
	eval {
	    $seq = Bio::PrimarySeq->new(
					-id => $et{$se},
					-seq => $string
					);
	};
	if ($@) {
	    print STDERR "Could not dump upstream sequence for ".$et{$se}." because of $@\n";
	    next;
	}
	if ($strand == -1) {
	    $seq = $seq->revcom;
	}
	my $seqio = new Bio::SeqIO('-format' => 'fasta',
				   '-fh' => \*STDOUT);
	$seqio->write_seq($seq);
    }
}
