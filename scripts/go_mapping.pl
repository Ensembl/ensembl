use strict;
use Bio::EnsEMBL::DBLoader;
use GD;

$| = 1;

my $iprgo = shift(@ARGV);
open (IPRGO,"<$iprgo");
while (<IPRGO>) {
    my ($ipr,$go) = split (/\t/);
    print STDERR "INTERPRO $ipr maps to GO $go\n";
}
exit;
#DB parameters
my $dbtype = 'rdb';
my $host   = 'ecs1c';
my $port   = '';
my $dbname = 'ensembl080';
my $dbuser = 'ensadmin';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::Obj';
my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";

print STDERR "Using $locator for db\n";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);
my $sp = $db->get_StaticGoldenPathAdaptor;

my $sth=$db->prepare("select distinct interpro_ac from interpro");
$sth->execute();
print STDERR "Getting all interpro domains...\n";
my @domains;
while (my $row= $sth->fetchrow_array) {
#    print STDERR "Getting $row\n";
    push (@domains,$row);
}

foreach my $domain (@domains) {
    print STDERR "Processing domain $domain\n";
    print STDERR "Finding genes...\n";
    my @genes = find_Genes($db,$domain);
    
    print STDERR "Processing genes...\n";
    my %bychr=process_Genes($domain,@genes);

    #print STDERR "Drawing genes...\n";
    #draw_Genes(%bychr);

    print STDERR "Dumping genes...\n";
    dump_genes($domain,@genes);
    
    #print STDERR "Checking next lowest...\n";
    #next_lowest(%bychr);
}

sub find_Genes {
    my ($db,$domain) = @_;

    my @genes;
    
    #Find all genes that have an interpro domain, 
    #and count how many domains it has in total
    my $sth = $db->prepare("select count(*),tr.gene from interpro i,protein_feature pf,transcript tr where i.interpro_ac = '$domain' and i.id = pf.hid and pf.translation = tr.translation group by tr.gene");

    $sth->execute();
    my $exnum;
    while (my $rowhash = $sth->fetchrow_hashref) {

	my $geneid = $rowhash->{'gene'};
        my $exnum  = $rowhash->{'count(*)'};
	my $gene;

        #Reference to a hash, containing gene id and domain number
	$gene->{id}    = $geneid;
        $gene->{exnum} = $exnum;
	push (@genes,$gene);
    }
    return @genes;
}

sub find_all_Genes {
    my ($db) = @_;

    my @genes;
    
    #Find all genes and get chr.location 
    my $sth = $db->prepare("select id from gene");

    $sth->execute();
    my $exnum;
    while (my $rowhash = $sth->fetchrow_hashref) {

	my $geneid = $rowhash->{'id'};
        #my $exnum  = $rowhash->{'count(*)'};
	my $gene;

        #Reference to a hash, containing gene id and domain number
	$gene->{id}    = $geneid;
        #$gene->{exnum} = $exnum;
	push (@genes,$gene);
    }
    return @genes;
}

sub process_Genes {
    my ($domain,@genes) = @_;
    
    my @located;
    my @sorted;
    my %bychr;

    if ($#genes > 0) {
		
	foreach my $gene (@genes) {

	    #Find chromosomal location
	    #Find approx.start, i.e. start of start_exon raw_contig
	    my $sth = $db->prepare("select STRAIGHT_JOIN p.chr_name,(e.seq_start-p.raw_start+p.chr_start),p.raw_ori,p.chr_start,p.chr_end from transcript tr,translation t,exon e,static_golden_path p where tr.gene = '".$gene->{id}."' and t.id = tr.translation and t.start_exon = e.id and e.contig = p.raw_id");
	    $sth->execute() || die("shit!");
	    my ($chr,$start,$raw_ori,$chr_start,$chr_end) = $sth->fetchrow_array;
	    $gene->{chr} = $chr;

	    #print STDERR "Got chr: $chr, start: $start, raw ori: $raw_ori, chr.start: $chr_start, chr_end: $chr_end\n";
	    if ($raw_ori == -1) {
		($start)=&_flip_coordinates ($start,$chr_start,$chr_end);
	    }
	    $gene->{start} = $start;

	    #Find approx.end, i.e. end of end_exon raw_contig
	    my $sth = $db->prepare("select STRAIGHT_JOIN (e.seq_end-p.raw_start+p.chr_start),p.raw_ori,p.chr_start,p.chr_end from transcript tr,translation t,exon e,static_golden_path p where tr.gene = '".$gene->{id}."' and t.id = tr.translation and t.end_exon = e.id and e.contig = p.raw_id");
	    $sth->execute() || die("shit!");
	    my ($end,$raw_ori,$chr_start,$chr_end) = $sth->fetchrow_array;
	    #print STDERR "and chr:$chr, end: $end, raw_ori: $raw_ori,chr_start: $chr_start,chr_end: $chr_end\n";
	    if ($raw_ori == -1) {
		($end)=&_flip_coordinates ($end,$chr_start,$chr_end);
	    }
	    $gene->{end} = $end;
	    #print STDERR "Gene ".$gene->{id}." has start $start and end $end\n";
	    push @located,$gene;
	}

	foreach my $gene (@located) {
	    push @{$bychr{$gene->{chr}}},$gene;
	}
    } 
    else {
	print STDERR "NO genes found for domain $domain\n";
    }
    return %bychr;
}

sub next_lowest {
    my (%bychr)=@_;
    my %freq;
    my $total=0;
    foreach my $chr (keys (%bychr)){
	my @genes = @{$bychr{$chr}};
	my @sorted = sort { $a->{start} <=> $b->{start}} @genes;
	for (my $i=0; $i <= $#sorted; $i++) {
	    my $gene=$sorted[$i];
	    #print STDERR "Gene ".$gene->{id}." is on chr ".$gene->{chr}." start ".$gene->{start}." and end ".$gene->{end}."\n";
	    if ($sorted[$i+1]) {
		my $next=$sorted[$i+1];
		my $dist=$next->{start}-$gene->{start};
		my $bin;
		my $min=10000;
		my $max= 100000000;
		my $increase=5000000;
		#binning
		if ($dist <= $min) {
		    next;
		}
		my $j;
		for ($j=0;$j<=$max;$j+=$increase) {
		    if (($dist > $j) && ($dist <= ($j+$increase))) {
			$bin= $j+$increase;
			$freq{$bin}++;
		    }
		}
		if ($dist > $max) {
		    $bin= $max;
		    $freq{$max+$increase}++;
		}
		if (defined $bin) {
		    #print $gene->{id}."\t$bin\n";
		}
		else {
		    #print $gene->{id}."\t$dist\n"; 
		}
		$total++;
	    }
	    else {
		#$freq{only}++;
		#print $gene->{id}."\tlast/only\n";
	    }
	}
	
    }
    foreach my $dist (sort {$a <=> $b} (keys(%freq))) {
	my $perc=($freq{$dist}/$total)*100;
	print "$dist\t$perc\n";
    }
}



sub draw_Genes {
    my (%bychr)=@_;
    
    foreach my $chr (keys (%bychr)){
	if (!$chr_length{$chr}) {
	    print STDERR "Not drawing gif for chromosome $chr\n";
	    next;
	}
	my ($im,$blue) = &create_gif;
	foreach my $gene (@{$bychr{$chr}}) {
	    #print STDERR "      gene ".$gene->{id}."\n";
	    my $mb=($y_length*($gene->{start}))/$chr_length{$chr};
	    my $mb2=($y_length*($gene->{end}))/$chr_length{$chr}+1;
	    my $x2 = 60 + ($gene->{exnum}*5);
	    $im->filledRectangle(60,$mb,$x2,$mb2,$blue);
	}
	my $gif = $chr.".gif";
	print STDERR "Writing $chr gif to $gif\n";
	open (GIF,">$gif") || die ("Could not open gif file!\n");
	print GIF $im->gif; 
    }
}

sub create_gif {
    my $im = new GD::Image(100,$y_length); 
    my $white = $im->colorAllocate(255,255,255); 
    my $black = $im->colorAllocate(0,0,0); 
    my $blue = $im->colorAllocate(0,0,255); 
    $im->line(50,10,50,590,$black);
    return ($im,$blue);
}

sub dump_genes {
   my ($domain,@genes)=@_;

    my $dna=$domain.".fa";
    open (DNA,">$dna");
    my $pep=$domain.".pep";
    open (PEP,">$pep");
    my $dna_seqio = Bio::SeqIO->new('-format' => 'Fasta' , -fh => \*DNA );
    my $pep_seqio = Bio::SeqIO->new('-format' => 'Fasta' , -fh => \*PEP );
    my @real_genes;
    foreach my $gene (@genes) {
	eval {
	    push @real_genes,$db->gene_Obj->get($gene->{id});
	};
	if ($@) {
	    print STDERR "Couldn't get gene ".$gene->{id}.", skipping\n";
	    next;
	}

    }

    foreach my $gene (@real_genes) {
	foreach my $trans ( $gene->each_Transcript ) {
	    # get out first exon. Tag it to clone and gene on this basis
	    my @exon = $trans->each_Exon;
	    my $fe = $exon[0];
	    my $tseq;
	    eval {
		$tseq = $trans->translate();
	    };
	    if ($@) {
		print STDERR "Can't translate ".$trans->id.", skipping\n";
		next;
	    }
	    $tseq->id($trans->id);
	    my $string;
	    foreach my $link ($gene->each_DBLink) {
		$string .= $link->database.":".$link->primary_id." ";
	    }
	    my $chr;
	    my $start;
	    my $end;
	    foreach my $hashgene (@genes) {
		if ($hashgene->{id} eq $gene->id) {
		    $chr=$hashgene->{chr};
		    $start=$hashgene->{start};
		    $end=$hashgene->{end};
		}
	    }
	    if ( $tseq->seq =~ /\*/ ) {
		print STDERR "translation has stop codons. Skipping! (in clone". $fe->clone_id .")\n";
		next;
	    }
	    $tseq->desc("Gene:".$gene->id." Clone:".$fe->clone_id . " Contig:" . $fe->contig_id." chr: ".$chr." start: $start end: $end dblinks: ".$string);
	    $pep_seqio->write_seq($tseq);
	    
	    my $seq = $trans->dna_seq();
	    $seq->id($trans->id);
	    $seq->desc("Gene:".$gene->id." Clone:".$fe->clone_id." chr: ".$chr." start: ".$start." end: ".$end." dblinks: ".$string);
	    $dna_seqio->write_seq($seq);
	}
    }
}

sub _flip_coordinates {
    my ($pos,$chr_start,$chr_end)=@_;
    
    die("need a pos") unless $pos;
    die("need a chromosome start") unless $chr_start;
    die("need a chromosome end") unless $chr_end;

    my $vc_pos=$chr_end+$chr_start-$pos;
    return ($vc_pos);
}
