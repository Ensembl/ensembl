#!/usr/local/bin/perl






use DBI;
use Getopt::Long;

my $host   = "ecs1d.sanger.ac.uk";
my $dbname = "mus_musculus_core_5_3";
my $dbuser = "ensro";
my $dbpass = undef;

&GetOptions( 
	     'host:s'   => \$host,
	     'dbname:s' => \$dbname,
	     'dbuser:s' => \$dbuser,
	     'dbpass:s' => \$dbpass );


$dsn = "dbi:mysql:host=$host;database=$dbname";
$db = DBI->connect( $dsn, $dbuser , $dbpass );

if( !defined $db ) {
    die "Unable to connect to database host $host database $dbname $dbpass";
}



#$sth = $db->prepare( "
#            SELECT gene_stable_id, transcript_stable_id, exon_stable_id, 
#	           exon_chrom_start, exon_chrom_end, exon_chrom_strand, 
#                   rank, start_rank, end_rank, seq_start, seq_end,
#		   chr_name 
#              FROM gene_structure
#          ORDER BY gene_stable_id, transcript_stable_id, rank
#" );
#$sth->execute();


$sth = $db->prepare( "
 SELECT gsi.stable_id as gene_stable_id,
         tsi.stable_id as transcript_stable_id,
         esi.stable_id as exon_stable_id,
         MIN(IF(sgp.raw_ori=1,
             ( e.seq_start+sgp.chr_start-sgp.raw_start ),
             ( sgp.chr_start+sgp.raw_end-e.seq_end ))) 
          as exon_chrom_start,
         MAX(IF(sgp.raw_ori=1,
             ( e.seq_end+sgp.chr_start-sgp.raw_start ),
             ( sgp.chr_start+sgp.raw_end-e.seq_start ))) 
          as exon_chrom_end, 
         e.strand * sgp.raw_ori as exon_chrom_strand,
         et.rank as rank,
         tl.start_exon_id,
         tl.seq_start as seq_start,
         tl.end_exon_id,
         tl.seq_end as seq_end,
         sgp.chr_name as chr_name,
         e.exon_id as exon_id
   FROM  exon e, exon_stable_id esi,
         exon_transcript et, transcript t,
         transcript_stable_id tsi,
         static_golden_path sgp,
         gene_stable_id gsi,
         translation tl
   WHERE e.contig_id = sgp.raw_id
     AND esi.exon_id = e.exon_id
     AND et.exon_id = e.exon_id
     AND t.transcript_id = et.transcript_id
     AND t.translation_id = tl.translation_id
     AND tsi.transcript_id = et.transcript_id
     AND gsi.gene_id = t.gene_id
   GROUP BY t.transcript_id, e.exon_id
   ORDER BY gene_stable_id, transcript_stable_id, rank
" );

$sth->execute();

my %gh;

while( my $h = $sth->fetchrow_hashref() ) {

  
  if( !exists $gh{$h->{gene_stable_id}} ) {
    # print_genes();
    # delete the other genes
    # %gh = ();
    $gh{$h->{gene_stable_id}} = {};
    $tr = 1; 
  }

  $g = $gh{$h->{gene_stable_id}};
  $g->{name} = $h->{gene_stable_id};
  $g->{chrom} = $h->{chr_name};	

  if( ! exists $g->{transcript}{$h->{transcript_stable_id}} ) {
    $g->{transcript}{$h->{transcript_stable_id}} = { 'rank' => $tr,
				 'name' => $h->{transcript_stable_id}};
    $tr++;
    # first exon comes as rank 1
    $er = 1;
  }

  $t = $g->{transcript}{$h->{transcript_stable_id}};

  if( ! exists $t->{exon}{$h->{exon_stable_id}} ) {
    
    $e = { 'name' => $h->{exon_stable_id},
	      'rank' => $h->{rank}};
    $t->{exon}{$h->{exon_stable_id}} = $e;

    # exon rank up
    $er++;

    if( $h->{start_exon_id} == $h->{exon_id} ) {
      $e->{startcodon} = $h->{seq_start};
    }
    
    if( $h->{end_exon_id} == $h->{exon_id}) {
      $e->{stopcodon} = $h->{seq_end};
    }

    $e->{strand} = $h->{exon_chrom_strand};

    $e->{seqstart} = $h->{exon_chrom_start};
    $e->{seqend} = $h->{exon_chrom_end};

    	
    if( $er ==2 && $tr == 2 ) {
      $g->{start} = $e->{seqstart};
      $g->{end} = $e->{seqend};
    } else {
      if( $g->{start} > $e->{seqstart} ) { $g->{start} = $e->{seqstart} }
      if( $g->{end} < $e->{seqend} ) { $g->{end} = $e->{seqend} }
    }
  }
}

&print_genes();

  
sub print_genes {
  my @sortedGenes;
  my @genes = values %gh;;
  
  local *FH;
  my $chrom;
  
  @sortedGenes = sort { $a->{chrom} cmp $b->{chrom} || $a->{start} <=> $b->{start} } @genes;
  
  print STDERR scalar( @sortedGenes )," Genes\n";

  foreach $g ( @sortedGenes ) {

    if( $g->{chrom} ne $chrom ) {
      close FH;
      open( FH, ">".$g->{chrom}.".gtf" ) or die( "Couldnt write result" );
    }

    @transcripts = values %{$g->{transcript}};
    @transcripts = sort { $a->{rank} <=> $b->{rank} } @transcripts; 
    foreach $t ( @transcripts ) {
      my @exons = values %{$t->{exon}};
      @exons = sort { $a->{rank} <=> $b->{rank} } @exons;
      my $coding = 0;	
      my $frame = 0;
      
      foreach $e ( @exons ) {
	my $gtfLineTail = 
	  "gene_id \"".$g->{name}."\"; ".
	  "transcript_id \"". $t->{name}. "\"; ".
	  "exon_number \"". $e->{rank}."\"; ".
	  "exon_id \"". $e->{name}."\"".
	  "\n";

	if( exists $e->{startcodon} ) {
	  my ( $start, $end );
	  if( $e->{strand} == 1 ) {
	    $start = $e->{seqstart} + $e->{startcodon} -1;
	    $end = $start + 2;
	  } else {
	    $end = $e->{seqend} - $e->{startcodon} + 1;
	    $start = $end -2;
	  }

	  print FH $g->{chrom},"\t",
	  "ENSEMBL\tstart_codon\t", 
	  $start,"\t",
	  $end,"\t",
	  ".\t",
	  ($e->{strand}==1?"+":"-"),"\t.\t",
	  $gtfLineTail;
	}
	print FH $g->{chrom},"\t",
	"ENSEMBL\texon\t", 
	$e->{seqstart},"\t",
	$e->{seqend},"\t",
	".\t",
	($e->{strand}==1?"+":"-"),"\t.\t",
	$gtfLineTail;

	if( exists $e->{stopcodon} ) {
	  my ( $start, $end );
	  if( $e->{strand} == 1 ) {
	    $end = $e->{seqstart} + $e->{stopcodon} - 1;
	    $start = $end - 2;
	  } else {
	    $start = $e->{seqend} - $e->{stopcodon} + 1;
	    $end = $start + 2;
	  }
	  print FH $g->{chrom},"\t",
	  "ENSEMBL\tstop_codon\t", 
	  $start,"\t",
	  $end,"\t",
	  ".\t",
	  ($e->{strand}==1?"+":"-"),"\t.\t",
	  $gtfLineTail;
	}

        # 0 not coding, 1 contains just startcodon
	# 2 contains only coding
	# 3 just stopcodon, 4 start and stopcodon
	
	if( $coding == 0 ) {
	  if( exists $e->{startcodon} ) {
	    if( exists $e->{stopcodon} ) {
	      $coding = 4;
	    } else {
	      $coding = 1;
	    }
	  }
	} elsif( $coding == 1 ) {
	  if( exists $e->{stopcodon} ) {
	    $coding = 4;
	  } else {
	    $coding = 2;
	  }
	} elsif( $coding == 2 ) {
          if( exists $e->{stopcodon} ) {
	    $coding = 3;
	  }
	} elsif( $coding == 3 ) {
	  $coding = 0;
	} elsif( $coding == 4 ) {
	  $coding = 0;
	}
	
	# construct CDS line
	my ( $cdsstart, $cdsend );
	
	if( $coding == 0 ) {
	  next;
	  # no CDS line
	} elsif( $coding == 1 ) {
	  if( $e->{strand} == 1 ) {
	    $cdsstart = $e->{startcodon} + $e->{seqstart} - 1;
	    $cdsend = $e->{seqend};
	  } else {
	    $cdsstart = $e->{seqstart};
	    $cdsend = $e->{seqend} - $e->{startcodon} + 1 ;
	  }
	} elsif( $coding == 2 ) {
	  $cdsstart = $e->{seqstart};
	  $cdsend = $e->{seqend};
	} elsif( $coding == 3 ) {
	  if( $e->{strand} == 1 ) {
	    $cdsstart = $e->{seqstart};
	    $cdsend = $e->{seqstart}+$e->{stopcodon} - 4;
	  } else {
	    $cdsstart = $e->{seqend} - $e->{stopcodon} + 4;
	    $cdsend = $e->{seqend};
	  }
	} elsif( $coding == 4 ) {
	  if( $e->{strand} == 1 ) {
	    $cdsstart = $e->{startcodon} + $e->{seqstart} - 1;
	    $cdsend = $e->{seqstart}+$e->{stopcodon} - 4;
	  } else {
	    $cdsstart = $e->{seqend} - $e->{stopcodon} + 4;
	    $cdsend = $e->{seqend} - $e->{startcodon} + 1 ;
	  }  
	}
	
	if( $cdsend >= $cdsstart ) {
	  print FH $g->{chrom},"\t",
	  "ENSEMBL\tCDS\t", 
	  $cdsstart,"\t",
	  $cdsend,"\t",
	  ".\t",
	  ($e->{strand}==1?"+":"-"),"\t",$frame,"\t",
	  $gtfLineTail;

	  $frame = ( $cdsend - $cdsstart + 1 - $frame ) % 3;
	} 
      }
    }
    # print "GENE:",$g->{name}, " ", $g->{fpc},"\n";
    $chrom = $g->{chrom};
  }
  close FH;
}
  
