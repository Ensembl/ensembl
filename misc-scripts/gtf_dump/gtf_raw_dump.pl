# this script uses the hardcoded db connection right in the beginning to
# an ensembl database and makes a gene gtf dump into genes_raw_coords.gtf
# The coord system used is raw contigs. Exons on more than 1 contig
# (sticky) get the contig name of theier first part and the size of all 
# partly exons together.
 
use DBI;

$dsn = "dbi:mysql:host=ensrv4.sanger.ac.uk;database=ensembl100";
$db = DBI->connect( $dsn, "ensro" );

$sth = $db->prepare( "
            SELECT t.gene, t.id, e.id, e.seq_start, 
                   e.seq_end, e.strand, e.phase, 
                   e.sticky_rank, sp.fpcctg_name, 
                   sp.fpcctg_start, sp.fpcctg_end, 
                   sp.raw_start, sp.raw_ori, tl.start_exon, 
                   tl.seq_start, tl.end_exon, tl.seq_end,
                   c.id
              FROM transcript t, translation tl, 
                   exon_transcript et , exon e, 
                   static_golden_path sp,
                   contig c
             WHERE et.transcript = t.id 
               AND tl.id = t.translation 
               AND et.exon = e.id 
               AND e.contig = sp.raw_id
               AND e.contig = c.internal_id 
          ORDER BY t.gene, t.id, et.rank, e.sticky_rank
" );
$sth->execute();



my %gh;

while( $colsRef = $sth->fetchrow_arrayref() ) {

  if( ! defined $_first ) {
    $_first = 1;
    next;
  }

  @cols = @$colsRef;
  
  if( !exists $gh{$cols[0]} ) {
    # print_genes();
    # delete the other genes
    # %gh = ();
    $gh{$cols[0]} = {};
    $tr = 1; 
  }

  $g = $gh{$cols[0]};

  if( $g->{buggy} ) {
    next;
  }

  $g->{name} = $cols[0];

  if( ! exists $g->{transcript}{$cols[1]} ) {
    $g->{transcript}{$cols[1]} = { 'rank' => $tr,
				 'name' => $cols[1]};
    $tr++;
    # first exon comes as rank 1
    $er = 1;
  }

  $t = $g->{transcript}{$cols[1]};

  if( exists $t->{exon}{$cols[2]} ) {
    # sticky exon
    # print "Sticky Exon ",$cols[7],".\n";
    if( $cols[7] == 1 ) {
      $g->{buggy} = 1;
      next;
      # buggy gene
    }

    my ( $addLength );
    $e = $t->{exon}{$cols[2]};

    $addLength = $cols[4]-$cols[3]+1;

    if( $e->{strand} == 1 ) {
      $e->{seqend} += $addLength;
    } else {
      $e->{seqstart} -= $addLength;
    }

  } else {
    
    $e = { 'name' => $cols[2],
	      'rank' => $er,
	  'contig' => $cols[17],
	 'phase' => $cols[6]};
          
    $t->{exon}{$cols[2]} = $e;

    # exon rank up
    $er++;

    if( $cols[13] eq $e->{name} ) {
      $e->{startcodon} = $cols[14];
    }
    if( $cols[15] eq $e->{name} ) {
      $e->{stopcodon} = $cols[16];
    }

    $e->{strand} = $cols[5];

    $e->{seqstart} = $cols[3];
    $e->{seqend} = $cols[4];

  }
}

&print_genes();

  
sub print_genes {
  my @genes = values %gh;;
  
  local *FH;
  
  open( FH, ">genes_raw_coords.gtf" ) or die( "Couldnt open output file" );
  print STDERR scalar( @genes )," Genes\n";

  foreach $g ( @genes ) {
    if( $g->{buggy} ) {
      next;
    }

    @transcripts = values %{$g->{transcript}};
    @transcripts = sort { $a->{rank} <=> $b->{rank} } @transcripts; 
    foreach $t ( @transcripts ) {
      my @exons = values %{$t->{exon}};
      @exons = sort { $a->{rank} <=> $b->{rank} } @exons;
	
      foreach $e ( @exons ) {
	if( exists $e->{startcodon} ) {
	  my ( $start, $end );
	  if( $e->{strand} == 1 ) {
	    $start = $e->{seqstart} + $e->{startcodon} -1;
	    $end = $start + 2;
	  } else {
	    $end = $e->{seqend} - $e->{startcodon} + 1;
	    $start = $end -2;
	  }

	  print FH $e->{contig},"\t",
	  "ENSEMBL\tstart_codon\t", 
	  $start,"\t",
	  $end,"\t",
	  "0\t",
	  ($e->{strand}==1?"+":"-"),"\t.\t",
	  "gene_id \"",$g->{name},"\"; ",
	  "transcript_id \"", $t->{name}, "\"; ",
	  "exon_number ", $e->{rank},"; ",
	  "exon_id \"", $e->{name},"\"",
	  "\n";
	}
	print FH $e->{contig},"\t",
	"ENSEMBL\tCDS\t", 
	$e->{seqstart},"\t",
	$e->{seqend},"\t",
	".\t",
	($e->{strand}==1?"+":"-"),"\t",$e->{phase},"\t",
	"gene_id \"",$g->{name},"\"; ",
	"transcript_id \"", $t->{name}, "\"; ",
	"exon_number ", $e->{rank},"; ",
	"exon_id \"", $e->{name},"\"",
	"\n";

	if( exists $e->{stopcodon} ) {
	  my ( $start, $end );
	  if( $e->{strand} == 1 ) {
	    $end = $e->{seqstart} + $e->{stopcodon} - 1;
	    $start = $end - 2;
	  } else {
	    $start = $e->{seqend} - $e->{stopcodon} + 1;
	    $end = $start + 2;
	  }
	  print FH $e->{contig},"\t",
	  "ENSEMBL\tstop_codon\t", 
	  $start,"\t",
	  $end,"\t",
	  "0\t",
	  ($e->{strand}==1?"+":"-"),"\t.\t",
	  "gene_id \"",$g->{name},"\"; ",
	  "transcript_id \"", $t->{name}, "\"; ",
	  "exon_number ", $e->{rank},"; ",
	  "exon_id \"", $e->{name},"\"",
	  "\n";
	}
      }
    }
    # print "GENE:",$g->{name}, " ", $g->{fpc},"\n";
  }
  close FH;
}
  
