#!/usr/local/bin/perl






use DBI;
use Getopt::Long;

my $host   = "kaka.sanger.ac.uk";
my $dbname = "current";
my $dbuser = "anonymous";
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


$sth = $db->prepare( "
            SELECT t.gene, t.id, e.id, e.seq_start, 
                   e.seq_end, e.strand, e.phase, 
                   e.sticky_rank, sp.fpcctg_name, 
                   sp.fpcctg_start, sp.fpcctg_end, 
                   sp.raw_start, sp.raw_ori, tl.start_exon, 
                   tl.seq_start, tl.end_exon, tl.seq_end 
              FROM transcript t, translation tl, 
                   exon_transcript et , exon e, 
                   static_golden_path sp 
             WHERE et.transcript = t.id 
               AND tl.id = t.translation 
               AND et.exon = e.id 
               AND e.contig = sp.raw_id 
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

  $g->{fpc} = $cols[8];
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

    my ( $seqstart, $seqend, $strand );
    $e = $t->{exon}{$cols[2]};

    if( $cols[12] == 1 ) {
      $seqstart = $cols[3]-$cols[11]+$cols[9];
      $seqend = $cols[4]-$cols[11]+$cols[9];
    } else {
      $seqend = $cols[10] - $cols[3]+$cols[11];
      $seqstart = $cols[10] - $cols[4]+$cols[11];
    }

    $strand = $cols[5]*$cols[12];
    
    if( $e->{strand} == 1 ) {
      if( $e->{seqend} +1 != $seqstart ) {
	print "COMPLAIN: ",$e->{name}, " ",$e->{seqend}," ",$e->{seqstart}," $seqstart $seqend\n";
      } else {
	# print "Gottit\n";
	$e->{seqend} = $seqend;
      }
    } else {
      if( $e->{seqstart} -1 != $seqend ) {
	print "COMPLAIN: ",$e->{seqend}," ",$e->{seqstart}," $seqstart $seqend\n";
      } else {
	# print "Gottit\n";
	$e->{seqstart} = $seqstart;
      }
    }
    if( $er == 2  && $tr == 2 ) {
      $g->{start} = $e->{seqstart};
    }
  } else {
    
    $e = { 'name' => $cols[2],
	      'rank' => $er};
    $t->{exon}{$cols[2]} = $e;

    # exon rank up
    $er++;

    if( $cols[13] eq $e->{name} ) {
      $e->{startcodon} = $cols[14];
    }
    if( $cols[15] eq $e->{name} ) {
      $e->{stopcodon} = $cols[16];
    }

    $e->{strand} = $cols[5]*$cols[12];

    if( $cols[12] == 1 ) {
      $e->{seqstart} = $cols[3]-$cols[11]+$cols[9];
      $e->{seqend} = $cols[4]-$cols[11]+$cols[9];
    } else {
      $e->{seqend} = $cols[10] - $cols[3]+$cols[11];
      $e->{seqstart} = $cols[10] - $cols[4]+$cols[11];
    }

    if( $er ==2 && $tr == 2 ) {
      $g->{start} = $e->{seqstart};
    }
  }
}

&print_genes();

  
sub print_genes {
  my @sortedGenes;
  my @genes = values %gh;;
  
  local *FH;
  
  @sortedGenes = sort { $a->{fpc} cmp $b->{fpc} ||
			  $a->{start} <=> $b->{start} } @genes;
  
  $fpc = undef;
  print STDERR scalar( @sortedGenes )," Genes\n";

  foreach $g ( @sortedGenes ) {
    if( $g->{buggy} ) {
      next;
    }

    if( $g->{fpc} ne $fpc ) {
      close FH;
      open( FH, ">".$g->{fpc}.".gtf" ) or die( "Couldnt write result" );
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

	  print FH $g->{fpc},"\t",
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
	print FH $g->{fpc},"\t",
	"ENSEMBL\texon\t", 
	$e->{seqstart},"\t",
	$e->{seqend},"\t",
	"0\t",
	($e->{strand}==1?"+":"-"),"\t.\t",
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
	  print FH $g->{fpc},"\t",
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
    $fpc = $g->{fpc};
  }
  close FH;
}
  
