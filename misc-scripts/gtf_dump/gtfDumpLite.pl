#!/usr/local/bin/perl






use DBI;
use Getopt::Long;

my $host   = "kaka.sanger.ac.uk";
my $dbname = "homo_sapiens_lite_4_28";
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
            SELECT gene_stable_id, transcript_stable_id, exon_stable_id, 
	           exon_chrom_start, exon_chrom_end, exon_chrom_strand, 
                   rank, start_rank, end_rank, seq_start, seq_end,
		   chr_name 
              FROM gene_structure
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

    if( $h->{rank} == $h->{start_rank} ) {
      $e->{startcodon} = $h->{seq_start};
    }
    
    if( $h->{rank} == $h->{end_rank}) {
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

	  print FH $g->{chrom},"\t",
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
	print FH $g->{chrom},"\t",
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
	  print FH $g->{chrom},"\t",
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
    $chrom = $g->{chrom};
  }
  close FH;
}
  
