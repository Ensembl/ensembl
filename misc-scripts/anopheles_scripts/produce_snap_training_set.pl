#!/usr/local/ensembl/bin/perl -w

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;


my $dbhost;
my $dbuser    = 'ensro';
my $dbname;
my $dbpass    = undef;
my $dnadbhost;
my $dnadbname;

my $genetype; # default genetype


$dbuser = "ensro";
GetOptions(
	   'dbname:s'    => \$dbname,
	   'dbhost:s'    => \$dbhost,
	   'dnadbhost:s' => \$dnadbhost,
	   'dnadbname:s' => \$dnadbname,
	   'genetype:s'  => \$genetype,
	  );

unless ( $dbname && $dbhost && $dnadbname && $dnadbhost ){
  print STDERR "Usage: $0 -dbname -dbhost (-genetype)\n";
  exit(0);
}

my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					       '-host'   => $dnadbhost,
					       '-user'   => $dbuser,
					       '-dbname' => $dnadbname,
					       '-pass'   => $dbpass,
					      );

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $dbhost,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					    '-dnadb'  => $dnadb,
					   );



print STDERR "connected to $dbname : $dbhost\n";

print STDERR "checking genes of type $genetype\n";

my $ok = 0;
my $nope = 0;


open (EX,">/acari/work1/mongin/training_set/exon_coordinates_filter.txt") || die;
open (GEN,">/acari/work1/mongin/training_set/genomic_filter.fa") || die;
#open (SPLICEUP,">/acari/work1/mongin/training_set/splice_gt.fa") || die;
#open (SPLICEDOWN,">/acari/work1/mongin/training_set/splice_ag.fa") || die;


my @gene_ids      = @{$db->get_GeneAdaptor->list_geneIds};
my $slice_adaptor = $db->get_SliceAdaptor;

GENE:
foreach my $gene_id ( @gene_ids){
  my $gene = $db->get_GeneAdaptor->fetch_by_dbID($gene_id,1);
  
  if ($genetype){
    next GENE unless ( $gene->type eq $genetype );
  }

   my @transcripts = @{$gene->get_all_Transcripts};
  
  foreach my $tr(@transcripts) {
      my $utr = &check_utr($tr);
      my $start = &check_start($tr);
      my $splice  = &check_splice_sites($tr);

      

      if (($utr)&&($start)&&($splice)) {
	  $ok++;
	  my $stable = $tr->stable_id;
	  
	  my $gen_slice = $slice_adaptor->fetch_by_transcript_stable_id($stable);
	  my $gen = $gen_slice->seq;
	  
	  my @genes = (@{$gen_slice->get_all_Genes});
	  
	  foreach my $g(@genes) {
	      my @transcs = (@{$g->get_all_Transcripts});

	      foreach my $trs(@transcs) {
		  if ($trs->stable_id eq $stable) {
		      
		
		      my @exons = (@{$trs->get_all_translateable_Exons});
		      
		      print GEN ">$stable\n$gen\n";
		      print EX ">$stable\n";
		      foreach my $ex(@exons) {
			  
			  if ($ex->strand == 1)  {
			      print EX "Exon\t".$ex->start."\t".$ex->end."\t".$stable."\n";
			  }
			  else {
			      print EX "Exon\t".$ex->end."\t".$ex->start."\t".$stable."\n";
			  }
		      }
		      
		  }
		  else {
		      $nope++;
		  }
		  
	      }
	  }
      }
      
  }
#  print STDERR "OK: $ok\nNope: $nope\n"
  
}


sub check_utr {
    my ($tr) = @_;
    my $cdna_length = length($tr->seq->seq);
    my $coding_end = $tr->cdna_coding_end;
    my $tl_start = $tr->translation->start;
    
    my $diff = abs($cdna_length - $coding_end);

    if (($tl_start > 3) && ($diff > 3)) {
	return 1;
    }
}

sub check_start {
    my ($tr) = @_;
    my $seq = $tr->translate->seq;
    
    if ($seq =~ /^M/) {
	return 1;
    }
}

sub check_splice_sites {
    my ($transcript) = @_;
    $transcript->sort;
    
    my $strand = $transcript->start_Exon->strand;
    my @exons  = @{$transcript->get_all_Exons};
  
    my $introns  = scalar(@exons) - 1 ; 
    if ( $introns <= 0 ){
	return (0,0,0,0);
    }
  
    my $correct        = 0;
  
    # all exons in the transcripts are in the same seqname coordinate system:
    my $slice = $transcript->start_Exon->contig;
  
    if ($strand == 1 ){
	
      INTRON:
	for (my $i=0; $i<$#exons; $i++ ){
	   
	    my $upstream_exon   = $exons[$i];
	   
	    my $downstream_exon = $exons[$i+1];

	    my $upstream_site;
	    my $downstream_site;

	    my $upstream_site1;
	    my $downstream_site1;


	     eval{
		 $upstream_site1 = 
		     $slice->subseq( ($upstream_exon->end     + 1), ($upstream_exon->end     + 13 ) );

		$downstream_site1 = 
		    $slice->subseq( ($downstream_exon->start - 13), ($downstream_exon->start - 1 ) );
		 
	    };
	    unless ( $upstream_site1 && $downstream_site1 ){
		print STDERR "problems retrieving sequence for splice sites\n$@";
		next INTRON;
	    }

	    eval{
		$upstream_site = 
		    $slice->subseq( ($upstream_exon->end     + 1), ($upstream_exon->end     + 2 ) );
		$downstream_site = 
		    $slice->subseq( ($downstream_exon->start - 2), ($downstream_exon->start - 1 ) );
	    };
	    unless ( $upstream_site && $downstream_site ){
		print STDERR "problems retrieving sequence for splice sites\n$@";
		next INTRON;
	    }
      
	    #print STDERR "upstream $upstream_site, downstream: $downstream_site\n";
	    ## good pairs of upstream-downstream intron sites:
	    ## ..###GT...AG###...   ...###AT...AC###...   ...###GC...AG###.
	    
	    ## bad  pairs of upstream-downstream intron sites (they imply wrong strand)
	    ##...###CT...AC###...   ...###GT...AT###...   ...###CT...GC###...
      
	    if (  ($upstream_site eq 'GT' && $downstream_site eq 'AG') ||
		  ($upstream_site eq 'AT' && $downstream_site eq 'AC') ||
		  ($upstream_site eq 'GC' && $downstream_site eq 'AG') ){
	#	print SPLICEUP "$upstream_site1\n";
	#	print SPLICEDOWN "$downstream_site1\n";
		$correct++;
	    }
	} # end of INTRON
  }
  elsif ( $strand == -1 ){
    
    #  example:
    #                                  ------CT...AC---... 
    #  transcript in reverse strand -> ######GA...TG###... 
    # we calculate AC in the slice and the revcomp to get GT == good site
    
  INTRON:
    for (my $i=0; $i<$#exons; $i++ ){
      my $upstream_exon   = $exons[$i];
      my $downstream_exon = $exons[$i+1];
      my $upstream_site;
      my $downstream_site;
      my $up_site;
      my $down_site;

      my $upstream_site1;
      my $downstream_site1;
      my $up_site1;
      my $down_site1;
      
      eval{
	$up_site1 = 
	  $slice->subseq( ($upstream_exon->start - 13), ($upstream_exon->start - 1) );
	$down_site1 = 
	  $slice->subseq( ($downstream_exon->end + 1), ($downstream_exon->end + 13 ) );
      };
      unless ( $up_site1 && $down_site1 ){
	print STDERR "problems retrieving sequence for splice sites\n$@";
	next INTRON;
    }
      ( $upstream_site1   = reverse(  $up_site1  ) ) =~ tr/ACGTacgt/TGCAtgca/;
      ( $downstream_site1 = reverse( $down_site1 ) ) =~ tr/ACGTacgt/TGCAtgca/;
      

      eval{
	$up_site = 
	  $slice->subseq( ($upstream_exon->start - 2), ($upstream_exon->start - 1) );
	$down_site = 
	  $slice->subseq( ($downstream_exon->end + 1), ($downstream_exon->end + 2 ) );
      };
      unless ( $up_site && $down_site ){
	print STDERR "problems retrieving sequence for splice sites\n$@";
	next INTRON;
      }
      ( $upstream_site   = reverse(  $up_site  ) ) =~ tr/ACGTacgt/TGCAtgca/;
      ( $downstream_site = reverse( $down_site ) ) =~ tr/ACGTacgt/TGCAtgca/;
      
      #print STDERR "upstream $upstream_site, downstream: $downstream_site\n";
      if (  ($upstream_site eq 'GT' && $downstream_site eq 'AG') ||
	    ($upstream_site eq 'AT' && $downstream_site eq 'AC') ||
	    ($upstream_site eq 'GC' && $downstream_site eq 'AG') ){
	$correct++;
	#print SPLICEUP "$upstream_site1\n";
	#print SPLICEDOWN "$downstream_site1\n";
      }
     
    } # end of INTRON
  }

    if ($correct == $introns) {
	return 1;
    }
}
