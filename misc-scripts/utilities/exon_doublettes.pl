# This script uses the hardcoded db connection
# It checks for exons with the same coordinate and 
# different id and produces SQL to repair the db


use DBI;

$dsn = "dbi:mysql:host=ecs1c.sanger.ac.uk;database=idmap_apr01";
$db = DBI->connect( $dsn, "ensadmin" );

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
      print STDERR "Buggy Gene on exon ",$cols[2],"\n"; 
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

# which exons have the same coords ...

print STDERR "Genes read!\n";

@exlist = ();
foreach $g ( values %gh ) {
  foreach $t (  values %{$g->{transcript}} ) {
    foreach $e ( values %{$t->{exon}} ) {
      $e->{fpc} = $g->{fpc};
      $e->{transcript} = $t;
      $e->{gene} = $g;
      
      push( @exlist, $e );
    }
  }
}

print STDERR "Read ",scalar( @exlist), " exons.\n";
@sortedExons = sort { ( $a->{fpc} cmp $b->{fpc} ) ||
			( $a->{strand} <=> $b->{strand} ) ||
			  ( $a->{seqstart} <=> $b->{seqstart} ) ||
			    ( $a->{seqend} <=> $b->{seqend} ) } @exlist;

print STDERR "Sorted ",scalar( @sortedExons), " exons.\n";


$lastExon = shift( @sortedExons );
while ( @sortedExons ) {
  $ex = shift( @sortedExons );
  if( equal_Exon( $lastExon, $ex ) ) {
    if( $lastExon->{transcript} == $ex->{transcript} ) {
      print STDERR "Same Transcript error ",$ex->{name}, " ", $lastExon->{name},"\n";
    } elsif( $ex->{gene} == $lastExon->{gene} ) {
      rename_Exon( $ex, $lastExon );
      next;
    } else {
      rename_Exon( $ex, $lastExon );
      if( gene_rep( $ex ) ne gene_rep( $lastExon ) ) {
	merge_Gen( gene_rep( $ex ), gene_rep( $lastExon) );
	gene_rep( $ex, $lastExon->{gene}{name} );
      }
      next;
    }
  }
  $lastExon = $ex;
}


sub equal_Exon {
  my ( $e1, $e2 ) = @_;
  if(( $e1->{fpc} == $e2->{fpc} ) &&
     ( $e1->{strand} == $e2->{strand} ) &&
     ( $e1->{seqstart} == $e2->{seqstart} ) &&
     ( $e1->{seqend} == $e2->{seqend} )) {
    return 1;
  } else {
    return 0;
  }
}

sub gene_rep {
  my ( $ex, $gen ) = @_;
  my $rep = $ex->{gene}{name};
  
  while( exists $geneRep{$rep} ) {
    $rep = $geneRep{$rep};
  }

  if( defined $gen ) {
    $geneRep{$rep} = $gen;
  }
  return $rep;
}



sub rename_Exon {
  my ( $e1, $e2 ) = @_;
  print "DELETE FROM exon WHERE id=\"",$e1->{name},"\";\n";
  print "UPDATE exon_transcript SET exon=\"",$e2->{name},"\" WHERE exon=\"",
  $e1->{name},"\";\n";
  print "UPDATE translation SET start_exon=\"",$e2->{name},"\" WHERE start_exon=\"",
  $e1->{name},"\";\n";
  print "UPDATE translation SET end_exon=\"",$e2->{name},"\" WHERE end_exon=\"",
  $e1->{name},"\";\n";
  print "UPDATE supporting_feature SET exon=\"",$e2->{name},"\" WHERE exon=\"",
  $e1->{name},"\";\n";
}

sub merge_Gen {
  my ( $g1, $g2 ) = @_;

  
  print "UPDATE transcript SET gene=\"$g2\" WHERE gene=\"$g1\";\n";
  print "DELETE FROM gene WHERE id = \"$g1\";\n";
  print "DELETE FROM genetype WHERE gene_id = \"$g1\";\n";
}


exit;



