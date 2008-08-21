package XrefMapper::homo_sapiens;

use  XrefMapper::BasicMapper;
use strict;
use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };

sub get_set_lists {

  return [["ExonerateGappedBest1", ["homo_sapiens","*"]]];

}

sub gene_description_filter_regexps {

  return ('^BA\S+\s+\(NOVEL PROTEIN\)\.?',
	  '^DJ\S+\s+\(NOVEL PROTEIN\)\.?',
	  '^LOC\d+\s*(PROTEIN)?\.?',
	  '^ORF.*',
	  '^PROTEIN C\d+ORF\d+\.*',
	  '\(CLONE \S+\)\s+',
	  '^BC\d+\_\d+\.?',
	  '^CGI\-\d+ PROTEIN\.?\;?',
	  '[0-9A-Z]{10}RIK PROTEIN[ \.]',
	  'R\d{5}_\d[ \.,].*',
	  'PROTEIN KIAA\d+[ \.].*',
	  'RIKEN CDNA [0-9A-Z]{10}[ \.]',
	  '^\(*HYPOTHETICAL\s+.*',
	  '^UNKNOWN\s+.*',
	  '^DKFZP[A-Z0-9]+\s+PROTEIN[\.;]?.*',
	  '^CHROMOSOME\s+\d+\s+OPEN\s+READING\s+FRAME\s+\d+\.?.*',
	  '^FKSG\d+\.?.*',
	  '^HSPC\d+\s+PROTEIN\.?.*',
	  '^KIAA\d+\s+PROTEIN\.?.*',
	  '^KIAA\d+\s+GENE\s+PRODUCT\.?.*',
	  '^HSPC\d+.*',
	  '^PRO\d+\s+PROTEIN\.?.*',
	  '^PRO\d+\.?.*',
	  '^FLJ\d+\s+PROTEIN.*',
	  '^PRED\d+\s+PROTEIN.*',
	  '^WUGSC:.*\s+PROTEIN\.?.*',
	  '^SIMILAR TO GENE.*',
	  '^SIMILAR TO PUTATIVE[ \.]',
	  '^SIMILAR TO HYPOTHETICAL.*',
	  '^SIMILAR TO (KIAA|LOC).*',
	  '^SIMILAR TO\s+$',
          '^WUGSC:H_.*',
          '^\s*\(?PROTEIN\)?\.?\s*$',
	  '^\s*\(?FRAGMENT\)?\.?\s*$',
          '^\s*\(?GENE\)?\.?\s*$',
	  '^\s*\(\s*\)\s*$',
          '^\s*\(\d*\)\s*[ \.]$');

}


sub get_canonical_name{
   return "HGNC";
}

sub species_specific_cleanup{
  my $self = shift;
  my $dbname = $self->get_canonical_name;

  print "Removing all $dbname from object_xref not on a Gene\n";
  my $remove_old_ones = (<<JSQL);
delete ox 
  from object_xref ox, xref x, external_db e
    where e.db_name like "$dbname" and 
          ox.ensembl_object_type != "Gene" and
          ox.xref_id = x.xref_id and
	  x.external_db_id = e.external_db_id;
JSQL

  #
  # First Delete all the hgnc object_xrefs not on a gene. (i.e these are copys).
  #

  my $sth = $self->core->dbc->prepare($remove_old_ones);

  $sth->execute() || die "Could not execute: \n$remove_old_ones \n";

  $sth->finish;

}


# For human we want to make a copy of the HGNC references on the genes and put them on 
# the "canonical" transcripts

sub species_specific_pre_attributes_set{
  my $self = shift;
  my $dbname = $self->get_canonical_name();

  my ($max_object_xref_id, $max_xref_id);

  my $sth = $self->core->dbc->prepare("SELECT MAX(object_xref_id) FROM object_xref");
  $sth->execute();
  $sth->bind_columns(\$max_object_xref_id);
  $sth->fetch;


  $sth = $self->core->dbc->prepare("SELECT MAX(xref_id) FROM xref");
  $sth->execute();
  $sth->bind_columns(\$max_xref_id);
  $sth->fetch;

  my $object_xref_id = $max_object_xref_id + 1;

  if($object_xref_id == 1){
    die "max_object_xref_id should not be 1\n";
  }

  my $sql = "select gene_id, canonical_transcript_id from gene";

  my $object_sql = (<<FSQL);
select x.xref_id, o.ensembl_id
  from xref x, external_db e, object_xref o
    where x.external_db_id = e.external_db_id and 
      e.db_name like "$dbname" and 
      o.xref_id = x.xref_id  and
      o.ensembl_object_type = "Gene";
FSQL


  my $file = $self->core->dir()."/hgnc_transcript_data.sql";

  open(HGNC_TRAN,">$file") || die "Could not open $file";

  my $sth = $self->core->dbc->prepare($sql);
  
  $sth->execute();
  my ($gene_id, $tran_id);
  $sth->bind_columns(\$gene_id,\$tran_id);
  my %gene_to_tran;
  while ($sth->fetch){
    $gene_to_tran{$gene_id} = $tran_id;
  }
  $sth->finish;

  $sth = $self->core->dbc->prepare($object_sql);

  $sth->execute();
  my ($xref_id);
  $sth->bind_columns(\$xref_id, \$gene_id);

  while ($sth->fetch){
    if(defined($gene_to_tran{$gene_id})){
      print HGNC_TRAN $object_xref_id++,"\t",
	$gene_to_tran{$gene_id},
	  "\tTranscript\t",
	    $xref_id,"\\N\n";
    }
    else{
      print STDERR "Could not find canonical for gene $gene_id\n";
    }
  }
  $sth->finish;

  close (HGNC_TRAN);

#
# import data
#
  $sth = $self->core->dbc->prepare("LOAD DATA LOCAL INFILE \'$file\' IGNORE INTO TABLE object_xref");
  print "Uploading data in $file to object_xref\n";
  $sth->execute() || die "error loading file $file\n";


##########################################################################################
# HGNC synonyms are special as they are only on some prioritys so copy all the synonyms
# from the xref database to the core database. Granted some will already be there but use
# update ignore just in case.
###########################################################################################


#get the xref synonyms store as a hash of arrays.


my $sql = 'select distinct(s.synonym), x.accession from xref x, synonym s, source so where so.source_id = x.source_id and x.xref_id = s.xref_id and so.name like "HGNC"';

#+--------------------------+-----------+
#| synonym                  | accession |
#+--------------------------+-----------+
#| FWP007                   | 7         |
#| S863-7                   | 7         |
#| CPAMD5                   | 7         |
#| AACT                     | 16        |
#| ACT                      | 16        |


  $sth = $self->xref->dbc->prepare($sql);
  
  $sth->execute();
  my ($syn, $acc);
  $sth->bind_columns(\$syn,\$acc);
  my %acc_to_syn;
  while ($sth->fetch){
    push @{$acc_to_syn{$acc}}, $syn;
  }
  $sth->finish;

# get the HGNC xrefs in the core and add the synonyms, plus cretae hash to get id from the display name
  
  my %display_label_to_id;

  my $hgnc_file  =  $self->core->dir()."/hgnc_syn.txt";
  open(SYN, ">$hgnc_file" ) || die "Could not open file $hgnc_file";

  $sql = 'select x.dbprimary_acc, x.xref_id, x.display_label from xref x, external_db e where e.external_db_id = x.external_db_id and e.db_name like "HGNC"';


  $sth = $self->core->dbc->prepare($sql);
  
  $sth->execute();
  my ($xref_id, $display_label);
  $sth->bind_columns(\$acc,\$xref_id,\$display_label);
  while($sth->fetch){
    $display_label_to_id{$display_label} = $acc;
    if(defined($acc_to_syn{$acc})){
      foreach my $a ( @{$acc_to_syn{$acc}} ){
	print SYN "$xref_id\t$a\n";
      }
    }
  }
  $sth->finish;

  close SYN;



#
# import hgnc synonyms
#
  $sth = $self->core->dbc->prepare("LOAD DATA LOCAL INFILE \'$hgnc_file\' IGNORE INTO TABLE external_synonym");
  print "Uploading data in $hgnc_file to external_synonym\n";
  $sth->execute() || die "error loading file $hgnc_file\n";
  $sth->finish;


#######################
#Do the naming bit now.
#######################

# get the vega external_sources


  my $xref_file        =  $self->core->dir()."/xref_vega_names.txt";
  my $object_xref_file =  $self->core->dir()."/object_xref_vega_names.txt";

  open(XREF, ">$xref_file" ) || die "Could not open file $xref_file";
  open(OBJECT_XREF, ">$object_xref_file") || die "Could not open file $object_xref_file";


#  my ($vega_gene_id, $vega_transcript_id, $vega_gene_like_id, $vega_transcript_like_id);

  my ($hgnc_curated_gene_id, $hgnc_automatic_gene_id, $clone_based_vega_gene_id, $clone_based_ensembl_gene_id);
  my ($hgnc_curated_tran_id, $hgnc_automatic_tran_id, $clone_based_vega_tran_id, $clone_based_ensembl_tran_id);

  $sth = $self->core->dbc->prepare("select external_db_id from external_db where db_name like ?");
  
  $sth->execute("HGNC_curated_gene");
  $sth->bind_columns(\$hgnc_curated_gene_id);
  $sth->fetch;
  
  $sth->execute("HGNC_automatic_gene");
  $sth->bind_columns(\$hgnc_automatic_gene_id);
  $sth->fetch;

  $sth->execute("Clone_based_vega_gene");
  $sth->bind_columns(\$clone_based_vega_gene_id);
  $sth->fetch;
  
  $sth->execute("Clone_based_ensembl_gene");
  $sth->bind_columns(\$clone_based_ensembl_gene_id);
  $sth->fetch;
  
  if(!defined($hgnc_curated_gene_id)){
    die "Could not find external database name HGNC_curated_gene\n";
  }
  if(!defined($hgnc_automatic_gene_id)){
    die "Could not find external database name HGNC_automatic_gene\n";
  }
  if(!defined($clone_based_vega_gene_id)){
    die "Could not find external database name Clone_based_vega_gene\n";
  }
  if(!defined($clone_based_ensembl_gene_id)){
    die "Could not find external database name Clone_based_ensembl_gene\n";
  }




  $sth->execute("HGNC_curated_transcript");
  $sth->bind_columns(\$hgnc_curated_tran_id);
  $sth->fetch;
  
  $sth->execute("HGNC_automatic_transcript");
  $sth->bind_columns(\$hgnc_automatic_tran_id);
  $sth->fetch;

  $sth->execute("Clone_based_vega_transcript");
  $sth->bind_columns(\$clone_based_vega_tran_id);
  $sth->fetch;
  
  $sth->execute("Clone_based_ensembl_transcript");
  $sth->bind_columns(\$clone_based_ensembl_tran_id);
  $sth->fetch;
  
  if(!defined($hgnc_curated_tran_id)){
    die "Could not find external database name HGNC_curated_transcript\n";
  }
  if(!defined($hgnc_automatic_tran_id)){
    die "Could not find external database name HGNC_automatic_transcript\n";
  }
  if(!defined($clone_based_vega_tran_id)){
    die "Could not find external database name Clone_based_vega_transcript\n";
  }
  if(!defined($clone_based_ensembl_tran_id)){
    die "Could not find external database name Clone_based_ensembl_transcript\n";
  }



  ###########################
  # Delete the old ones.
  ###########################


  my $del_vega_sql = "delete o from object_xref o, xref x where x.xref_id = o.xref_id and x.external_db_id in ( $hgnc_curated_gene_id, $hgnc_automatic_gene_id, $clone_based_vega_gene_id, $clone_based_ensembl_gene_id,$hgnc_automatic_tran_id, $clone_based_ensembl_tran_id)";

  $sth = $self->core->dbc->prepare($del_vega_sql);
  $sth->execute();
 
  $del_vega_sql = "delete x from xref x where x.external_db_id in ($hgnc_curated_gene_id, $hgnc_automatic_gene_id, $clone_based_vega_gene_id, $clone_based_ensembl_gene_id,$hgnc_automatic_tran_id, $clone_based_ensembl_tran_id)";

  $sth = $self->core->dbc->prepare($del_vega_sql);
  $sth->execute();



  ######################################################
  # Get the current max values for xref and object_xref
  ######################################################


  my $ensembl = $self->core;
  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-dbconn => $ensembl->dbc);
  my $gene_adaptor = $db->get_GeneAdaptor();

  if(!defined($max_xref_id) or $max_xref_id == 0){
    die "Sorry dont belives there are no xrefs  max = 0\n";
  }
  if(!defined($max_object_xref_id) or $max_object_xref_id == 0){
    die "Sorry dont belives there are no object_xrefs  max = 0\n";
  }


  #####################################
  # get synonyms to test NONE problems
  #####################################

  my %synonym;
  $sth = $self->core->dbc->prepare('select es.synonym, x.display_label from external_synonym es, xref x, external_db e where x.xref_id = es.xref_id and x.external_db_id = e.external_db_id and e.db_name = "'.$dbname.'"' );
  $sth->execute();
  my ($syn, $name);
  $sth->bind_columns(\$syn,\$name);
  while($sth->fetch){
    $synonym{$syn} = $name;
  }
  $sth->finish;

  #Some hgnc synonyms are missing so use entrez gene instead a hack i know but will be fixed later.
  $sth = $self->core->dbc->prepare('select es.synonym, x.display_label from external_synonym es, xref x, external_db e where x.xref_id = es.xref_id and x.external_db_id = e.external_db_id and e.db_name = "EntrezGene"' );
  $sth->execute();
  $sth->bind_columns(\$syn,\$name);
  while($sth->fetch){
    $synonym{$syn} = $name;
  }
  $sth->finish;

  ###########################
  # Process each Gene
  ###########################

  my @genes = @{$gene_adaptor->fetch_all()};
  
  my $total_gene_vega = 0;
  my $total_gene = 0;
  my $total_clone_name = 0;
  my $hgnc_count = 0;
  
  while (my $gene = shift @genes){
    #    if($gene->biotype ne "protein_coding"){
    #    $types{$gene->biotype}++;
    #      next;
    #  }
    my @dbentries = @{$gene->get_all_DBEntries()};
    my @HGNC=();
    my @VEGA_NAME=();
    my $CLONE_NAME = undef;
    foreach my $dbe (@dbentries){
      if($dbe->dbname eq "HGNC"){
	push @HGNC, $dbe->display_id;
      }
    }
    if(scalar(@HGNC)){
      $hgnc_count++;
    }
#    my @has_vega;
    my %no_vega; # hash now as we want to sort by stable id and hence need a key value pair
#    my @no_vega;
    my $vega_count = 0;
    foreach my $tr (@{$gene->get_all_Transcripts}){
      my $VEGA = undef;
      
      my $count = 0;
      foreach my $dbe (@{$tr->get_all_DBEntries}){
	if($dbe->dbname eq "HGNC_curated_transcript" or $dbe->dbname eq "Clone_based_vega_transcript"){
	  my($hgnc_bit, $num) = split(/-0\d\d/,$dbe->display_id);
	  if($hgnc_bit =~ /[.]/){
	    $CLONE_NAME = $hgnc_bit;
	  }
	  else{
	    $VEGA = $hgnc_bit;
	    $vega_count++;
	    my $multiple = 1;
	    if(scalar(@VEGA_NAME)){
	      foreach my $name (@VEGA_NAME){
		if($name ne $VEGA){
		  if(defined($synonym{$name}) and uc($synonym{$name}) eq $VEGA){
		    $multiple = 0;
		  }
		  if(defined($synonym{$VEGA}) and uc($synonym{$VEGA}) eq uc($name)){
		    $multiple = 0;
		  }
		}
		else{
		  $multiple = 0;
		}
	      }
	      if($multiple){
		push @VEGA_NAME , $hgnc_bit;
	      }
	    }
	    else{
	      push @VEGA_NAME , $hgnc_bit;
	    }
	  }	
	  $count++;
	}      
      } # end of dbentries fro this transcript
      if($count == 0){
	$no_vega{$tr->stable_id} = $tr->dbID;
      }
      if($count > 1){
	print "Problem: ".$tr->stable_id." has more than one vega_transcript\n";
      }
      if($count == 1){
	if(defined($VEGA) and scalar(@HGNC)){
	  my $found = 0;
	  foreach my $hgnc_name (@HGNC){
	    if(uc($VEGA) eq uc($hgnc_name)){
	      $found = 1;
	    }
	    elsif(defined($synonym{$VEGA}) and uc($synonym{$VEGA}) eq uc($hgnc_name)){
	      $found = 1;
	    }
	  }
	  if(!$found){
	    print "Problem: ".$gene->stable_id." linked to hgnc (".join(', ',@HGNC).")   BUT ".$tr->stable_id." linked to vega_transcript $VEGA????\n";	
	  }
	}
#	push @has_vega, $tr->dbID;
      }
    } # end for each transcript
    ####################################################################################
    # if there is at least one transcript
    # set vega_gene
    # loop through no_vega array and set vega_transcript_like for each starting at 101
    ####################################################################################
    if(scalar(@VEGA_NAME) > 1){
      print "Warning: gene ".$gene->stable_id." has more than one vega_transcript these are (".join(', ',@VEGA_NAME).")\n";
    }	
    if($vega_count){
      foreach my $name (@VEGA_NAME){
	$max_xref_id++;
	$max_object_xref_id++;
	my $id = $display_label_to_id{$name};
	if(!defined($id)){
	  $id = $name;
	  print "Warning Could not find id for $name\n";
	}
	print XREF  $max_xref_id . "\t"
	  . $hgnc_curated_gene_id. "\t" . $id . "\t".$name."\t" . "0" . "\t". "\n" ;
	print OBJECT_XREF "$max_object_xref_id\t".$gene->dbID."\tGene\t" .$max_xref_id . "\t\\N\n";
      }

      my $name = $VEGA_NAME[0];
      my $tran_name_ext = 201;
      foreach my $tran (sort keys %no_vega){
	$max_xref_id++;
	$max_object_xref_id++;
	print XREF  $max_xref_id . "\t"
	  . $hgnc_automatic_tran_id. "\t" . $name."-".$tran_name_ext . "\t".$name."-".$tran_name_ext."\t" . "0" . "\t". "\n" ;
	print OBJECT_XREF "$max_object_xref_id\t".$no_vega{$tran}."\tTranscript\t" .$max_xref_id . "\t\\N\n";
	$tran_name_ext++;
      }
    }
    
    ####################################################################################
    # if no vega_transcript but hgnc
    # set vega_gene_like to hgnc
    # loop through both arrays and set vega_transcript_like to hgnc-101 etc
    ####################################################################################
    elsif(scalar(@HGNC)){
      foreach my $name (@HGNC){
	$max_xref_id++;
	$max_object_xref_id++;
	my $id = $display_label_to_id{$name};
	if(!defined($id)){
	  $id = $name;
	  print "Warning Could not find id for $name\n";
	}
	print XREF  $max_xref_id . "\t"
	  . $hgnc_automatic_gene_id. "\t" . $id . "\t".$name."\t" . "0" . "\t". "\n" ;
	print OBJECT_XREF "$max_object_xref_id\t".$gene->dbID."\tGene\t" .$max_xref_id . "\t\\N\n";
      }

      my $name = $HGNC[0];
      my $tran_name_ext = 201;
      foreach my $tran (sort keys %no_vega){
	$max_xref_id++;
	$max_object_xref_id++;
	print XREF  $max_xref_id . "\t"
	  . $hgnc_automatic_tran_id. "\t" . $name."-".$tran_name_ext . "\t".$name."-".$tran_name_ext."\t" . "0" . "\t". "\n" ;
	print OBJECT_XREF "$max_object_xref_id\t".$no_vega{$tran}."\tTranscript\t" .$max_xref_id . "\t\\N\n";
	$tran_name_ext++;
      }
    }	
    
    
    ####################################################################################
    # if no vega_transcript and no hgnc use clone name
    # set vega_gene_like to clone name
    # loop through both arrays and set vega_transcript_like to clone_name-101 etc
    ####################################################################################
    else{
      if(defined($CLONE_NAME)){
	$max_xref_id++;
	$max_object_xref_id++;
	print XREF  $max_xref_id . "\t"
	  . $clone_based_vega_gene_id. "\t" . $CLONE_NAME . "\t".$CLONE_NAME."\t" . "0" . "\t". "\n" ;
	print OBJECT_XREF "$max_object_xref_id\t".$gene->dbID."\tGene\t" .$max_xref_id . "\t\\N\n";
	
	my $tran_name_ext = 201;
	foreach my $tran (sort keys %no_vega){
	  $max_xref_id++;
	  $max_object_xref_id++;
	  print XREF  $max_xref_id . "\t"
	    . $clone_based_vega_tran_id. "\t" . $CLONE_NAME."-".$tran_name_ext . "\t".$CLONE_NAME."-".$tran_name_ext."\t" . "0" . "\t". "\n" ;
	  print OBJECT_XREF "$max_object_xref_id\t".$no_vega{$tran}."\tTranscript\t" .$max_xref_id . "\t\\N\n";
	  $tran_name_ext++;
	}
	
      }
      else{
	#get the clone name
	my $new_gene = $gene->transform('clone');
	my $new_clone_name = undef;
	if(defined($new_gene)){
	  $new_clone_name = $new_gene->slice->seq_region_name;
	}
	else{
	  # on more than one clone?? try  
	  my $slice = $gene->slice->sub_Slice($gene->start,$gene->end,$gene->strand);
	  my $clone_projection = $slice->project('clone');
	  foreach my $seg (@$clone_projection) {
	    my $clone = $seg->to_Slice();
	    $new_clone_name = $clone->seq_region_name;
	  }
	  if(!defined($new_clone_name)){
	    print "PROJECT failed for ".$gene->stable_id."\n";
	    next;
	  }
	}
	# store the data
	$max_xref_id++;
	$max_object_xref_id++;
	print XREF  $max_xref_id . "\t"
	  . $clone_based_ensembl_gene_id. "\t" . $new_clone_name . "\t".$new_clone_name."\t" . "0" . "\t". "\n" ;
	print OBJECT_XREF "$max_object_xref_id\t".$gene->dbID."\tGene\t" .$max_xref_id . "\t\\N\n";
	
	my $tran_name_ext = 201;
	foreach my $tran (sort keys %no_vega){
	  $max_xref_id++;
	  $max_object_xref_id++;
	  print XREF  $max_xref_id . "\t"
	    . $clone_based_ensembl_tran_id. "\t" . $new_clone_name."-".$tran_name_ext . "\t".$new_clone_name."-".$tran_name_ext."\t" . "0" . "\t". "\n" ;
	  print OBJECT_XREF "$max_object_xref_id\t".$no_vega{$tran}."\tTranscript\t" .$max_xref_id . "\t\\N\n";
	  $tran_name_ext++;
	}
	

	
      }


    }
    
    if($vega_count){
      $total_gene_vega++;
    }
    $total_gene++;
    #  print "Finished Gene ".$gene->stable_id."\t$HGNC\t $vega_count\n"
  }

  close XREF;
  close OBJECT_XREF;



#
# import data
#
  $sth = $self->core->dbc->prepare("LOAD DATA LOCAL INFILE \'$xref_file\' IGNORE INTO TABLE xref");
  print "Uploading data in $file to xref\n";
  $sth->execute() || die "error loading file $file\n";

  $sth = $self->core->dbc->prepare("LOAD DATA LOCAL INFILE \'$object_xref_file\' IGNORE INTO TABLE object_xref");
  print "Uploading data in $file to object_xref\n";
  $sth->execute() || die "error loading file $file\n";
  
}






1;
