package XrefMapper::homo_sapiens;

use  XrefMapper::BasicMapper;

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
# the "canonical" transcripts. For now this means, go through all the HGNC xrefs which are all on genes
# and if they are DIRECT ENSG then copy to the longest transcript of that gene. If they are DEPENDENT then 
# copy xref to the transcript the dependent one is on.

sub species_specific_pre_attributes_set{
  my $self = shift;
  my $dbname = $self->get_canonical_name();

  my $object_xref_id = $max_object_xref_id + 1;


  my $direct_sql = (<<ESQL);
select x.xref_id, x.info_text
  from xref x, external_db e
    where x.external_db_id = e.external_db_id and 
      e.db_name like "$dbname" and
      x.info_type = "DIRECT"
ESQL


  my $dependent_sql = (<<FSQL);
select x.xref_id, x.info_text
  from xref x, external_db e
    where x.external_db_id = e.external_db_id and 
      e.db_name like "$dbname" and
      x.info_type = "DEPENDENT"
FSQL

  my $get_gene_stable_id = (<<GSQL);
select gsi.stable_id, t.transcript_id, tr.translation_id
  from gene_stable_id gsi, transcript t 
    left join translation tr on tr.transcript_id = t.transcript_id
       where  gsi.gene_id     = t.gene_id 
GSQL

  my $longest_tran_sql = (<<HSQL);
select gsi.gene_id , gsi.stable_id, t.transcript_id, (t.seq_region_end - t.seq_region_start) as size
  from gene_stable_id gsi, transcript t 
    where t.gene_id = gsi.gene_id order by gsi.stable_id, size 
HSQL


  my $object_xref_sql = (<<ISQL);
select ox.xref_id, ox.ensembl_object_type, ox.ensembl_id
  from object_xref ox, xref x
    where x.dbprimary_acc like ? and x.xref_id = ox.xref_id 
ISQL


  my $translation_sql = "select transcript_id, translation_id from translation";
 

 
  my $tsi_2_id_sql = "select stable_id, transcript_id from transcript_stable_id";



#
# find the longest transcript for each gene.
#

  my %stable_id_2_transcript_id;
  my %gene_id_2_transcript_id;
  my $translation_id_2_transcript_id;

  my $sth = $self->core->dbc->prepare($longest_tran_sql);

  $sth->execute();
  my ($g_id, $gsi, $t_id, $size);
  $sth->bind_columns(\$g_id, \$gsi, \$t_id, \$size);


  
  while ($sth->fetch){
    $stable_id_2_transcript_id{$gsi} = $t_id;
    $gene_id_2_transcript_id{$g_id} = $t_id;
  }
  $sth->finish;


# now add the transcript_stable_ids

  $sth = $self->core->dbc->prepare($tsi_2_id_sql);

  my ($tsi, $tran_id);
  $sth->execute();
  $sth->bind_columns(\$tsi, \$tran_id);

  while ($sth->fetch){
    $stable_id_2_transcript_id{$tsi} = $tran_id;
  }
  $sth->finish;


# what about translations ??? Add if needed.
# but they should never be on the translations.

#
# open file to store new xref matches
#
  my $file = $self->core->dir()."/hgnc_transcript_data.sql";

  open(HGNC_TRAN,">$file") || die "Could not open $file";


#
# Process the DIRECT ones first.
#
  $sth = $self->core->dbc->prepare($direct_sql);

  $sth->execute();

  my ($xref_id, $text);
  $sth->bind_columns(\$xref_id, \$text);

  while($sth->fetch){

    if($text =~ /between (\S+) and /){
      my $stable_id =  $1;
      if(defined($stable_id_2_transcript_id{$stable_id})){
	print HGNC_TRAN $object_xref_id++,"\t",
	  $stable_id_2_transcript_id{$stable_id},
	  "\tTranscript\t",
	  $xref_id,"\\N\n";
      }
      else{
	print STDERR "Could not find transcript id for stable id $stable_id\n";
      }
    }
    else{
      print STDERR "$text\n failed regular expression test\n";
    }
  }


#
# Process the dependent ones.
#

  ( $gsi, $transcript_id, $translation_id ) = ( undef, undef, undef );
  $sth = $self->core->dbc->prepare($get_gene_stable_id);
  $sth->execute();

  $sth->bind_columns(\$gsi, \$transcript_id, \$translation_id);
  my %id_2_gsi;
  while($sth->fetch){
    $id_2_gsi{'Transcript'.$transcript_id} = $gsi;
#    print STDERR "LOADING: Transcript $transcript_id -> $gsi\n";
    if(defined($translation_id)){
      $id_2_gsi{'Translation'.$translation_id} = $gsi;
    }
  }
  $sth->finish;


  $sth = $self->core->dbc->prepare($dependent_sql);
  $sth->execute();

  $sth->bind_columns(\$xref_id, \$text);

  while($sth->fetch){

    if($text =~ /Generated via (\S+)/){
      $acc_2_xref_id{$1} = $xref_id;
    }
  }
  $sth->finish;



  $sth = $self->core->dbc->prepare($translation_sql);
  $sth->execute();

  my ($transcript_id, $translation_id);

  $sth->bind_columns(\$transcript_id, \$translation_id);

  while($sth->fetch){

      $translation_id_2_transcript_id{$translation_id} = $transcript_id;
  }
  $sth->finish;




  ( $xref_id, $xref_type, $ensembl_id ) = ( undef, undef, undef );

#
# Need to change this so that only get object_xrefs fro this genes transcripts and translations
#
#
  $sth = $self->core->dbc->prepare($object_xref_sql);


  foreach my $acc (keys %acc_2_xref_id){
    $sth->execute($acc) || die "Could not execute $!\n";
    $sth->bind_columns(\$xref_id, \$xref_type, \$ensembl_id);
    my $i = 0;
    while($sth->fetch){
      $i++
    }

    # matched to more than one transcript/translation so put later
    # for canonicals "There can only be one".  
    if($i > 0){ 
      my $gsi = $id_2_gsi{$xref_type.$ensembl_id};
      if(defined($gsi)){
	if(defined( $stable_id_2_transcript_id{$gsi})){
	  print HGNC_TRAN $object_xref_id++,"\t",
	    $stable_id_2_transcript_id{$gsi},
	      "\tTranscript\t",
		$acc_2_xref_id{$acc},"\\N\n";
	}
	else{
	  print STDERR "Could not get transcript_id for gene stable id $gsi\n";
	}
      }
      else{
	print STDERR "Could not find gene stable id for $xref_type $ensembl_id\n";
      }
      next;
    }
    
 

    if ($xref_type =~/Transcript/){
      print HGNC_TRAN $object_xref_id++,"\t",
	$ensembl_id,
	  "\tTranscript\t",
	    $acc_2_xref_id{$acc},"\\N\n";
    }
    elsif($xref_type=~/Translation/){
      if(defined( $translation_id_2_transcript_id{$ensembl_id})){
	print HGNC_TRAN $object_xref_id++,"\t",
	  $translation_id_2_transcript_id{$ensembl_id},
	    "\tTranscript\t",
	      $acc_2_xref_id{$acc},"\\N\n";
      }
      else{
	print STDERR "Could not find transcript id for translation_id $ensembl_id\n";
      }
    }
    elsif($xref_type=~/Gene/){
      if(defined( $gene_id_2_transcript_id{$ensembl_id})){
        print HGNC_TRAN $object_xref_id++,"\t",
	  $gene_id_2_transcript_id{$ensembl_id},
	    "\tTranscript\t",
	      $acc_2_xref_id{$acc},"\\N\n";
      }
      else{
	print STDERR "Could not find transcript id for gene_id $ensembl_id\n";
      }
    }
    else{
      print STDERR "ensembl type NOT Gene,Transcript or translation but $xref_type\n";
    }
  }





  close (HGNC_TRAN);

#
# import data
#
  $sth = $self->core->dbc->prepare("LOAD DATA LOCAL INFILE \'$file\' IGNORE INTO TABLE object_xref");
  print "Uploading data in $file to object_xref\n";
  $sth->execute() || die "error loading file $file\n";

          
#
# Some dependent accs match to more than one transcript which may not necessarily be in the same
# gene. In this case we get genes with object_xref to XXXXs but no Transcripts for that gene due to the
# moving to the largest transcript idea.
# 

# 1) get  all the object_xrefs linking genes to hugo.
# 2)  for each of these genes see if the transcripts have the same hugo
# 3) if not delete the object_xref for the gene one.

my $get_all_object_xref_genes =(<<KSQL);
   select ox.ensembl_id, ox.object_xref_id, ox.xref_id
     from object_xref ox, xref x, external_db e
       where ox.xref_id = x.xref_id and 
             x.external_db_id = e.external_db_id and
             e.db_name like '$dbname' and
             ox.ensembl_object_type = 'Gene'
KSQL

my $get_all_object_xref_transcripts_for_a_gene =(<<LSQL);
   select ox.xref_id
     from object_xref ox, xref x, external_db e, transcript t
       where ox.xref_id = x.xref_id and 
             x.external_db_id = e.external_db_id and
             t.transcript_id = ox.ensembl_id and
             e.db_name like '$dbname' and
             ox.ensembl_object_type = 'Transcript' and
             t.gene_id = ? 
LSQL




  my ($gene_id, $object_xref_id, $xref_id);
  $sth = $self->core->dbc->prepare($get_all_object_xref_genes);
  $sth->execute();

  $sth->bind_columns(\$gene_id, \$object_xref_id, \$xref_id);

  my %gene_id_2_object_xref;
#  my %gene_id_2_xref_id;
  while($sth->fetch){
    $gene_id_2_object_xref{$gene_id."|".$xref_id} = $object_xref_id;
#    $gene_id_2_xref_id{$gene_id}     = $xref_id;
     
  }
  $sth->finish;


  $file = $self->core->dir()."/hgnc_transcript_data_cleanup.sql";

  open(HGNC_TRAN,">$file") || die "Could not open $file";

  foreach my $key (keys %gene_id_2_object_xref){
    my ($gene_id, $gene_xref_id)  = split(/\|/,$key) ;  
    $sth = $self->core->dbc->prepare($get_all_object_xref_transcripts_for_a_gene);
    $sth->execute($gene_id);
    $sth->bind_columns(\$xref_id);

    my $found = 0;
    while($sth->fetch){
      if($xref_id == $gene_xref_id){
         $found = 1;
      } 
    } 
    if(!$found){
      if(defined($gene_id_2_transcript_id{$gene_id})){
	print HGNC_TRAN $object_xref_id++,"\t",
	  $gene_id_2_transcript_id{$gene_id},
	  "\tTranscript\t",
	  $gene_xref_id,"\\N\n";
      }
      else{
	print STDERR "Could not find transcript id for gene id $gene_id\n";
        die "horribly";
      }
    }
  }


  close (HGNC_TRAN);

#
# import data
#

  $sth = $self->core->dbc->prepare("LOAD DATA LOCAL INFILE \'$file\' IGNORE INTO TABLE object_xref");
  print "Uploading data in $file to object_xref\n";
  $sth->execute() || die "error loading file $file\n";



}

1;
