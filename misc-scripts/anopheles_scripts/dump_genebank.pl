#!/usr/local/ensembl/bin/perl

use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Bio::EnsEMBL::Exon;

my $host      = 'ecs1a';
my $dbuser    = 'ensro';
my $dbname    = 'anopheles_gambiae_core_19_2b';
my $dbpass    = '';

GetOptions(
	   'dbname:s'    => \$dbname,
	   'dbhost:s'    => \$host,
	   'dbuser:s'    => \$dbuser,
	   'dbpass:s'    => \$dbpass,
	   );
print STDERR "Connecting to $host, $dbname\n";
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $host,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					   );

my $clone_adapt = $db->get_CloneAdaptor();
my $gene_adapt = $db->get_GeneAdaptor();
my $slice_adapt = $db->get_SliceAdaptor();

#get handle for an appropriate mart database
$host = 'ensembldb';
$dbuser = 'anonymous';
$dbname = 'ensembl_mart_19_3';
print STDERR "Connecting to $host, $dbname\n";
my $martdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $host,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					   );

my %scafmap;
my %old_prot_ann;
my %kill;
my %gene_name;
my %utr;

# get the old (Celera) scaffold names
open (SCAFMAP,"/ecs2/work6/mh4/gbdump_19_2b/input/accessions") || die "can't open SCAFMAP file";
while(<SCAFMAP>) {
    chomp;
    my ($new,$a,$old) = split;
    $scafmap{$new} = $old;
}
close SCAFMAP;

open (OLD_PROT_ANN, "/ecs2/work6/mh4/gbdump_19_2b/input/protein_annotation_rev") || die "can't find old protein annotation file";
while (<OLD_PROT_ANN>) {
 chomp;
 my ($scaffold, $psi, $prot_id, $old_id, $line) = split;
 $old_prot_ann{$psi}->{'scaffold'} = $scaffold if $scaffold ne 'x';
 $old_prot_ann{$psi}->{'prot_id'} = $prot_id  if $prot_id ne 'x';
 $old_prot_ann{$psi}->{'old_id'} = $old_id if $old_id ne 'x';
 $old_prot_ann{$psi}->{'line'} = $line if $line ne 'x';
}
close OLD_PROT_ANN;

# read kill list of proteins that have 'moved' within a scaffold
open (KILL, "/ecs2/work6/mh4/gbdump_19_2b/input/jump_within_scaff_2v2a") || die "can't find kill list";
while (<KILL>) {
  my ($kill_tsi, $kill_psi) = split;
  $kill{$kill_psi}=1;
}
close KILL;

# read file dumped from gene name database linking sequence_ref (sequence names as used in core xrefs)
# to official gene name and description (and gene_id)
# making all the keys lc facilitates case insensitive matching later
open (GENE_NAME, "/ecs2/work6/mh4/gbdump_19_2b/input/gene_names_rev.txt") || die "can't find gene name file";
while (<GENE_NAME>) {
 chomp;
 my ($seq_ref, $gene_name_id, $symbol, $description) = split /\t/;
 $gene_name{lc $seq_ref}->{'real_seq_ref'} = $seq_ref;
 $gene_name{lc $seq_ref}->{'gene_name_id'} = $gene_name_id;
 $gene_name{lc $seq_ref}->{'symbol'} = $symbol;
 $gene_name{lc $seq_ref}->{'description'} = $description if (($description ne '')&&($description ne 'NULL'));
}
close GENE_NAME;


#Get all of the distinct scaffolds (could get just from assembly table??)
my $query1 = "select distinct(c.name), a.superctg_ori from clone c, assembly a where a.superctg_name = c.name";
my $sth1 = $db->prepare($query1);
$sth1->execute();

# temporary counters for gene stuff
my $genes_multiscaff=0;
my $genes_bacterial=0;
my $genes_kept=0;
my $genes_named=0;
my $orths=0;
my $gene_orth=0;

# temporary counters for protein stuff
my ($fetch_count, $count_old_id, $count_prot_id, $count_new,$proteins_named,$moved,$proteins_described);
$fetch_count=$count_old_id=$count_prot_id=$count_new=$proteins_named=$moved=$proteins_described=0;


while (my ($clone_name,$ori) = $sth1->fetchrow_array) {

  open (OUT,">/ecs2/work6/mh4/gbdump_19_2b/output/$clone_name.tbl") || die "can't open tbl file";
  open (SEQ,">/ecs2/work6/mh4/gbdump_19_2b/output/$clone_name.fsa") || die "can't open fsa file";

  my $slice = $slice_adapt->fetch_by_clone_accession($clone_name);
  if ($ori == -1) {
    $slice = $slice->invert();
  }
  print STDERR "CLONE: $clone_name\n";

  my $chr_name = $slice->chr_name;
  my $clone_seq = $slice->seq;
  my $length_clone = length($clone_seq);
  my $old_clone_name = $scafmap{$clone_name};

# different header for seq file if UNKN v. on a real chromosome

  if ($chr_name !~ /UNKN/) {
    print SEQ ">gnl|WGS:AAAB|$old_clone_name|gb|$clone_name [organism=Anopheles gambiae str. PEST] [tech=wgs] [chromosome=$chr_name]\n$clone_seq\n";
  }
  else {
    print SEQ ">gnl|WGS:AAAB|$old_clone_name|gb|$clone_name [organism=Anopheles gambiae str. PEST] [tech=wgs]\n$clone_seq\n";
  }
  print OUT ">Feature gnl|WGS:AAAB|$old_clone_name|gb|$clone_name\n";

  my @genes = @{$slice->get_all_Genes};

  foreach my $gene(@genes) {
# these are the 'original' start and end
    my $start = $gene->start;
    my $end = $gene->end;
#    print STDERR "ID: ".$gene->stable_id."\tSTRAND: ".$gene->strand."\tSTART: ".$start."\tEND: ".$end."\n";

# skip gene if it runs off the edge of the scaffold
    if (($start < 0) || ($start > $length_clone) || ($end < 0) || ($end > $length_clone)) {
      $genes_multiscaff++;
      next;
    }

# tag gene if bacterial contaminant, print as misc feature and skip to next gene
    elsif ($gene->type eq 'bacterial_contaminant') {
      $genes_bacterial++;
      if ($gene->strand == 1) {
	print OUT "$start\t$end\tmisc_feature\n";
      }
      else {
	print OUT "$end\t$start\tmisc_feature\n";
      }
      print OUT "\t\t\tnote\tpossible bacterial gene - sequence may be bacterial contaminant\n";
      next;
    }

    else {
      $genes_kept++;

      my $new_gene = &checks($gene);

# &checks strips UTR from any gene with UTR adjacent to non-Met start or non-stop end - also sets %utr  if are remainning (valid) UTRs
# note that $gene and $new_gene are the same thing after this point?
# I think both are the result of running &checks subroutine

      $start = $new_gene->start;
      $end = $new_gene->end;
      my $new_gene_dbID = $new_gene->dbID;
      my $gene_id = $gene->stable_id;

      my $cdna_start;
      my $cdna_end;
      my $coding_start;
      my $coding_end;

# print the gene coordinates after checking for UTRs
# where a UTR not present in the munged gene add < or > as appropriate

      if (($utr{$new_gene_dbID}->{'up'}==1)&&($utr{$new_gene_dbID}->{'down'}==1)) {
	if ($gene->strand == 1) {
	  print OUT "$start\t$end\tgene\n";
	}
	else {
	  print OUT "$end\t$start\tgene\n";
	}
      }

      elsif (($utr{$new_gene_dbID}->{'up'} == 1)&&($utr{$new_gene_dbID}->{'down'}!=1)) {
	if ($gene->strand == 1) {
	  print OUT "$start\t>$end\tgene\n";
	}
	else {
	  print OUT "$end\t>$start\tgene\n";
	}
      }

      elsif (($utr{$new_gene_dbID}->{'up'} != 1)&&($utr{$new_gene_dbID}->{'down'}==1)) {
	if ($gene->strand == 1) {
	  print OUT "<$start\t$end\tgene\n";
	}
	else {
	  print OUT "<$end\t$start\tgene\n";
	}
      }

      elsif (($utr{$new_gene_dbID}->{'up'} != 1)&&($utr{$new_gene_dbID}->{'down'}!=1)) {
	if ($gene->strand == 1) {
	  print OUT "<$start\t>$end\tgene\n";
	}
	else {
	  print OUT "<$end\t>$start\tgene\n";
	}
      }

# fetch official gene name(s), if any, but use only if a single one
# uses xrefs that are sequence_refs to get gene names from %gene_name hash
# note use of lc when comparing to the previously-made-lc keys of %gene_name
      my (@seq_names,%symbols,$gene_name);
      foreach my $db_entry_g (@{$gene->get_all_DBLinks()}) {
	if ($db_entry_g->dbname eq 'Anopheles_symbol') {
	  push (@seq_names, $db_entry_g->display_id);
	}
      }
      foreach my $name (@seq_names) {
	$symbols{$gene_name{lc $name}->{'symbol'}}+=1;
      }
      if ((scalar keys %symbols) == 1) {
	foreach my $symbol_name (keys %symbols) {
	  $gene_name = $symbol_name;
	  print OUT "\t\t\tgene\t$gene_name\n";
#	  print STDERR "$gene_id has single name $gene_name\n";
	}
	$genes_named++;
      }
#      else {	
#	foreach (keys %symbols) {
#	  print STDERR "$gene_id has multiple names $_\n";
#	}
#      }

# always print ensembl stable id as locus_tag
      print OUT "\t\t\tlocus_tag\t$gene_id\n";

# query Mart to get the display id for any drosophila orthologues
      my @homologues;
      my $query2 = "select display_id from agambiae_ensemblgene_homologs_dmelanogaster_dm where gene_stable_id like ?";
      my $sth2 = $martdb->prepare($query2);
      $sth2->execute($gene_id);
      while ((my $dros_name) = $sth2->fetchrow_array) {
	push @homologues,$dros_name;
	$orths++;
      }
      if (scalar @homologues >0) {
	$gene_orth++;
	print OUT "\t\t\tnote\tsimilar to Drosophila ";
	print OUT join(', ',@homologues),"\n";
      }

# get the transcript(s) and print mRNA and protein features

      my @new_transcripts = @{$new_gene->get_all_Transcripts};	
      foreach my $new_tr(@new_transcripts) {
	$fetch_count++;	
	&print_transcript_coordinates($new_tr);
	&print_translation_coordinates($new_tr,$db,$clone_name,$gene_name);
      }
# temporary progress report
    if ($fetch_count % 100 == 0) {
      print STDERR "Fetched $fetch_count proteins so far\n";
    }
# end of else (i.e. this real gene)
    }
# end of gene loop
  }
# end of this scaffold
  close(OUT);
  close(SEQ);
}

print STDERR "Finished\nRemoved as multiscaff: $genes_multiscaff\tRemoved as bacterial: $genes_bacterial\tKept: $genes_kept\nWith names: $genes_named\tGenes with Dros orths: $gene_orth\tTotal dros orthologues: $orths\n";
print STDERR "Fetched proteins: $fetch_count, have 'moved' $moved, named $proteins_named, described $proteins_described\n";
print STDERR "Assigned $count_old_id old & prot ids, $count_prot_id prot id only, $count_new with neither\n";


## checks subroutine ##
# takes gene and looks at its transcript(s)
# if the coding sequence has UTR next to non-Met start or non-stop end does something complicated!!
# maybe ... makes new set of exons and sets translation appropriately so that a revised transcript is made without those dodgy UTRs
# then removes transcripts from gene and replaces with new ones

sub checks {

  my ($gene) = @_;
  my $gene_dbID = $gene->dbID; #probably redundant
  my @transcripts = @{$gene->get_all_Transcripts};
  my @new_transcripts;

  foreach my $tr(@transcripts) {

    my $tr_dbID = $tr->dbID;
    my $c_start = $tr->cdna_coding_start;
    my $c_end = $tr->cdna_coding_end;
    my $spl_seq = $tr->spliced_seq;

    my $tl_start = $tr->translation->start;
    my $tl_end = $tr->translation->end;
    my $new_cdna_start;
    my $new_cdna_end;

#  my $pep = $tr->translate->seq;
# ensembl translate method now strips off any final stop *
# so hack to make this script work is to get translation via simple bio::seq object

  my $able_seq = $tr->translateable_seq;
  my $able = Bio::Seq->new( -display_id => '$translation_name',
                               -seq => $able_seq);
  my $pep = $able->translate()->seq;

# next part looks for non-M start or non-* end and
# if found makes $new_cDNA_start/end that have any UTR sliced off
# ?? saves any valid UTR offset as $tl_start or sets to 1 if none
# ?? bizarrely sets $tl_end to the position of the CDS end in transcript coordiantes - used later on

    if ($pep =~ /^M/) {
      $tl_start = $c_start;
      $new_cdna_start = 1;
    } else {
      $tl_start = 1;
      $new_cdna_start = $c_start;
    }

    $tl_end = $c_end - $c_start + $tl_start;

    if ($pep =~ /\*$/) {
      $new_cdna_end = length( $spl_seq );
    } else {
      $new_cdna_end = $c_end;
    }

# get a set of (genomic) exon coordinates, but only for the span defined above
    my @exon_coords = $tr->cdna2genomic($new_cdna_start,$new_cdna_end);
    my @new_exons;

# make new exons for each pair of coords
    for my $exon_coord ( @exon_coords ) {
      if ($exon_coord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
	# print STDERR $exon_coord->start."\t".$exon_coord->end."\t".$exon_coord->strand."\t".$exon_coord->id."\n";
		
	my $new_exon = new Bio::EnsEMBL::Exon($exon_coord->start,$exon_coord->end,$exon_coord->strand);
	$new_exon->contig($exon_coord->id);
	$new_exon->phase(0);
	$new_exon->end_phase(0);
	push (@new_exons,$new_exon);
	#print STDERR "NEW EX: $new_exon\n";
      }
    }

    my $tl_start_exon;
    my $tl_end_exon;
    my $seen;
	
# print STDERR "TRANSLATION: $tl_start\t$tl_end\n";

# run through new exons (should be in order??), find the ones with translation start and end, and make a copy

    foreach my $exon  (@new_exons ) {
      if( ($tl_start > $exon->length) && (! defined $seen) ) { #eh????
	#	print STDERR "Translation_start: $tl_start\n";
	$tl_start -= $exon->length;                            #eh????
      } elsif (! defined $seen) {
	$seen = 1;
	$tl_start_exon = $exon;
      }
      if(  $tl_end > $exon->length) {
	$tl_end -= $exon->length;
      } else {
	$tl_end_exon = $exon;
	last;
      }
# i.e. subtract exon length from the end-CDS position until reach last exon which will have coding end in or at end of it
    }

# remove existing exons from transcript and add the new ones
    $tr->flush_Exons();		
    foreach my $ne (@new_exons) {
      $tr->add_Exon($ne);
    }

# unsets transcript cdna_coding_start/end
# then sets properties of the transcript's translation using copied start/end exons, $tl_start from way before and $tr_end as munged above
	
    $tr->{'cdna_coding_start'} =undef;
    $tr->{'cdna_coding_end'} =undef;

# ?? not sure what this used to do - was already commented out
#    $tr->cdna_coding_start($new_cdna_start);
#    $tr->cdna_coding_end($new_cdna_end);	
#    print STDERR "CDNA start: ".$tr->cdna_coding_end."\t".$tl_start."\n";
#    print STDERR "EX: ".$tl_end_exon."\n";

    my $translation = $tr->translation;
    $translation->start_Exon($tl_start_exon);
    $translation->end_Exon($tl_end_exon);
    $translation->start($tl_start);
    $translation->end($tl_end);

    $tr->{'translation'} = [];
    $tr->translation($translation);
#    print STDERR "TR: .".$tr->translation->start."\n";

    push(@new_transcripts,$tr);	
  }

# removes existing transcripts from gene, and replaces with new one(s)
# also stores existance of UTRs in %utr

  $gene->{'_transcript_array'} =[];
	
  foreach my $new_tr(@new_transcripts) {
    my $cdna_length = length($new_tr->seq->seq);
    my $coding_start = $new_tr->cdna_coding_start;
    my $coding_end = $new_tr->cdna_coding_end;
    my $gene_dbid = $gene->dbID;
    my $tl_start = $new_tr->translation->start;
    my $coding_length = $coding_end - $coding_start + 1;

#    print STDERR $new_tr->stable_id."\tCoding start: $coding_start\tcoding end: $coding_end\tCoding length: $coding_length\tcDNA length: $cdna_length\n ";

    if ($tl_start > 1) {
      $utr{$gene_dbid}->{'up'} = 1;
    }
    if ($coding_end != $cdna_length) {
      $utr{$gene_dbid}->{'down'} = 1;
    }

    $gene->add_Transcript($new_tr);
  }
  return $gene;
}



sub print_transcript_coordinates {
# prints mRNA locations and mRNA product tag

  my ($tr) = @_;

  my $coding_start = $tr->cdna_coding_start;
  my $coding_end = $tr->cdna_coding_end;
  my $cdna_length = length($tr->seq->seq);
  my $count_ex = 0;
  my $tr_name = $tr->stable_id;
  my $tr_dbID = $tr->dbID;

  my @exons = @{$tr->get_all_Exons};
  my $nb_ex = scalar(@exons);
  my $strand = $exons[$nb_ex-1]->strand;

  foreach my $ex(@exons) {
    $count_ex++;	
    my $ex_start = $ex->start;
    my $ex_end = $ex->end;
#    print STDERR "$coding_start\t$coding_end\t$cdna_length\n";

#One exon case
    if ($nb_ex == 1) {
      if (($coding_start > 1)&& ($coding_end != $cdna_length)) {
	if ($strand == 1) {
	  print OUT "$ex_start\t$ex_end\tmRNA\n";
	}
	else {
	  print OUT "$ex_end\t$ex_start\tmRNA\n";
	}
      }

      elsif (($coding_start == 1)&& ($coding_end != $cdna_length)) {
	if ($strand == 1) {
	  print OUT "<$ex_start\t$ex_end\tmRNA\n";
	}
	else {
	  print OUT "<$ex_end\t$ex_start\tmRNA\n";
	}
      }

      elsif (($coding_start > 1)&& ($coding_end == $cdna_length)) {
	if ($strand == 1) {
	  print OUT "$ex_start\t>$ex_end\tmRNA\n";
	}
	else {
	  print OUT "$ex_end\t>$ex_start\tmRNA\n";
	}
      }

      elsif (($coding_start == 1)&& ($coding_end == $cdna_length)) {
	if ($strand == 1) {
	  print OUT "<$ex_start\t>$ex_end\tmRNA\n";
	}
	else {
	  print OUT "<$ex_end\t>$ex_start\tmRNA\n";
	}
      }
    }
	
#Multiple exons case	
    elsif ($nb_ex > 1) {

#First exon
      if ($count_ex == 1) {
#5' UTR
	if ($coding_start > 1) { 
	  if ($strand == 1) {
	    print OUT "$ex_start\t$ex_end\tmRNA\n";
	  }
	  else {
	    print OUT "$ex_end\t$ex_start\tmRNA\n";
	  }
	}
#no 5' UTR partial sequence
	elsif  ($coding_start == 1) {
	  if ($strand == 1) {
	    print OUT "<$ex_start\t$ex_end\tmRNA\n";
	  }
	  else {
	    print OUT "<$ex_end\t$ex_start\tmRNA\n";
	  }
	}
      }
#End first exon

#Last exon
      elsif ($count_ex == $nb_ex) {
	if ($coding_end == $cdna_length)  {
	  if ($strand == 1) {
	    print OUT "$ex_start\t>$ex_end\n";
	  }
	  else {
	    print OUT "$ex_end\t>$ex_start\n";
	  }
	}
	elsif($coding_end != $cdna_length)  {
	  if ($strand == 1) {
	    print OUT "$ex_start\t$ex_end\n";
	  }
	  else {
	    print OUT "$ex_end\t$ex_start\n";
	  }
	}
      }
#End last exon

#other exons
      else {
	if ($strand == 1) {
	  print OUT "$ex_start\t$ex_end\n";
	}
	else {
	  print OUT "$ex_end\t$ex_start\n";
	}
      }
    }
  }
  print OUT "\t\t\tproduct\t$tr_name\n";
}


sub print_translation_coordinates {
# prints CDS locations and protein tags

  my ($tr,$db,$clone_name,$gene_name) = @_;

  my @exons = @{$tr->get_all_Exons};
  my $nb_ex = scalar(@exons);
  my $strand = $exons[$nb_ex-1]->strand;

  my $count_lation;
  my $translation = $tr->translation;
  my $translation_name = $translation->stable_id;
  my $tr_dbID = $translation->dbID;
#  my $seq = $tr->translate->seq;
# ensembl translate method now strips off any final stop *
# so hack to make this script work is to get translation via simple bio::seq object

  my $able_seq = $tr->translateable_seq;
  my $able = Bio::Seq->new( -display_id => '$translation_name',
                               -seq => $able_seq);
  my $seq = $able->translate()->seq;


# print STDERR "SEQ: $seq\n";
  my ($first) = $seq =~ /(^\S)/;
  my ($last) = $seq =~ /(\S$)/;
# print STDERR "FIRST: $first\tLAST: $last\n";

  my @translateable = @{$tr->get_all_translateable_Exons};
  my $nb1 = scalar(@translateable);

  foreach my $lation (@translateable) {
    $count_lation++;	
    my $tr_start = $lation->start;
    my $tr_end = $lation->end;
	
#only one exon
    if ($nb1 == 1) {
#both methionine and stop codon are included
      if (($first =~ /M/)&&($last =~ /\*/)) {
#takes in account the strand			
	if ($strand == 1) {
	  print OUT "$tr_start\t$tr_end\tCDS\n";
	}
	else {
	  print OUT "$tr_end\t$tr_start\tCDS\n";
	}
      }
#Only methionine included
      elsif (($first =~ /M/) && ($last !~ /\*/)) {
	if ($strand == 1) {
	  print OUT "$tr_start\t>$tr_end\tCDS\n";	
	}
	else {
	  print OUT "$tr_end\t>$tr_start\tCDS\n";
	}
      }
#Only stop codon included
      elsif (($first !~ /M/) && ($last =~ /\*/)) {
	if ($strand == 1) {
	  print OUT "<$tr_start\t$tr_end\tCDS\n";	
	}
	else {
	  print OUT "<$tr_end\t$tr_start\tCDS\n";
	}
      }
#Neither feature is included
      elsif (($first !~ /M/) && ($last !~ /\*/)) {
	if ($strand == 1) {
	  print OUT "<$tr_start\t>$tr_end\tCDS\n";
	}
	else {
	  print OUT "<$tr_end\t>$tr_start\tCDS\n";
	}
      }
      else {
	die;
      }
    }
#End one exon case
	
    elsif ($nb1 > 1) {
#Multiple exons case
#First exon case
      if ($count_lation == 1) {
#Has a methionine
	if ($first =~ /M/) { 
	  if ($strand == 1) {
	    print OUT "$tr_start\t$tr_end\tCDS\n";
	  }
	  else {
	    print OUT "$tr_end\t$tr_start\tCDS\n";
	  }
	}
#Partial feature
	else {
	  if ($strand == 1) {
	    print OUT "<$tr_start\t$tr_end\tCDS\n";
	  }
	  else {
	    print OUT "<$tr_end\t$tr_start\tCDS\n";
	  }
	}
      }
#end first exon case

#last exon case
      elsif ($count_lation == $nb1) {
#last exon partial feature (no stop)
	if ($last !~ /\*/) {
	  if ($strand == 1) {
	    print OUT "$tr_start\t>$tr_end\n";
	  }
	  else {
	    print OUT "$tr_end\t>$tr_start\n";
	  }
	}
#last exon has a stop
	else {
#takes in account the strand
	  if ($strand == 1) {
	    print OUT "$tr_start\t$tr_end\n";
	  }
	  else {
	    print OUT "$tr_end\t$tr_start\n";
	  }
	}
      }
#end of last exon case

#other exons
      else {
	if ($strand == 1) {
	  print OUT "$tr_start\t$tr_end\n";
	}
	else {
	  print OUT "$tr_end\t$tr_start\n";
	}
      }
#end other exons
    }
  }
#end printing coordinates


#start printing tags

# print protein stable id as product
  print OUT "\t\t\tproduct\t$translation_name\n";

# construct new id from the old_prot_ann hash only
# print the official protein_id string
# and print the old id (if any) as note, to make it visible
# *but* do not get any old id if the scaffolds do not match
# or if on moved-on-scaffold kill list
  my ($old_id2print, $protein_id2print);

  if ( ((defined $old_prot_ann{$translation_name}->{'scaffold'})&&($clone_name ne $old_prot_ann{$translation_name}->{'scaffold'}))||(defined $kill{$translation_name}) ) {
    $moved++;
#    print STDERR "$translation_name identified as moved\n";
  }
  else {
# assign old ids (if any) if not moved
    $old_id2print = $old_prot_ann{$translation_name}->{'old_id'} if defined($old_prot_ann{$translation_name}->{'old_id'});
    $protein_id2print = $old_prot_ann{$translation_name}->{'prot_id'} if defined($old_prot_ann{$translation_name}->{'prot_id'});
  }

# now print id line and note depending on whether anything assigned
  if ((defined($old_id2print))&&(defined($protein_id2print))) {
    $count_old_id++;
    print OUT "\t\t\tprotein_id\tgnl|WGS:AAAB|$old_id2print|gb|$protein_id2print\n";
    print OUT "\t\t\tnote\t$old_id2print\n";
  }
  elsif (defined($protein_id2print)) {
    $count_prot_id++;
    print OUT "\t\t\tprotein_id\tgnl|WGS:AAAB|$translation_name|gb|$protein_id2print\n";
  }
  else {
    $count_new++;
    print OUT "\t\t\tprotein_id\tgnl|WGS:AAAB|$translation_name\n";
  }

# if find single protein name ...
# print as note (but not if case-insensitive-same as gene name), and
# print description, if any, as a 2nd product tag (not as prot_desc which is invisible to TrEMBL folk)
  my @protein_names;
  foreach my $db_entry (@{$translation->get_all_DBEntries()}) {
    if ($db_entry->dbname eq 'Anopheles_symbol') {
      push (@protein_names, $db_entry->display_id);
    }
  }
  if (scalar @protein_names == 1) {
    if (lc $protein_names[0] ne lc $gene_name) {
      print OUT "\t\t\tnote\t$protein_names[0]\n";
      $proteins_named++;
    }
    if (defined  $gene_name{lc $protein_names[0]}->{'description'}) {
      my $description = $gene_name{lc $protein_names[0]}->{'description'};
      print OUT "\t\t\tproduct\t$description\n";
      $proteins_described++;
    }
  }

# get from sub and print all ensembl-mapped interpro ids
    my @interpro = &get_protein_annotation($tr_dbID,$db);
    foreach my $ipr(@interpro) {
	print OUT "\t\t\tdb_xref\tInterPro:$ipr\n";
    }

# important note re evidence!
    print OUT "\t\t\tevidence\tnot_experimental\n";
}


sub get_protein_annotation {
  my ($tr_dbID,$db) = @_;
  my @interpro;

  my $query = "select distinct(i.interpro_ac) from interpro i, protein_feature pf where pf.hit_id = i.id and translation_id = $tr_dbID";

  my $sth = $db->prepare($query);
  $sth->execute;
  while (my $ipr = $sth->fetchrow) {
    push (@interpro,$ipr);
  }
  return @interpro;
}
