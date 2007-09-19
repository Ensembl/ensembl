#!/usr/local/ensembl/bin/perl

use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Bio::EnsEMBL::Exon;

my $host      = 'ens-livemirror';
my $dbuser    = 'ensro';
my $dbname    = 'anopheles_gambiae_core_45_3h';
my $dbport    = 3306;
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
					    '-port'   => $dbport,
					   );

my $gene_adapt = $db->get_GeneAdaptor();
my $slice_adapt = $db->get_SliceAdaptor();
my $axfa =  $db->get_AssemblyExceptionFeatureAdaptor();

#get handle for an appropriate mart database
$dbname = 'ensembl_mart_45';
print STDERR "Connecting to $host, $dbname\n";
my $martdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $host,
					    '-user'   => $dbuser,
					    '-dbname' => $dbname,
					    '-pass'   => $dbpass,
					    '-port'   => $dbport,
					   );


#directory with input files and for output
my $inputdir = '/lustre/work1/ensembl/mh4/Genbank_dump_45/Input';
my $outdir = '/lustre/work1/ensembl/mh4/Genbank_dump_45/Output/Outfiles_2007-08-29';


###############################################
##  load some data from files
###############################################

my %scafmap;
my %old_prot_ann;
# my %moved; # not using in 45 dump
my %gene_name;
my %autoname;
my %mir_evidence;
my %old_ensangg;
my %idmap_gene;
my %idmap_prot;
my %idmap_mrna;
my %aedes;
my %dros;

# get the old (Celera) scaffold names
open (SCAFMAP,"$inputdir/accessions") || die "can't open SCAFMAP file";
while(<SCAFMAP>) {
    chomp;
    my ($new,$a,$old) = split;
    $scafmap{$new} = $old;
}
close SCAFMAP;

# get the existing protein tags - use a file derived from NCBI-generated list AAAB.modified.2006.p2g
# AAAB01000034    ENSANGP00000027317      EAL42448        old_id
# where old_id may be celera/ebi, ENSANGPnnn, ENSANGPnnn_37
# different to before, where only celera/ebi ids were stored here
# multiscaff cases have been excluded
open (OLD_PROT_ANN, "$inputdir/protein_annotation") || die "can't find old protein annotation file";
while (<OLD_PROT_ANN>) {
 chomp;
 my ($scaffold, $psi, $prot_id, $old_id) = split;
 $old_prot_ann{$psi}->{'scaffold'} = $scaffold ;
 $old_prot_ann{$psi}->{'prot_id'} = $prot_id ;
 $old_prot_ann{$psi}->{'old_id'} = $old_id ;
}
close OLD_PROT_ANN;

#################################################################################
#  Deleted reading file re proteins that have 'moved' within or between scaffolds
#  Not using in in v45 dump - all mapped ids will have location overlap
#################################################################################

# read file derived from gene name database and xref table
# links sequence names as used in core xrefs (= symbol unless have alt trans)
# to official gene name for use in gene tag
# and to edited description suitable for use in product tag of CDS
# making all the keys lc facilitates case insensitive matching later
# optional 4th column with 'correct' v45 gene stable id
# manually added for cases where a gene symbol is used for >1 gene in v45
# ACE1    ACE1    acetylcholinesterase1
# AgGr32a GPRGR32 gustatory receptor
# CPFL6   CPFL6   cuticular protein 6 from CPFL family   AGAP010907

open (GENE_NAME, "$inputdir/symbol_desc_lookup") || die "can't find gene name file";
while (<GENE_NAME>) {
 chomp;
 my ($seq_ref, $symbol, $description, $designated) = split /\t/;
 $gene_name{lc $seq_ref}->{'real_seq_ref'} = $seq_ref;
 $gene_name{lc $seq_ref}->{'symbol'} = $symbol;
 $gene_name{lc $seq_ref}->{'description'} = $description if (($description ne '')&&($description ne 'NULL'));
 $gene_name{lc $seq_ref}->{'designated'} = $designated if ($designated);

# foreach my $key (keys  %{$gene_name{lc $seq_ref}}) {
#   print sTDERR "$key\t$gene_name{lc $seq_ref}{$key}\n"
# }
# print "\n";

}
close GENE_NAME;


#*# changes pending for v45
#*# want to have Interpro/Panther/auotannotation descr look-up for automated second-choice product naming
#*# add reading file with asociations between AGAMxxxxx-PA and descriptions/names
#*# and storing in hash %autoname


# read file that associates mir (miRNA) names / descriptions with the source of evidence that miRBase (aka Sam G-J) use to identify them
# these were all Dmel miRNas, and I have their miRBase accessions in the file
# this file was created semi-manually using the data at miRBase
# the aga names are derived from the dros ones, but need to check some entries to see exactly which one is which
#miRNA aga-mir-210       MI0000376
#miRNA aga-mir-7         MI0000127

open (MIR, "$inputdir/mirna_list_evidence") || die "can't find miRNA evidence file";
while (<MIR>) {
 chomp;
 my ($name, $evidence) = split /\t/;
 $mir_evidence{$name} = $evidence;
}
close MIR;


# read a list of ENSANGGxxx identifiers in the previous release
# to be printed as old_locus_tag with AgaP_ if have id mapping
# also record (column2 = 1)  if the previous entry had an old locus tag itself
# if so, will also print that as old_locus_tag as plain ENSANGG

open (OLD_ENSANGG,  "$inputdir/ENSANGG_all_in_GB_locustag_status") || die "can't find file with status of old ids";
while (<OLD_ENSANGG>) {
  my ($gid, $status) = split;
  $old_ensangg{$gid} = $status;
}
close OLD_ENSANGG;
#print STDERR "Number of old gene ids parsed: ".scalar keys(%old_ensangg)."\n";


# read lists of mappings from old identifiers to new identifiers
# only needed for dump of 45 because of change from ENSANGG/T/P to AGAP style
# decided to make list where each old and each new id appear only once
# and where minimum score is 0.2
# new_id\told_id\tscore
# AGAP005596-PA   ENSANGP00000031652      0.469468

open IDMAPG, "$inputdir/best_old_for_new_final_gene" || die "can't find file with gene id mapping";
while (<IDMAPG>) {
  my ($newid, $old_id) = split;
  $idmap_gene{$newid} = $old_id;
}
close IDMAPG;

open IDMAPP, "$inputdir/best_old_for_new_final_translation" || die "can't find file with protein id mapping";
while (<IDMAPP>) {
  my ($newid, $old_id) = split;
  $idmap_prot{$newid} = $old_id;
}
close IDMAPP;

open IDMAPM, "$inputdir/best_old_for_new_final_transcript" || die "can't find file with transcript id mapping";
while (<IDMAPM>) {
  my ($newid, $old_id) = split;
  $idmap_mrna{$newid} = $old_id;
}
close IDMAPM;

# read files to lookup protein_id from aedes VB stable id
# these are for use in inference (evidence) tags
# AAEL006014-PA   EAT42428.1
# and uniprot acc.sv from drosphila CGxxxx-Px
# CG8340-PA       P32234.2
# in both cases, neither column is unique
# but decided ok just to take any (ie the last one)!

open AEDES, "$inputdir/aedes_VB_2_proteinid_lookup" || die "can't find file with aedes lookup";
while (<AEDES>)  {
  my ($stable_id, $protein) = split;
  $aedes{$stable_id} = $protein;
}
close AEDES;

open DROS, "$inputdir/drosophila_ensembl_2_uniprot_lookup" || die "can't find file with dros lookup";
while (<DROS>)  {
  my ($stable_id, $protein) = split;
  $dros{$stable_id} = $protein;
}
close DROS;


##############################
##  other preliminary stuff
##############################

# set official locus_tag prefix for this project
my $locus_tag_prefix = "AgaP";

# hash needed later in gene stuff
my %utr;

######################## comment out while testing ################
# get total gene count, and counts by biotype, so can check later to see if any were not retrieved
my %gene_counts;
my $gene_total;
foreach my $gene_id (@{$gene_adapt->list_stable_ids}) {
  $gene_total ++;
  my $gene = $gene_adapt->fetch_by_stable_id($gene_id);
  my $biotype = $gene->biotype;
  $gene_counts{$biotype} ++;
}
print "Gene total: $gene_total\n";
foreach my $type (keys %gene_counts) {
  print "Biotype $type total: ".$gene_counts{$type}."\n";
}
######################## comment out while testing ################



#############################
##  scaffold info section
#############################
# find whether chromosomal (and which) v bacterial contaminant v haplotype scaffold
# for chr ones, get details of any part(s) that are considered a haplotype segment
# get details of subregions represented by alternative haplotype scaff of segment
# include a tag in the file name to help ncbi distinguish bacterial & hap scaff ones for manual intervention


#Get all of the scaffolds
my @scaffolds = @{$slice_adapt->fetch_all('scaffold',undef,1)};
print STDERR "Fetched all ".scalar @scaffolds." scaffolds\n";

# global counters for gene stuff
my $genes_off_scaff=0;
my $genes_ncRNA=0;
my $genes_found=0;
my $genes_named=0;
my $undesignated=0;
my $orths = 0;
my $gene_orth=0;

# global counters etc for mRNA stuff
my ($transcripts_described,$transcripts_autonamed,$transcripts_plain);
$transcripts_described=$transcripts_autonamed=$transcripts_plain=0;

# global counters etc for CDS stuff
my %pid_count;
my ($fetch_count,$proteins_named,$proteins_described,$proteins_autonamed,$proteins_plain);
$fetch_count=$proteins_named=$proteins_described=$proteins_autonamed=$proteins_plain=0;
my ($man_c, $snap_c, $comm_c, $inf_c, $default_c, $error_c);
$man_c=$snap_c= $comm_c= $inf_c=$default_c=$error_c= 0;


# global counters etc for scaff stuff
my %scaffprop;
my $scaffs_seen = 0;
my ($bacc,$hapc,$realc,$fakec);
$bacc=$hapc=$realc=$fakec=0;

# real chr name will only be set where scaff is mapped to a chr - need it because some scaff have 2 chr mappings and want to know the non-hap one for printing.
my $real_chr_name;

# run through scaffolds in name order
foreach my $scaff (sort {$a->seq_region_name cmp $b->seq_region_name} @scaffolds) {
  last if $scaffs_seen >9000;
  my $scaff_name = $scaff -> seq_region_name;
#  next unless ($scaff_name eq "AAAB01008905");  # limit to one scaff for testing
  print STDERR "Processing scaffold $scaff_name\n" if $scaff_name =~/00$/;
  $scaffs_seen ++;

# get basic properties needed for all scaffolds
  my $scaff_length = $scaff -> length;
  my $scaff_old_name = $scafmap{$scaff_name};
# make scaff_name for printing - this can get modified later to tag the different scaffold categories
  my $scaff_name_printable = $scaff_name;
  my $scaff_seq = $scaff -> seq;

# Distinguish the different kinds of scaffold

  if (scalar @{$scaff->get_all_Attributes('bacterial')} == 1) {
    print STDERR "Bacterial attrib for $scaff_name\n";
    $scaffprop{$scaff_name}{'bacterial'} = 1;
    $scaff_name_printable = $scaff_name.".bac";
  }
  else {
#    print STDERR "No bacterial attrib for $scaff_name\n";
    my @segments = @{$scaff -> project('chromosome')};

# If there is only one projection segment
# then scaff is *either a haplotype scaffold, *or entirely on a chromosome
# is there a gotcha in that I have aligned a few haplotype scaffolds along only part of their length?? Probably not as the other bit will not have a chr projection
    if (scalar @segments == 1) {
      my $chr_slice = $segments[0] -> to_Slice;
      my $chr_name = $chr_slice -> seq_region_name;
      if ($chr_name =~ /hap/) {
	print STDERR "Hap scaff for $scaff_name\n";
	$scaffprop{$scaff_name}{'hapscaff'} = 1;
# tag the scaffold name with alt for printing purposes
	$scaff_name_printable = $scaff_name.".alt";
# need to know what the name of the scaffold that makes up the real chr at this point
# first fetch the alternative chr slice then project it down to its scaffold
	my @alternative = @{$axfa -> fetch_all_by_Slice($segments[0] -> to_Slice)};
	print STDERR "ERROR haplotype scaff $scaff_name goes to more than 1 place\n" if (scalar @alternative != 1);
	my $altslice = $alternative[0] -> alternate_slice();
# project it back down to a scaffold
	my @backsegments = @{$altslice -> project('scaffold')};
	print STDERR "ERROR $scaff_name alt slice projects to more than 1 scaffold\n" if (scalar @backsegments != 1);
	my $alt_scaff_slice = $backsegments[0] -> to_Slice;
	my $alt_scaff_name = $alt_scaff_slice -> seq_region_name;
	$scaffprop{$scaff_name}{'hapscaff_altscaff'} = $alt_scaff_name;	
      }
      elsif ($chr_name =~ /2R|2L|3L|3R|X/) {
	$scaffprop{$scaff_name}{'real_chr'} = $chr_name;
	$scaffprop{$scaff_name}{'real_chr_slice'} = $chr_slice;
	$real_chr_name = $chr_name;
	print STDERR "Real chr scaff for $scaff_name\n";
      }

# if want to print info re being a Y chr scaff, have this
# if not add Y_unplaced to UNKN elsif and include in fake with no name print in fsa header
      elsif ($chr_name =~ /Y_unplaced/) {
	$real_chr_name = "Y";
	$scaffprop{$scaff_name}{'real_chr'} = $chr_name;
	$scaffprop{$scaff_name}{'real_chr_slice'} = $chr_slice;
	print STDERR "Y chr scaff for $scaff_name\n";
#	print STDERR $scaffprop{$scaff_name}{'real_chr'}."\n";
      }

      elsif ($chr_name =~ /UNKN/) {
	$scaffprop{$scaff_name}{'fake_chr'} = $chr_name;
	print STDERR "Unkn chr scaff for $scaff_name\n";
      }
      else {
	print STDERR "ERROR - not matched chr name <<$chr_name>> in single seg projection\n";
      }
    }

# if there are 2 or more projection segments it should be part a single real chr and 1 or more part(s) haplotype segment
    else {
      print STDERR "Complex for $scaff_name - ".scalar @segments." segments\n";
      my $count_real=0;
      foreach my $seg (@segments) {
	my $chr_slice = $seg -> to_Slice;
	my $chr_name = $chr_slice -> seq_region_name;
# if hap seg, will need to know the name of the scaffold that makes up the real chr at this point
# record details and leave till later as only going in feature table
# should be only 1 real chromosome
	if ($chr_name =~ /hap/) {
	  print STDERR "$chr_name  \n";
	  $scaffprop{$scaff_name}{'hapseg'}{$chr_name}{'scaff_start'} = $seg->from_start();
	  $scaffprop{$scaff_name}{'hapseg'}{$chr_name}{'scaff_end'} = $seg->from_end();
	  $scaffprop{$scaff_name}{'hapseg'}{$chr_name}{'hap_slice'} = $chr_slice;
	}
	elsif ($chr_name =~ /2R|2L|3L|3R|X/) {
	  $count_real ++;
	  $real_chr_name=$chr_name;
	  print STDERR "$chr_name  \n";
	  $scaffprop{$scaff_name}{'real_chr'}{$chr_name}{'scaff_start'} = $seg->from_start();
	  $scaffprop{$scaff_name}{'real_chr'}{$chr_name}{'scaff_end'} = $seg->from_end();
	  $scaffprop{$scaff_name}{'real_chr'}{$chr_name}{'real_slice'} = $chr_slice;
	}
	else {
	  print STDERR "ERROR - not matched chr name <<$chr_name>> in multiseg ".scalar @segments." projection\n";
	}
      }
      print STDERR "ERROR - $count_real instead of a single real chr mapping for $scaff_name\n" if ($count_real !=1);
    }
  }

#########################################################################
# by this point should have the info needed to print out scaffold def lines
# do all the .fsa printing here to keep things tidy
# also include in this section printing of misc_features in .tbl file that relate to alternative assembly stuff
# for a bac scaff - nothing more to print - go to next scaff
# for a hap scaff - nothing more to print - *not* specifying its location in the chr scaff that it matches, as GenBank say they don't like references to coordinates in other entries - go to next scaff
# for a scaff on unkn - no scaff info to print - proceed to gene info
#########################################################################

  open (OUT,">$outdir/$scaff_name_printable.tbl") || die "can't open tbl file";
  open (SEQ,">$outdir/$scaff_name_printable.fsa") || die "can't open fsa file";

# if it is a bacterial scaff
  if ($scaffprop{$scaff_name}{'bacterial'}) {
    $bacc++;
    print SEQ ">gnl|WGS:AAAB|$scaff_old_name|gb|$scaff_name [organism=Anopheles gambiae str. PEST] [tech=wgs] [note=probable bacterial contamination]\n";
    print SEQ "$scaff_seq\n";
    print OUT ">Feature gnl|WGS:AAAB|$scaff_old_name|gb|$scaff_name\n";
    next;  # no other annotation - go to next scaff
  }

# if it is a simple haplotype scaff
  elsif ($scaffprop{$scaff_name}{'hapscaff'}) {
    $hapc++;
    print SEQ ">gnl|WGS:AAAB|$scaff_old_name|gb|$scaff_name [organism=Anopheles gambiae str. PEST] [tech=wgs] [note=probable alternative assembly of part of scaffold ".$scaffprop{$scaff_name}{'hapscaff_altscaff'}."]\n";
    print SEQ "$scaff_seq\n";
    print OUT ">Feature gnl|WGS:AAAB|$scaff_old_name|gb|$scaff_name\n";
    next;  # no other annotation - go to next scaff
  }

# if it is a simple assembly scaff on a fake chr - as currently set that means UNKN only
  elsif ($scaffprop{$scaff_name}{'fake_chr'}) {
    $fakec++;
    print SEQ ">gnl|WGS:AAAB|$scaff_old_name|gb|$scaff_name [organism=Anopheles gambiae str. PEST] [tech=wgs] [note=component of genomic assembly AgamP3]\n";
    print SEQ "$scaff_seq\n";
    print OUT ">Feature gnl|WGS:AAAB|$scaff_old_name|gb|$scaff_name\n";
  }

# if it is a assembly scaff on a real chr - or, as currently set on Y_unplaced
# note that the defline is going to look the same whether or not it has part(s) that are haplotype segments
# details of those will appear only in the .tbl file
  elsif ($scaffprop{$scaff_name}{'real_chr'}) {
    $realc++;
    print SEQ ">gnl|WGS:AAAB|$scaff_old_name|gb|$scaff_name [organism=Anopheles gambiae str. PEST] [tech=wgs] [chromosome=$real_chr_name] [note=component of assembly AgamP3]\n";
    print SEQ "$scaff_seq\n";
    print OUT ">Feature gnl|WGS:AAAB|$scaff_old_name|gb|$scaff_name\n";

#####################################################################################
# printing of info about hapsegs and any target regions for alt assemblies
# For Y - no scaff info to print - proceed to gene info
# (but unless it causes problems take it through the other'real' stuff for simplicity)
# For a scaff on 'real' chr
# 1: print already part-retrieved info re any overlap hap seg it contains
# 2: retrieve info re any segments that are the chromosomal target of hap scaffs or segs & print that

# declare chr slice because will be making within 2 separate 'if' blocks depending on whether it is simple or has hap overlap segs
    my $chr_slice_for_alts;

# 1 - does it have segment(s) that are overlap hap segs?
    if ($scaffprop{$scaff_name}{'hapseg'}) {
      foreach my $hapsegchr (keys %{$scaffprop{$scaff_name}{'hapseg'}}) {
# get the the 'real' scaffold corresponding to the stored hap chr slice
# same code as used above for hap scaffs
	my $seg_start = $scaffprop{$scaff_name}{'hapseg'}{$hapsegchr}{'scaff_start'};
	my $seg_end = $scaffprop{$scaff_name}{'hapseg'}{$hapsegchr}{'scaff_end'};
	my @alternative = @{$axfa -> fetch_all_by_Slice($scaffprop{$scaff_name}{'hapseg'}{$hapsegchr}{'hap_slice'})};
	print STDERR "ERROR hap seg of $scaff_name $hapsegchr goes to more than 1 place\n" if (scalar @alternative != 1);
	my $altslice = $alternative[0] -> alternate_slice();
	my @backsegments = @{$altslice -> project('scaffold')};
	print STDERR "ERROR $scaff_name hap seg projection goes to more than 1 scaffold\n" if (scalar @backsegments != 1);
	my $alt_scaff_slice = $backsegments[0] -> to_Slice;
	my $alt_scaff_name = $alt_scaff_slice -> seq_region_name;
	print OUT "$seg_start\t$seg_end\tmisc_feature\n";
	print OUT "\t\t\tnote\tprobable alternative assembly of part of scaffold $alt_scaff_name\n";
      }

# now get the real chromosome projection in the complicated way for use in 2. below
      $chr_slice_for_alts = $scaffprop{$scaff_name}{'real_chr'}{$real_chr_name}{'real_slice'}
    }
# if no hap segs then the single projection slice is stored differently
    else {
      $chr_slice_for_alts = $scaffprop{$scaff_name}{'real_chr_slice'};
    }

# 2 - does it have any regions that are targets of alternative assemblies?
# use its real chromosomal projection slice (defined above) and fetch all assembly exception features
# slightly surpringly, get an empty array rather than undef if there are no ass exception features - this is useful here as it saves me another 'if' block!
# they seem to come back in reverse order of position on the target scaffolds, so add a sort

    my @alternatives = @{$axfa -> fetch_all_by_Slice($chr_slice_for_alts)};
    print STDERR "$scaff_name has ".scalar @alternatives." alt assemblies on it\n";
    foreach my $alt (sort {$a->start <=> $b->start} @alternatives) {

# first find out about the real chr region of this part of the target scaffold - need to project back to scaffold
      my $target_chr_slice = ($alt -> slice) -> sub_Slice($alt -> start, $alt -> end);
      my @target_segments = @{$target_chr_slice -> project('scaffold')};
      print STDERR "ERROR target seg of $scaff_name is not a single projection\n" if (scalar @target_segments != 1);
      my $target_scaff_slice = $target_segments[0] -> to_Slice;
      print STDERR "ERROR $scaff_name target seg did not come back as $scaff_name\n" if ($scaff_name ne $target_scaff_slice -> seq_region_name);
      my $target_start =  $target_scaff_slice -> start;
      my $target_end =  $target_scaff_slice -> end;
      my $altslice = $alt -> alternate_slice();

# now find the name of the alternative scaffold that aligns here
# and whether it is a simple hap scaff, or a segment of something that also has real bits
      my @backsegments = @{$altslice -> project('scaffold')};
      print STDERR "ERROR $scaff_name alt assembly projection goes to more than 1 scaffold\n" if (scalar @backsegments != 1);
      my $alt_scaff_slice = $backsegments[0] -> to_Slice;
      my $alt_scaff_name = $alt_scaff_slice -> seq_region_name;
# take the whole scaffold and see how many chr projections it has: if 1 it is a simple hap scaff
      my $alt_scaff_whole = $slice_adapt->fetch_by_region('scaffold',$alt_scaff_name);
      my $extra_text = "";
      if (scalar @{$alt_scaff_whole -> project('chromosome')} >1) {
	print STDERR "$scaff_name is target for $alt_scaff_name which projects to real chr as well\n";
	$extra_text = "part of ";
      }
# and finally print out the details
      print OUT "$target_start\t$target_end\tmisc_feature\n";
      print OUT "\t\t\tnote\ta probable alternative assembly of this region is represented by ".$extra_text."scaffold $alt_scaff_name\n";
    }

# end of printing info about hapsegs and any target regions for alt assemblies
####################################################################################

  }
  # should have printed all the scaffolds at this point
  else {
    print STDERR "ERROR: unclassified scaffold $scaff_name will not be printed\n"
  }


###########################################################################
##  Gene annotation starts here
###########################################################################

# first check if there are any genes that hang off the end of the slice (if I just call fetch all genes on each scaffold these will be ignored)
# refer to separate script find_multiscaff_genes to get details of the 'real' multiscaff genes
# can then enter those details manually in the tbl files

  my @segments = @{$scaff -> project('chromosome')};
  foreach my $seg (@segments) {
    my $chr_slice = $seg -> to_Slice;
    my $chr_name = $chr_slice -> seq_region_name;
    my $chr_slice_length = $chr_slice -> length;
    my $chr_slice_strand = $chr_slice -> strand;
    if ($chr_name =~ /2R|2L|3L|3R|X|UNKN/) {
# keep track of scaffold coordinates for this segment
      my $proj_start = $seg -> from_start;
      my $proj_end = $seg -> from_end;
      my @chr_genes = @{$chr_slice -> get_all_Genes};
      foreach my $chr_gene(@chr_genes) {
	my $chr_gene_start = $chr_gene->start;
	my $chr_gene_end = $chr_gene->end;
	if (($chr_gene_start < 1) || ($chr_gene_end > $chr_slice_length)) {
	  print STDERR "TOO BIG: ".$chr_gene->stable_id." extends beyond scaffold $scaff_name:$proj_start:$proj_end on chromosome $chr_name and will not be fetched.  Gene $chr_gene_start:$chr_gene_end\n";
	  $genes_off_scaff++;
	}
      }
    }
  }

# now get all genes contained on each scaffold

  my @genes = @{$scaff->get_all_Genes};

# come back in F order-on-scaff for scaffs that are +1 ori on csome and R order-on-scaff for scaffs that are -1 ori
# rather than jsut reverse the latter, don't rely on this behaviour and do a proper sort
  @genes = sort {$a->start <=> $b->start} @genes;

  foreach my $gene(@genes) {
# these are the 'original' start and end
    my $start = $gene->start;
    my $end = $gene->end;
#    print STDERR "BEFORE ".$gene->stable_id."\t".$start."\t".$end."\t".$gene->strand."\n";
    $genes_found++;


###############################
#  treat ncRNA genes separately
###############################

# in release 37/38 and 45 names are not stored as xrefs but in gene description field
# this will likely change in the future
# also likely to be other ncRNA biotypes to deal with
    my $gene_biotype = $gene->biotype;
    if ($gene_biotype eq 'miRNA' || $gene_biotype eq 'tRNA') {
      $genes_ncRNA ++;
      my $gene_id = $gene->stable_id;
      $gene_biotype ='misc_RNA' if ($gene_biotype eq 'miRNA');
# first print a gene feature with no gene name, but location, locus tag, xref and evidence
      if ($gene->strand == 1) {
	print OUT "$start\t$end\tgene\n";
      }
      else {
	print OUT "$end\t$start\tgene\n";
      }

# locus tag , and old locus tag(s) if any, and note re old gene stable id
      print OUT "\t\t\tlocus_tag\t$locus_tag_prefix"."_$gene_id\n";
      my $old_gene_id;
      if ($idmap_gene{$gene_id}) {
	$old_gene_id = $idmap_gene{$gene_id};
	print OUT "\t\t\told_locus_tag\t$locus_tag_prefix"."_$old_gene_id\n";
	if ($old_ensangg{$old_gene_id} == 1) {
	  print OUT "\t\t\told_locus_tag\t$old_gene_id\n";
	}
	else {
	  print OUT "\t\t\tnote\t$old_gene_id\n";
	}
      }

# db_xref to VB
      print OUT "\t\t\tdb_xref\tVectorBase:$gene_id\n";
# then add an evidence line, depending on kind of gene (and pseudo tag if needed)
      if ($gene_biotype eq 'tRNA') {
	print OUT "\t\t\tinference\tprofile:tRNAscan-SE:1.23\n";
	print OUT "\t\t\tpseudo\n" if ($gene->description =~ /Pseudo/i);
      }
# note that the evidence tag below will only be valid when the misc_RNA is miRNA from miRBase
      if ($gene_biotype eq 'misc_RNA') {
	if (my $mirbase_evidence = $mir_evidence{$gene->description}) {
	  print OUT "\t\t\tinference\tsimilar to RNA sequence, other RNA:miRBase:$mirbase_evidence\n";
	}
	else {
	  print STDERR "ERROR: miRNA ".$gene->description." has no evidence\n";
	}
      }

#second, print a tRNA or misc_RNA feature, with product only
      if ($gene->strand == 1) {
	print OUT "$start\t$end\t$gene_biotype\n";
      }
      else {
	print OUT "$end\t$start\t$gene_biotype\n";
      }
# tRNA genes with no proper product need special handling
      if (($gene->description =~ /tRNA-Pseudo/i) || ($gene->description =~ /tRNA-Undet/i)) {
	print OUT "\t\t\tproduct\t"."tRNA-Xxx"."\n";
      }
      else {
	print OUT "\t\t\tproduct\t".$gene->description."\n";
      }
      next;  # next gene
    }
    elsif ($gene_biotype ne 'protein_coding') {
      print STDERR "ERROR: Gene ".$gene->stable_id." biotype is not protein_coding, miRNA or tRNA but $gene_biotype - will not be handled properly\n";
    }

##############################################
#  protein-coding genes
#############################################

# gene coord details
####################

    my $new_gene = &checks($gene);

# &checks suboutine strips UTR from any gene with UTR adjacent to non-Met start or non-stop end
# also sets %utr if are remaining (valid) UTRs
# note that $gene and $new_gene are the same thing after this point as return $gene and set $new_gene to it.   doh!
# i.e. both are the result of running &checks subroutine

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
# end gene coord details
########################################

# gene tags
###########
#*#  changed for v45:  get gene display xref directly (Ano symbol)- don't use DBlinks
#*#  changed for v45:  the stuff dealing with multiple fetches is then not needed
#*#  changed for v45:  but use symbol_desc_lookup file with manual addition of 'designated' column to avoid using same symbol twice

# gene symbol
# note use of lc when comparing to the previously-made-lc keys of %gene_name
    my $gene_name;
    if ($gene -> display_xref) {
      my $dbentry = $gene -> display_xref;
      if ($dbentry -> dbname eq "Anopheles_symbol") {
	$gene_name = $dbentry -> display_id;
	if ( ($gene_name{lc $gene_name}{'designated'}) && ($gene_name{lc $gene_name}{'designated'} ne $gene_id) ) {
	  $undesignated++;
# undef $gene_name so can carry info over to CDS sub-routine
	  $gene_name = undef;
	}
	else {
	  print OUT "\t\t\tgene\t$gene_name\n";
	  $genes_named++;
	}
      }
    }

# locus tag(s)
# always print ensembl stable id as locus_tag - now need a locus_tag prefix
    my $locus_tag_gene = $locus_tag_prefix.'_'.$gene_id;
    print OUT "\t\t\tlocus_tag\t$locus_tag_gene\n";

# print old_locus_tag if the gene id was present in the previous release
#*# changed for v45 - added a second old locus tag if previous entry had an old locus tag.  Format of old_ensangg will be changed
#*# changed for v45 - added another look-up to see which if any old-style ensembl id is mapped here!!
    my $old_gene_id;
    if ($idmap_gene{$gene_id}) {
      $old_gene_id = $idmap_gene{$gene_id};
      print OUT "\t\t\told_locus_tag\t$locus_tag_prefix"."_$old_gene_id\n";
#      print STDERR "Printing old locus tag $locus_tag_prefix"."_$old_gene_id for $gene_id\n";
      if ($old_ensangg{$old_gene_id} == 1) {
	print OUT "\t\t\told_locus_tag\t$old_gene_id\n";
	#  print STDERR "Printing old old locus tag $old_gene_id for $gene_id\n";
      }
      else {
	print OUT "\t\t\tnote\t$old_gene_id\n" ;
	#  print STDERR "NOT printing old old locus tag $old_gene_id for $gene_id but printing note\n";
	
      }
    }
    else {
#      print STDERR "NOT printing old locus tag or ENSANGG note for $gene_id\n";
    }

# Dros orthologue(s)
# query Mart to get the display id for any drosophila orthologues - easier than using compara
# note that most dros orth in this table have both display id and stable id but a few only have stable id!
    my (%homologues,@homologue_list,$dros_id,$dros_symbol);
    my $query2 = "select homol_stable_id,display_id from agambiae_gene_ensembl__homologs_dmelanogaster__dm  where gene_stable_id like ?";
    my $sth2 = $martdb->dbc->prepare($query2);
    $sth2->execute($gene_id);
    while (($dros_id, $dros_symbol) = $sth2->fetchrow_array) {
      if ($dros_symbol) {
	$orths ++;
	$homologues{$dros_symbol}+=1;
      }
      else {
	$orths ++;
	$homologues{$dros_id}+=1;
      }
    }
    @homologue_list = (keys %homologues);
    if (scalar @homologue_list >0) {
      $gene_orth++;
      print OUT "\t\t\tnote\tsimilar to Drosophila ";
      print OUT join(', ', sort @homologue_list),"\n";
    }


# db_xref to VB
    print OUT "\t\t\tdb_xref\tVectorBase:$gene_id\n";


# get the transcript(s) and print mRNA and protein features
###########################################################
# details are in sub routines

    my @new_transcripts = @{$new_gene->get_all_Transcripts};	
    foreach my $new_tr(@new_transcripts) {
      $fetch_count++;	
      &print_transcript_coordinates($new_tr, $gene_name,$locus_tag_gene, $gene_id);
      &print_translation_coordinates($new_tr,$db,$scaff_name,$gene_name,$locus_tag_gene, $gene_id);

# temporary progress report
      if ($fetch_count % 100 == 0) {
	print STDERR "Up to $scaff_name - Fetched $fetch_count proteins so far\n";
      }
# end of trancript loop
    }

# end of gene loop
  }

# end of this scaffold
  close(OUT);
  close(SEQ);
}

#######################
## print overall stats
#######################
print STDERR "Finished\n";

print "Scaffolds seen:  $scaffs_seen\n";
print "\tBac scaff: $bacc\n\tHap scaff: $hapc\n\tReal chr: $realc\n\tFake chr: $fakec\n";

print  "Genes that run off end of a scaffold: $genes_off_scaff\n";
print  "\tIf there are genes running off scaffold, multiscaff genes may exist and be lost.\n\tNB Genes will be counted twice if they are on 2 scaffolds.\n\tUse script find_multiscaff_genes to generate coordinates for manual addition to tbl file.\n";

print  "Genes found on scaffolds: $genes_found\n\tof which ncRNA genes: $genes_ncRNA\nGenes with names: $genes_named\tUn-named because of multiply used gene name: $undesignated\nGenes with Dros orths: $gene_orth\tTotal dros orthologues: $orths\n";

print  "Fetched proteins: $fetch_count, note with distinct name: $proteins_named\n";

print "Product details - CDS:\n";
print "\tfrom gene_name $proteins_described\n\tfrom autoname $proteins_autonamed\n\tstable_id only $proteins_plain\n";
print "Product details - mRNA:\n";
print "\tfrom gene_name $transcripts_described\n\tfrom autoname $transcripts_autonamed\n\tstable_id only $transcripts_plain\n";

print  "Protein id asssignment details:\n";
foreach my $pid_type (sort keys %pid_count) {
  print  "\t$pid_type: ".$pid_count{$pid_type}."\n";
}

print  "Inference details:\n";
print "\tManual note: $man_c\n\tDefault note: $default_c\n\tSnap: $snap_c\n\tInference tag(s): $inf_c\n\tComm or genename only no note: $comm_c\n\tError (other):  $error_c\n";




#######################
## checks subroutine ##
#######################
# takes gene and looks at its transcript(s)
# if the coding sequence has UTR next to non-Met start or non-stop end does something complicated!!
# maybe ... makes new set of exons and sets translation appropriately so that a revised transcript is made without those dodgy UTRs
# then removes transcripts from gene and replaces with new ones

sub checks {

  my ($gene) = @_;
  my $gene_dbID = $gene->dbID;  #probably redundant
  my $gene_slice = $gene->slice;  #need access to the slice that is the current scaffold
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
# ?? bizarrely sets $tl_end to the position of the CDS end in transcript coordinates - used later on

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

# make new exons for each pair of coords - different to previous way of doing this
    for my $exon_coord ( @exon_coords ) {
      if ($exon_coord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
# print STDERR $exon_coord->start."\t".$exon_coord->end."\t".$exon_coord->strand."\t".$exon_coord->id."\n";
	my $new_exon = new Bio::EnsEMBL::Exon(-START	=> $exon_coord->start,
					      -END          => $exon_coord->end,
					      -STRAND       => $exon_coord->strand,
					      -SLICE        => $gene_slice,
					      -PHASE        => 0,
					      -END_PHASE    => 0,
					     );
	push (@new_exons,$new_exon);
#  print STDERR $new_exon->start."\t".$new_exon->end."\t".$new_exon->strand."\tnew exon\n";
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


###############################################
##  sub to print transcript info
##  prints mRNA locations and mRNA product tag
###############################################
sub print_transcript_coordinates {
  my ($tr, $gene_name, $locus_tag_gene, $gene_id) = @_;

  my $coding_start = $tr->cdna_coding_start;
  my $coding_end = $tr->cdna_coding_end;
  my $cdna_length = length($tr->seq->seq);
  my $count_ex = 0;
  my $tr_name = $tr->stable_id;
  my $tr_dbID = $tr->dbID;

#*# changed for v45 - also need to get info about CDS to print in mRNA section
  my $translation = $tr->translation;
  my $translation_name = $translation->stable_id;

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


# start printing tags
######################

#*# changed for v45 - adding gene and locus_tag
  if (defined $gene_name) {
    print OUT "\t\t\tgene\t$gene_name\n";
  }
  print OUT "\t\t\tlocus_tag\t$locus_tag_gene\n";

#*# changed for v45 - mRNA product tag name should match CDS product
#*# VB identifier in brackets is however -RA not -PA form
#*#  ugly - uses near-duplication of code in the CDS section

  my $description;
  if ( (defined $gene_name) && ($tr ->display_xref) ) {
    my $dbentry = $tr -> display_xref;
    if ($dbentry -> dbname eq "Anopheles_symbol") {
      my $prot_name = $dbentry -> display_id;

## omitting note re distinct name - used in CDS tag sub-routine only
#      if (lc $prot_name ne lc $gene_name) {
#	print OUT "\t\t\tnote\t$prot_name\n";
#	$proteins_named++;	
#      }

      if ( $gene_name{lc $prot_name}{'description'} ) {
	$description = $gene_name{lc $prot_name}{'description'};
	print OUT "\t\t\tproduct\t$description ($tr_name)\n";
	$transcripts_described++;
      }
    }
  }

# look for autoname description if no gene-name-description-derived product
# may need to add not-taking autoname if 'hypothetical protein' etc
## note deliberately using translation_name as the look-up, not tr(anscript) name
  unless (defined $description) {
    if ($autoname{$translation_name}) {
      print OUT "\t\t\tproduct\t$autoname{$translation_name} ($tr_name)\n";
      $transcripts_autonamed++;
    }
    else {
      print OUT "\t\t\tproduct\t$tr_name\n";
      $transcripts_plain ++;
    }
  }


#*# changed for v45 - print transcript_id (for NCBI internal purposes) with prefix
#*# changed for v45 - add protein_id tag to mRNA
#*# protein_id generation is ugly - simplification/duplication of code in the CDS section
#*# qv for comments and checks and counts

  my $prefix = 'gnl|WGS:AAAB|';

  print OUT "\t\t\ttranscript_id\t$prefix$tr_name\n";

  my $pid;
  my $old_translation_name;

  if ($idmap_prot{$translation_name}) {
    $old_translation_name = $idmap_prot{$translation_name};

    if (defined $old_prot_ann{$old_translation_name}) {
      $pid = $prefix.$old_prot_ann{$old_translation_name}{'old_id'}."|gb|".$old_prot_ann{$old_translation_name}{'prot_id'};
      $pid_count{'old-id-found-mRNA-section'} ++;
    }
    else {
#  id mapped ENSANGP found, but this was not in previous release
#  will generate a warning in CDS section but proceed to assign 'new' protein id
      $pid =  $prefix.$translation_name;
      $pid_count{'old-id-found-but-not-in-gb-mRNA-section'} ++;
    }
  }

  else {
    $pid = $prefix.$translation_name;
    $pid_count{'new-protein-mRNA-section'} ++;
  }

  print OUT "\t\t\tprotein_id\t$pid\n";


#*# changed for v45
#*# - explicitly linked mRNA and its CDS
#*# - added note re old-style ENSANGT if id mapped
#*# - added db_xrefs to VB

  print OUT "\t\t\tnote\t$tr_name coding for $translation_name\n";

  if ($idmap_mrna{$tr_name}) {
    print OUT "\t\t\tnote\t".$idmap_mrna{$tr_name}."\n";
  }

  print OUT "\t\t\tdb_xref\tVectorBase:$gene_id\n";
  print OUT "\t\t\tdb_xref\tVectorBase:$tr_name\n";

}


#########################################
##  sub to print CDS info
##  prints CDS locations and protein tags
#########################################

sub print_translation_coordinates {
  my ($tr,$db,$scaff_name,$gene_name,$locus_tag_gene, $gene_id) = @_;

# First print coodinates of cds
#################################
  my @exons = @{$tr->get_all_Exons};
  my $nb_ex = scalar(@exons);
  my $strand = $exons[$nb_ex-1]->strand;

  my $count_lation;
  my $transcript_name = $tr->stable_id;
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

  my ($first) = $seq =~ /^(\S)/;
  my ($last) = $seq =~ /(\S)$/;
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


# Start printing CDS tags
#########################

#*# changed for v45 - added print of locus_tag tag and gene name (if defined)
  if (defined $gene_name) {
    print OUT "\t\t\tgene\t$gene_name\n";
  }
  print OUT "\t\t\tlocus_tag\t$locus_tag_gene\n";

#*# changed for v45 - added print of (internal use) transcript_id - use same prefix

  my $prefix = 'gnl|WGS:AAAB|';

  print OUT "\t\t\ttranscript_id\t$prefix$transcript_name\n";

# construct ncbi protein id string
# from the old_prot_ann hash if stable_id existed in previous release
# and print the old id (s) (if any) as notes, to make it visible
#*# changed for v45
#*# - no stable id should have moved between or within scaffolds, so no longer need hacking to cope with this
#*# - new step needed to find if have an ENSANGP that matches AGAPxxx-PA
#*# - old_protein_ann look-up now has details of old protein id in all cases
#*# - print old ENSANGP stable id as note
#*# - print additional note if the protein_id includes an old style ebiP or Celera id

  my $pid;
  my $old_translation_name;

# find if have old ENSANGP stable id mapped to this AGAPxxx-PA
  if ($idmap_prot{$translation_name}) {
    $old_translation_name = $idmap_prot{$translation_name};
      print OUT "\t\t\tnote\t$old_translation_name\n";

# find if this ENSANGP was used in previous release
    if (defined $old_prot_ann{$old_translation_name}) {

    # double-check that this ENSANGP has not moved scaffolds
      print STDERR "WARNING - scaffold mismatch for protein $translation_name old name $old_translation_name - was ".$old_prot_ann{$old_translation_name}{'scaffold'}." now $scaff_name\n" if ($old_prot_ann{$old_translation_name}{'scaffold'} ne $scaff_name);

#################################################################################
#  Deleted section re proteins that have 'moved' within or between scaffolds
#  Not using in  in 45_3h dump - all mapped ids will have location overlap
#################################################################################

      if (defined $old_prot_ann{$old_translation_name}{'old_id'}) {
	# construct protein id to match that used last time
	$pid_count{'old-id-found'} ++;
	$pid = $prefix.$old_prot_ann{$old_translation_name}{'old_id'}."|gb|".$old_prot_ann{$old_translation_name}{'prot_id'};
        # also print note with the old protein identifier if Celera or ebiP style
	if ( ( $old_prot_ann{$old_translation_name}{'old_id'} =~ /^ebiP/ ) || ( $old_prot_ann{$old_translation_name}{'old_id'} =~ /^agCP/ ) ) {
	  print OUT "\t\t\tnote\t".$old_prot_ann{$old_translation_name}{'old_id'}."\n";
	}
      }
      else {
	print STDERR "ERROR: failed to retrieve details for old stable id $old_translation_name that mapped to new $translation_name despite presence in previous release\n";
      }
    }

    else {
    # possible error - check these
    # id mapped ENSANGP found for AGAPxxx-PA but that ENSANGP not in previous genebank release??
    # probably means previous release's protein did not make it into Genbank
    # print with 'new' id anyway.
      print STDERR "WARNING:  $old_translation_name match for $translation_name but $old_translation_name not in previous release???? - protein _id printed\n";
      $pid_count{'old-id-found-but-not-in-gb'} ++;
      $pid =  $prefix.$translation_name;
    }
  }

  else {
    # new protein - or at least not matched to any in previous genbank release
    # construct protein id with AGAPxxx-PA  (but no EAA)
    # no note to print
    $pid_count{'new-protein'} ++;
    $pid = $prefix.$translation_name;
  }
  print OUT "\t\t\tprotein_id\t$pid\n";

#*# changed for v45 - explicitly link mRNA and its CDS
  print OUT "\t\t\tnote\t$translation_name encoded by $transcript_name\n";


#*# changed for v45 - protein stable id and description if any (from genename or autoname) will be used within a single product tag, rather than printing separately
#*# - use transcript display xref to retrieve single possibly-different name, not DBentries
#*# - require that there is an defined gene name (this handles the cases where gene name is used multiple times and I have found this gene to be undesignated)
#*# - description look-up from autoname-based hash if no genename one

  my $description;
  if ( (defined $gene_name) && ($tr ->display_xref) ) {
    my $dbentry = $tr -> display_xref;
    if ($dbentry -> dbname eq "Anopheles_symbol") {
      my $prot_name = $dbentry -> display_id;
      if (lc $prot_name ne lc $gene_name) {
	print OUT "\t\t\tnote\t$prot_name\n";
	$proteins_named++;	
      }
      if ( $gene_name{lc $prot_name}{'description'} ) {
	$description = $gene_name{lc $prot_name}{'description'};
	print OUT "\t\t\tproduct\t$description ($translation_name)\n";
	$proteins_described++;
      }
    }
  }

# look for autoname description if no gene-name-description-derived product
# may need to add not-taking autoname if 'hypothetical protein' etc
  unless (defined $description) {
    if ($autoname{$translation_name}) {
      print OUT "\t\t\tproduct\t$autoname{$translation_name} ($translation_name)\n";
      $proteins_autonamed++;
    }
    else {
      print OUT "\t\t\tproduct\t$translation_name\n";
      $proteins_plain ++;
    }
  }




# db_xref to VB
#*# changed for v45 - add gene db_xref
  print OUT "\t\t\tdb_xref\tVectorBase:$gene_id\n";
  print OUT "\t\t\tdb_xref\tVectorBase:$translation_name\n";

# db_xref(s) to InterPro
# get from sub and print all ensembl-mapped interpro ids
  my @interpro = &get_protein_annotation($tr_dbID,$db);
  foreach my $ipr(@interpro) {
    print OUT "\t\t\tdb_xref\tInterPro:$ipr\n";
  }

# source of prediction and inference re evidence
# would like to use inference tags but strict rules
# find if it has any supporting evidence - if not then it is a snap and can use inference tag.
# snap 'version' required - I use date c/o Ian Holm as SNAP has no formal versions.
#*# changed for v45 - include some inference tags
#*# changed for v45 - add note if manual

# first see if there is any evidence
#*# unchanged for v45 - added section to give count of distinct hit name
#*# but laborious and not used so commented out
  my $evidence = &get_evidence($tr,$db);
  my ($manual, $comm, @infered_from);

  if ($evidence == 0) {
    print OUT "\t\t\tinference\tab initio prediction:SNAP:2003-07-30\n";
    $snap_c ++;
  }

# otherwise, proceed to find the tsfs
# use presence of .sv to spot those that can be use in inference tags
  else {
    my @tsfs = @{ $tr->get_all_supporting_features };
    foreach my $tsf (@tsfs) {
      my $logic = $tsf->analysis->logic_name;
      my $hname = $tsf -> hseqname;
      if ($logic eq 'manual_models') {
	$manual = 1;
      }
      elsif ($logic eq 'community_models') {
	$comm = 1;
      }
      elsif ($logic eq 'anopheles_cdna_est') {
	if ($hname =~ /\.\d$/)  {
	  push @infered_from, "similar to RNA sequence:INSD:$hname";
	}
	else {
	  $comm = 1;
	}
      }
      elsif ($logic eq 'anopheles_protein') {
	if ($hname =~ /\.\d$/)  {
	  push @infered_from, "similar to AA sequence:UniProtKB:$hname";
	}
	else {
	  $comm = 1;
	}
      }
      elsif (($logic eq 'other_protein') && ($hname =~ /\.\d$/) ) {
	push @infered_from, "similar to AA sequence:UniProtKB:$hname";
      }
      elsif ($logic eq 'aedes_protein') {
	if ($aedes{$hname})  {
	  push @infered_from, "similar to AA sequence:INSD:$aedes{$hname}";
	}
	else {
	  print STDERR "Failed Aedes 2 protein_id inference lookup for $hname\n";
	}
      }
      elsif ($logic eq 'drosophila_protein') {
	if ($dros{$hname})  {
	  push @infered_from, "similar to AA sequence:UniProtKB:$dros{$hname}";
	}
	else {
	  print STDERR "Failed Drosophila 2 Uniprot inference lookup for $hname\n";
	}
      }

    }

# have now processed all the tsfs
# if found a manual, just print that as a note (no inference tag, even in rare case where might have a ano protein one)
    if ($manual == 1) {
      print OUT "\t\t\tnote\tmanually reviewed by VectorBase\n";
      $man_c ++;
    }

# if comm or genename-based with no other inference print nothing at all
    elsif ( ($comm == 1) &&  (scalar @infered_from == 0) ) {
      $comm_c ++;
    }

# print default note if not manual or comm, and have no inference tags
# as currently implemented, these will be mostly est gene ones
    elsif ( (scalar @infered_from == 0) && ($comm == 0) )  {
      print OUT "\t\t\tnote\tidentified by similarity to sequences in INSD and/or UniProtKB databases\n";
      $default_c ++;
    }

# otherwise print the inference tags (should be 1 or 2 only)
    elsif (scalar @infered_from > 0) {
      $inf_c ++;
      foreach my $inference (@infered_from) {
	print  OUT "\t\t\tinference\t$inference\n";
      }
    }
    else {
      print STDERR "ERROR: unexpected - no note or inference or comm for CDS $translation_name\n";
      $error_c ++;
    }
  }

}


#########################################################
# sub routine to get interpro ids associated with protein
#########################################################
sub get_protein_annotation {
  my ($tr_dbID,$db) = @_;
  my @interpro;

  my $query = "select distinct(i.interpro_ac) from interpro i, protein_feature pf where pf.hit_id = i.id and translation_id = $tr_dbID";

  my $sth = $db->dbc->prepare($query);
  $sth->execute;
  while (my $ipr = $sth->fetchrow) {
    push (@interpro,$ipr);
  }
  return @interpro;
}


################################################
# sub routine to find if any supporting evidence
################################################

#*# change for v45 - count number of distinct hit names in exon sfs
# problem - $tr at this point has new exons that don't have sfs
# can ony get evidence by direct sql querying using $tr dbID
# commented out section in long version below will get proper count
# but decided to just find 0 v non-zero
# returning $evidence = 0 will produce SNAP inference

sub get_evidence {
  my ($tr,$db) = @_;
  my $evidence;

  my $query = "select feature_id from exon_transcript, supporting_feature where exon_transcript.exon_id=supporting_feature.exon_id and transcript_id=? limit 1";

  my $sth = $db->dbc->prepare($query);
  $sth->execute($tr->dbID);
  if ($sth->rows == 0) {
    $evidence = 0;
  }
  else {
    $evidence = 1;
  }

#########  long version #################################
#  my %hits;
#  my $query1 = "select feature_id, feature_type from exon_transcript, supporting_feature where exon_transcript.exon_id=supporting_feature.exon_id and transcript_id=?";
#  my $sth1 = $db->dbc->prepare($query1);
#
#  my $query2 = "select hit_name from dna_align_feature where dna_align_feature_id = ?";
#  my $sth2 = $db->dbc->prepare($query2);
#  my $query3 = "select hit_name from protein_align_feature where protein_align_feature_id = ?";
#  my $sth3 = $db->dbc->prepare($query3);
#
#  $sth1->execute($tr->dbID);
#  while ( my ($fid, $ftype) = $sth1 -> fetchrow_array ) {
#    if ($ftype eq 'dna_align_feature') {
#      $sth2->execute($fid);
#      while ( my $hit = $sth2 -> fetchrow_array ) {
#	$hits{$hit} +=1;
#      }
#    }
#    elsif ($ftype eq 'protein_align_feature') {
#      $sth3->execute($fid);
#      while ( my $hit = $sth3 -> fetchrow_array ) {
#	$hits{$hit} +=1;
#      }
#    }
#  }
#
#  $evidence = scalar(keys %hits);
#
#  print STDERR "sub evidence for ".$tr->dbID." has $evidence supporting hit names\t";
#  my @hit_array = keys %hits;
#  print STDERR "@hit_array\n";
#
#########  long version end #################################

  return $evidence;
}
