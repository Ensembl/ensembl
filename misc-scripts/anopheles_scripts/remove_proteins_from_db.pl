#!/usr/local/ensembl/bin/perl -w

# script to remove proteins from core database
# leaves transcript and tags gene as bacterial contaminant
# input is a list of protein stable ids
# need to change the following tables:
# mapping_session, stable_id_event, gene_archive, peptide_archive
# object_xref, protein_feature (no identity_xref so no need to change table in this case)
# gene, transcript (note no display_ids to remove in this case)
# translation & translation_stable_id (remove lines), gene_stable_id (version+)

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor (-host=>'ecs1a',
					     -user=>'ensadmin',
					     -dbname=>'anopheles_gambiae_core_19_2b',
					     -pass =>'ensembl'
					    );

my $pa = $db -> get_ProteinAdaptor();
my $ta = $db -> get_TranscriptAdaptor();
my $ga = $db -> get_GeneAdaptor();

# list of ENSANGPnnn.v  for putative bacterial contaminants on UNKN chr
open PROTEIN_LIST, "/nfs/acari/mh4/mozzy/data/bacterial_list_unkn" or die "Can't open list of protein stable ids";

# make array of these translation_stable_ids & make non-redundant array of corresponding gene_stable_ids

my (@protein_list);
while (<PROTEIN_LIST>) {
  chomp;
  /^(\w+)\.(\d)/;
  push @protein_list, $1;
}
print STDERR "Protein list has ".scalar @protein_list." entries\n";
close PROTEIN_LIST;

my (@gene_list, %gene_list_temp);
foreach my $protein_from_list (@protein_list) {
  my $protein = $pa -> fetch_by_translation_stable_id($protein_from_list);
  my $gene_from_list = $protein -> gene;
  my $gene_stable_id = $gene_from_list -> stable_id;
  $gene_list_temp{"$gene_stable_id"} +=1;
}
@gene_list = keys %gene_list_temp;
print STDERR "Gene list has ".scalar @gene_list." entries\n";


## id mapping section ##
# choose new mapping_session_id
my $session = 3;

# make mapping session
my $sth_session = $db -> prepare ("insert mapping_session values ($session,'anopheles_gambiae_core_19_2a','anopheles_gambiae_core_19_2b',NOW())");
$sth_session -> execute();

# take opportunity to do a fix on the dud db names in the previous mapping session
my $sth_fix = $db -> prepare ("update mapping_session set old_db_name='anopheles_gambiae_core_16_2', new_db_name='anopheles_gambiae_core_17_2a', created=20030903153726 where mapping_session_id=1");
$sth_fix -> execute();

# handle to use for all stable id event insertions
my $sth_event = $db->prepare("insert stable_id_event values(?,?,?,?,?,?)");

my ($stable_id, $version);
my $pep_changes=0;
my $gene_changes=0;
my $transcript_copies=0;
my $pep_copies=0;
my $gene_copies=0;

# get each protein from stable_id table
my $sth_prot = $db->prepare("select stable_id, version from translation");
$sth_prot -> execute();
$sth_prot -> bind_columns(\$stable_id, \$version);
while ($sth_prot -> fetch()) {
  my $found_protein = 0;
# check if on 'kill' list
  foreach my $list (@protein_list) {
    if ($list eq $stable_id) {
      $sth_event -> execute($stable_id,$version,'NULL',0,$session,'translation');
      $pep_changes +=1;
      $found_protein =1;
      last;
    }
  }
# if not on kill list ....
  next if ($found_protein);
  $sth_event -> execute($stable_id,$version,$stable_id,$version,$session,'translation');
  $pep_copies+=1;
}
print STDERR "Peptides - lines changed $pep_changes, lines copied $pep_copies\n";

# get each transcript
my $sth_transcript = $db->prepare("select stable_id, version from transcript");
$sth_transcript -> execute();
$sth_transcript -> bind_columns(\$stable_id, \$version);
while ($sth_transcript -> fetch()) {
  $sth_event -> execute($stable_id,$version,$stable_id,$version,$session,'transcript');
  $transcript_copies+=1;
}
print STDERR "Transcripts - lines copied $transcript_copies\n";

# get each gene
my $sth_gene = $db->prepare("select stable_id, version from gene");
$sth_gene -> execute();
$sth_gene -> bind_columns(\$stable_id, \$version);
while ($sth_gene -> fetch()) {
  my $found_gene = 0;
  foreach my $list (@gene_list) {
    if ($list eq $stable_id) {
      $sth_event -> execute($stable_id,$version,$stable_id,$version+1,$session,'gene');
      $gene_changes +=1;
      $found_gene =1;
      last;
    }
  }
  next if ($found_gene);
  $sth_event -> execute($stable_id,$version,$stable_id,$version,$session,'gene');
  $gene_copies+=1;
}
print STDERR "Genes - lines changed $gene_changes, lines copied $gene_copies\n";


# file to store record of pep archive additions
open PEP_ARCH, ">/nfs/acari/mh4/mozzy/data/bacterial_pep_arch" or die "Can't open pep archive output file";
#file to store record of gene archive additions
open GENE_ARCH, ">/nfs/acari/mh4/mozzy/data/bacterial_gene_arch" or die "Can't open gene archive output file";

# prepare peptide & gene archive table change
my $sth_parch = $db->prepare("insert peptide_archive values(?,?,?)");
my $sth_garch = $db->prepare("insert gene_archive values(?,?,?,?,?,?,?)");

open PROTEIN_LIST, "/nfs/acari/mh4/mozzy/data/bacterial_list_unkn" or die "Can't re-open list of protein stable ids";
my $protein_id_count =0;

while (<PROTEIN_LIST>) {
  chomp;
  /^(\w+)\.(\d)/;
  my $protein_stable_id = $1;
  my $version = $2;
  $protein_id_count+=1;
  my $transcript = $ta -> fetch_by_translation_stable_id($protein_stable_id);
  my $seq = $transcript -> translate() -> seq();
  print PEP_ARCH ">$protein_stable_id.$version\n$seq\n";
  $sth_parch -> execute($protein_stable_id,$version,$seq);

  my $transcript_stable_id = $transcript -> stable_id;
  my $transcript_v = $transcript -> version();
  my $protein = $pa -> fetch_by_translation_stable_id($protein_stable_id);
  my $gene = $protein -> gene;
  my $gene_stable_id = $gene -> stable_id;
  my $gene_v = $gene -> version();
  print GENE_ARCH ">$gene_stable_id\t$gene_v\t$transcript_stable_id\t$transcript_v\t$protein_stable_id\t$version\t$session\n";
  $sth_garch -> execute($gene_stable_id,$gene_v,$transcript_stable_id,$transcript_v,$protein_stable_id,$version,$session);
}
print STDERR "found $protein_id_count protein ids and archived\n";
close PROTEIN_LIST;


## section to update gene, transcipt and translation tables plus delete xrefs and protein_features ##
#  prepare to get internal gene, transcript, translation id's for each protein
my $sth_ids = $db->prepare ("select g.gene_id, t.transcript_id, tl.translation_id from gene g, transcript t, translation tl where tl.stable_id = ? and tl.translation_id = t.translation_id and t.gene_id = g.gene_id");

open IDLIST, ">/nfs/acari/mh4/mozzy/data/bacterial_list_unkn_ids" or die "Can't open id output file";
print IDLIST "protein_stable_id\tgene_id\ttranscript_id\ttranslation_id\n";

#  prepare table changes
my $sth_feature = $db->prepare ("delete from protein_feature where translation_id = ?");
my $sth_object = $db->prepare ("delete from object_xref where ensembl_id = ?");
my $sth_transcript = $db->prepare("update transcript set translation_id=0 where transcript_id = ?");
my $sth_translation = $db->prepare("delete from translation where translation_id= ?");
my $sth_gene = $db->prepare("update gene set type='bacterial_contaminant' where gene_id = ?");

# get internal id's for each protein
my ($gene_id, $transcript_id, $translation_id);
my $protein_id_count =0;
foreach my $protein_stable_id (@protein_list) {
  $protein_id_count +=1;
  $sth_ids -> execute ($protein_stable_id);
  $sth_ids -> bind_columns(\$gene_id, \$transcript_id, \$translation_id);
  while ($sth_ids -> fetch()) {
    print IDLIST "$protein_stable_id\t$gene_id\t$transcript_id\t$translation_id\n";
    # delete any entries from protein_feature table
    $sth_feature -> execute($translation_id);
    # delete any rows from object_xref table
    $sth_object -> execute($translation_id);
    # set transcript so has no translation
    $sth_transcript -> execute($transcript_id);
    # delete entries from translation table
    $sth_translation -> execute($translation_id);
    # change gene type
    $sth_gene -> execute($gene_id);
  }
}
print STDERR "found $protein_id_count protein ids and made changes in 6 tables\n";

# finally, update the gsi version using the non-redundant list
my $sth_gsi = $db ->prepare ("update gene set version= ? where stable_id= ?");

my $gene_id_count =0;
foreach my $gene_stable_id (@gene_list) {
  $gene_id_count +=1;
  my $gene = $ga->fetch_by_stable_id($gene_stable_id);
  my $gene_v = $gene ->version();
  my $gene_v_new = $gene_v +1;
  $sth_gsi -> execute($gene_v_new,$gene_stable_id);
}
print STDERR "updated $gene_id_count gsi versions\n";
