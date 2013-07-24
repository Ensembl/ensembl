#!/usr/local/ensembl/bin/perl -w

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Getopt::Long;

use strict;

my $pass;
my $patchtype_file = "./data/patch_type";
my $dbname;
my $host;
my $user;
my $port = 3306;
my $central_coord_system = 'supercontig';
my $toplevel_coord_system = 'chromosome';

&GetOptions(
            'pass=s'         => \$pass,
            'patchtype_file=s' => \$patchtype_file,
            'host=s'         => \$host,
            'dbname=s'       => \$dbname,
            'user=s'         => \$user,
            'port=n'         => \$port,
           );

open(TYPE,"<".$patchtype_file)         || die "Could not open file $patchtype_file";
open(SQL,">delete_patch.sql") || die "Could not open delete_patch.sql for writing\n";
#connect to the database

my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    '-host'    => $host,
    '-user'    => $user,
    '-pass'    => $pass,
    '-dbname'  => $dbname,
    '-species' => "load"
    );

my $get_all_synonyms_sth = $dba->dbc->prepare("select synonym from seq_region, seq_region_synonym, coord_system where seq_region.name in (select seq_region.name from seq_region, seq_region_attrib, attrib_type where
attrib_type.code in ('patch_fix','patch_novel') and attrib_type.attrib_type_id = seq_region_attrib.attrib_type_id and seq_region_attrib.seq_region_id = seq_region.seq_region_id) and seq_region.coord_system_id=coord_system.coord_system_id and coord_system.name = '".$central_coord_system."' and seq_region.seq_region_id=seq_region_synonym.seq_region_id")
|| die "Could not prepare to get synonyms";

if ($get_all_synonyms_sth == 0) {
  throw("Could not prepare to get synonyms");
}
my $get_name_sth = $dba->dbc->prepare("select name from seq_region,seq_region_synonym where synonym = ? and seq_region_synonym.seq_region_id = seq_region.seq_region_id")
|| die "Could not prepare to get name";

if ($get_name_sth == 0) {
  throw("Could not prepare to get name");
}
my $get_synonym_sth = $dba->dbc->prepare("select synonym from seq_region_synonym, seq_region, coord_system where seq_region. name = ? and seq_region.coord_system_id=coord_system.coord_system_id and coord_system.name = ? and seq_region.seq_region_id=seq_region_synonym.seq_region_id")
  || die "Could not prepare to get synonym";

if ($get_synonym_sth == 0) {
  throw("Could not prepare to get synonym");
}
my $get_seq_region_id_sth = $dba->dbc->prepare("select seq_region_id from seq_region, coord_system where seq_region. name = ? and seq_region.coord_system_id=coord_system.coord_system_id and
coord_system.name = ?")
  || die "Could not prepare to get seq_region_id";

if ($get_seq_region_id_sth == 0) {
  throw("Could not prepare to get seq_region_id");
}
#get the full list of existing synonyms/accessions
my $existing_synonym;
my %existing_synonyms;

$get_all_synonyms_sth->execute() || die "problem executing seq_region_id";
$get_all_synonyms_sth->bind_columns(\$existing_synonym) || die "problem binding";
while($get_all_synonyms_sth->fetch()){
  print "Synonym ".$existing_synonym."\n";
  $existing_synonyms{$existing_synonym} = 1;
}


#for the new files
while (<TYPE>) {
  chomp;
  next if (/^#/);
  my ($alt_scaf_name,$alt_scaf_acc,$type) = split(/\t/,$_);
  if (!$alt_scaf_name || !$alt_scaf_acc || !$type) {
    throw("Unable to find name, accession or type");
  }

  #check to see if this is an updated patch
  my $synonym;
  $get_synonym_sth->execute($alt_scaf_name, $central_coord_system) || die "problem executing get synonym";
  $get_synonym_sth->bind_columns(\$synonym) || die "problem binding";
  $get_synonym_sth->fetch();

  if(defined($synonym)){
    #if patch still exists (even if modified) record by setting to 0
    #modified patches are dealt with here already
    $existing_synonyms{$synonym} = 0;
    #unchanged
    if($synonym eq $alt_scaf_acc){
      print $alt_scaf_name." exists and is unchanged\n";
    }
    #modified
    else{
      print $alt_scaf_name." has been modified acc ".$synonym." to ".$alt_scaf_acc."\n";
      remove_patch($alt_scaf_name);
    }
  }
  #new
  else{
    print $alt_scaf_name." is a new patch\n";
  }

}
#removed
my @exist_accs = keys  %existing_synonyms;
foreach my $acc (@exist_accs){
  if($existing_synonyms{$acc}){
    my $name;
    $get_name_sth->execute($acc) || die "problem executing get name";
    $get_name_sth->bind_columns(\$name) || die "problem binding";
    $get_name_sth->fetch();
    print $name." with acc ".$acc." has been removed from the new set\n";
    remove_patch($name);
  }
}


sub remove_patch{

  my $alt_scaf_name = shift;
  #get the patch seq_region_ids (chrom and scaffold)
  my($scaf_id, $chrom_id);
  $get_seq_region_id_sth->execute($alt_scaf_name, $central_coord_system) || die "problem executing seq_region_id";
  $get_seq_region_id_sth->bind_columns(\$scaf_id) || die "problem binding";
  $get_seq_region_id_sth->fetch();

  $get_seq_region_id_sth->execute($alt_scaf_name, $toplevel_coord_system) || die "problem executing seq_region_id";
  $get_seq_region_id_sth->bind_columns(\$chrom_id) || die "problem binding";
  $get_seq_region_id_sth->fetch();

  #check for components only used in patch
  my $scaf_comp;
  my $get_region_components_sth = $dba->dbc->prepare("select distinct cmp_seq_region_id from assembly where asm_seq_region_id = ?") || die "Couldn't get components";
  $get_region_components_sth->execute($scaf_id) || die "problem executing get components";
  $get_region_components_sth->bind_columns(\$scaf_comp) || die "problem binding";
  COMP: while($get_region_components_sth->fetch()){
    #check each component
    my $asm_id;
    my $get_asm_ids_sth = $dba->dbc->prepare("select distinct asm_seq_region_id from assembly where cmp_seq_region_id = ?") || die "Couldn't get asm ids";
    $get_asm_ids_sth->execute($scaf_comp) || die "problem executing get components";
    $get_asm_ids_sth->bind_columns(\$asm_id) || die "problem binding";
    while($get_asm_ids_sth->fetch()){
      if($asm_id != $chrom_id and $asm_id != $scaf_id){
        #print "Used outside of the patch\n";
        next COMP;
      }
    }
    #component only used in patch
    #print $scaf_comp." only used in patch\n";
    print SQL "delete from dna where seq_region_id = ".$scaf_comp.";\n";
    print SQL "delete from seq_region_attrib where seq_region_id = ".$scaf_comp.";\n";
    print SQL "delete from seq_region where seq_region_id = ".$scaf_comp.";\n";
  }

  print SQL "delete from seq_region_synonym where seq_region_id = ".$scaf_id.";\n";

  print SQL "delete from seq_region_attrib where seq_region_id = ".$chrom_id.";\n";

  print SQL "delete from assembly_exception where seq_region_id = ".$chrom_id.";\n";

  print SQL "delete from assembly where asm_seq_region_id = ".$chrom_id.";\n";
  print SQL "delete from assembly where asm_seq_region_id = ".$scaf_id.";\n";

  print SQL "delete from seq_region where seq_region_id = ".$chrom_id.";\n";
  print SQL "delete from seq_region where seq_region_id = ".$scaf_id.";\n";

  print SQL "delete from splicing_event where seq_region_id = ".$chrom_id.";\n";

  print SQL "delete from marker_feature where seq_region_id = ".$chrom_id.";\n";
}

