use strict;
use warnings;

use SeqStoreConverter::BasicConverter;

package SeqStoreConverter::DanioRerio;

use vars qw(@ISA);

@ISA = qw(SeqStoreConverter::BasicConverter);

sub create_coord_systems {
  my $self = shift;

  $self->debug("DanioRerio Specific: creating chromosome, supercontig, clone "
               . " and chunk coordinate systems");

  my $target = $self->target();
  my $dbh    = $self->dbh();

  my $ass_def = $self->get_default_assembly();

  my @coords =
    (["chromosome" , $ass_def, "top_level,default_version"],
     ["supercontig", undef, "default_version"],
     ["clone"      , undef, "default_version"],
     ["chunk"      , undef, "default_version,sequence_level"]);

  my @assembly_mappings =  ("chromosome:$ass_def|chunk",
                            "clone|chunk",
                            "supercontig|chunk");

  my %cs = (gene                   => ['supercontig','chromosome'],
             transcript             => ['supercontig','chromosome'],
             exon                   => ['supercontig','chromosome'],
             dna_align_feature      => ['chunk'],
             protein_align_feature  => ['chunk'],
             marker_feature         => ['chunk'],
             simple_feature         => ['chunk'],
             repeat_feature         => ['chunk'],
             qtl_feature            => ['chunk'],
             misc_feature           => ['chunk'],
             prediction_transcript  => ['chunk'],
             karyotype              => ['chromosome']);

  $self->debug("Building coord_system table");

  my $sth = $dbh->prepare("INSERT INTO $target.coord_system " .
                           "(name, version, attrib) VALUES (?,?,?)");

  my %coord_system_ids;

  foreach my $cs (@coords) {
    $sth->execute(@$cs);
    $coord_system_ids{$cs->[0]} = $sth->{'mysql_insertid'};
  }
  $sth->finish();

  $self->debug("Building meta_coord table");
  $sth = $dbh->prepare("INSERT INTO $target.meta_coord VALUES (?, ?)");
  foreach my $feature_type (keys %cs) {
    foreach my $coord_sys (@{$cs{$feature_type}}) {
      $sth->execute($feature_type, $coord_system_ids{$coord_sys});
    }
  }
  $sth->finish();



  $self->debug("Adding assembly.mapping entries to meta table");

  $sth = $dbh->prepare("INSERT INTO $target.meta(meta_key, meta_value) " .
                       "VALUES ('assembly.mapping', ?)");

  foreach my $mapping (@assembly_mappings) {
    $sth->execute($mapping);
  }

  $sth->finish();


  return;
}


sub create_seq_regions {
  my $self = shift;

  my $source = $self->source();
  my $target = $self->target();
  my $dbh    = $self->dbh();


  #
  # Turn all of the contents of the contig table into 'chunks' and 
  # give them arbitrary names like chunk1, chunk2. Keep old internal
  # ids for conveneience.
  #

  $self->debug("DanioRerio Specific: creating chunk seq_regions");

  my $sth = $dbh->prepare
    ("INSERT INTO $target.seq_region (seq_region_id, name, coord_system_id, " .
     "                                length) ".
     "SELECT ctg.contig_id, concat('chunk', ctg.contig_id), " .
     "       cs.coord_system_id, ctg.length " .
     "FROM   $source.contig ctg, $target.coord_system cs " .
     "WHERE  cs.name = 'chunk'");

  $sth->execute();

  $sth->finish();


  my $insert_sth = $dbh->prepare
    ("INSERT INTO $target.seq_region (name, coord_system_id, length) " .
     "VALUES (?,?,?)");

  my $tmp_chr_insert_sth = $dbh->prepare
    ("INSERT INTO $target.tmp_chr_map (old_id, new_id) VALUES (?, ?)");

  my $tmp_supercontig_insert_sth = $dbh->prepare
    ("INSERT INTO $target.tmp_superctg_map (name, new_id) VALUES (?,?)");

  my $tmp_clone_insert_sth = $dbh->prepare
    ("INSERT INTO $target.tmp_cln_map (old_id, new_id) VALUES (?,?)");


  #
  # Turn real clones into clones
  #
  $self->debug("DanioRerio Specific: creating clone seq_regions");

  my $select_sth = $dbh->prepare
    ("SELECT ctg.contig_id, ctg.name, ctg.length " .
     "FROM   $source.contig ctg " .
     "WHERE  ctg.name not like 'ctg%' and ctg.name not like 'NA%'");

  my $cs_id = $self->get_coord_system_id('clone');

  $select_sth->execute();

  my ($old_id, $name, $length);
  $select_sth->bind_columns(\$old_id, \$name, \$length);

  while ($select_sth->fetch()) {
    #insert into seq_region table
    $insert_sth->execute($name, $cs_id, $length);
    #copy old/new mapping into temporary table
    $tmp_clone_insert_sth->execute($old_id, $insert_sth->{'mysql_insertid'});
  }

  $select_sth->finish();

  #
  # Turn real chromosomes into chromosomes
  #
  $self->debug("DanioRerio Specific: creating chromosome seq_regions");

  $select_sth = $dbh->prepare
    ("SELECT chr.chromosome_id, chr.name, chr.length " .
     "FROM   $source.chromosome chr " .
     "WHERE  length(chr.name) <= 2");

  $cs_id = $self->get_coord_system_id('chromosome');

  $select_sth->execute();

  $select_sth->bind_columns(\$old_id, \$name, \$length);

  my %chr_id_added;

  while ($select_sth->fetch()) {
    #insert into seq_region table
    $insert_sth->execute($name, $cs_id, $length);
    #copy old/new mapping into temporary table
    $tmp_chr_insert_sth->execute($old_id, $insert_sth->{'mysql_insertid'});
    $chr_id_added{$old_id} = 1;
  }

  $select_sth->finish();

  #
  # Turn supercontigs into supercontigs
  #
  $self->debug("DanioRerio Specific: creating supercontig seq_regions");

  $select_sth = $dbh->prepare
    ("SELECT a.chromosome_id, a.superctg_name, " .
     "       MAX(a.chr_end) - MIN(a.chr_start) + 1 " .
     "FROM   $source.assembly a, $target.coord_system cs " .
     "GROUP BY a.superctg_name");

  $select_sth->execute();
  $select_sth->bind_columns(\$old_id, \$name, \$length);

  $cs_id = $self->get_coord_system_id('supercontig');

  while ($select_sth->fetch()) {
    #insert into seq_region table
    $insert_sth->execute($name, $cs_id, $length);
    #copy old/new mapping into temporary table
   $tmp_supercontig_insert_sth->execute($name,$insert_sth->{'mysql_insertid'});
  }

  $select_sth->finish();
  $tmp_chr_insert_sth->finish();
  $tmp_supercontig_insert_sth->finish();
  $tmp_clone_insert_sth->finish();
  $insert_sth->finish();
}



sub create_assembly {
  my $self = shift;

  #chromosomes are made of chunks
  $self->assembly_contig_chromosome();

  #supercontigs are made of chunks
  $self->assembly_contig_supercontig();

  #clones are made of chunks
  $self->assembly_contig_clone();

  return;
}


sub assembly_contig_clone {
  my $self = shift;


  $self->debug("DanioRerio Specific: building assembly table - chunk/clone");
  #this is easy, there is simply one entire chunk for a given clone

  my $source = $self->source();
  my $target = $self->target();
  my $dbh    = $self->dbh();

  $dbh->do
    ("INSERT INTO $target.assembly (asm_seq_region_id, cmp_seq_region_id, " .
     "                         asm_start, asm_end, cmp_start, cmp_end, ori) " .
     "SELECT tcm.new_id, tcm.old_id, 1, sr.length, 1, sr.length, 1 " .
     "FROM $target.tmp_cln_map tcm, $target.seq_region sr " .
     "WHERE sr.seq_region_id = tcm.new_id");
}



1;
