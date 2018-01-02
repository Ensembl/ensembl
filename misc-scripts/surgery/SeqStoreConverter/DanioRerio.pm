=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

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
    (["chromosome" , $ass_def, "default_version", 1],
     ["supercontig", $ass_def, "default_version", 2],
     ["clone"      , undef, "default_version", 3],
     ["chunk"      , undef, "default_version,sequence_level", 4]);

  my @assembly_mappings =  ("chromosome:$ass_def|chunk",
                            "clone|chunk",
                            "supercontig:$ass_def|chunk",
                            "chromosome:$ass_def|chunk|clone",
                            "supercontig:$ass_def|chunk|clone",
                            "chromosome:$ass_def|chunk|supercontig");

  $self->debug("Building coord_system table");

  my $sth = $dbh->prepare
    ("INSERT INTO $target.coord_system (name, version, attrib, rank) " .
     "VALUES (?,?,?,?)");

  my %coord_system_ids;

  foreach my $cs (@coords) {
    $sth->execute(@$cs);
    $coord_system_ids{$cs->[0]} = $sth->{'mysql_insertid'};
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
  # create a temporary table to hold the ids of all 'toplevel'
  # seq_regions.  Keep the old chromosome_id, and the new seq_region_id
  #
  $dbh->do
    ("CREATE TEMPORARY TABLE $target.tmp_toplevel_map " .
     "(old_id INT, new_id INT, INDEX new_idx(new_id), INDEX old_idx(old_id))");

  my $tmp_toplevel_insert_sth = $dbh->prepare
    ("INSERT INTO $target.tmp_toplevel_map (old_id, new_id) VALUES (?,?)");


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
    my $new_id = $insert_sth->{'mysql_insertid'};
    $tmp_chr_insert_sth->execute($old_id, $new_id);
    $tmp_toplevel_insert_sth->execute($old_id, $new_id);
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
    my $new_id = $insert_sth->{'mysql_insertid'};
    $tmp_supercontig_insert_sth->execute($name,$new_id);

    if(!$chr_id_added{$old_id}) {
      $chr_id_added{$old_id} = 1;
      $tmp_toplevel_insert_sth->execute($old_id, $new_id);
    }
  }

  $select_sth->finish();
  $tmp_chr_insert_sth->finish();
  $tmp_supercontig_insert_sth->finish();
  $tmp_clone_insert_sth->finish();
  $tmp_toplevel_insert_sth->finish();
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



# we need to override the transfer of the genes since danio genes can be on
# supercontigs and on chromosomes
sub transfer_genes {
  my $self = shift;

  my $target = $self->target();
  my $source = $self->source();
  my $dbh    = $self->dbh();

  #
  # Transfer the gene table
  #

  $self->debug("DanioRerio Specific: Building gene table");

  # first transfer genes on chromosomes

  $dbh->do
    ("INSERT INTO $target.gene " .
     "SELECT g.gene_id, g.type, g.analysis_id, toplev.new_id, " .
     "MIN(IF (a.contig_ori=1,(e.contig_start+a.chr_start-a.contig_start)," .
     "       (a.chr_start+a.contig_end-e.contig_end ))) as start, " .
     "MAX(IF (a.contig_ori=1,(e.contig_end+a.chr_start-a.contig_start), " .
     "       (a.chr_start+a.contig_end-e.contig_start))) as end, " .
     "       a.contig_ori*e.contig_strand as strand, " .
     "       g.display_xref_id " .
     "FROM   $source.transcript t, $source.exon_transcript et, " .
     "       $source.exon e, $source.assembly a, $source.gene g, " .
     "       $target.tmp_toplevel_map toplev " .
     "WHERE  t.transcript_id = et.transcript_id " .
     "AND    et.exon_id = e.exon_id " .
     "AND    e.contig_id = a.contig_id " .
     "AND    g.gene_id = t.gene_id " .
     "AND    a.chromosome_id = toplev.old_id " .
     "GROUP BY g.gene_id");


  #
  # Transfer the transcript table
  #

  $self->debug("DanioRerio Specific: Building transcript table ");
  $dbh->do
    ("INSERT INTO $target.transcript " .
     "SELECT t.transcript_id, t.gene_id, toplev.new_id, " .
     "MIN(IF (a.contig_ori=1,(e.contig_start+a.chr_start-a.contig_start)," .
     "       (a.chr_start+a.contig_end-e.contig_end ))) as start, " .
     "MAX(IF (a.contig_ori=1,(e.contig_end+a.chr_start-a.contig_start), " .
     "       (a.chr_start+a.contig_end-e.contig_start))) as end, " .
     "       a.contig_ori*e.contig_strand as strand, " .
     "       t.display_xref_id " .
     "FROM   $source.transcript t, $source.exon_transcript et, " .
     "       $source.exon e, $source.assembly a, " .
     "       $target.tmp_toplevel_map toplev " .
     "WHERE  t.transcript_id = et.transcript_id " .
     "AND    et.exon_id = e.exon_id " .
     "AND    e.contig_id = a.contig_id " .
     "AND    a.chromosome_id = toplev.old_id " .
     "GROUP BY t.transcript_id");

  #
  # Transfer the exon table
  #

  $self->debug("DanioRerio Specific: Building exon table ");

  $dbh->do
    ("INSERT INTO $target.exon " .
     "SELECT e.exon_id, toplev.new_id, " .
     "MIN(IF (a.contig_ori=1,(e.contig_start+a.chr_start-a.contig_start)," .
     "       (a.chr_start+a.contig_end-e.contig_end ))) as start, " .
     "MAX(IF (a.contig_ori=1,(e.contig_end+a.chr_start-a.contig_start), " .
     "       (a.chr_start+a.contig_end-e.contig_start))) as end, " .
     "       a.contig_ori*e.contig_strand as strand, " .
     "       e.phase, e.end_phase " .
     "FROM   $source.transcript t, $source.exon_transcript et, " .
     "       $source.exon e, $source.assembly a, $source.gene g, " .
     "       $target.tmp_toplevel_map toplev " .
     "WHERE  t.transcript_id = et.transcript_id " .
     "AND    et.exon_id = e.exon_id " .
     "AND    e.contig_id = a.contig_id " .
     "AND    g.gene_id = t.gene_id " .
     "AND    a.chromosome_id = toplev.old_id " .
     "GROUP BY e.exon_id");

  #
  # Transfer translation table
  #

  $self->debug("Building translation table");

  $dbh->do
    ("INSERT INTO $target.translation " .
     "SELECT tl.translation_id, ts.transcript_id, tl.seq_start, " .
     "       tl.start_exon_id, tl.seq_end, tl.end_exon_id " .
     "FROM $source.transcript ts, $source.translation tl " .
     "WHERE ts.translation_id = tl.translation_id");

  return;
}



sub set_top_level {
  my $self = shift;

  my $target = $self->target();
  my $dbh = $self->dbh();

  my $attrib_type_id = $self->add_attrib_code();

  $self->debug("DanioRerio Specific: Setting toplevel attributes of " .
               "seq_regions");

  my $sth = $dbh->prepare("DELETE FROM $target.seq_region_attrib " .
                          "WHERE attrib_type_id = ?");
  $sth->execute($attrib_type_id);
  $sth->finish();

  $sth = $dbh->prepare("INSERT INTO $target.seq_region_attrib " .
                       '            (seq_region_id, attrib_type_id, value) ' .
                       "SELECT toplev.new_id, $attrib_type_id, 1 " .
                       "FROM   $target.tmp_toplevel_map toplev ");

  $sth->execute();
  $sth->finish();
}



1;
