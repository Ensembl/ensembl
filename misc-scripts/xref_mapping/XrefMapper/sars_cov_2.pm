=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

package XrefMapper::sars_cov_2;

use  XrefMapper::BasicMapper;

use vars '@ISA';

@ISA = qw{ XrefMapper::BasicMapper };
use strict;


sub transcript_names_from_gene {
  my $self = shift;
  my $core_dbi = $self->core->dbc;

  my $reset_sth = $core_dbi->prepare("UPDATE transcript SET display_xref_id = null");
  $reset_sth->execute();
  $reset_sth->finish;

  my $del_xref_sth = $core_dbi->prepare("DELETE x FROM xref x, object_xref ox, external_db e WHERE x.xref_id = ox.xref_id AND ensembl_object_type = 'Transcript' AND e.external_db_id = x.external_db_id AND e.db_name like '%_trans_name'");
  my $del_ox_sth = $core_dbi->prepare("DELETE ox FROM object_xref ox LEFT JOIN xref x ON x.xref_id = ox.xref_id WHERE isnull(x.xref_id)");

  my $xref_id_sth = $core_dbi->prepare("SELECT max(xref_id) FROM xref");
  my $ox_id_sth = $core_dbi->prepare("SELECT max(object_xref_id) FROM object_xref");
  my $ins_xref_sth = $core_dbi->prepare("INSERT IGNORE into xref (xref_id, external_db_id, dbprimary_acc, display_label, version, description, info_type, info_text) values(?, ?, ?, ?, 0, ?, 'MISC', ?)");
  my $ins_ox_sth = $core_dbi->prepare("INSERT into object_xref (object_xref_id, ensembl_id, ensembl_object_type, xref_id) values(?, ?, 'Transcript', ?)");
  my $update_tran_sth = $core_dbi->prepare("UPDATE transcript t SET t.display_xref_id= ? WHERE t.transcript_id=?");
  my $get_genes = $core_dbi->prepare("SELECT g.gene_id, e.db_name, x.dbprimary_acc, x.display_label, x.description FROM gene g, xref x, external_db e where g.display_xref_id = x.xref_id and e.external_db_id = x.external_db_id");
  my $get_transcripts = $core_dbi->prepare("SELECT transcript_id FROM transcript WHERE gene_id = ? ORDER BY seq_region_start, seq_region_end");
  my $get_source_id = $core_dbi->prepare("SELECT external_db_id FROM external_db WHERE db_name like ?");

  $get_genes->execute();
  my ($gene_id, $external_db, $external_db_id, $acc, $label, $description, $transcript_id, $xref_id, $ox_id, $ext, $reuse_xref_id, $info_text);
  $get_genes->bind_columns(\$gene_id, \$external_db, \$acc, \$label, \$description);
  $xref_id_sth->execute();
  $xref_id_sth->bind_columns(\$xref_id);
  $xref_id_sth->fetch();
  $ox_id_sth->execute();
  $ox_id_sth->bind_columns(\$ox_id);
  $ox_id_sth->fetch();

  $del_xref_sth->execute();
  while ($get_genes->fetch()) {
    $get_source_id->execute($external_db . "_trans_name");
    $get_source_id->bind_columns(\$external_db_id);
    $get_source_id->fetch();
    $get_transcripts->execute($gene_id);
    $get_transcripts->bind_columns(\$transcript_id);
    while ($get_transcripts->fetch) {
      $xref_id++;
      $ox_id++;
      $info_text = 'via gene ' . $acc;
      $ins_xref_sth->execute($xref_id, $external_db_id, $label, $label, $description, $info_text);
      $ins_ox_sth->execute($ox_id, $transcript_id, $xref_id);
      $update_tran_sth->execute($xref_id, $transcript_id);
    }
  }

  $del_ox_sth->execute();
  $del_ox_sth->finish();

  $del_xref_sth->finish();
  $xref_id_sth->finish();
  $ox_id_sth->finish();
  $get_genes->finish();
  $get_source_id->finish();
  $get_transcripts->finish();
  $ins_xref_sth->finish();
  $ins_ox_sth->finish();
  $update_tran_sth->finish();
}

1;
