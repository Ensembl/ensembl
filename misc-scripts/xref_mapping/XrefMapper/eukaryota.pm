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


package XrefMapper::eukaryota;

use  XrefMapper::BasicMapper;

use vars qw(@ISA);

@ISA = qw(XrefMapper::BasicMapper);


=head2 set_methods

 Overrides the default exonerate method and non default methods which should be used for 
 one or more sources.

=cut

sub set_methods{
 
  my $default_method = 'ExonerateGappedBest1';
  my %override_method_for_source = (  
	   ExonerateGappedBest5 => ['RefSeq_mRNA','RefSeq_mRNA_predicted', 'RefSeq_ncRNA', 'RefSeq_ncRNA_predicted' ],
         );

  return $default_method, \%override_method_for_source;
}


=head2 gene_display_xref_sources

 Overrides the list of sources to use for assigning gene names

=cut

sub gene_display_xref_sources {
    my $self     = shift;

    print STDERR "getting the list of external_dbs for assigning gene names from eukaryota.pm\n";

    my @list = qw(
                 TAIR_SYMBOL
                 RFAM
                 RNAMMER
                 TRNASCAN_SE
                 Uniprot_gn
                 ENA_GENE
                 BROAD_U_maydis
                 BROAD_F_oxysporum
                 BROAD_G_zeae
                 BROAD_G_moniliformis
                 BROAD_P_infestans
                 phyra_jgi_v1.1
                 physo1_jgi_v1.1
	   	 phatr_jgi_v2
		 phatr_jgi_v2_bd
                 PGD_GENE
                 Mycgr3_jgi_v2.0_gene
                 BROAD_Magnaporthe_DB
                 PHYTOZOME_GMAX_GENE
               );
    
    my %ignore;

    
    #don't use EntrezGene labels dependent on predicted RefSeqs

    $ignore{'EntrezGene'} =<<IEG;
SELECT DISTINCT ox.object_xref_id
  FROM object_xref ox, dependent_xref dx, 
       xref xmas, xref xdep, 
       source smas, source sdep
    WHERE ox.xref_id = dx.dependent_xref_id AND
          dx.dependent_xref_id = xdep.xref_id AND
          dx.master_xref_id = xmas.xref_id AND
          xmas.source_id = smas.source_id AND
          xdep.source_id = sdep.source_id AND
          smas.name like "Refseq%predicted" AND
          sdep.name like "EntrezGene" AND
          ox.ox_status = "DUMP_OUT" AND
          ox.master_xref_id = dx.master_xref_id
IEG

    #don't use labels starting with LOC

    $ignore{'LOC_prefix'} =<<LOCP;
SELECT object_xref_id
  FROM object_xref JOIN xref USING(xref_id) JOIN source USING(source_id)
   WHERE ox_status = 'DUMP_OUT' AND label REGEXP '^LOC[[:digit:]]+'
LOCP

    return [\@list,\%ignore];
}


=head2 transcript_display_xref_sources

 Overrides the list of sources to use for assigning transcript names

=cut

sub transcript_display_xref_sources {
    my $self     = shift;

    print STDERR "getting the list of external_dbs for assigning transcript names from eukaryota.pm\n";

    my @list = qw(
                 RFAM
                 RNAMMER
                 TRNASCAN_SE
                 Uniprot_gn_trans_name
                 ENA_GENE
                 BROAD_U_maydis
                 BROAD_F_oxysporum
                 BROAD_G_zeae
                 BROAD_G_moniliformis
                 BROAD_P_infestans
                 phyra_jgi_v1.1
                 physo1_jgi_v1.1
	   	 phatr_jgi_v2
		 phatr_jgi_v2_bd
                 PGD_GENE
                 Mycgr3_jgi_v2.0_gene
                 BROAD_Magnaporthe_DB
                 PHYTOZOME_GMAX_GENE
               );
    
    my %ignore;

    
    #don't use EntrezGene labels dependent on predicted RefSeqs

    $ignore{'EntrezGene'} =<<IEG;
SELECT DISTINCT ox.object_xref_id
  FROM object_xref ox, dependent_xref dx, 
       xref xmas, xref xdep, 
       source smas, source sdep
    WHERE ox.xref_id = dx.dependent_xref_id AND
          dx.dependent_xref_id = xdep.xref_id AND
          dx.master_xref_id = xmas.xref_id AND
          xmas.source_id = smas.source_id AND
          xdep.source_id = sdep.source_id AND
          smas.name like "Refseq%predicted" AND
          sdep.name like "EntrezGene" AND
          ox.ox_status = "DUMP_OUT" AND
          ox.master_xref_id = dx.master_xref_id
IEG

    #don't use labels starting with LOC

    $ignore{'LOC_prefix'} =<<LOCP;
SELECT object_xref_id
  FROM object_xref JOIN xref USING(xref_id) JOIN source USING(source_id)
   WHERE ox_status = 'DUMP_OUT' AND label REGEXP '^LOC[[:digit:]]+'
LOCP

    return [\@list,\%ignore];
}


=head2 gene_description_sources

 Overrides the list of external_db entries to use for assigning gene descriptions

=cut

sub gene_description_sources {
  return (
          "TAIR_LOCUS",
          "PomBase_GENE",
          "PomBase_TRANSCRIPT",
          "Uniprot/SWISSPROT",
          "Uniprot/SPTREMBL",
          "BROAD_U_maydis",
          "BROAD_F_oxysporum",
          "BROAD_G_zeae",
          "BROAD_G_moniliformis",
          "BROAD_P_infestans",
          "phyra_jgi_v1.1",
          "physo1_jgi_v1.1",
	  "phatr_jgi_v2",
	  "phatr_jgi_v2_bd",
          "PGD_GENE",
          "BROAD_Magnaporthe_DB",
          "PGSC_GENE",
          "PHYTOZOME_GMAX_GENE",
          "RFAM",
          "TRNASCAN_SE",
          "RNAMMER",
         );
}


=head2 transcript_names_from_gene

 Overrides the transcript names logic assignment from gene names
 Avoid adding '-\d+' suffix to any of them

=cut


sub transcript_names_from_gene {
  my $self = shift;

  print "Assigning transcript names from gene names\n" if ($self->verbose);

  my $reset_sth = $self->core->dbc->prepare("UPDATE transcript SET display_xref_id = null");
  $reset_sth->execute();
  $reset_sth->finish;

  my $xref_id_sth = $self->core->dbc->prepare("SELECT max(xref_id) FROM xref");
  my $ox_id_sth = $self->core->dbc->prepare("SELECT max(object_xref_id) FROM object_xref");
  my $del_xref_sth = $self->core->dbc->prepare("DELETE x FROM xref x, object_xref ox WHERE x.xref_id = ox.xref_id AND ensembl_object_type = 'Transcript' AND display_label REGEXP '-2[0-9]{2}\$'");
  my $reuse_xref_sth = $self->core->dbc->prepare("SELECT xref_id FROM xref x WHERE external_db_id = ? AND display_label = ? AND version = 0 AND description = ? AND info_type = 'MISC' AND info_text = 'via gene name'");
  my $del_ox_sth = $self->core->dbc->prepare("DELETE ox FROM object_xref ox LEFT JOIN xref x ON x.xref_id = ox.xref_id WHERE isnull(x.xref_id)");
  my $ins_xref_sth = $self->core->dbc->prepare("INSERT IGNORE into xref (xref_id, external_db_id, dbprimary_acc, display_label, version, description, info_type, info_text) values(?, ?, ?, ?, 0, ?, 'MISC', 'via gene name')");
  my $ins_ox_sth = $self->core->dbc->prepare("INSERT into object_xref (object_xref_id, ensembl_id, ensembl_object_type, xref_id) values(?, ?, 'Transcript', ?)");
  my $update_tran_sth = $self->core->dbc->prepare("UPDATE transcript t SET t.display_xref_id= ? WHERE t.transcript_id=?");

  my $get_genes = $self->core->dbc->prepare("SELECT g.gene_id, e.db_name, x.dbprimary_acc, x.display_label, x.description FROM gene g, xref x, external_db e where g.display_xref_id = x.xref_id and e.external_db_id = x.external_db_id");
  my $get_transcripts = $self->core->dbc->prepare("SELECT transcript_id FROM transcript WHERE gene_id = ? ORDER BY seq_region_start, seq_region_end");
  my $get_source_id = $self->core->dbc->prepare("SELECT external_db_id FROM external_db WHERE db_name like ?");

  $get_genes->execute();
  my ($gene_id, $external_db, $external_db_id, $acc, $label, $description, $transcript_id, $xref_id, $ox_id, $ext, $reuse_xref_id);
  $get_genes->bind_columns(\$gene_id, \$external_db, \$acc, \$label, \$description);
  $xref_id_sth->execute();
  $xref_id_sth->bind_columns(\$xref_id);
  $xref_id_sth->fetch();
  $ox_id_sth->execute();
  $ox_id_sth->bind_columns(\$ox_id);
  $ox_id_sth->fetch();
  $del_xref_sth->execute();
  while ($get_genes->fetch()) {
    my $ext = '';
    my $index=0;
    $get_source_id->execute($external_db . "_trans_name");
    $get_source_id->bind_columns(\$external_db_id);
    $get_source_id->fetch();
    $get_transcripts->execute($gene_id);
    $get_transcripts->bind_columns(\$transcript_id);
    while ($get_transcripts->fetch) {
      $xref_id++;
      $ox_id++;
      if ($ext ne '') {
	  $reuse_xref_sth->execute($external_db_id, $label . '-' . $ext, $description);
      }
      else {
	  $reuse_xref_sth->execute($external_db_id, $label, $description);
      }
      $reuse_xref_sth->bind_columns(\$reuse_xref_id);
      if ($reuse_xref_sth->fetch()) {
        $ins_ox_sth->execute($ox_id, $transcript_id, $reuse_xref_id);
        $update_tran_sth->execute($reuse_xref_id, $transcript_id);
      } else {
	  if ($ext ne '') {
	      $ins_xref_sth->execute($xref_id, $external_db_id, $label. "-" . $ext, $label . "-" . $ext, $description);
	  }
	  else {
	      $ins_xref_sth->execute($xref_id, $external_db_id, $label, $label, $description);
	  }
        $ins_ox_sth->execute($ox_id, $transcript_id, $xref_id);
        $update_tran_sth->execute($xref_id, $transcript_id);
      }
      $index++;
    }
  }

  $del_xref_sth->finish();
  $del_ox_sth->execute();
  $del_ox_sth->finish();
  $reuse_xref_sth->finish();
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
