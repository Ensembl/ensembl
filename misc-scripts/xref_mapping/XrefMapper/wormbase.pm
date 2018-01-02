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

package XrefMapper::wormbase;

use  XrefMapper::BasicMapper;

use vars qw(@ISA);

@ISA = qw(XrefMapper::BasicMapper);


# This module is activated by specifying "taxon=wormbase" in the mapping input file
# It contains some common config for worms maintained by WormBase (i.e. having genes
# with WBGene ids etc)

sub set_methods{
 
  my $default_method = 'ExonerateGappedBest1';
  my %override_method_for_source = (  
    ExonerateGappedBest5_subtran => ['RefSeq_mRNA',
                                     'RefSeq_mRNA_predicted', 
                                     'RefSeq_ncRNA', 
                                     'RefSeq_ncRNA_predicted' ],
      );

  return $default_method, \%override_method_for_source;
}



sub set_display_xrefs {
  my ($self) = @_;

  # wormbase_gene
  # wormbase_transcript
  # wormbase_locus
  # wormpep_id

  print "Building Transcript and Gene display_xrefs using WormBase direct xrefs\n" if ($self->verbose);


  # strategy:
  # - populate transcript display xref with wormbase_transcript
  # - populate gene display xref with wormbase_locus 

  my (%external_dbs, %gene_display_xrefs, %tran_display_xrefs);

  #
  # Get external_db ids for the sources we are interested in
  #
  my $edb_sth = $self->core->dbc->prepare("SELECT external_db_id, db_name from external_db WHERE db_name like 'wormbase%'");
  $edb_sth->execute;
  while( my ($edb_id, $edb_name) = $edb_sth->fetchrow_array ) {
    $external_dbs{$edb_name} = $edb_id;
  }
  $edb_sth->finish;

  if (not exists $external_dbs{wormbase_transcript} or
      not exists $external_dbs{wormbase_locus}) {
    print "Could not find wormbase_transcript and wormbase_locus in external_db table, so doing nothing\n" if $self->verbose;
    return;
  }

  my $query_gseq_sth =  $self->core->dbc->prepare("SELECT ensembl_id, x.xref_id " .
                                                 "FROM object_xref ox, xref x ". 
                                                 "WHERE ox.xref_id = x.xref_id AND external_db_id = " . $external_dbs{wormbase_gseqname});
  $query_gseq_sth->execute();
  while( my ($gid, $xid) = $query_gseq_sth->fetchrow_array) {
    $gene_display_xrefs{$gid} = $xid;
  }
  $query_gseq_sth->finish;
  


  #
  # Some genes will have a locus name. Over-write display xrefs for those that do
  #
  my $query_gene_sth = $self->core->dbc->prepare("SELECT ensembl_id, x.xref_id " .
                                                 "FROM object_xref ox, xref x ". 
                                                 "WHERE ox.xref_id = x.xref_id AND external_db_id = " . $external_dbs{wormbase_locus});

  $query_gene_sth->execute();
  while( my ($gid, $xid) = $query_gene_sth->fetchrow_array) {
    $gene_display_xrefs{$gid} = $xid;
  }
  $query_gene_sth->finish;


  #
  # Get the wormbase_transcript xrefs for the genes
  #
  my $query_tran_sth = $self->core->dbc->prepare("SELECT ensembl_id, x.xref_id " . 
                                                 "FROM object_xref ox, xref x " . 
                                                 "WHERE ox.xref_id = x.xref_id and external_db_id = " . $external_dbs{wormbase_transcript});
  $query_tran_sth->execute();
  while( my ($gid, $xid) = $query_tran_sth->fetchrow_array) {
    $tran_display_xrefs{$gid} = $xid;
  }
  $query_tran_sth->finish;


  #
  # finally, update
  #
  my $reset_sth = $self->core->dbc->prepare("UPDATE gene SET display_xref_id = null");
  $reset_sth->execute();
  $reset_sth->finish;
 
  $reset_sth = $self->core->dbc->prepare("UPDATE transcript SET display_xref_id = null");
  $reset_sth->execute();
  $reset_sth->finish;

  my $update_gene_sth = $self->core->dbc->prepare("UPDATE gene g SET g.display_xref_id= ? WHERE g.gene_id=?");
  my $update_tran_sth = $self->core->dbc->prepare("UPDATE transcript t SET t.display_xref_id= ? WHERE t.transcript_id=?");

  foreach my $gid (keys %gene_display_xrefs) {
    $update_gene_sth->execute( $gene_display_xrefs{$gid}, $gid );
  }
  $update_gene_sth->finish;

  foreach my $tid (keys %tran_display_xrefs) {
    $update_tran_sth->execute( $tran_display_xrefs{$tid}, $tid );
  }
  $update_tran_sth->finish;

  print "Updated display xrefs in core for genes and transcripts\n" if $self->verbose;
}


# over-ride the following, to ensure that our carefully constructed transcript
# display ids are not stamped over by the default behaviour (propagation from
# gene)
sub transcript_names_from_gene {
  return;
}


sub gene_description_sources {

  return ("RFAM",
          "RNAMMER",
          "TRNASCAN_SE",
	  "miRBase",
          "HGNC",
          "IMGT/GENE_DB",
	  "Uniprot/SWISSPROT",
	  "RefSeq_peptide",
          "Uniprot/SPTREMBL",
      );

}


sub gene_description_filter_regexps {

  return ( '^(Protein \S+\s*)+$',
           '^Uncharacterized protein\s*\S+\s*',
           '^Uncharacterized protein\s*',
           '^Putative uncharacterized protein\s*\S+\s*',
           '^Putative uncharacterized protein\s*',
           '^Hypothetical protein\s*\S+\s*',
   );

}

1;
