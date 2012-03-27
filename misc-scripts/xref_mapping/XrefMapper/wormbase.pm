package XrefMapper::wormbase;

use  XrefMapper::BasicMapper;

use vars qw(@ISA);

@ISA = qw(XrefMapper::BasicMapper);


# Abstact superclass, not intended to be created directly. Contains
# common methods for the "tierII" nematodes in WormBase


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

  #
  # Get the wormbase_locus xrefs for the genes
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


sub gene_description_filter_regexps {

  return ('[0-9A-Z]+\.\d*[A-Z]* PROTEIN[ \.]',
	  '\(\d[A-Z]\d+\)\.',
	  '\([0-9A-Z]+\.\d*[A-Z]* PROTEIN\)[ \.]',
	  '^\(*HYPOTHETICAL\s+.*',
	  '^\s*\(FRAGMENT\)\.?\s*$' );

}

1;
