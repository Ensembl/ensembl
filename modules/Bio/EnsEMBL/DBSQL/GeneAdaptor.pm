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



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME
Bio::EnsEMBL::DBSQL::GeneAdaptor - Database adaptor for the retrieval and
storage of Gene objects

=head1 SYNOPSIS

  use Bio::EnsEMBL::Registry;

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
  );

  $gene_adaptor =
    Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "gene" );

  $gene = $gene_adaptor->fetch_by_dbID(1234);

  $gene = $gene_adaptor->fetch_by_stable_id('ENSG00000184129');

  @genes = @{ $gene_adaptor->fetch_all_by_external_name('BRCA2') };

  $slice_adaptor =
    Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "slice" );

  $slice =
    $slice_adaptor->fetch_by_region( 'chromosome', '1', 1, 1000000 );

  @genes = @{ $gene_adaptor->fetch_all_by_Slice($slice) };

=head1 DESCRIPTION

This is a database aware adaptor for the retrieval and storage of gene
objects.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::GeneAdaptor;

use strict;

use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning );
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;

use vars '@ISA';
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);

# _tables
#  Arg [1]    : none
#  Description: PROTECTED implementation of superclass abstract method.
#               Returns the names, aliases of the tables to use for queries.
#  Returntype : list of listrefs of strings
#  Exceptions : none
#  Caller     : internal
#  Status     : Stable

sub _tables {
  return (['gene', 'g'], ['xref', 'x'], ['external_db', 'exdb']);
}

# _columns
#  Arg [1]    : none
#  Example    : none
#  Description: PROTECTED implementation of superclass abstract method.
#               Returns a list of columns to use for queries.
#  Returntype : list of strings
#  Exceptions : none
#  Caller     : internal
#  Status     : Stable

sub _columns {
  my ($self) = @_;

  my $created_date  = $self->db()->dbc()->from_date_to_seconds("g.created_date");
  my $modified_date = $self->db()->dbc()->from_date_to_seconds("g.modified_date");

  return ('g.gene_id', 'g.seq_region_id', 'g.seq_region_start', 'g.seq_region_end', 'g.seq_region_strand', 'g.analysis_id', 'g.biotype', 'g.display_xref_id', 'g.description', 'g.source', 'g.is_current', 'g.canonical_transcript_id', 'g.stable_id', 'g.version', $created_date, $modified_date, 'x.display_label', 'x.dbprimary_acc', 'x.description', 'x.version', 'exdb.db_name', 'exdb.status', 'exdb.db_release', 'exdb.db_display_name', 'x.info_type', 'x.info_text');
}

sub _left_join {
  return (['xref', "x.xref_id = g.display_xref_id"], ['external_db', "exdb.external_db_id = x.external_db_id"]);
}

=head2 list_dbIDs

  Example    : @gene_ids = @{$gene_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all genes in the current db
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.
  Returntype : Listref of Ints
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub list_dbIDs {
  my ($self, $ordered) = @_;

  return $self->_list_dbIDs("gene", undef, $ordered);
}

=head2 list_stable_ids

  Example    : @stable_gene_ids = @{$gene_adaptor->list_stable_ids()};
  Description: Gets an listref of stable ids for all genes in the current db
  Returntype : reference to a list of strings
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub list_stable_ids {
  my ($self) = @_;

  return $self->_list_dbIDs("gene", "stable_id");
}

sub list_seq_region_ids {
  my $self = shift;

  return $self->_list_seq_region_ids('gene');
}

=head2 fetch_by_display_label

  Arg [1]    : String $label - display label of gene to fetch
  Example    : my $gene = $geneAdaptor->fetch_by_display_label("BRCA2");
  Description: Returns the gene which has the given display label or undef if
               there is none. If there are more than 1, the gene on the 
               reference slice is reported or if none are on the reference,
               the first one is reported.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_display_label {
  my $self  = shift;
  my $label = shift;

  my $constraint = "x.display_label = ? AND g.is_current = 1";
  $self->bind_param_generic_fetch($label, SQL_VARCHAR);
  my @genes = @{$self->generic_fetch($constraint)};
  my $gene;
  if (scalar(@genes) > 1) {
    foreach my $gene_tmp (@genes) {
      if ($gene_tmp->slice->is_reference) {
        $gene = $gene_tmp;
      }
      last if ($gene);
    }
    if (!$gene) {
      $gene = $genes[0];
    }

  } elsif (scalar(@genes) == 1) {
    $gene = $genes[0];
  }

  return $gene;
} ## end sub fetch_by_display_label

=head2 fetch_all_by_display_label

  Arg [1]    : String $label - display label of genes to fetch
  Example    : my @genes = @{$geneAdaptor->fetch_all_by_display_label("PPP1R2P1")};
  Description: Returns all genes which have the given display label or undef if
               there are none. 
  Returntype : listref of Bio::EnsEMBL::Gene objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_display_label {
  my $self  = shift;
  my $label = shift;

  my $constraint = "x.display_label = ? AND g.is_current = 1";
  $self->bind_param_generic_fetch($label, SQL_VARCHAR);
  my $genes = $self->generic_fetch($constraint);

  return $genes;
}

=head2 fetch_by_stable_id

  Arg [1]    : String $id 
               The stable ID of the gene to retrieve
  Example    : $gene = $gene_adaptor->fetch_by_stable_id('ENSG00000148944');
  Description: Retrieves a gene object from the database via its stable id.
               The gene will be retrieved in its native coordinate system (i.e.
               in the coordinate system it is stored in the database). It may
               be converted to a different coordinate system through a call to
               transform() or transfer(). If the gene or exon is not found
               undef is returned instead.
  Returntype : Bio::EnsEMBL::Gene or undef
  Exceptions : if we cant get the gene in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_stable_id {
  my ($self, $stable_id) = @_;

  my $constraint = "g.stable_id = ? AND g.is_current = 1";
  $self->bind_param_generic_fetch($stable_id, SQL_VARCHAR);
  my ($gene) = @{$self->generic_fetch($constraint)};

  # If we didn't get anything back, desperately try to see if there's
  # a version number in the stable_id
  if(!defined($gene) && (my $vindex = rindex($stable_id, '.'))) {
      $gene = $self->fetch_by_stable_id_version(substr($stable_id,0,$vindex),
						substr($stable_id,$vindex+1));
  }

  return $gene;
}

=head2 fetch_by_stable_id_version

  Arg [1]    : String $id 
               The stable ID of the gene to retrieve
  Arg [2]    : Integer $version
               The version of the stable_id to retrieve
  Example    : $gene = $gene_adaptor->fetch_by_stable_id('ENSG00000148944', 14);
  Description: Retrieves a gene object from the database via its stable id and version.
               The gene will be retrieved in its native coordinate system (i.e.
               in the coordinate system it is stored in the database). It may
               be converted to a different coordinate system through a call to
               transform() or transfer(). If the gene or exon is not found
               undef is returned instead.
  Returntype : Bio::EnsEMBL::Gene or undef
  Exceptions : if we cant get the gene in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_stable_id_version {
    my ($self, $stable_id, $version) = @_;

    # Enforce that version be numeric
    return unless($version =~ /^\d+$/);

    my $constraint = "g.stable_id = ? AND g.version = ? AND g.is_current = 1";
    $self->bind_param_generic_fetch($stable_id, SQL_VARCHAR);
    $self->bind_param_generic_fetch($version, SQL_INTEGER);
    my ($gene) = @{$self->generic_fetch($constraint)};

    return $gene;
}

=head2 fetch_all_by_source

  Arg [1]    : String $source
               listref of $sources
               The source of the gene to retrieve. You can have as an argument a reference
               to a list of sources
  Example    : $genes = $gene_adaptor->fetch_all_by_source('havana'); 
               $genes = $gene_adaptor->fetch_all_by_source(['ensembl', 'vega']);
  Description: Retrieves an array reference of gene objects from the database via its source or sources.
               The gene will be retrieved in its native coordinate system (i.e.
               in the coordinate system it is stored in the database). It may
               be converted to a different coordinate system through a call to
               transform() or transfer(). If the gene or exon is not found
               undef is returned instead.
  Returntype  : listref of Bio::EnsEMBL::Gene
  Exceptions : if we cant get the gene in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_source {
  my ($self, $source) = @_;
  my @genes = @{$self->generic_fetch($self->source_constraint($source))};
  return \@genes;
}

=head2 source_constraint 

  Arg [1]    : String $source
               listref of $sources
               The source of the gene to retrieve. You can have as an argument a reference
               to a list of sources
  Description: Used internally to generate a SQL constraint to restrict a transcript query by source
  Returntype  : String
  Exceptions : If source is not supplied
  Caller     : general
  Status     : Stable

=cut

sub source_constraint {
  my ($self, $sources, $inline_variables) = @_;
  my $constraint = "g.is_current = 1";
  my $in_statement = $self->generate_in_constraint($sources, 'g.source', SQL_VARCHAR, $inline_variables);
  $constraint .= " and $in_statement";
  return $constraint;
}

=head2 count_all_by_source

  Arg [1]     : String $source
                listref of $source
                The source of the gene to retrieve. You can have as an argument a reference
                to a list of sources
  Example     : $cnt = $gene_adaptor->count_all_by_source('ensembl'); 
                $cnt = $gene_adaptor->count_all_by_source(['havana', 'vega']);
  Description : Retrieves count of gene objects from the database via its source or sources.
  Returntype  : integer
  Caller      : general
  Status      : Stable

=cut

sub count_all_by_source {
  my ($self, $source) = @_;
  return $self->generic_count($self->source_constraint($source));
}

=head2 fetch_all_by_biotype 

  Arg [1]    : String $biotype 
               listref of $biotypes
               The biotype of the gene to retrieve. You can have as an argument a reference
               to a list of biotypes
  Example    : $gene = $gene_adaptor->fetch_all_by_biotype('protein_coding'); 
               $gene = $gene_adaptor->fetch_all_by_biotypes(['protein_coding', 'sRNA', 'miRNA']);
  Description: Retrieves an array reference of gene objects from the database via its biotype or biotypes.
               The genes will be retrieved in its native coordinate system (i.e.
               in the coordinate system it is stored in the database). It may
               be converted to a different coordinate system through a call to
               transform() or transfer(). If the gene or exon is not found
               undef is returned instead.
  Returntype  : listref of Bio::EnsEMBL::Gene
  Exceptions : if we cant get the gene in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_biotype {
  my ($self, $biotype) = @_;
  my @genes = @{$self->generic_fetch($self->biotype_constraint($biotype))};
  return \@genes;
}

=head2 biotype_constraint 

  Arg [1]    : String $biotypes 
               listref of $biotypes
               The biotype of the gene to retrieve. You can have as an argument a reference
               to a list of biotypes
  Description: Used internally to generate a SQL constraint to restrict a gene query by biotype
  Returntype  : String
  Exceptions : If biotype is not supplied
  Caller     : general
  Status     : Stable

=cut

sub biotype_constraint {
  my ($self, $biotypes, $inline_variables) = @_;
  my $constraint = "g.is_current = 1";
  my $in_statement = $self->generate_in_constraint($biotypes, 'g.biotype', SQL_VARCHAR, $inline_variables);
  $constraint .= " and $in_statement";
  return $constraint;
}

=head2 count_all_by_biotype 

  Arg [1]     : String $biotype 
                listref of $biotypes
                The biotype of the gene to retrieve. You can have as an argument a reference
                to a list of biotypes
  Example     : $cnt = $gene_adaptor->count_all_by_biotype('protein_coding'); 
                $cnt = $gene_adaptor->count_all_by_biotypes(['protein_coding', 'sRNA', 'miRNA']);
  Description : Retrieves count of gene objects from the database via its biotype or biotypes.
  Returntype  : integer
  Caller      : general
  Status      : Stable

=cut

sub count_all_by_biotype {
  my ($self, $biotype) = @_;
  return $self->generic_count($self->biotype_constraint($biotype));
}

sub fetch_all {
  my ($self)     = @_;
  my $constraint = 'g.biotype != "LRG_gene" and g.is_current = 1';
  my @genes      = @{$self->generic_fetch($constraint)};
  return \@genes;
}

=head2 fetch_all_versions_by_stable_id 

  Arg [1]     : String $stable_id 
                The stable ID of the gene to retrieve
  Example     : $gene = $gene_adaptor->fetch_all_versions_by_stable_id
                  ('ENSG00000148944');
  Description : Similar to fetch_by_stable_id, but retrieves all versions of a
                gene stored in the database.
  Returntype  : listref of Bio::EnsEMBL::Gene
  Exceptions  : if we cant get the gene in given coord system
  Caller      : general
  Status      : At Risk

=cut

sub fetch_all_versions_by_stable_id {
  my ($self, $stable_id) = @_;

  my $constraint = "g.stable_id = ?";
  $self->bind_param_generic_fetch($stable_id, SQL_VARCHAR);
  return $self->generic_fetch($constraint);
}

=head2 fetch_by_exon_stable_id

  Arg [1]    : String $id
               The stable id of an exon of the gene to retrieve
  Example    : $gene = $gene_adptr->fetch_by_exon_stable_id('ENSE00000148944');
  Description: Retrieves a gene object from the database via an exon stable id.
               The gene will be retrieved in its native coordinate system (i.e.
               in the coordinate system it is stored in the database). It may
               be converted to a different coordinate system through a call to
               transform() or transfer(). If the gene or exon is not found
               undef is returned instead.
  Returntype : Bio::EnsEMBL::Gene or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_exon_stable_id {
  my ($self, $stable_id, $version) = @_;

  my $sql = qq(
      SELECT t.gene_id
        FROM transcript as t,
             exon_transcript as et,
             exon as e
       WHERE t.transcript_id = et.transcript_id 
         AND et.exon_id = e.exon_id
         AND e.stable_id = ?
         AND e.is_current = 1
  );

  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $stable_id, SQL_VARCHAR);
  $sth->execute();

  my ($dbID) = $sth->fetchrow_array();

  return undef if (!defined($dbID));

  my $gene = $self->fetch_by_dbID($dbID);

  return $gene;
} ## end sub fetch_by_exon_stable_id

=head2 fetch_all_by_domain

  Arg [1]    : String $domain
               The domain to fetch genes from
  Example    : my @genes = @{ $gene_adaptor->fetch_all_by_domain($domain) };
  Description: Retrieves a listref of genes whose translation contain interpro
               domain $domain. The genes are returned in their native coord
               system (i.e. the coord_system they are stored in). If the coord
               system needs to be changed, then tranform or transfer should be
               called on the individual objects returned.
  Returntype : list of Bio::EnsEMBL::Genes
  Exceptions : none
  Caller     : domainview
  Status     : Stable

=cut

sub fetch_all_by_domain {
  my ($self, $domain) = @_;

  throw("domain argument is required") unless ($domain);

  my $sth = $self->prepare(
  qq(
  SELECT    tr.gene_id
  FROM      interpro i,
            protein_feature pf,
            transcript tr,
            translation tl,
            seq_region sr,
            coord_system cs
  WHERE     cs.species_id = ?
    AND     cs.coord_system_id = sr.coord_system_id
    AND     sr.seq_region_id = tr.seq_region_id
    AND     tr.is_current = 1
    AND     tr.transcript_id = tl.transcript_id
    AND     tl.translation_id = pf.translation_id
    AND     pf.hit_name = i.id
    AND     i.interpro_ac = ?
  GROUP BY  tr.gene_id));

  $sth->bind_param(1, $self->species_id(), SQL_VARCHAR);
  $sth->bind_param(2, $domain,             SQL_VARCHAR);

  $sth->execute();

  my @array = @{$sth->fetchall_arrayref()};
  $sth->finish();

  my @gene_ids = map { $_->[0] } @array;

  return $self->fetch_all_by_dbID_list(\@gene_ids);
} ## end sub fetch_all_by_domain

=head2 fetch_all_by_Slice_and_external_dbname_link

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to fetch genes on.
  Arg [2]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Arg [3]    : (optional) boolean $load_transcripts
               if true, transcripts will be loaded immediately
               rather than lazy loaded later.
  Arg [4]    : String
               Name of the external database to fetch the Genes by
  Example    : @genes = @{
                 $ga->fetch_all_by_Slice_and_external_dbname_link(
                                          $slice, undef, undef, "HGNC" ) };
  Description: Overrides superclass method to optionally load
               transcripts immediately rather than lazy-loading them
               later.  This is more efficient when there are a lot
               of genes whose transcripts are going to be used. The
               genes are then filtered to return only those with
               external database links of the type specified
  Returntype : reference to list of genes
  Exceptions : thrown if exon cannot be placed on transcript slice
  Caller     : 
  Status     : Stable

=cut

sub fetch_all_by_Slice_and_external_dbname_link {
  my ($self, $slice, $logic_name, $load_transcripts, $db_name) = @_;

  # Get the external_db_id(s) from the name.
  my $dbentry_adaptor = $self->db()->get_DBEntryAdaptor();
  my $external_db_ids = $dbentry_adaptor->get_external_db_ids($db_name, undef, 'ignore release');

  if (scalar(@{$external_db_ids}) == 0) {
    my $external_db_names = $dbentry_adaptor->get_distinct_external_dbs();
    my $available = join("\n", map { "\t${_}"} @{$external_db_names});
    warning sprintf("Could not find external database " . "'%s' in the external_db table\n" . "Available are:\n%s", $db_name, $available);
    return [];
  }

  # Get the gene_ids for those with links.
  my %linked_genes;

  foreach my $local_external_db_id (@{$external_db_ids}) {
    my @linked_genes = $dbentry_adaptor->list_gene_ids_by_external_db_id($local_external_db_id);
    $linked_genes{$_} = 1 for @linked_genes;
  }
  
  # Get all the genes on the slice and filter by the gene ids list
  my $genes = $self->fetch_all_by_Slice($slice, $logic_name, $load_transcripts);
  my $genes_passed = [ grep { exists $linked_genes{$_->dbID()} } @{$genes} ];
  return $genes_passed;
} ## end sub fetch_all_by_Slice_and_external_dbname_link

=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to fetch genes on.
  Arg [2]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Arg [3]    : (optional) boolean $load_transcripts
               if true, transcripts will be loaded immediately rather than
               lazy loaded later.
  Arg [4]    : (optional) string $source
               the source name of the features to obtain.
  Arg [5]    : (optional) string biotype
                the biotype of the features to obtain.
  Example    : @genes = @{$gene_adaptor->fetch_all_by_Slice()};
  Description: Overrides superclass method to optionally load transcripts
               immediately rather than lazy-loading them later.  This
               is more efficient when there are a lot of genes whose
               transcripts are going to be used.
  Returntype : reference to list of genes 
  Exceptions : thrown if exon cannot be placed on transcript slice
  Caller     : Slice::get_all_Genes
  Status     : Stable

=cut

sub fetch_all_by_Slice {
  my ($self, $slice, $logic_name, $load_transcripts, $source, $biotype) = @_;

  my $constraint = 'g.is_current = 1';

  if (defined($source)) {
    $constraint .= " and g.source = '$source'";
  }
  if (defined($biotype)) {
    my $inline_variables = 1;
    $constraint .= " and ".$self->generate_in_constraint($biotype, 'g.biotype', SQL_VARCHAR, $inline_variables);
  }

  my $genes = $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint, $logic_name);

  # If there are 0 genes, still do lazy-loading.
  if (!$load_transcripts || @$genes < 1) {
    return $genes;
  }

  # Preload all of the transcripts now, instead of lazy loading later,
  # faster than one query per transcript.

  # First check if transcripts are already preloaded.
  # FIXME: Should check all transcripts.
  if (exists($genes->[0]->{'_transcript_array'})) {
    return $genes;
  }

  # Get extent of region spanned by transcripts.
  my ($min_start, $max_end);
  foreach my $g (@$genes) {
    if (!defined($min_start) || $g->seq_region_start() < $min_start) {
      $min_start = $g->seq_region_start();
    }
    if (!defined($max_end) || $g->seq_region_end() > $max_end) {
      $max_end = $g->seq_region_end();
    }
  }

  my $ext_slice;

  if ($min_start >= $slice->start() && $max_end <= $slice->end()) {
    $ext_slice = $slice;
  } else {
    my $sa = $self->db()->get_SliceAdaptor();
    $ext_slice = $sa->fetch_by_region($slice->coord_system->name(), $slice->seq_region_name(), $min_start, $max_end, $slice->strand(), $slice->coord_system->version());
  }

  # Associate transcript identifiers with genes.

  my %g_hash = map { $_->dbID => $_ } @{$genes};

  my $g_id_str = join(',', keys(%g_hash));

  my $sth = $self->prepare("SELECT gene_id, transcript_id " . "FROM   transcript " . "WHERE  gene_id IN ($g_id_str)");

  $sth->execute();

  my ($g_id, $tr_id);
  $sth->bind_columns(\($g_id, $tr_id));

  my %tr_g_hash;

  while ($sth->fetch()) {
    $tr_g_hash{$tr_id} = $g_hash{$g_id};
  }

  my $ta = $self->db()->get_TranscriptAdaptor();
  my $transcripts = $ta->fetch_all_by_Slice($ext_slice, 1, undef, sprintf("t.transcript_id IN (%s)", join(',', sort { $a <=> $b } keys(%tr_g_hash))));

  # Move transcripts onto gene slice, and add them to genes.
  foreach my $tr (@{$transcripts}) {
    if (!exists($tr_g_hash{$tr->dbID()})) { next }

    my $new_tr;
    if ($slice != $ext_slice) {
      $new_tr = $tr->transfer($slice);
      if (!defined($new_tr)) {
        throw("Unexpected. " . "Transcript could not be transfered onto Gene slice.");
      }
    } else {
      $new_tr = $tr;
    }

    $tr_g_hash{$tr->dbID()}->add_Transcript($new_tr);
  }

  return $genes;
} ## end sub fetch_all_by_Slice

=head2 count_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to count genes on.
  Arg [2]    : (optional) biotype(s) string or arrayref of strings 
                the biotype of the features to count.
  Arg [1]    : (optional) string $source
               the source name of the features to count.
  Example    : $cnt = $gene_adaptor->count_all_by_Slice();
  Description: Method to count genes on a given slice, filtering by biotype and source
  Returntype : integer
  Exceptions : thrown if exon cannot be placed on transcript slice
  Status     : Stable
  Caller     : general
=cut

sub count_all_by_Slice {
  my ($self, $slice, $biotype, $source) = @_;

  my $constraint = 'g.is_current = 1';
  if (defined($source)) {
    $constraint .= " and g.source = '$source'";
  }
  if (defined($biotype)) {
    $constraint .= " and " . $self->biotype_constraint($biotype);
  }

  return $self->count_by_Slice_constraint($slice, $constraint);
}

=head2 fetch_by_transcript_id

  Arg [1]    : Int $trans_id
               Unique database identifier for the transcript whose gene should
               be retrieved. The gene is returned in its native coord
               system (i.e. the coord_system it is stored in). If the coord
               system needs to be changed, then tranform or transfer should
               be called on the returned object. undef is returned if the
               gene or transcript is not found in the database.
  Example    : $gene = $gene_adaptor->fetch_by_transcript_id(1241);
  Description: Retrieves a gene from the database via the database identifier
               of one of its transcripts.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_transcript_id {
  my ($self, $trans_id) = @_;

  # this is a cheap SQL call
  my $sth = $self->prepare(
  qq(
      SELECT tr.gene_id
      FROM transcript tr
      WHERE tr.transcript_id = ?
  ));

  $sth->bind_param(1, $trans_id, SQL_INTEGER);
  $sth->execute();

  my ($geneid) = $sth->fetchrow_array();

  $sth->finish();

  return undef if (!defined $geneid);

  my $gene = $self->fetch_by_dbID($geneid);
  return $gene;
}

=head2 fetch_by_transcript_stable_id

  Arg [1]    : string $trans_stable_id
               transcript stable ID whose gene should be retrieved
  Example    : my $gene = $gene_adaptor->fetch_by_transcript_stable_id
                 ('ENST0000234');
  Description: Retrieves a gene from the database via the stable ID of one of
               its transcripts
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_transcript_stable_id {
  my ($self, $trans_stable_id) = @_;

  my $sth = $self->prepare(
  qq(
        SELECT  gene_id
  FROM	transcript
        WHERE   stable_id = ?
        AND     is_current = 1
    ));

  $sth->bind_param(1, $trans_stable_id, SQL_VARCHAR);
  $sth->execute();

  my ($geneid) = $sth->fetchrow_array();
  $sth->finish;

  return undef if (!defined $geneid);

  my $gene = $self->fetch_by_dbID($geneid);
  return $gene;
}

=head2 fetch_by_translation_stable_id

  Arg [1]    : String $translation_stable_id
               The stable id of a translation of the gene to be obtained
  Example    : my $gene = $gene_adaptor->fetch_by_translation_stable_id
                 ('ENSP00000278194');
  Description: Retrieves a gene via the stable id of one of its translations.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_translation_stable_id {
  my ($self, $translation_stable_id) = @_;

  my $sth = $self->prepare(
  qq(
        SELECT  tr.gene_id
  FROM    transcript tr,
                translation tl
  WHERE   tl.stable_id = ?
        AND     tr.transcript_id = tl.transcript_id
        AND     tr.is_current = 1
    ));

  $sth->bind_param(1, $translation_stable_id, SQL_VARCHAR);
  $sth->execute();

  my ($geneid) = $sth->fetchrow_array();
  $sth->finish;
  if (!defined $geneid) {
    return undef;
  }
  return $self->fetch_by_dbID($geneid);
}

=head2 fetch_all_by_external_name

  Arg [1]    : String $external_name
               The external identifier for the gene to be obtained
  Arg [2]    : (optional) String $external_db_name
               The name of the external database from which the
               identifier originates.
  Arg [3]    : Boolean override. Force SQL regex matching for users
               who really do want to find all 'NM%'
  Example    : @genes = @{$gene_adaptor->fetch_all_by_external_name('BRCA2')}
               @many_genes = @{$gene_adaptor->fetch_all_by_external_name('BRCA%')}
  Description: Retrieves a list of genes with an external database
               identifier $external_name. The genes returned are in
               their native coordinate system, i.e. in the coordinate
               system they are stored in the database in.  If another
               coordinate system is required then the Gene::transfer or
               Gene::transform method can be used.
               SQL wildcards % and _ are supported in the $external_name,
               but their use is somewhat restricted for performance reasons.
               Users that really do want % and _ in the first three characters
               should use argument 3 to prevent optimisations
  Returntype : listref of Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : goview, general
  Status     : Stable

=cut

sub fetch_all_by_external_name {
  my ($self, $external_name, $external_db_name, $override) = @_;

  my $entryAdaptor = $self->db->get_DBEntryAdaptor();

  my @ids = $entryAdaptor->list_gene_ids_by_extids($external_name, $external_db_name, $override);

  my %genes_by_dbIDs = map { $_->dbID(), $_ } @{$self->fetch_all_by_dbID_list(\@ids)};

  my @features = map { $genes_by_dbIDs{$_} } @ids;
  my @reference = grep { $_->slice()->is_reference() } @features;
  my @non_reference = grep { ! $_->slice()->is_reference() } @features;
  return [ @reference, @non_reference ];
}

=head2 fetch_all_by_description

  Arg [1]    : String of description
  Example    : $gene_list = $gene_adaptor->fetch_all_by_description('RNA%');
  Description: Fetches genes by their textual description. Fully supports SQL
               wildcards, since getting an exact hit is unlikely.
  Returntype : listref of Bio::EnsEMBL::Gene

=cut

sub fetch_all_by_description {
    my ($self,$description) = @_;
    
    my $constraint = "g.description LIKE ?";
    $self->bind_param_generic_fetch($description, SQL_VARCHAR);
    return $self->generic_fetch($constraint);
}

=head2 fetch_all_by_GOTerm

  Arg [1]   : Bio::EnsEMBL::OntologyTerm
              The GO term for which genes should be fetched.

  Example:  @genes = @{
              $gene_adaptor->fetch_all_by_GOTerm(
                $go_adaptor->fetch_by_accession('GO:0030326') ) };

  Description   : Retrieves a list of genes that are associated with
                  the given GO term, or with any of its descendent
                  GO terms.  The genes returned are in their native
                  coordinate system, i.e. in the coordinate system
                  in which they are stored in the database.  If
                  another coordinate system is required then the
                  Gene::transfer or Gene::transform method can be
                  used.

  Return type   : listref of Bio::EnsEMBL::Gene
  Exceptions    : Throws of argument is not a GO term
  Caller        : general
  Status        : Stable

=cut

sub fetch_all_by_GOTerm {
  my ($self, $term) = @_;

  assert_ref($term, 'Bio::EnsEMBL::OntologyTerm');
  if ($term->ontology() ne 'GO') {
    throw('Argument is not a GO term');
  }

  my $entryAdaptor = $self->db->get_DBEntryAdaptor();

  my %unique_dbIDs;
  foreach my $accession (map { $_->accession() } ($term, @{$term->descendants()})) {
    my @ids = $entryAdaptor->list_gene_ids_by_extids($accession, 'GO');
    foreach my $dbID (@ids) { $unique_dbIDs{$dbID} = 1 }
  }

  my @result = @{$self->fetch_all_by_dbID_list([sort { $a <=> $b } keys(%unique_dbIDs)])};

  return \@result;
}

=head2 fetch_all_by_ontology_linkage_type

  Arg [1]   : (optional) string $db_name
              The database name to search for. Defaults to GO
  Arg [2]   : string $linkage_type
              Linkage type to search for e.g. IMP

  Example:    my $genes = $gene_adaptor->fetch_all_by_ontology_linkage_type('GO', 'IMP');
              my $genes = $gene_adaptor->fetch_all_by_ontology_linkage_type(undef, 'IMP');

  Description   : Retrieves a list of genes that are associated with
                  the given ontology linkage type.  The genes returned 
                  are in their native coordinate system, i.e. in the 
                  coordinate system in which they are stored in the database.
  Return type   : listref of Bio::EnsEMBL::Gene
  Exceptions    : Throws if a linkage type is not given
  Caller        : general
  Status        : Stable

=cut

sub fetch_all_by_ontology_linkage_type {
  my ($self, $db_name, $linkage_type) = @_;
  $db_name = 'GO' if ! defined $db_name;
  throw "No linkage type given" if ! defined $linkage_type;

  my $dbentry_adaptor = $self->db->get_DBEntryAdaptor();
  my $external_db_ids = $dbentry_adaptor->get_external_db_ids($db_name, undef, 'ignore release');
  if (scalar(@{$external_db_ids}) == 0) {
    warning sprintf("Could not find external database '%s' in the external_db table", $db_name);
    return [];
  }

  # Get the gene_ids for those with links.
  my %unique_dbIDs;
  foreach my $local_external_db_id (@{$external_db_ids}) {
    my @gene_ids = $dbentry_adaptor->list_gene_ids_by_external_db_id($local_external_db_id, $linkage_type);
    $unique_dbIDs{$_} = 1 for @gene_ids;
  }

  # Get all the genes and return
  return $self->fetch_all_by_dbID_list([keys %unique_dbIDs]);
}

=head2 fetch_all_by_GOTerm_accession

  Arg [1]   : String
              The GO term accession for which genes should be
              fetched.

  Example   :

    @genes =
      @{ $gene_adaptor->fetch_all_by_GOTerm_accession(
        'GO:0030326') };

  Description   : Retrieves a list of genes that are associated with
                  the given GO term, or with any of its descendent
                  GO terms.  The genes returned are in their native
                  coordinate system, i.e. in the coordinate system
                  in which they are stored in the database.  If
                  another coordinate system is required then the
                  Gene::transfer or Gene::transform method can be
                  used.

  Return type   : listref of Bio::EnsEMBL::Gene
  Exceptions    : Throws of argument is not a GO term accession
  Caller        : general
  Status        : Stable

=cut

sub fetch_all_by_GOTerm_accession {
  my ($self, $accession) = @_;

  if ($accession !~ /^GO:/) {
    throw('Argument is not a GO term accession');
  }

  my $goAdaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'Ontology', 'OntologyTerm');

  my $term = $goAdaptor->fetch_by_accession($accession);

  return $self->fetch_all_by_GOTerm($term);
}

=head2 fetch_all_alt_alleles

  Arg [1]    : Bio::EnsEMBL::Gene $gene
               The gene to fetch alternative alleles for
  Arg [2]    : Boolean (optional)
               Ask the method to warn about any gene without an alt allele 
               group. Defaults to false
  Example    : my @alt_genes = @{ $gene_adaptor->fetch_all_alt_alleles($gene) };
               foreach my $alt_gene (@alt_genes) {
                 print "Alternate allele: " . $alt_gene->stable_id() . "\n" ;
               }
  Description: Retrieves genes which are alternate alleles to a provided gene.
               Alternate alleles in Ensembl are genes which are similar and are
               on an alternative haplotype of the same region. There are not 
               currently very many of these. This method will return a 
               reference to an empty list if no alternative alleles are found.
  Returntype : ArrayRef of Bio::EnsEMBL::Gene objects
  Exceptions : throw if incorrect arg provided
               warning if gene arg does not have an entry in an alt allele and if
               the warn flag is true
  Caller     : Gene::get_all_alt_alleles
  Status     : Stable

=cut

sub fetch_all_alt_alleles {
  my $self = shift;
  my $gene = shift;
  my $warn = shift;

  if (!ref($gene) || !$gene->isa('Bio::EnsEMBL::Gene')) {
    throw('Bio::EnsEMBL::Gene argument is required');
  }

  my $gene_id = $gene->dbID();

  if (!$gene_id) {
    warning('Cannot retrieve alternate alleles for gene without dbID');
    return [];
  }

  my $aaga = $self->db->get_adaptor('AltAlleleGroup');
  my $aag = $aaga->fetch_by_gene_id($gene->dbID);
  unless ($aag) {
    if ($warn) {
      warning("Supplied gene has no alternative alleles"); 
    }
    return [];
  }
  # query for all alternative genes. do not filter 
  # the representative but do filter this gene out
  return $aag->get_all_Genes(undef, [$gene]);
} ## end sub fetch_all_alt_alleles

=head2 is_ref

  Arg [1]    : Gene dbID
  Description: Used to determine whether a given Gene is the representative 
               Gene of an alt allele group. If it does not have an alternative
               allele that is more representative, then this ID will be said to
               be representative.
  Returntype : Boolean - True for yes or no alternatives  

=cut

sub is_ref {
  my ($self, $gene_id) = @_;
  my $aag = $self->db->get_adaptor('AltAlleleGroup')->fetch_by_gene_id($gene_id);
  if (defined($aag)) {
      if ($aag->rep_Gene_id == $gene_id) {
          return 1;
      } else {
          return 0;
      }
  } else {
      return 1;
  }
  throw("Unhandled circumstance in GeneAdaptor->is_ref");
}

=head2 store_alt_alleles


  Arg [1]    : reference to list of Bio::EnsEMBL::Genes $genes
  Example    : $gene_adaptor->store_alt_alleles([$gene1, $gene2, $gene3]);
  Description: DEPRECATED. Switch to using AltAlleleGroup and the 
               AltAlleleGroupAdaptor which supports more complex queries

               This method creates a group of alternative alleles (i.e. locus)
               from a set of genes. The genes should be genes from alternate
               haplotypes which are similar. The genes must already be stored
               in this database. WARNING - now that more fine-grained support
               for alt_alleles has been implemented, this method is rather coarse.
               Consider working directly with AltAlleleGroup and 
               AltAlleleGroupAdaptor.
  Returntype : int alt_allele_group_id or undef if no alt_alleles were stored
  Exceptions : throw on incorrect arguments
               throw on sql error (e.g. duplicate unique id)
  Caller     : general
  Status     : Stable

=cut

sub store_alt_alleles {
  my $self  = shift;
  my $genes = shift;

  warning "Unsupported. Switch to using AltAlleleGroupAdaptor::store() and AltAlleleGroups";

  if (!ref($genes) eq 'ARRAY') {
    throw('List reference of Bio::EnsEMBL::Gene argument expected.');
  }
  my @genes     = @$genes;
  my $num_genes = scalar(@genes);
  if ($num_genes < 2) {
    warning('At least 2 genes must be provided to construct alternative alleles (gene id: ' . $genes[0]->dbID() . '). Ignoring.');
    return;
  }

  my $allele_list;
  foreach my $gene (@$genes) {
      my $aa_record = [];
      push @$aa_record, $gene->dbID;
      my %type = {};
      if ($gene->slice->is_reference()) {
          $type{'IS_REPRESENTATIVE'} = 1;
      }
      push @$aa_record, \%type; 
      push @$allele_list, $aa_record;
  }

  my $aag = Bio::EnsEMBL::AltAlleleGroup->new(
    -MEMBERS => $allele_list,
  );
  if (scalar( @{$aag->get_all_members_with_type('IS_REPRESENTATIVE')} ) != 1) {
    warning('Inappropriate number of alternative alleles on the reference sequence. Ignoring.');
    return;
  }
  
  my $aaga = $self->db->get_adaptor('AltAlleleGroup');
  return $aaga->store($aag);
} ## end sub store_alt_alleles

=head2 store

  Arg [1]    : Bio::EnsEMBL::Gene $gene
               The gene to store in the database
  Arg [2]    : ignore_release in xrefs [default 1] set to 0 to use release info 
               in external database references
  Arg [3]    : prevent coordinate recalculation if you are persisting 
               transcripts with this gene
  Arg [4]    : prevent copying supporting features across exons
               increased speed for lost accuracy
  Example    : $gene_adaptor->store($gene);
  Description: Stores a gene in the database.
  Returntype : the database identifier (dbID) of the newly stored gene
  Exceptions : thrown if the $gene is not a Bio::EnsEMBL::Gene or if 
               $gene does not have an analysis object
  Caller     : general
  Status     : Stable

=cut

sub store {
  my ($self, $gene, $ignore_release, $skip_recalculating_coordinates, $skip_exon_sf) = @_;

  if (!ref $gene || !$gene->isa('Bio::EnsEMBL::Gene')) {
    throw("Must store a gene object, not a $gene");
  }
  if (!defined($ignore_release)) {
    $ignore_release = 1;
  }
  my $db = $self->db();

  if ($gene->is_stored($db)) {
    return $gene->dbID();
  }

  # ensure coords are correct before storing
  $gene->recalculate_coordinates();

  my $analysis = $gene->analysis();
  throw("Genes must have an analysis object.") if (!defined($analysis));

  my $analysis_id;
  if ($analysis->is_stored($db)) {
    $analysis_id = $analysis->dbID();
  } else {
    $analysis_id = $db->get_AnalysisAdaptor->store($analysis);
  }

  my $type = $gene->biotype || "";

  # default to is_current = 1 if this attribute is not set
  my $is_current = $gene->is_current;
  $is_current = 1 unless (defined($is_current));

  my $original             = $gene;
  my $original_transcripts = $gene->get_all_Transcripts();

  my $seq_region_id;

  ($gene, $seq_region_id) = $self->_pre_store($gene);

  my @columns = qw(
               biotype
               analysis_id
               seq_region_id
               seq_region_start
               seq_region_end
               seq_region_strand
               description
               source
               is_current
               canonical_transcript_id
  );

  my @canned_columns;
  my @canned_values;

  if (defined($gene->stable_id)) {
      push @columns, 'stable_id', 'version';

      my $created  = $self->db->dbc->from_seconds_to_date($gene->created_date());
      my $modified = $self->db->dbc->from_seconds_to_date($gene->modified_date());

      if ($created) {
	push @canned_columns, 'created_date';
	push @canned_values,  $created;
      }
      if ($modified) {
	push @canned_columns, 'modified_date';
	push @canned_values,  $modified;
      }
      
  }

  my $columns = join(', ', @columns, @canned_columns);
  my $values  = join(', ', ('?') x @columns, @canned_values);
  my $store_gene_sql = qq(
        INSERT INTO gene ( $columns ) VALUES ( $values )
  );

  my $sth = $self->prepare($store_gene_sql);
  $sth->bind_param(1,  $type,                SQL_VARCHAR);
  $sth->bind_param(2,  $analysis_id,         SQL_INTEGER);
  $sth->bind_param(3,  $seq_region_id,       SQL_INTEGER);
  $sth->bind_param(4,  $gene->start(),       SQL_INTEGER);
  $sth->bind_param(5,  $gene->end(),         SQL_INTEGER);
  $sth->bind_param(6,  $gene->strand(),      SQL_TINYINT);
  $sth->bind_param(7,  $gene->description(), SQL_LONGVARCHAR);
  $sth->bind_param(8,  $gene->source(),      SQL_VARCHAR);
  $sth->bind_param(9,  $is_current,          SQL_TINYINT);

  # Canonical transcript ID will be updated later.
  # Set it to zero for now.
  $sth->bind_param(10, 0, SQL_TINYINT);


  if (defined($gene->stable_id)) {

    $sth->bind_param(11, $gene->stable_id, SQL_VARCHAR);
    $sth->bind_param(12, $gene->version,   SQL_INTEGER);
  }

  $sth->execute();
  $sth->finish();

  my $gene_dbID = $self->last_insert_id('gene_id', undef, 'gene');

  # store the dbentries associated with this gene
  my $dbEntryAdaptor = $db->get_DBEntryAdaptor();

  foreach my $dbe (@{$gene->get_all_DBEntries}) {
    $dbEntryAdaptor->store($dbe, $gene_dbID, "Gene", $ignore_release);
  }

  # We allow transcripts not to share equal exons and instead have
  # copies.  For the database we still want sharing though, to have
  # easier time with stable ids. So we need to have a step to merge
  # exons together before store.
  my %exons;

  foreach my $trans (@{$gene->get_all_Transcripts}) {
    foreach my $e (@{$trans->get_all_Exons}) {
      my $key = $e->hashkey();
      if (exists $exons{$key}) {
        $trans->swap_exons($e, $exons{$key}, $skip_exon_sf);
      } else {
        $exons{$key} = $e;
      }
    }
  }

  my $transcript_adaptor = $db->get_TranscriptAdaptor();

  my $transcripts = $gene->get_all_Transcripts();

  my $new_canonical_transcript_id;
  for (my $i = 0; $i < @$transcripts; $i++) {
    my $new = $transcripts->[$i];
    my $old = $original_transcripts->[$i];

    $transcript_adaptor->store($new, $gene_dbID, $analysis_id, $skip_recalculating_coordinates);

    if (!defined($new_canonical_transcript_id) && $new->is_canonical()) {
      $new_canonical_transcript_id = $new->dbID();
    }

  # update the original transcripts since we may have made copies of
  # them by transforming the gene
    $old->dbID($new->dbID());
    $old->adaptor($new->adaptor());

    if ($new->translation) {
      $old->translation->dbID($new->translation()->dbID);
      $old->translation->adaptor($new->translation()->adaptor);
    }
  }

  if (defined($new_canonical_transcript_id)) {
  # Now the canonical transcript has been stored, so update the
  # canonical_transcript_id of this gene with the new dbID.
    my $sth = $self->prepare(
      q(
        UPDATE gene
        SET canonical_transcript_id = ?
        WHERE gene_id = ?)
    );

    $sth->bind_param(1, $new_canonical_transcript_id, SQL_INTEGER);
    $sth->bind_param(2, $gene_dbID, SQL_INTEGER);

    $sth->execute();
    $sth->finish();
  }

  # update gene to point to display xref if it is set
  if (my $display_xref = $gene->display_xref) {
    my $dxref_id;
    if ($display_xref->is_stored($db)) {
      $dxref_id = $display_xref->dbID();
    } else {
      $dxref_id = $dbEntryAdaptor->exists($display_xref);
    }

  if (defined($dxref_id)) {
    my $sth = $self->prepare("UPDATE gene SET display_xref_id = ? WHERE gene_id = ?");
      $sth->bind_param(1, $dxref_id,  SQL_INTEGER);
      $sth->bind_param(2, $gene_dbID, SQL_INTEGER);
      $sth->execute();
      $sth->finish();
      $display_xref->dbID($dxref_id);
      $display_xref->adaptor($dbEntryAdaptor);
      $display_xref->dbID($dxref_id);
      $display_xref->adaptor($dbEntryAdaptor);
    } else {
      warning("Display_xref " . $display_xref->dbname() . ":" . $display_xref->display_id() . " is not stored in database.\n" . "Not storing relationship to this gene.");
      $display_xref->dbID(undef);
      $display_xref->adaptor(undef);
    }
  }

  # store gene attributes if there are any
  my $attr_adaptor = $db->get_AttributeAdaptor();
  $attr_adaptor->store_on_Gene($gene_dbID, $gene->get_all_Attributes);

  # set the adaptor and dbID on the original passed in gene not the
  # transfered copy
  $original->adaptor($self);
  $original->dbID($gene_dbID);

  return $gene_dbID;
} ## end sub store

=head2 remove

  Arg [1]    : Bio::EnsEMBL::Gene $gene
               the gene to remove from the database
  Example    : $gene_adaptor->remove($gene);
  Description: Removes a gene completely from the database. All associated
               transcripts, exons, stable_identifiers, descriptions, etc.
               are removed as well. Use with caution!
  Returntype : none
  Exceptions : throw on incorrect arguments 
               warning if gene is not stored in this database
  Caller     : general
  Status     : Stable

=cut

sub remove {
  my $self = shift;
  my $gene = shift;

  if (!ref($gene) || !$gene->isa('Bio::EnsEMBL::Gene')) {
    throw("Bio::EnsEMBL::Gene argument expected.");
  }

  if (!$gene->is_stored($self->db())) {
    warning("Cannot remove gene " . $gene->dbID() . ". Is not stored in " . "this database.");
    return;
  }

  # remove all object xrefs associated with this gene

  my $dbe_adaptor = $self->db()->get_DBEntryAdaptor();
  foreach my $dbe (@{$gene->get_all_DBEntries()}) {
    $dbe_adaptor->remove_from_object($dbe, $gene, 'Gene');
  }

  # remove all alternative allele entries associated with this gene
  my $sth = $self->prepare("DELETE FROM alt_allele WHERE gene_id = ?");
  $sth->bind_param(1, $gene->dbID, SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  # remove the attributes associated with this transcript
  my $attrib_adaptor = $self->db->get_AttributeAdaptor;
  $attrib_adaptor->remove_from_Gene($gene);

  # remove all of the transcripts associated with this gene
  my $transcriptAdaptor = $self->db->get_TranscriptAdaptor();
  foreach my $trans (@{$gene->get_all_Transcripts()}) {
    $transcriptAdaptor->remove($trans);
  }

  # remove this gene from the database

  $sth = $self->prepare("DELETE FROM gene WHERE gene_id = ? ");
  $sth->bind_param(1, $gene->dbID, SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  # unset the gene identifier and adaptor thereby flagging it as unstored

  $gene->dbID(undef);
  $gene->adaptor(undef);

  return;
} ## end sub remove

=head2 get_Interpro_by_geneid

  Arg [1]    : String $gene_stable_id
               The stable ID of the gene to obtain
  Example    : @i = @{
                  $gene_adaptor->get_Interpro_by_geneid(
                    $gene->stable_id() ) };
  Description: Gets interpro accession numbers by gene stable id. A hack really
               - we should have a much more structured system than this.
  Returntype : listref of strings (Interpro_acc:description)
  Exceptions : none 
  Caller     : domainview
  Status     : Stable

=cut

sub get_Interpro_by_geneid {
  my ($self, $gene_stable_id) = @_;

  my $sql = qq(
  SELECT    i.interpro_ac,
            x.description
  FROM      transcript t,
            translation tl,
            protein_feature pf,
            interpro i,
            xref x,
            gene g
  WHERE     g.stable_id = ?
    AND     t.gene_id = g.gene_id
    AND     t.is_current = 1
    AND     tl.transcript_id = t.transcript_id
    AND     tl.translation_id = pf.translation_id
    AND     i.id = pf.hit_name
    AND     i.interpro_ac = x.dbprimary_acc);

  my $sth = $self->prepare($sql);

  $sth->bind_param(1, $gene_stable_id, SQL_VARCHAR);

  $sth->execute;

  my @out;
  my %h;
  while ((my $arr = $sth->fetchrow_arrayref())) {
    if ($h{$arr->[0]}) { next; }
    $h{$arr->[0]} = 1;
    my $string = $arr->[0] . ":" . $arr->[1];
    push(@out, $string);
  }

  return \@out;
} ## end sub get_Interpro_by_geneid

=head2 update

  Arg [1]    : Bio::EnsEMBL::Gene $gene
               The gene to update
  Example    : $gene_adaptor->update($gene);
  Description: Updates the type, analysis, display_xref, is_current and
               description of a gene in the database.
  Returntype : None
  Exceptions : thrown if the $gene is not a Bio::EnsEMBL::Gene
  Caller     : general
  Status     : Stable

=cut

sub update {
  my ($self, $gene) = @_;
  my $update = 0;

  if (!defined $gene || !ref $gene || !$gene->isa('Bio::EnsEMBL::Gene')) {
    throw("Must update a gene object, not a $gene");
  }

  my $update_gene_sql = qq(
       UPDATE gene
          SET biotype = ?,
              analysis_id = ?,
              display_xref_id = ?,
              description = ?,
              is_current = ?,
              canonical_transcript_id = ?
        WHERE gene_id = ?
  );

  my $display_xref = $gene->display_xref();
  my $display_xref_id;

  if ($display_xref && $display_xref->dbID()) {
    $display_xref_id = $display_xref->dbID();
  } else {
    $display_xref_id = undef;
  }

  my $sth = $self->prepare($update_gene_sql);

  $sth->bind_param(1, $gene->biotype(),        SQL_VARCHAR);
  $sth->bind_param(2, $gene->analysis->dbID(), SQL_INTEGER);
  $sth->bind_param(3, $display_xref_id,        SQL_INTEGER);
  $sth->bind_param(4, $gene->description(),    SQL_VARCHAR);
  $sth->bind_param(5, $gene->is_current(),     SQL_TINYINT);

  if (defined($gene->canonical_transcript())) {
    $sth->bind_param(6, $gene->canonical_transcript()->dbID(), SQL_INTEGER);
  } else {
    $sth->bind_param(6, 0, SQL_INTEGER);
  }

  $sth->bind_param(7, $gene->dbID(), SQL_INTEGER);

  $sth->execute();

  # maybe should update stable id ???
} ## end sub update


=head2 update_coords

  Arg [1]    : Bio::EnsEMBL::Gene $gene
               The gene to update
  Example    : $gene_adaptor->update_coords($gene);
  Description: In the event of a transcript being removed, coordinates for the Gene
               need to be reset, but update() does not do this. update_coords 
               fills this niche
  Returntype : None
  Exceptions : thrown if the $gene is not supplied
  Caller     : general

=cut

sub update_coords {
  my ($self, $gene) = @_;
  throw('Must have a gene to update in order to update it') unless ($gene);
  $gene->recalculate_coordinates;
  my $update_sql = qq(
    UPDATE gene
       SET seq_region_start = ?,
           seq_region_end = ?
       WHERE gene_id = ?
    );
  my $sth = $self->prepare($update_sql);
  $sth->bind_param(1, $gene->start);
  $sth->bind_param(2, $gene->end);
  $sth->bind_param(3, $gene->dbID);
  $sth->execute();
}

# _objs_from_sth

#  Arg [1]    : StatementHandle $sth
#  Arg [2]    : Bio::EnsEMBL::AssemblyMapper $mapper
#  Arg [3]    : Bio::EnsEMBL::Slice $dest_slice
#  Description: PROTECTED implementation of abstract superclass method.
#               responsible for the creation of Genes
#  Returntype : listref of Bio::EnsEMBL::Genes in target coordinate system
#  Exceptions : none
#  Caller     : internal
#  Status     : Stable

sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;

  #
  # This code is ugly because an attempt has been made to remove as many
  # function calls as possible for speed purposes.  Thus many caches and
  # a fair bit of gymnastics is used.
  #

  my $sa             = $self->db()->get_SliceAdaptor();
  my $aa             = $self->db()->get_AnalysisAdaptor();
  my $dbEntryAdaptor = $self->db()->get_DBEntryAdaptor();

  my @genes;
  my %analysis_hash;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;

  my (
    $gene_id,                 $seq_region_id,     $seq_region_start,
    $seq_region_end,          $seq_region_strand, $analysis_id,
    $biotype,                 $display_xref_id,   $gene_description,
    $source,                  $is_current,
    $canonical_transcript_id, $stable_id,         $version,
    $created_date,            $modified_date,     $xref_display_label,
    $xref_primary_acc,        $xref_description,  $xref_version,
    $external_db,             $external_status,   $external_release,
    $external_db_name,        $info_type,         $info_text
  );

  $sth->bind_columns(\(
                      $gene_id,                 $seq_region_id,     $seq_region_start,
                      $seq_region_end,          $seq_region_strand, $analysis_id,
                      $biotype,                 $display_xref_id,   $gene_description,
                      $source,                  $is_current,
                      $canonical_transcript_id, $stable_id,         $version,
                      $created_date,            $modified_date,     $xref_display_label,
                      $xref_primary_acc,        $xref_description,  $xref_version,
                      $external_db,             $external_status,   $external_release,
                      $external_db_name,        $info_type,         $info_text
                    ) );

  my $dest_slice_start;
  my $dest_slice_end;
  my $dest_slice_strand;
  my $dest_slice_length;
  my $dest_slice_cs;
  my $dest_slice_sr_name;
  my $dest_slice_sr_id;
  my $asma;

  if ($dest_slice) {
    $dest_slice_start   = $dest_slice->start();
    $dest_slice_end     = $dest_slice->end();
    $dest_slice_strand  = $dest_slice->strand();
    $dest_slice_length  = $dest_slice->length();
    $dest_slice_cs      = $dest_slice->coord_system();
    $dest_slice_sr_name = $dest_slice->seq_region_name();
    $dest_slice_sr_id   = $dest_slice->get_seq_region_id();
    $asma               = $self->db->get_AssemblyMapperAdaptor();
  }

  FEATURE: while($sth->fetch()) {

    #get the analysis object
    my $analysis = $analysis_hash{$analysis_id} ||= $aa->fetch_by_dbID($analysis_id);
    $analysis_hash{$analysis_id} = $analysis;

    #need to get the internal_seq_region, if present
    $seq_region_id = $self->get_seq_region_id_internal($seq_region_id);
    my $slice = $slice_hash{"ID:".$seq_region_id};

    if (!$slice) {
      $slice                            = $sa->fetch_by_seq_region_id($seq_region_id);
      $slice_hash{"ID:".$seq_region_id} = $slice;
      $sr_name_hash{$seq_region_id}     = $slice->seq_region_name();
      $sr_cs_hash{$seq_region_id}       = $slice->coord_system();
    }

    #obtain a mapper if none was defined, but a dest_seq_region was
    if(!$mapper && $dest_slice && !$dest_slice_cs->equals($slice->coord_system)) {
      $mapper = $asma->fetch_by_CoordSystems($dest_slice_cs, $slice->coord_system);
    }

    my $sr_name = $sr_name_hash{$seq_region_id};
    my $sr_cs   = $sr_cs_hash{$seq_region_id};

    #
    # remap the feature coordinates to another coord system
    # if a mapper was provided
    #

    if ($mapper) {

      if (defined $dest_slice && $mapper->isa('Bio::EnsEMBL::ChainedAssemblyMapper') ) {
        ($seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand) =
         $mapper->map($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand, $sr_cs, 1, $dest_slice);

      } else {
        ($seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand) =
         $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand, $sr_cs);
      }

      #skip features that map to gaps or coord system boundaries
      next FEATURE if (!defined($seq_region_id));

      #get a slice in the coord system we just mapped to
      $slice = $slice_hash{"ID:".$seq_region_id} ||= $sa->fetch_by_seq_region_id($seq_region_id);
    }

    #
    # If a destination slice was provided convert the coords.
    #
    if (defined($dest_slice)) {
      my $seq_region_len = $dest_slice->seq_region_length();

      if ( $dest_slice_strand == 1 ) {
        $seq_region_start = $seq_region_start - $dest_slice_start + 1;
        $seq_region_end   = $seq_region_end - $dest_slice_start + 1;

        if ( $dest_slice->is_circular ) {
        # Handle circular chromosomes.

          if ( $seq_region_start > $seq_region_end ) {
            # Looking at a feature overlapping the chromosome origin.

            if ( $seq_region_end > $dest_slice_start ) {
              # Looking at the region in the beginning of the chromosome
              $seq_region_start -= $seq_region_len;
            }
            if ( $seq_region_end < 0 ) {
              $seq_region_end += $seq_region_len;
            }
          } else {
            if ($dest_slice_start > $dest_slice_end && $seq_region_end < 0) {
              # Looking at the region overlapping the chromosome
              # origin and a feature which is at the beginning of the
              # chromosome.
              $seq_region_start += $seq_region_len;
              $seq_region_end   += $seq_region_len;
            }
          }
        }
      } else {

        my $start = $dest_slice_end - $seq_region_end + 1;
        my $end = $dest_slice_end - $seq_region_start + 1;

        if ($dest_slice->is_circular()) {

          if ($dest_slice_start > $dest_slice_end) {
            # slice spans origin or replication

            if ($seq_region_start >= $dest_slice_start) {
              $end += $seq_region_len;
              $start += $seq_region_len if $seq_region_end > $dest_slice_start;

            } elsif ($seq_region_start <= $dest_slice_end) {
              # do nothing
            } elsif ($seq_region_end >= $dest_slice_start) {
              $start += $seq_region_len;
              $end += $seq_region_len;

            } elsif ($seq_region_end <= $dest_slice_end) {
              $end += $seq_region_len if $end < 0;

            } elsif ($seq_region_start > $seq_region_end) {
              $end += $seq_region_len;
            }

          } else {

            if ($seq_region_start <= $dest_slice_end and $seq_region_end >= $dest_slice_start) {
              # do nothing
            } elsif ($seq_region_start > $seq_region_end) {
              if ($seq_region_start <= $dest_slice_end) {
                $start -= $seq_region_len;
              } elsif ($seq_region_end >= $dest_slice_start) {
                $end += $seq_region_len;
              }
            }
          }
        }

        $seq_region_start = $start;
        $seq_region_end = $end;
        $seq_region_strand *= -1;

      } ## end else [ if ( $dest_slice_strand...)]

      # Throw away features off the end of the requested slice or on
      # different seq_region.
      if ($seq_region_end < 1
          || $seq_region_start > $dest_slice_length
          || ($dest_slice_sr_id != $seq_region_id)) {
        next FEATURE;
      }
      $slice = $dest_slice;
    }

    my $display_xref;

    if ($display_xref_id) {
      $display_xref = Bio::EnsEMBL::DBEntry->new_fast({
        'dbID'            => $display_xref_id,
        'adaptor'         => $dbEntryAdaptor,
        'display_id'      => $xref_display_label,
        'primary_id'      => $xref_primary_acc,
        'version'         => $xref_version,
        'description'     => $xref_description,
        'release'         => $external_release,
        'dbname'          => $external_db,
        'db_display_name' => $external_db_name,
        'info_type'       => $info_type,
        'info_text'       => $info_text
      });
      $display_xref->status($external_status);
    }

    # Finally, create the new Gene.
    push(
      @genes,
      $self->_create_feature_fast(
      'Bio::EnsEMBL::Gene', {
       'analysis'                => $analysis,
       'biotype'                 => $biotype,
       'start'                   => $seq_region_start,
       'end'                     => $seq_region_end,
       'strand'                  => $seq_region_strand,
       'adaptor'                 => $self,
       'slice'                   => $slice,
       'dbID'                    => $gene_id,
       'stable_id'               => $stable_id,
       'version'                 => $version,
       'created_date'            => $created_date || undef,
       'modified_date'           => $modified_date || undef,
       'description'             => $gene_description,
       'external_name'           => undef,                      # will use display_id
                                                                # from display_xref
       'external_db'             => $external_db,
       'external_status'         => $external_status,
       'display_xref'            => $display_xref,
       'source'                  => $source,
       'is_current'              => $is_current,
       'canonical_transcript_id' => $canonical_transcript_id}));

  } ## end while ($sth->fetch())

  return \@genes;
} ## end sub _objs_from_sth

=head2 cache_gene_seq_mappings

  Example    : $gene_adaptor->cache_gene_seq_mappings();
  Description: caches all the assembly mappings needed for genes
  Returntype : None
  Exceptions : None
  Caller     : general
  Status     : At Risk
             : New experimental code

=cut

sub cache_gene_seq_mappings {
  my ($self) = @_;

  # get the sequence level to map too

  my $sql = 'SELECT name ' . 'FROM coord_system ' . 'WHERE attrib like "%%sequence_level%%"' . 'AND species_id = ?';

  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $self->species_id(), SQL_INTEGER);
  $sth->execute();

  my $sequence_level = $sth->fetchrow_array();

  $sth->finish();

  my $csa = $self->db->get_CoordSystemAdaptor();
  my $ama = $self->db->get_AssemblyMapperAdaptor();

  my $cs1 = $csa->fetch_by_name($sequence_level);

  # get level to map to two

  my $mcc   = $self->db->get_MetaCoordContainerAdaptor();
  my $csnew = $mcc->fetch_all_CoordSystems_by_feature_type('gene');

  foreach my $cs2 (@$csnew) {
    my $am = $ama->fetch_by_CoordSystems($cs1, $cs2);
    $am->register_all();
  }

} ## end sub cache_gene_seq_mappings

=head2 fetch_all_by_exon_supporting_evidence

  Arg [1]    : String $hit_name
               Name of supporting feature
  Arg [2]    : String $feature_type 
               one of "dna_align_feature" or "protein_align_feature"
  Arg [3]    : (optional) Bio::Ensembl::Analysis
  Example    : $genes = $gene_adaptor->fetch_all_by_exon_supporting_evidence(
                  'XYZ', 'dna_align_feature');
  Description: Gets all the genes with transcripts with exons which have a
               specified hit on a particular type of feature. Optionally filter
               by analysis.
  Returntype : Listref of Bio::EnsEMBL::Gene
  Exceptions : If feature_type is not of correct type.
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_exon_supporting_evidence {
  my ($self, $hit_name, $feature_type, $analysis) = @_;

  if ($feature_type !~ /(dna)|(protein)_align_feature/) {
    throw("feature type must be dna_align_feature or protein_align_feature");
  }

  my ($anal_from, $anal_where);
  if($analysis) {
    $anal_from = ", analysis a ";
    $anal_where = "AND a.analysis_id = f.analysis_id AND a.analysis_id=? ";
  }

  my $sql = qq(
      SELECT DISTINCT(g.gene_id)
        FROM gene g,
             transcript t,
             exon_transcript et,
             supporting_feature sf,
             $feature_type f
             $anal_from
       WHERE g.gene_id = t.gene_id
         AND g.is_current = 1
         AND t.transcript_id = et.transcript_id
         AND et.exon_id = sf.exon_id
         AND sf.feature_id = f.${feature_type}_id
         AND sf.feature_type = ?
         AND f.hit_name=?
         $anal_where
  );

  my $sth = $self->prepare($sql);

  $sth->bind_param(1, $feature_type, SQL_VARCHAR);
  $sth->bind_param(2, $hit_name,     SQL_VARCHAR);
  $sth->bind_param(3, $analysis->dbID(), SQL_INTEGER) if ($analysis);

  $sth->execute();

  my @genes;

  while (my $id = $sth->fetchrow_array) {
    my $gene = $self->fetch_by_dbID($id);
    push(@genes, $gene) if $gene;
  }

  return \@genes;
} ## end sub fetch_all_by_exon_supporting_evidence

=head2 fetch_all_by_transcript_supporting_evidence

  Arg [1]    : String $hit_name
               Name of supporting feature
  Arg [2]    : String $feature_type 
               one of "dna_align_feature" or "protein_align_feature"
  Arg [3]    : (optional) Bio::Ensembl::Analysis
  Example    : $genes = $gene_adaptor->fetch_all_by_transcript_supporting_evidence('XYZ', 'dna_align_feature');
  Description: Gets all the genes with transcripts with evidence for a
               specified hit on a particular type of feature. Optionally filter
               by analysis.
  Returntype : Listref of Bio::EnsEMBL::Gene.
  Exceptions : If feature_type is not of correct type.
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_transcript_supporting_evidence {
  my ($self, $hit_name, $feature_type, $analysis) = @_;

  if ($feature_type !~ /(dna)|(protein)_align_feature/) {
    throw("feature type must be dna_align_feature or protein_align_feature");
  }

  my ($anal_from, $anal_where);
  if($analysis) {
    $anal_from = ", analysis a ";
    $anal_where = "AND a.analysis_id = f.analysis_id AND a.analysis_id=? ";
  }

  my $sql = qq(
      SELECT DISTINCT(g.gene_id)
        FROM gene g,
             transcript t,
             transcript_supporting_feature sf,
             $feature_type f
             $anal_from
       WHERE g.gene_id = t.gene_id
         AND g.is_current = 1
         AND t.transcript_id = sf.transcript_id
         AND sf.feature_id = f.${feature_type}_id
         AND sf.feature_type = ?
         AND f.hit_name=?
         $anal_where
  );

  my $sth = $self->prepare($sql);

  $sth->bind_param(1, $feature_type, SQL_VARCHAR);
  $sth->bind_param(2, $hit_name,     SQL_VARCHAR);
  $sth->bind_param(3, $analysis->dbID(), SQL_INTEGER) if ($analysis);

  $sth->execute();

  my @genes;

  while (my $id = $sth->fetchrow_array) {
    my $gene = $self->fetch_by_dbID($id);
    push(@genes, $gene) if $gene;
  }

  return \@genes;
} ## end sub fetch_all_by_transcript_supporting_evidence

1;

