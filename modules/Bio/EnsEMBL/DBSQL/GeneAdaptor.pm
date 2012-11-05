
=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

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

  return ('g.gene_id', 'g.seq_region_id', 'g.seq_region_start', 'g.seq_region_end', 'g.seq_region_strand', 'g.analysis_id', 'g.biotype', 'g.display_xref_id', 'g.description', 'g.status', 'g.source', 'g.is_current', 'g.canonical_transcript_id', 'g.canonical_annotation', 'g.stable_id', 'g.version', $created_date, $modified_date, 'x.display_label', 'x.dbprimary_acc', 'x.description', 'x.version', 'exdb.db_name', 'exdb.status', 'exdb.db_release', 'exdb.db_display_name', 'x.info_type', 'x.info_text');
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

  return $gene;
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

  Arg [1]    : String $biotype 
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
  my ($self, $biotype) = @_;
  if (!defined $biotype) {
	throw("Biotype or listref of biotypes expected");
  }
  my $constraint;
  if (ref($biotype) eq 'ARRAY') {
	$constraint = "g.biotype IN (";
	foreach my $b (@{$biotype}) {
	  $constraint .= "?,";
	  $self->bind_param_generic_fetch($b, SQL_VARCHAR);
	}
	chop($constraint);    #remove last , from expression
	$constraint .= ") and g.is_current = 1";

  } else {
	$constraint = "g.biotype = ? and g.is_current = 1";
	$self->bind_param_generic_fetch($biotype, SQL_VARCHAR);
  }
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
  Arg [4]    : Name of the external database
  Example    : @genes = @{
                 $ga->fetch_all_by_Slice_and_external_dbname_link(
                                          $slice, undef, undef, "HUGO" ) };
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
  my $sth = $self->prepare("SELECT external_db_id FROM external_db WHERE db_name = ?");

  $sth->bind_param(1, $db_name, SQL_VARCHAR);
  $sth->execute();

  my $external_db_id;
  $sth->bind_columns(\$external_db_id);

  my @external_db_ids;
  while ($sth->fetch()) {
	push(@external_db_ids, $external_db_id);
  }

  if (scalar(@external_db_ids) == 0) {
	warn sprintf("Could not find external database " . "'%s' in the external_db table\n" . "Available are:\n", $db_name);

	$sth = $self->prepare("SELECT DISTINCT db_name FROM external_db");

	$sth->execute();
	$sth->bind_columns(\$external_db_id);

	while ($sth->fetch()) {
	  warn "\t$external_db_id\n";
	}
	return [];
  }

  # Get the gene_ids for those with links.
  my $dbe_adaptor = $self->db()->get_DBEntryAdaptor();

  my %linked_genes;
  foreach $external_db_id (@external_db_ids) {
	my @linked_genes = $dbe_adaptor->list_gene_ids_by_external_db_id($external_db_id);

	foreach my $gene_id (@linked_genes) {
	  $linked_genes{$gene_id} = 1;
	}
  }

  # Get all the genes on the slice.
  my $genes = $self->SUPER::fetch_all_by_Slice_constraint($slice, 'g.is_current = 1', $logic_name);

  # Create a list of those that are in the gene_ids list.
  my @genes_passed;
  foreach my $gene (@$genes) {
	if (exists($linked_genes{$gene->dbID()})) {
	  push(@genes_passed, $gene);
	}
  }

  # Return the list of those that passed.
  return \@genes_passed;
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
	$constraint .= " and g.biotype = '$biotype'";
  }

  my $genes = $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint, $logic_name);

  # If there are less than two genes, still do lazy-loading.
  if (!$load_transcripts || @$genes < 2) {
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

  my %genes_by_dbIDs =
	map { $_->dbID(), $_ } @{$self->fetch_all_by_dbID_list(\@ids)};

  my @result = map { $genes_by_dbIDs{$_} } @ids;

  return \@result;
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
  Example    : my @alt_genes = @{ $gene_adaptor->fetch_all_alt_alleles($gene) };
               foreach my $alt_gene (@alt_genes) {
                 print "Alternate allele: " . $alt_gene->stable_id() . "\n" ;
               }
  Description: Retrieves genes which are alternate alleles to a provided gene.
               Alternate alleles in Ensembl are genes which are similar and are
               on an alternative haplotype of the same region. There are not 
               currently very many of these. This method will return a 
               reference to an empty list if no alternative alleles are found.
  Returntype : listref of Bio::EnsEMBL::Genes
  Exceptions : throw if incorrect arg provided
               warning if gene arg does not have dbID
  Caller     : Gene::get_all_alt_alleles
  Status     : Stable

=cut

sub fetch_all_alt_alleles {
  my $self = shift;
  my $gene = shift;

  if (!ref($gene) || !$gene->isa('Bio::EnsEMBL::Gene')) {
	throw('Bio::EnsEMBL::Gene argument is required');
  }

  my $gene_id = $gene->dbID();

  if (!$gene_id) {
	warning('Cannot retrieve alternate alleles for gene without dbID');
	return [];
  }

  my $sth = $self->prepare("SELECT aa1.gene_id " . "FROM   alt_allele aa1, alt_allele aa2 " . "WHERE  aa1.alt_allele_id = aa2.alt_allele_id " . "AND    aa2.gene_id = ? " . "AND    aa1.gene_id <> ?");

  $sth->bind_param(1, $gene_id, SQL_INTEGER);
  $sth->bind_param(2, $gene_id, SQL_INTEGER);
  $sth->execute();

  my @alt_ids;
  my $row;
  while ($row = $sth->fetchrow_arrayref()) {
	push @alt_ids, $row->[0];
  }
  $sth->finish();

  if (@alt_ids) {
	return $self->fetch_all_by_dbID_list(\@alt_ids);
  }

  return [];
} ## end sub fetch_all_alt_alleles

sub is_ref {
  my ($self, $gene_id) = @_;
  my $is_not_ref;

  # easier to find if it is not an alt_Allele do this and then negate it.
  my $sth = $self->prepare("select count(1) from alt_allele where gene_id = $gene_id and is_ref = 0");
  $sth->execute();
  $sth->bind_columns(\$is_not_ref);
  $sth->fetch;

  if (defined($is_not_ref) and $is_not_ref) {
	return 0;
  }

  return 1;
}

=head2 store_alt_alleles


  Arg [1]    : reference to list of Bio::EnsEMBL::Genes $genes
  Example    : $gene_adaptor->store_alt_alleles([$gene1, $gene2, $gene3]);
  Description: This method creates a group of alternative alleles (i.e. locus)
               from a set of genes. The genes should be genes from alternate
               haplotypes which are similar. The genes must already be stored
               in this database. 
  Returntype : int alt_allele_id or undef if no alt_alleles were stored
  Exceptions : throw on incorrect arguments
               throw on sql error (e.g. duplicate unique id)
  Caller     : general
  Status     : Stable

=cut

sub store_alt_alleles {
  my $self  = shift;
  my $genes = shift;

  if (!ref($genes) eq 'ARRAY') {
	throw('List reference of Bio::EnsEMBL::Gene argument expected.');
  }

  my @genes     = @$genes;
  my $num_genes = scalar(@genes);

  if ($num_genes < 2) {
	warning('At least 2 genes must be provided to construct alternative alleles (gene id: ' . $genes[0]->dbID() . '). Ignoring.');
	return;
  }

  my @is_ref;
  my @ref_genes     = ();
  my @non_ref_genes = ();
  my @gene_ids      = ();

  foreach my $gene (@genes) {

	if (!ref($gene) || !$gene->isa('Bio::EnsEMBL::Gene')) {
	  throw('List reference of Bio::EnsEMBL::Gene argument expected.');
	}

	my $gene_id = $gene->dbID();

	if (!$gene_id) {
	  throw('Genes must have dbIDs in order to construct alternative alleles.');
	} else {
	  push @gene_ids, $gene_id;
	}

	my $is_ref = $gene->slice->is_reference();

	push @is_ref, $is_ref;

	if ($is_ref) {
	  push @ref_genes, $gene->dbID();
	} else {
	  push @non_ref_genes, $gene->dbID();
	}
  }
  if (scalar(@ref_genes) > 1) {
	warning('More than one alternative allele on the reference sequence (gene ids: ' . join(',', @ref_genes) . '). Ignoring.');
	return;
  }

  #
  #insert the first gene seperately in order to get a unique identifier for
  #the set of alleles
  #

  my $sth = $self->prepare("INSERT INTO alt_allele (gene_id, is_ref) VALUES (?,?)");
  $sth->bind_param(1, $gene_ids[0], SQL_INTEGER);
  $sth->bind_param(2, $is_ref[0],   SQL_INTEGER);
  eval { $sth->execute(); };
  my $alt_allele_id = $sth->{'mysql_insertid'};

  if (!$alt_allele_id || $@) {
	throw("An SQL error occured inserting alternative alleles:\n$@");
  }
  $sth->finish();
  #
  # Insert all subsequent alt alleles using the alt_allele identifier
  # from the first insert
  #

  $sth = $self->prepare("INSERT INTO alt_allele (alt_allele_id, gene_id, is_ref) " . "VALUES (?,?,?)");

  for (my $i = 1; $i < $num_genes; $i++) {

	$sth->bind_param(1, $alt_allele_id, SQL_INTEGER);
	$sth->bind_param(2, $gene_ids[$i],  SQL_INTEGER);
	$sth->bind_param(3, $is_ref[$i],    SQL_INTEGER);
	eval { $sth->execute(); };

	if ($@) {
	  # an error occured, revert the db to the previous state
	  $sth = $self->prepare("DELETE FROM alt_allele WHERE alt_allele_id = ?");
	  $sth->bind_param(1, $alt_allele_id, SQL_INTEGER);
	  $sth->execute();
	  $sth->finish();
	  throw("An SQL error occured inserting alternative alleles:\n$@");
	}
  }

  $sth->finish();

  return $alt_allele_id;
} ## end sub store_alt_alleles

=head2 store

  Arg [1]    : Bio::EnsEMBL::Gene $gene
               The gene to store in the database
  Arg [2]    : ignore_release in xrefs [default 1] set to 0 to use release info 
               in external database references
  Example    : $gene_adaptor->store($gene);
  Description: Stores a gene in the database.
  Returntype : the database identifier (dbID) of the newly stored gene
  Exceptions : thrown if the $gene is not a Bio::EnsEMBL::Gene or if 
               $gene does not have an analysis object
  Caller     : general
  Status     : Stable

=cut

sub store {
  my ($self, $gene, $ignore_release) = @_;

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

  my $store_gene_sql = qq(
        INSERT INTO gene
           SET biotype = ?,
               analysis_id = ?,
               seq_region_id = ?,
               seq_region_start = ?,
               seq_region_end = ?,
               seq_region_strand = ?,
               description = ?,
               source = ?,
               status = ?,
               is_current = ?,
               canonical_transcript_id = ?,
               canonical_annotation = ?
  );

  if (defined($gene->stable_id)) {
	my $created  = $self->db->dbc->from_seconds_to_date($gene->created_date());
	my $modified = $self->db->dbc->from_seconds_to_date($gene->modified_date());
	$store_gene_sql .= ", stable_id = ?, version = ?, created_date = " . $created . " , modified_date = " . $modified;

  }

  # column status is used from schema version 34 onwards (before it was
  # confidence)

  my $sth = $self->prepare($store_gene_sql);
  $sth->bind_param(1,  $type,                SQL_VARCHAR);
  $sth->bind_param(2,  $analysis_id,         SQL_INTEGER);
  $sth->bind_param(3,  $seq_region_id,       SQL_INTEGER);
  $sth->bind_param(4,  $gene->start(),       SQL_INTEGER);
  $sth->bind_param(5,  $gene->end(),         SQL_INTEGER);
  $sth->bind_param(6,  $gene->strand(),      SQL_TINYINT);
  $sth->bind_param(7,  $gene->description(), SQL_LONGVARCHAR);
  $sth->bind_param(8,  $gene->source(),      SQL_VARCHAR);
  $sth->bind_param(9,  $gene->status(),      SQL_VARCHAR);
  $sth->bind_param(10, $is_current,          SQL_TINYINT);

  # Canonical transcript ID will be updated later.
  # Set it to zero for now.
  $sth->bind_param(11, 0, SQL_TINYINT);

  $sth->bind_param(12, $gene->canonical_annotation(), SQL_VARCHAR);

  if (defined($gene->stable_id)) {

	$sth->bind_param(13, $gene->stable_id, SQL_VARCHAR);
	my $version = ($gene->version()) ? $gene->version() : 1;
	$sth->bind_param(14, $version, SQL_INTEGER);
  }

  $sth->execute();
  $sth->finish();

  my $gene_dbID = $sth->{'mysql_insertid'};

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
		$trans->swap_exons($e, $exons{$key});
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

	$transcript_adaptor->store($new, $gene_dbID, $analysis_id);

	if (!defined($new_canonical_transcript_id)
		&& $new->is_canonical())
	{
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

  # store unconventional transcript associations if there are any
  my $utaa = $db->get_UnconventionalTranscriptAssociationAdaptor();
  foreach my $uta (@{$gene->get_all_unconventional_transcript_associations()}) {
	$utaa->store($uta);
  }

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

  # remove any unconventional transcript associations involving this gene

  $sth = $self->prepare("DELETE FROM unconventional_transcript_association " . "WHERE gene_id = ? ");
  $sth->bind_param(1, $gene->dbID, SQL_INTEGER);
  $sth->execute();
  $sth->finish();

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
  Description: Updates the type, analysis, display_xref, status, is_current and
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
              status = ?,
              description = ?,
              is_current = ?,
              canonical_transcript_id = ?,
              canonical_annotation = ?
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
  $sth->bind_param(4, $gene->status(),         SQL_VARCHAR);
  $sth->bind_param(5, $gene->description(),    SQL_VARCHAR);
  $sth->bind_param(6, $gene->is_current(),     SQL_TINYINT);

  if (defined($gene->canonical_transcript())) {
	$sth->bind_param(7, $gene->canonical_transcript()->dbID(), SQL_INTEGER);
  } else {
	$sth->bind_param(7, 0, SQL_INTEGER);
  }

  $sth->bind_param(8, $gene->canonical_annotation(), SQL_VARCHAR);
  $sth->bind_param(9, $gene->dbID(), SQL_INTEGER);

  $sth->execute();

  # maybe should update stable id ???
} ## end sub update

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

  my ($gene_id, $seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand, $analysis_id, $biotype, $display_xref_id, $gene_description, $status, $source, $is_current, $canonical_transcript_id, $canonical_annotation, $stable_id, $version, $created_date, $modified_date, $xref_display_id, $xref_primary_acc, $xref_desc, $xref_version, $external_db, $external_status, $external_release, $external_db_name, $info_type, $info_text);

  $sth->bind_columns(\($gene_id, $seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand, $analysis_id, $biotype, $display_xref_id, $gene_description, $status, $source, $is_current, $canonical_transcript_id, $canonical_annotation, $stable_id, $version, $created_date, $modified_date, $xref_display_id, $xref_primary_acc, $xref_desc, $xref_version, $external_db, $external_status, $external_release, $external_db_name, $info_type, $info_text));

  my $asm_cs;
  my $cmp_cs;
  my $asm_cs_vers;
  my $asm_cs_name;
  my $cmp_cs_vers;
  my $cmp_cs_name;

  if ($mapper) {
	$asm_cs      = $mapper->assembled_CoordSystem();
	$cmp_cs      = $mapper->component_CoordSystem();
	$asm_cs_name = $asm_cs->name();
	$asm_cs_vers = $asm_cs->version();
	$cmp_cs_name = $cmp_cs->name();
	$cmp_cs_vers = $cmp_cs->version();
  }

  my $dest_slice_start;
  my $dest_slice_end;
  my $dest_slice_strand;
  my $dest_slice_length;
  my $dest_slice_sr_name;
  my $dest_slice_sr_id;

  if ($dest_slice) {
	$dest_slice_start   = $dest_slice->start();
	$dest_slice_end     = $dest_slice->end();
	$dest_slice_strand  = $dest_slice->strand();
	$dest_slice_length  = $dest_slice->length();
	$dest_slice_sr_name = $dest_slice->seq_region_name();
	$dest_slice_sr_id   = $dest_slice->get_seq_region_id();
  }

FEATURE: while ($sth->fetch()) {
	#get the analysis object
	my $analysis = $analysis_hash{$analysis_id} ||= $aa->fetch_by_dbID($analysis_id);

	#need to get the internal_seq_region, if present
	$seq_region_id = $self->get_seq_region_id_internal($seq_region_id);
	my $slice = $slice_hash{"ID:" . $seq_region_id};

	if (!$slice) {
	  $slice                              = $sa->fetch_by_seq_region_id($seq_region_id);
	  $slice_hash{"ID:" . $seq_region_id} = $slice;
	  $sr_name_hash{$seq_region_id}       = $slice->seq_region_name();
	  $sr_cs_hash{$seq_region_id}         = $slice->coord_system();
	}

	my $sr_name = $sr_name_hash{$seq_region_id};
	my $sr_cs   = $sr_cs_hash{$seq_region_id};

	#
	# remap the feature coordinates to another coord system
	# if a mapper was provided
	#
	if ($mapper) {

	  if (defined $dest_slice
		  && $mapper->isa('Bio::EnsEMBL::ChainedAssemblyMapper'))
	  {
		($seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand) = $mapper->map($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand, $sr_cs, 1, $dest_slice);

	  } else {

		($seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand) = $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand, $sr_cs);
	  }

	  #skip features that map to gaps or coord system boundaries
	  next FEATURE if (!defined($seq_region_id));

	  #get a slice in the coord system we just mapped to
	  #      if($asm_cs == $sr_cs || ($cmp_cs != $sr_cs && $asm_cs->equals($sr_cs))) {
	  $slice = $slice_hash{"ID:" . $seq_region_id} ||= $sa->fetch_by_seq_region_id($seq_region_id);
	  #      } else {
	  #        $slice = $slice_hash{"NAME:$sr_name:$asm_cs_name:$asm_cs_vers"} ||=
	  #          $sa->fetch_by_region($asm_cs_name, $sr_name, undef, undef, undef,
	  #                               $asm_cs_vers);
	  #      }
	}

	#
	# If a destination slice was provided convert the coords.
	#
	if (defined($dest_slice)) {
	  if ($dest_slice_strand == 1) {
		# Positive strand.

		$seq_region_start = $seq_region_start - $dest_slice_start + 1;
		$seq_region_end   = $seq_region_end - $dest_slice_start + 1;

		if ($dest_slice->is_circular()) {
		  # Handle cicular chromosomes.

		  if ($seq_region_start > $seq_region_end) {
			# Looking at a feature overlapping the chromsome origin.

			if ($seq_region_end > $dest_slice_start) {
			  # Looking at the region in the beginning of the
			  # chromosome.
			  $seq_region_start -= $dest_slice->seq_region_length();
			}

			if ($seq_region_end < 0) {
			  $seq_region_end += $dest_slice->seq_region_length();
			}

		  } else {

			if (   $dest_slice_start > $dest_slice_end
				&& $seq_region_end < 0)
			{
			  # Looking at the region overlapping the chromosome
			  # origin and a feature which is at the beginning of the
			  # chromosome.
			  $seq_region_start += $dest_slice->seq_region_length();
			  $seq_region_end   += $dest_slice->seq_region_length();
			}
		  }

		} ## end if ($dest_slice->is_circular...)

	  } else {
		# Negative strand.

		if (   $dest_slice->is_circular()
			&& $seq_region_start > $seq_region_end)
		{
		  # Handle cicular chromosomes.

		  if ($seq_region_end > $dest_slice_start) {
			# Looking at the region in the beginning of the
			# chromosome.
			$seq_region_start = $dest_slice_end - $seq_region_end + 1;
			$seq_region_end   = $seq_region_end - $dest_slice->seq_region_length - $dest_slice_start + 1;
		  } else {
			my $tmp_seq_region_start = $seq_region_start;
			$seq_region_start = $dest_slice_end - $seq_region_end - $dest_slice->seq_region_length + 1;
			$seq_region_end   = $dest_slice_end - $tmp_seq_region_start + 1;
		  }

		} else {
		  # Non-circular chromosome.

		  my $tmp_seq_region_start = $seq_region_start;
		  $seq_region_start = $dest_slice_end - $seq_region_end + 1;
		  $seq_region_end   = $dest_slice_end - $tmp_seq_region_start + 1;
		}

		$seq_region_strand = -$seq_region_strand;

	  } ## end else [ if ($dest_slice_strand...)]

	  # Throw away features off the end of the requested slice or on
	  # different seq_region.

	  if (   $seq_region_end < 1
		  || $seq_region_start > $dest_slice_length
		  || ($dest_slice_sr_id ne $seq_region_id))
	  {
		next FEATURE;
	  }

	  $slice = $dest_slice;
	} ## end if (defined($dest_slice...))

	my $display_xref;

	if ($display_xref_id) {
	  $display_xref = Bio::EnsEMBL::DBEntry->new_fast({'dbID'            => $display_xref_id,
													   'adaptor'         => $dbEntryAdaptor,
													   'display_id'      => $xref_display_id,
													   'primary_id'      => $xref_primary_acc,
													   'version'         => $xref_version,
													   'description'     => $xref_desc,
													   'release'         => $external_release,
													   'dbname'          => $external_db,
													   'db_display_name' => $external_db_name,
													   'info_type'       => $info_type,
													   'info_text'       => $info_text});
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
		 'status'                  => $status,
		 'source'                  => $source,
		 'is_current'              => $is_current,
		 'canonical_transcript_id' => $canonical_transcript_id,
		 'canonical_annotation'    => $canonical_annotation}));

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

  my $anal_from = ", analysis a " if ($analysis);
  my $anal_where = "AND a.analysis_id = f.analysis_id AND a.analysis_id=? "
	if ($analysis);

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

  my $anal_from = ", analysis a " if ($analysis);
  my $anal_where = "AND a.analysis_id = f.analysis_id AND a.analysis_id=? "
	if ($analysis);

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

=head2 fetch_nearest_Gene_by_Feature

  Arg [1]    : Feature object
  Example    : $genes = $gene_adaptor->fetch_nearest_Gene_by_Feature($feat);
  Description: Gets the nearest gene to the feature 
  Returntype : Listref of Bio::EnsEMBL::Gene, EMPTY list if no nearest
  Caller     : general
  Status     : UnStable

=cut

sub fetch_nearest_Gene_by_Feature {
  my $self = shift;
  my $feat = shift;

  my $stranded = shift;
  my $stream   = shift;    # 1 up stream -1 downstream
  my @genes;

  my $strand = $feat->strand;
  if (defined($stream) and !$strand) {
	warn("stream specified but feature has no strand so +ve strand will be used");
	$strand = 1;
  }
  my $min_dist = 999;
  my $gene_id  = 0;

  my $overlapping = $feat->get_overlapping_Genes();

  return $overlapping if (defined(@{$overlapping}[0]));

  my $seq_region_id = $feat->slice->adaptor->get_seq_region_id($feat->slice);
  my $start         = ($feat->start + $feat->slice->start) - 1;
  my $end           = ($feat->end + $feat->slice->start) - 1;

  my @gene_ids;
  if (!defined($stream) or $stream == 0) {

	my $sql1 = "select g.gene_id, (? - g.seq_region_end)  as 'dist' from gene g where ";
	if ($stranded) {
	  $sql1 .= "g.seq_region_strand = " . $strand . " and ";
	}
	$sql1 .= "seq_region_id = ? and g.seq_region_end < ? order by dist limit 10";

	#
	# MAYBE set the result of prepare to be static in case lots of calls.
	#
	my $sql1_sth = $self->prepare($sql1) || die "Could not prepare $sql1";
	$sql1_sth->execute($start, $seq_region_id, $start)
	  || die "Could not execute sql";
	$sql1_sth->bind_columns(\$gene_id, \$min_dist)
	  || die "Could mot bin columns";

	my $last_dist = 99999999999999999;
	while ($sql1_sth->fetch()) {
	  if ($min_dist <= $last_dist) {
		push @gene_ids, $gene_id;
		$last_dist = $min_dist;
	  }
	}
	$sql1_sth->finish();

	my $sql2 = "select g.gene_id, (g.seq_region_start - ?)  as 'dist' from gene g  where ";
	if ($stranded) {
	  $sql2 .= "g.seq_region_strand = " . $feat->strand . " and ";
	}
	$sql2 .= "seq_region_id = ? and g.seq_region_start > ? order by dist limit 10";

	my $sql2_sth = $self->prepare($sql2) || die "could not prepare $sql2";

	my ($tmp_min_dist, $tmp_gene_id);
	$sql2_sth->execute($end, $seq_region_id, $end)
	  || die "Could not execute sql";
	$sql2_sth->bind_columns(\$tmp_gene_id, \$tmp_min_dist)
	  || die "Could mot bin columns";
	my $first = 1;
	while ($sql2_sth->fetch()) {
	  if ($tmp_min_dist <= $last_dist) {
		if ($first) {
		  $first = 0;
		  if ($tmp_min_dist < $last_dist) {
			@gene_ids = ();    #reset
		  }
		}
		push @gene_ids, $tmp_gene_id;
		$last_dist = $tmp_min_dist;
	  }
	}
	$sql2_sth->finish();

  } elsif (($stream*$strand) == 1) {
	my $sql1 = "select g.gene_id, (? - g.seq_region_end)  as 'dist' from gene g where ";
	if ($stranded) {
	  $sql1 .= "g.seq_region_strand = " . $strand . " and ";
	}
	$sql1 .= "seq_region_id = ? and g.seq_region_end < ? order by dist limit 10";

	#
	# MAYBE set the result of prepare to be static in case lots of calls.
	#
	my $sql1_sth = $self->prepare($sql1) || die "Could not prepare $sql1";
	$sql1_sth->execute($start, $seq_region_id, $start)
	  || die "Could not execute sql";
	$sql1_sth->bind_columns(\$gene_id, \$min_dist)
	  || die "Could mot bin columns";

	my $last_dist;
	my $first = 1;
	while ($sql1_sth->fetch()) {
	  if ($first) {
		$first = 0;
	  } else {
		next if ($min_dist > $last_dist);
	  }
	  push @gene_ids, $gene_id;
	  $last_dist = $min_dist;
	}
	$sql1_sth->finish();
  } elsif (($stream*$strand) == -1) {

	my $sql2 = "select g.gene_id, (g.seq_region_start - ?)  as 'dist' from gene g  where ";
	if ($stranded) {
	  $sql2 .= "g.seq_region_strand = " . $feat->strand . " and ";
	}
	$sql2 .= "seq_region_id = ? and g.seq_region_start > ? order by dist limit 10";

	my $sql2_sth = $self->prepare($sql2) || die "could not prepare $sql2";

	my ($tmp_min_dist, $tmp_gene_id);
	$sql2_sth->execute($end, $seq_region_id, $end)
	  || die "Could not execute sql";
	$sql2_sth->bind_columns(\$tmp_gene_id, \$tmp_min_dist)
	  || die "Could mot bin columns";
	my $first = 1;
	my $last_dist;
	while ($sql2_sth->fetch()) {
	  if ($first) {
		$first = 0;
	  } else {
		next if ($tmp_min_dist > $last_dist);
	  }
	  push @gene_ids, $tmp_gene_id;
	  $last_dist = $tmp_min_dist;
	}
	$sql2_sth->finish();
  } else {
	die "Invalid stream or strand must be -1, 0 or 1\n";
  }

  foreach my $gene_id (@gene_ids) {
	push @genes, $self->fetch_by_dbID($gene_id);
  }
  return \@genes;

} ## end sub fetch_nearest_Gene_by_Feature

##########################
#                        #
#  DEPRECATED METHODS    #
#                        #
##########################

=head2 fetch_by_maximum_DBLink

 DEPRECATED - use fetch_all_by_external_name instead

=cut

sub fetch_by_maximum_DBLink {
  my ($self, $external_id) = @_;

  deprecate("use fetch_all_by_external_name instead");

  my $genes = $self->fetch_all_by_external_name($external_id);

  my $biggest;
  my $max  = 0;
  my $size = scalar(@$genes);
  if ($size > 0) {
	foreach my $gene (@$genes) {
	  my $size = scalar(@{$gene->get_all_Exons});
	  if ($size > $max) {
		$biggest = $gene;
		$max     = $size;
	  }
	}
	return $biggest;
  }
  return;
}

=head2 get_display_xref

  DEPRECATED use $gene->display_xref

=cut

sub get_display_xref {
  my ($self, $gene) = @_;

  deprecate("display xref should retrieved from Gene object directly");

  if (!defined $gene) {
	throw("Must call with a Gene object");
  }

  my $sth = $self->prepare(
	qq(
      SELECT e.db_name,
             x.display_label,
             x.xref_id
      FROM   gene g, 
             xref x, 
             external_db e
      WHERE  g.gene_id = ?
        AND  g.display_xref_id = x.xref_id
        AND  x.external_db_id = e.external_db_id
  ));

  $sth->bind_param(1, $gene->dbID, SQL_INTEGER);
  $sth->execute();

  my ($db_name, $display_label, $xref_id) = $sth->fetchrow_array();
  if (!defined $xref_id) {
	return undef;
  }

  my $db_entry = Bio::EnsEMBL::DBEntry->new(-dbid       => $xref_id,
											-adaptor    => $self->db->get_DBEntryAdaptor(),
											-dbname     => $db_name,
											-display_id => $display_label);

  return $db_entry;
} ## end sub get_display_xref

=head2 get_description

  DEPRECATED, use gene->get_description

=cut

sub get_description {
  my ($self, $dbID) = @_;

  deprecate("Gene description should be loaded on gene retrieval. Use gene->get_description()");

  if (!defined $dbID) {
	throw("must call with dbID");
  }

  my $sth = $self->prepare(
	"SELECT description 
                            FROM   gene_description 
                            WHERE  gene_id = ?");

  $sth->bind_param(1, $dbID, SQL_INTEGER);
  $sth->execute();

  my @array = $sth->fetchrow_array();
  return $array[0];
}

=head2 fetch_by_Peptide_id

  DEPRECATED, use fetch_by_translation_stable_id()

=cut

sub fetch_by_Peptide_id {
  my ($self, $translation_stable_id) = @_;

  deprecate("Please use better named fetch_by_translation_stable_id \n" . caller(2));

  $self->fetch_by_translation_stable_id($translation_stable_id);
}

=head2 get_stable_entry_info

  DEPRECATED use $gene->stable_id instead

=cut

sub get_stable_entry_info {
  my ($self, $gene) = @_;

  deprecated("stable id info is loaded on default, no lazy loading necessary");

  if (!defined $gene || !ref $gene || !$gene->isa('Bio::EnsEMBL::Gene')) {
	throw("Needs a gene object, not a $gene");
  }

  my $created_date  = $self->db->dbc->from_date_to_seconds("created_date");
  my $modified_date = $self->db->dbc->from_date_to_seconds("modified_date");

  my $sth = $self->prepare("SELECT stable_id, " . $created_date . "," . $modified_date . ", version FROM gene WHERE gene_id = ?");

  $sth->bind_param(1, $gene->dbID, SQL_INTEGER);
  $sth->execute();

  my @array = $sth->fetchrow_array();
  $gene->{'stable_id'} = $array[0];
  $gene->{'created'}   = $array[1];
  $gene->{'modified'}  = $array[2];
  $gene->{'version'}   = $array[3];

  return 1;
} ## end sub get_stable_entry_info

=head2 fetch_all_by_DBEntry

  DEPRECATED - Use fetch_all_by_external_name instead

=cut

sub fetch_all_by_DBEntry {
  my $self = shift;

  deprecate('Use fetch_all_by_external_name instead.');

  return $self->fetch_all_by_external_name(@_);
}

1;

