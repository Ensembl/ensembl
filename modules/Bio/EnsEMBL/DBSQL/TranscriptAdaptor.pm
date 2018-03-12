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

Bio::EnsEMBL::DBSQL::TranscriptAdaptor - An adaptor which performs database
interaction relating to the storage and retrieval of Transcripts

=head1 SYNOPSIS

  use Bio::EnsEMBL::Registry;

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
  );

  $transcript_adaptor =
    Bio::EnsEMBL::Registry->get_adaptor( 'Human', 'Core',
    'Transcript' );

  $transcript = $transcript_adaptor->fetch_by_dbID(1234);

  $transcript =
    $transcript_adaptor->fetch_by_stable_id('ENST00000201961');

  $slice =
    $slice_adaptor->fetch_by_region( 'Chromosome', '3', 1, 1000000 );
  @transcripts = @{ $transcript_adaptor->fetch_all_by_Slice($slice) };

  ($transcript) =
    @{ $transcript_adaptor->fetch_all_by_external_name('NP_065811.1') };

=head1 DESCRIPTION

This adaptor provides a means to retrieve and store information related
to Transcripts.  Primarily this involves the retrieval or storage of
Bio::EnsEMBL::Transcript objects from a database.

See Bio::EnsEMBL::Transcript for details of the Transcript class.

=cut

package Bio::EnsEMBL::DBSQL::TranscriptAdaptor;

use strict;

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning );
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );

use vars qw(@ISA);
@ISA = qw( Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor );


# _tables
#
#  Description: PROTECTED implementation of superclass abstract method.
#               Returns the names, aliases of the tables to use for queries.
#  Returntype : list of listrefs of strings
#  Exceptions : none
#  Caller     : internal
#  Status     : Stable

sub _tables {
  return (
    [ 'transcript',           't' ],
    [ 'xref',                 'x' ],
    [ 'external_db',          'exdb' ] );
}


#_columns
#
#  Description: PROTECTED implementation of superclass abstract method.
#               Returns a list of columns to use for queries.
#  Returntype : list of strings
#  Exceptions : none
#  Caller     : internal
#  Status     : Stable

sub _columns {
  my ($self) = @_;

  my $created_date =
    $self->db()->dbc()->from_date_to_seconds("created_date");
  my $modified_date =
    $self->db()->dbc()->from_date_to_seconds("modified_date");

  my @columns = 
    (
     't.transcript_id',     't.seq_region_id',
     't.seq_region_start',  't.seq_region_end',
     't.seq_region_strand', 't.analysis_id',
     't.gene_id',           't.is_current',
     't.stable_id',         't.version',
     $created_date,         $modified_date,
     't.description',       't.biotype',
     'exdb.db_name',
     'exdb.status',         'exdb.db_display_name',
     'x.xref_id',           'x.display_label',
     'x.dbprimary_acc',     'x.version',
     'x.description',       'x.info_type',
     'x.info_text',         'exdb.db_release'
    );

  $self->schema_version > 74 and push @columns, 't.source';

  return @columns;
}

sub _left_join {
  return (
    [ 'xref',                 "x.xref_id = t.display_xref_id" ],
    [ 'external_db',          "exdb.external_db_id = x.external_db_id" ]
  );
}


=head2 fetch_by_stable_id

  Arg [1]    : String $stable_id 
               The stable id of the transcript to retrieve
  Example    : my $tr = $tr_adaptor->fetch_by_stable_id('ENST00000309301');
  Description: Retrieves a transcript via its stable id.
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_stable_id {
  my ($self, $stable_id) = @_;

  my $constraint = "t.stable_id = ? AND t.is_current = 1";

  $self->bind_param_generic_fetch($stable_id,SQL_VARCHAR);

  my ($transcript) = @{ $self->generic_fetch($constraint) };

  # If we didn't get anything back, desperately try to see if there's
  # a version number in the stable_id
  if(!defined($transcript) && (my $vindex = rindex($stable_id, '.'))) {
      $transcript = $self->fetch_by_stable_id_version(substr($stable_id,0,$vindex),
						substr($stable_id,$vindex+1));
  }

  return $transcript;
}


=head2 fetch_by_stable_id_version

  Arg [1]    : String $id 
               The stable ID of the transcript to retrieve
  Arg [2]    : Integer $version
               The version of the stable_id to retrieve
  Example    : $tr = $tr_adaptor->fetch_by_stable_id('ENST00000309301', 3);
  Description: Retrieves a transcript object from the database via its 
               stable id and version.
               The transcript will be retrieved in its native coordinate system (i.e.
               in the coordinate system it is stored in the database). It may
               be converted to a different coordinate system through a call to
               transform() or transfer(). If the transcript is not found
               undef is returned instead.
  Returntype : Bio::EnsEMBL::Transcript or undef
  Exceptions : if we cant get the transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_stable_id_version {
    my ($self, $stable_id, $version) = @_;

    # Enforce that version be numeric
    return unless($version =~ /^\d+$/);

    my $constraint = "t.stable_id = ? AND t.version = ? AND t.is_current = 1";
    $self->bind_param_generic_fetch($stable_id, SQL_VARCHAR);
    $self->bind_param_generic_fetch($version, SQL_INTEGER);
    my ($transcript) = @{$self->generic_fetch($constraint)};

    return $transcript;
}

sub fetch_all {
  my ($self) = @_;

  my $constraint = 't.biotype != "LRG_gene" and t.is_current = 1';
  my @trans  = @{ $self->generic_fetch($constraint) };
  return \@trans ;
}

=head2 fetch_all_versions_by_stable_id 

  Arg [1]     : String $stable_id 
                The stable ID of the transcript to retrieve
  Example     : my $tr = $tr_adaptor->fetch_all_version_by_stable_id
                  ('ENST00000309301');
  Description : Similar to fetch_by_stable_id, but retrieves all versions of a
                transcript stored in the database.
  Returntype  : listref of Bio::EnsEMBL::Transcript objects
  Exceptions  : if we cant get the gene in given coord system
  Caller      : general
  Status      : At Risk

=cut

sub fetch_all_versions_by_stable_id {
  my ($self, $stable_id) = @_;

  my $constraint = "t.stable_id = ?";

  $self->bind_param_generic_fetch($stable_id,SQL_VARCHAR);

  return $self->generic_fetch($constraint);
}


=head2 fetch_by_translation_stable_id

  Arg [1]    : String $transl_stable_id
               The stable identifier of the translation of the transcript to 
               retrieve
  Example    : my $tr = $tr_adaptor->fetch_by_translation_stable_id
                  ('ENSP00000311007');
  Description: Retrieves a Transcript object using the stable identifier of
               its translation.
  Returntype : Bio::EnsEMBL::Transcript or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_translation_stable_id {
  my ($self, $transl_stable_id ) = @_;

  my $sth = $self->prepare(qq(
      SELECT t.transcript_id
      FROM   translation tl,
             transcript t
      WHERE  tl.stable_id = ?
      AND    tl.transcript_id = t.transcript_id
      AND    t.is_current = 1
  ));

  $sth->bind_param(1, $transl_stable_id, SQL_VARCHAR);
  $sth->execute();

  my ($id) = $sth->fetchrow_array;
  $sth->finish;
  if ($id){
    return $self->fetch_by_dbID($id);
  } elsif(my $vindex = rindex($transl_stable_id, '.')) {
    return $self->fetch_by_translation_stable_id_version(substr($transl_stable_id,0,$vindex),
							 substr($transl_stable_id,$vindex+1));
  } else {
      return undef;
  }
}

=head2 fetch_by_translation_stable_id_version

  Arg [1]    : String $transl_stable_id
               The stable identifier of the translation of the transcript to 
               retrieve
  Arg [2]    : Integer $version
               The version of the translation of the transcript to retrieve
  Example    : my $tr = $tr_adaptor->fetch_by_translation_stable_id_version
                  ('ENSP00000311007', 2);
  Description: Retrieves a Transcript object using the stable identifier and
               version of its translation.
  Returntype : Bio::EnsEMBL::Transcript or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_translation_stable_id_version {
  my ($self, $transl_stable_id, $transl_version ) = @_;

  # Enforce that version be numeric
  return unless($transl_version =~ /^\d+$/);

  my $sth = $self->prepare(qq(
      SELECT t.transcript_id
      FROM   translation tl,
             transcript t
      WHERE  tl.stable_id = ?
      AND    tl.version = ?
      AND    tl.transcript_id = t.transcript_id
      AND    t.is_current = 1
  ));

  $sth->bind_param(1, $transl_stable_id, SQL_VARCHAR);
  $sth->bind_param(2, $transl_version, SQL_INTEGER);
  $sth->execute();

  my ($id) = $sth->fetchrow_array;
  $sth->finish;
  if ($id){
    return $self->fetch_by_dbID($id);
  } else {
    return undef;
  }
}

=head2 fetch_by_translation_id

  Arg [1]    : Int $id
               The internal identifier of the translation whose transcript
               is to be retrieved
  Example    : my $tr = $tr_adaptor->fetch_by_translation_id($transl->dbID);
  Description: Given the internal identifier of a translation this method 
               retrieves the transcript associated with that translation.
               If the transcript cannot be found undef is returned instead.
  Returntype : Bio::EnsEMBL::Transcript or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_translation_id {
  my ( $self, $p_dbID ) = @_;

  if ( !defined($p_dbID) ) {
    throw("dbID argument is required");
  }

  my $sth =
    $self->prepare(   "SELECT transcript_id "
                    . "FROM   translation "
                    . "WHERE  translation_id = ?" );

  $sth->bind_param( 1, $p_dbID, SQL_INTEGER );
  $sth->execute();

  my ($dbID) = $sth->fetchrow_array();
  $sth->finish();

  if ($dbID) {
    return $self->fetch_by_dbID($dbID);
  }

  return undef;
}

=head2 fetch_all_by_Gene

  Arg [1]    : Bio::EnsEMBL::Gene $gene
               The gene to fetch transcripts of
  Example    : my $gene = $gene_adaptor->fetch_by_stable_id('ENSG0000123');
               my @transcripts = { $tr_adaptor->fetch_all_by_Gene($gene) };
  Description: Retrieves Transcript objects for given gene. Puts Genes slice
               in each Transcript. 
  Returntype : Listref of Bio::EnsEMBL::Transcript objects
  Exceptions : none
  Caller     : Gene->get_all_Transcripts()
  Status     : Stable

=cut

sub fetch_all_by_Gene {
  my ( $self, $gene ) = @_;

  my $constraint = "t.gene_id = " . $gene->dbID();

  # Use the fetch_all_by_Slice_constraint method because it handles the
  # difficult Haps/PARs and coordinate remapping.

  # Get a slice that entirely overlaps the gene.  This is because we
  # want all transcripts to be retrieved, not just ones overlapping
  # the slice the gene is on (the gene may only partially overlap the
  # slice).  For speed reasons, only use a different slice if necessary
  # though.

  my $gslice = $gene->slice();

  if ( !defined($gslice) ) {
    throw("Gene must have attached slice to retrieve transcripts.");
  }

  my $slice;

  if ( $gene->start() < 1 || $gene->end() > $gslice->length() ) {
    if ( $gslice->is_circular() ) {
      $slice = $gslice;
    } else {
      $slice = $self->db->get_SliceAdaptor->fetch_by_Feature($gene);
    }
  } else {
    $slice = $gslice;
  }

  my $transcripts =
    $self->fetch_all_by_Slice_constraint( $slice, $constraint );

  if ( $slice != $gslice ) {
    my @out;
    foreach my $tr ( @{$transcripts} ) {
      push( @out, $tr->transfer($gslice) );
    }
    $transcripts = \@out;
  }

  my $canonical_t = $gene->canonical_transcript();

  foreach my $t ( @{$transcripts} ) {
    if ( $t->equals($canonical_t) ) {
      $t->is_canonical(1);
      last;
    }
  }

  return $transcripts;
} ## end sub fetch_all_by_Gene


=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to fetch transcripts on
  Arg [2]    : (optional) Boolean $load_exons
               If true, exons will be loaded immediately rather than
               lazy loaded later
  Arg [3]    : (optional) String $logic_name
               The logic name of the type of features to obtain
  ARG [4]    : (optional) String $constraint
               An extra contraint.
  Example    : my @transcripts = @{ $tr_adaptor->fetch_all_by_Slice($slice) };
  Description: Overrides superclass method to optionally load exons
               immediately rather than lazy-loading them later. This
               is more efficient when there are a lot of transcripts whose
               exons are going to be used.
  Returntype : Listref of Bio::EnsEMBL::Transcript objects
  Exceptions : thrown if exon cannot be placed on transcript slice
  Caller     : Slice::get_all_Transcripts
  Status     : Stable

=cut

sub fetch_all_by_Slice {
  my ( $self, $slice, $load_exons, $logic_name, $constraint, $source, $biotype ) = @_;

  if (defined $constraint and $constraint ne '') {
    $constraint .= ' AND t.is_current = 1';
  } else {
    $constraint .= 't.is_current = 1';
  }
  if (defined($source)) {
    $constraint .= " and t.source = '$source'";
  }
  if (defined($biotype)) {
    my $inline_variables = 1;
    $constraint .= " and ".$self->generate_in_constraint($biotype, 't.biotype', SQL_VARCHAR, $inline_variables);
  }

  my $transcripts = $self->SUPER::fetch_all_by_Slice_constraint( $slice, $constraint, $logic_name);

  # if there are 0 transcripts still do lazy-loading
  if ( !$load_exons || @$transcripts < 1 ) {
    return $transcripts;
  }

  # preload all of the exons now, instead of lazy loading later
  # faster than 1 query per transcript

  # first check if the exons are already preloaded
  # @todo FIXME: Should test all exons.
  if ( exists( $transcripts->[0]->{'_trans_exon_array'} ) ) {
    return $transcripts;
  }

  # get extent of region spanned by transcripts
  my ($min_start, $max_end);
  my $ext_slice;

  unless ($slice->is_circular()) {
    foreach my $t (@$transcripts) {
      if (!defined($min_start) || $t->seq_region_start() < $min_start) {
	$min_start = $t->seq_region_start();
      }
      if (!defined($max_end) || $t->seq_region_end() > $max_end) {
	$max_end = $t->seq_region_end();
      }
    }

    if ($min_start >= $slice->start() && $max_end <= $slice->end()) {
      $ext_slice = $slice;
    } else {
      my $sa = $self->db()->get_SliceAdaptor();
      $ext_slice = $sa->fetch_by_region($slice->coord_system->name(), $slice->seq_region_name(), $min_start, $max_end, $slice->strand(), $slice->coord_system->version());
    }

  } else {
    # feature might be crossing the origin of replication (i.e. seq_region_start > seq_region_end)
    # the computation of min_start|end based on seq_region_start|end is not safe
    # use feature start/end relative to the slice instead
    my ($min_start_feature, $max_end_feature);
    foreach my $t (@$transcripts) {
      if (!defined($min_start) || ($t->start >= 0 && $t->start() < $min_start)) {
  	$min_start = $t->start();
  	$min_start_feature = $t;
      }
      if (!defined($max_end) || ($t->end() >= 0 && $t->end() > $max_end)) {
  	$max_end = $t->end();
  	$max_end_feature = $t;
      }
    }
    
    # now we can reassign min_start|end to seq_region_start|end of
    # the feature which spans the largest region
    $min_start = $min_start_feature->seq_region_start();
    $max_end = $max_end_feature->seq_region_end();

    my $sa = $self->db()->get_SliceAdaptor();
    $ext_slice = 
      $sa->fetch_by_region($slice->coord_system->name(), 
  			   $slice->seq_region_name(), 
  			   $min_start, 
  			   $max_end, 
  			   $slice->strand(), 
  			   $slice->coord_system->version());
  }



  # associate exon identifiers with transcripts

  my %tr_hash = map { $_->dbID => $_ } @{$transcripts};

  my $tr_id_str = join( ',', keys(%tr_hash) );

  my $sth =
    $self->prepare( "SELECT transcript_id, exon_id, rank "
      . "FROM exon_transcript "
      . "WHERE transcript_id IN ($tr_id_str)" );

  $sth->execute();

  my ( $tr_id, $ex_id, $rank );
  $sth->bind_columns( \( $tr_id, $ex_id, $rank ) );

  my %ex_tr_hash;

  while ( $sth->fetch() ) {
    $ex_tr_hash{$ex_id} ||= [];
    push( @{ $ex_tr_hash{$ex_id} }, [ $tr_hash{$tr_id}, $rank ] );
  }

  my $ea    = $self->db()->get_ExonAdaptor();
  my $exons = $ea->fetch_all_by_Slice_constraint(
    $ext_slice,
    sprintf( "e.exon_id IN (%s)",
      join( ',', sort { $a <=> $b } keys(%ex_tr_hash) ) ) );

  # move exons onto transcript slice, and add them to transcripts
  foreach my $ex ( @{$exons} ) {
    my $new_ex;
    if ( $slice != $ext_slice ) {
      $new_ex = $ex->transfer($slice);
      if ( !defined($new_ex) ) {
        throw("Unexpected. "
            . "Exon could not be transfered onto Transcript slice." );
      }
    } else {
      $new_ex = $ex;
    }

    foreach my $row ( @{ $ex_tr_hash{ $new_ex->dbID() } } ) {
      my ( $tr, $rank ) = @{$row};
      $tr->add_Exon( $new_ex, $rank );
    }
  }

  my $tla = $self->db()->get_TranslationAdaptor();

  # load all of the translations at once
  $tla->fetch_all_by_Transcript_list($transcripts);

  return $transcripts;
} ## end sub fetch_all_by_Slice


=head2 fetch_all_by_external_name

  Arg [1]    : String $external_name
               An external identifier of the transcript to be obtained
  Arg [2]    : (optional) String $external_db_name
               The name of the external database from which the
               identifier originates.
  Arg [3]    : Boolean override. Force SQL regex matching for users
               who really do want to find all 'NM%'
  Example    : my @transcripts =
                  @{ $tr_adaptor->fetch_all_by_external_name( 'NP_065811.1') };
               my @more_transcripts = 
                  @{$tr_adaptor->fetch_all_by_external_name( 'NP_0658__._')};
  Description: Retrieves all transcripts which are associated with
               an external identifier such as a GO term, Swissprot
               identifer, etc.  Usually there will only be a single
               transcript returned in the list reference, but not
               always.  Transcripts are returned in their native
               coordinate system, i.e. the coordinate system in which
               they are stored in the database.  If they are required
               in another coordinate system the Transcript::transfer or
               Transcript::transform method can be used to convert them.
               If no transcripts with the external identifier are found,
               a reference to an empty list is returned.
               SQL wildcards % and _ are supported in the $external_name
               but their use is somewhat restricted for performance reasons.
               Users that really do want % and _ in the first three characters
               should use argument 3 to prevent optimisations
  Returntype : listref of Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_external_name {
  my ( $self, $external_name, $external_db_name, $override) = @_;

  my $entryAdaptor = $self->db->get_DBEntryAdaptor();

  my @ids =
    $entryAdaptor->list_transcript_ids_by_extids( $external_name,
                                                  $external_db_name, $override );

  my @features = @{ $self->fetch_all_by_dbID_list( \@ids ) };
  my @reference = grep { $_->slice()->is_reference() } @features;
  my @non_reference = grep { ! $_->slice()->is_reference() } @features;
  return [ @reference, @non_reference ];
}

=head2 fetch_all_by_GOTerm

  Arg [1]   : Bio::EnsEMBL::OntologyTerm
              The GO term for which transcripts should be fetched.

  Example:  @transcripts = @{
              $transcript_adaptor->fetch_all_by_GOTerm(
                $go_adaptor->fetch_by_accession('GO:0030326') ) };

  Description   : Retrieves a list of transcripts that are
                  associated with the given GO term, or with any of
                  its descendent GO terms.  The transcripts returned
                  are in their native coordinate system, i.e. in
                  the coordinate system in which they are stored
                  in the database.  If another coordinate system
                  is required then the Transcript::transfer or
                  Transcript::transform method can be used.

  Return type   : listref of Bio::EnsEMBL::Transcript
  Exceptions    : Throws of argument is not a GO term
  Caller        : general
  Status        : Stable

=cut

sub fetch_all_by_GOTerm {
  my ( $self, $term ) = @_;

  assert_ref( $term, 'Bio::EnsEMBL::OntologyTerm' );
  if ( $term->ontology() ne 'GO' ) {
    throw('Argument is not a GO term');
  }

  my $entryAdaptor = $self->db->get_DBEntryAdaptor();

  my %unique_dbIDs;
  foreach my $accession ( map { $_->accession() }
                          ( $term, @{ $term->descendants() } ) )
  {
    my @ids =
      $entryAdaptor->list_transcript_ids_by_extids( $accession, 'GO' );
    foreach my $dbID (@ids) { $unique_dbIDs{$dbID} = 1 }
  }

  my @result = @{
    $self->fetch_all_by_dbID_list(
                              [ sort { $a <=> $b } keys(%unique_dbIDs) ]
    ) };

  return \@result;
} ## end sub fetch_all_by_GOTerm

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
  my ( $self, $accession ) = @_;

  if ( $accession !~ /^GO:/ ) {
    throw('Argument is not a GO term accession');
  }

  my $goAdaptor =
    Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Ontology',
                                         'OntologyTerm' );

  my $term = $goAdaptor->fetch_by_accession($accession);

  return $self->fetch_all_by_GOTerm($term);
}

=head2 fetch_by_display_label

  Arg [1]    : String $label - display label of transcript to fetch
  Example    : my $tr = $tr_adaptor->fetch_by_display_label("BRCA2");
  Description: Returns the transcript which has the given display label or
               undef if there is none. If there are more than 1, only the first
               is reported.
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_display_label {
  my $self = shift;
  my $label = shift;

  my $constraint = "x.display_label = ? AND t.is_current = 1";

  $self->bind_param_generic_fetch($label,SQL_VARCHAR);

  my ($transcript) = @{ $self->generic_fetch($constraint) };

  return $transcript;
}


=head2 fetch_all_by_exon_stable_id

  Arg [1]    : String $stable_id 
               The stable id of an exon in a transcript
  Example    : my $tr = $tr_adaptor->fetch_all_by_exon_stable_id
                  ('ENSE00000309301');
  Description: Retrieves a list of transcripts via an exon stable id.
  Returntype : Listref of Bio::EnsEMBL::Transcript objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_exon_stable_id {
  my ($self, $stable_id) = @_;

  my @trans ;

  my $sth = $self->prepare(qq(
      SELECT t.transcript_id 
      FROM exon_transcript et, exon e, transcript t
      WHERE e.exon_id = et.exon_id
      AND et.transcript_id = t.transcript_id
      AND e.stable_id = ?
      AND t.is_current = 1
  ));
  
  $sth->bind_param(1, $stable_id, SQL_VARCHAR);
  $sth->execute();

  while( my $id = $sth->fetchrow_array ) {
    my $transcript = $self->fetch_by_dbID($id);
    push(@trans, $transcript) if $transcript;
  }

  if (!@trans) {
    return undef;
  }

  return \@trans;
}

=head2 fetch_all_by_source

  Arg [1]    : String $source
               listref of $sources
               The source of the transcript to retrieve. You can have as an argument a reference
               to a list of sources
  Example    : $transcripts = $transcript_adaptor->fetch_all_by_source('havana'); 
               $transcripts = $transcript_adaptor->fetch_all_by_source(['ensembl', 'vega']);
  Description: Retrieves an array reference of transcript objects from the database via its source or sources.
               The transcript will be retrieved in its native coordinate system (i.e.
               in the coordinate system it is stored in the database). It may
               be converted to a different coordinate system through a call to
               transform() or transfer(). If the gene or exon is not found
               undef is returned instead.
  Returntype  : listref of Bio::EnsEMBL::Transcript
  Exceptions : if we cant get the gene in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_source {
  my ($self, $source) = @_;
  my @transcripts = @{$self->generic_fetch($self->source_constraint($source))};
  return \@transcripts;
}

=head2 source_constraint 

  Arg [1]    : String $source
               listref of $sources
               The source of the transcript to retrieve. You can have as an argument a reference
               to a list of sources
  Description: Used internally to generate a SQL constraint to restrict a transcript query by source
  Returntype  : String
  Exceptions : If source is not supplied
  Caller     : general
  Status     : Stable

=cut

sub source_constraint {
  my ($self, $sources, $inline_variables) = @_;
  my $constraint = "t.is_current = 1";
  my $in_statement = $self->generate_in_constraint($sources, 't.source', SQL_VARCHAR, $inline_variables);
  $constraint .= " and $in_statement";
  return $constraint;
}

=head2 count_all_by_source

  Arg [1]     : String $source
                listref of $source
                The source of the transcript to retrieve. You can have as an argument a reference
                to a list of sources
  Example     : $cnt = $transcript_adaptor->count_all_by_source('ensembl'); 
                $cnt = $transcript_adaptor->count_all_by_source(['havana', 'vega']);
  Description : Retrieves count of transcript objects from the database via its source or sources.
  Returntype  : integer
  Caller      : general
  Status      : Stable

=cut

sub count_all_by_source {
  my ($self, $source) = @_;
  return $self->generic_count($self->source_constraint($source));
}

=head2 count_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to count transcripts on.
  Arg [2]    : (optional) biotype(s) string or arrayref of strings
                the biotype of the features to count.
  Arg [1]    : (optional) string $source
               the source name of the features to count.
  Example    : $cnt = $transcript_adaptor->count_all_by_Slice();
  Description: Method to count transcripts on a given slice, filtering by biotype and source
  Returntype : integer
  Exceptions : thrown if exon cannot be placed on transcript slice
  Status     : Stable
  Caller     : general
=cut

sub count_all_by_Slice {
  my ($self, $slice, $biotype, $source) = @_;

  my $constraint = 't.is_current = 1';
  if (defined($source)) {
        $constraint .= " and t.source = '$source'";
  }
  if (defined($biotype)) {
        $constraint .= " and " . $self->biotype_constraint($biotype);
  }

  return $self->count_by_Slice_constraint($slice, $constraint);
}

=head2 fetch_all_by_biotype 

  Arg [1]    : String $biotype 
               listref of $biotypes
               The biotype of the transcript to retrieve. You can have as an argument a reference
               to a list of biotypes
  Example    : $gene = $transcript_adaptor->fetch_all_by_biotype('protein_coding'); 
               $gene = $transcript_adaptor->fetch_all_by_biotypes(['protein_coding', 'sRNA', 'miRNA']);
  Description: Retrieves an array reference of transcript objects from the database via its biotype or biotypes.
               The transcript will be retrieved in its native coordinate system (i.e.
               in the coordinate system it is stored in the database). It may
               be converted to a different coordinate system through a call to
               transform() or transfer(). If the gene or exon is not found
               undef is returned instead.
  Returntype  : listref of Bio::EnsEMBL::Transcript
  Exceptions : if we cant get the gene in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_biotype {
  my ($self, $biotype) = @_;
  my @transcripts = @{$self->generic_fetch($self->biotype_constraint($biotype))};
  return \@transcripts;
}

=head2 biotype_constraint 

  Arg [1]    : String $biotypes 
               listref of $biotypes
               The biotype of the transcript to retrieve. You can have as an argument a reference
               to a list of biotypes
  Description: Used internally to generate a SQL constraint to restrict a transcript query by biotype
  Returntype  : String
  Exceptions : If biotype is not supplied
  Caller     : general
  Status     : Stable

=cut

sub biotype_constraint {
  my ($self, $biotypes, $inline_variables) = @_;
  my $constraint = "t.is_current = 1";
  my $in_statement = $self->generate_in_constraint($biotypes, 't.biotype', SQL_VARCHAR, $inline_variables);
  $constraint .= " and $in_statement";
  return $constraint;
}

=head2 count_all_by_biotype 

  Arg [1]     : String $biotype 
                listref of $biotypes
                The biotype of the transcript to retrieve. You can have as an argument a reference
                to a list of biotypes
  Example     : $cnt = $transcript_adaptor->count_all_by_biotype('protein_coding'); 
                $cnt = $transcript_adaptor->count_all_by_biotypes(['protein_coding', 'sRNA', 'miRNA']);
  Description : Retrieves count of transcript objects from the database via its biotype or biotypes.
  Returntype  : integer
  Caller      : general
  Status      : Stable

=cut

sub count_all_by_biotype {
  my ($self, $biotype) = @_;
  return $self->generic_count($self->biotype_constraint($biotype));
}

=head2 store

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript
               The transcript to be written to the database
  Arg [2]    : Int $gene_dbID
               The identifier of the gene that this transcript is associated 
               with
  Arg [3]    : DEPRECATED (optional) Int $analysis_id
               The analysis_id to use when storing this gene. This is for 
               backward compatibility only and used to fall back to the gene
               analysis_id if no analysis object is attached to the transcript
               (which you should do for new code).
  Arg [4]    : prevent coordinate recalculation if you are persisting 
               transcripts with this gene
  Example    : $transID = $tr_adaptor->store($transcript, $gene->dbID);
  Description: Stores a transcript in the database and returns the new
               internal identifier for the stored transcript.
  Returntype : Int 
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub store {
  my ( $self, $transcript, $gene_dbID, $analysis_id, $skip_recalculating_coordinates ) = @_;

  if (    !ref($transcript)
       || !$transcript->isa('Bio::EnsEMBL::Transcript') )
  {
    throw("$transcript is not a EnsEMBL transcript - not storing");
  }

  my $db = $self->db();

  if ( $transcript->is_stored($db) ) {
    return $transcript->dbID();
  }

  # Force lazy-loading of exons and ensure coords are correct.
  # If we have been told not to do this then skip doing this
  # and we assume the user knows what they are doing. You have been
  # warned
  if(! $skip_recalculating_coordinates) {
    $transcript->recalculate_coordinates();
  }

  my $is_current = ( defined( $transcript->is_current() )
                     ? $transcript->is_current()
                     : 1 );

  # store analysis
  my $analysis = $transcript->analysis();
  my $new_analysis_id;

  if ($analysis) {
    if ( $analysis->is_stored($db) ) {
      $new_analysis_id = $analysis->dbID;
    } else {
      $new_analysis_id = $db->get_AnalysisAdaptor->store($analysis);
    }
  } elsif ($analysis_id) {
    # Fall back to analysis passed in (usually from gene) if analysis
    # wasn't set explicitely for the transcript. This is deprectated
    # though.
    warning(   "You should explicitely attach "
             . "an analysis object to the Transcript. "
             . "Will fall back to Gene analysis, "
             . "but this behaviour is deprecated." );
    $new_analysis_id = $analysis_id;
  } else {
    throw("Need an analysis_id to store the Transcript.");
  }

  #
  # Store exons - this needs to be done before the possible transfer
  # of the transcript to another slice (in _prestore()).  Transfering
  # results in copies being made of the exons and we need to preserve
  # the object identity of the exons so that they are not stored twice
  # by different transcripts.
  #
  my $exons       = $transcript->get_all_Exons();
  my $exonAdaptor = $db->get_ExonAdaptor();
  foreach my $exon ( @{$exons} ) {
    $exonAdaptor->store($exon);
  }

  my $original_translation = $transcript->translation();
  my $original             = $transcript;
  my $seq_region_id;
  ( $transcript, $seq_region_id ) = $self->_pre_store($transcript);

  # First store the transcript without a display xref.  The display xref
  # needs to be set after xrefs are stored which needs to happen after
  # transcript is stored.

  #
  # Store transcript
  #

#  my $store_transcript_sql = 
#    sprintf "INSERT INTO transcript SET gene_id = ?, analysis_id = ?, seq_region_id = ?, seq_region_start = ?, seq_region_end = ?, seq_region_strand = ?,%s biotype = ?, description = ?, is_current = ?, canonical_translation_id = ?", ($self->schema_version > 74)?" source = ?,":'';

  my @columns = qw(
            gene_id
            analysis_id
            seq_region_id
            seq_region_start
            seq_region_end
            seq_region_strand
  );

  push @columns, 'source' if ($self->schema_version > 74);

  push @columns, qw(
            biotype
            description
            is_current
            canonical_translation_id
  );

  my @canned_columns;
  my @canned_values;

  if ( defined( $transcript->stable_id() ) ) {
      push @columns, 'stable_id', 'version';

      my $created = $self->db->dbc->from_seconds_to_date($transcript->created_date());
      my $modified = $self->db->dbc->from_seconds_to_date($transcript->modified_date());

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
  my $store_transcript_sql = qq(
        INSERT INTO transcript ( $columns ) VALUES ( $values )
  );

  my $tst = $self->prepare($store_transcript_sql);
  my $i = 0;
  $tst->bind_param( ++$i,  $gene_dbID,                 SQL_INTEGER );
  $tst->bind_param( ++$i,  $new_analysis_id,           SQL_INTEGER );
  $tst->bind_param( ++$i,  $seq_region_id,             SQL_INTEGER );
  $tst->bind_param( ++$i,  $transcript->start(),       SQL_INTEGER );
  $tst->bind_param( ++$i,  $transcript->end(),         SQL_INTEGER );
  $tst->bind_param( ++$i,  $transcript->strand(),      SQL_TINYINT );

  $self->schema_version > 74 and 
    $tst->bind_param( ++$i,  $transcript->source(),      SQL_VARCHAR );

  $tst->bind_param( ++$i,  $transcript->biotype(),     SQL_VARCHAR );
  $tst->bind_param( ++$i,  $transcript->description(), SQL_LONGVARCHAR );
  $tst->bind_param( ++$i, $is_current,                SQL_TINYINT );

  # If the transcript has a translation, this is updated later:
  $tst->bind_param( ++$i, undef, SQL_INTEGER );

  if ( defined( $transcript->stable_id() ) ) {

    $tst->bind_param( ++$i, $transcript->stable_id(), SQL_VARCHAR );
    $tst->bind_param( ++$i, $transcript->version(),   SQL_INTEGER );
  }

  $tst->execute();
  $tst->finish();

  my $transc_dbID = $self->last_insert_id('transcript_id', undef, 'transcript');

  #
  # Store translation
  #

  my $alt_translations =
    $transcript->get_all_alternative_translations();
  my $translation = $transcript->translation();

  if ( defined($translation) ) {
    # Make sure that the start and end exon are set correctly.
    my $start_exon = $translation->start_Exon();
    my $end_exon   = $translation->end_Exon();

    if ( !defined($start_exon) ) {
      throw("Translation does not define a start exon.");
    }

    if ( !defined($end_exon) ) {
      throw("Translation does not defined an end exon.");
    }

    # If the dbID is not set, this means the exon must have been a
    # different object in memory than the the exons of the transcript.
    # Try to find the matching exon in all of the exons we just stored.
    if ( !defined( $start_exon->dbID() ) ) {
      my $key = $start_exon->hashkey();
      ($start_exon) = grep { $_->hashkey() eq $key } @$exons;

      if ( defined($start_exon) ) {
        $translation->start_Exon($start_exon);
      } else {
        throw(   "Translation's start_Exon does not appear "
               . "to be one of the exons in "
               . "its associated Transcript" );
      }
    }

    if ( !defined( $end_exon->dbID() ) ) {
      my $key = $end_exon->hashkey();
      ($end_exon) = grep { $_->hashkey() eq $key } @$exons;

      if ( defined($end_exon) ) {
        $translation->end_Exon($end_exon);
      } else {
        throw(   "Translation's end_Exon does not appear "
               . "to be one of the exons in "
               . "its associated Transcript." );
      }
    }

    my $old_dbid = $translation->dbID();
    $db->get_TranslationAdaptor()->store( $translation, $transc_dbID );

    # Need to update the canonical_translation_id for this transcript.

    my $sth = $self->prepare(
      q(
      UPDATE transcript
      SET canonical_translation_id = ?
      WHERE transcript_id = ?)
    );

    $sth->bind_param( 1, $translation->dbID(), SQL_INTEGER );
    $sth->bind_param( 2, $transc_dbID,         SQL_INTEGER );

    $sth->execute();

    # Set values of the original translation, we may have copied it when
    # we transformed the transcript.
    $original_translation->dbID( $translation->dbID() );
    $original_translation->adaptor( $translation->adaptor() );
  } ## end if ( defined($translation...))

  #
  # Store the alternative translations, if there are any.
  #

  if ( defined($alt_translations)
       && scalar( @{$alt_translations} ) > 0 )
  {
    foreach my $alt_translation ( @{$alt_translations} ) {
      my $start_exon = $alt_translation->start_Exon();
      my $end_exon   = $alt_translation->end_Exon();

      if ( !defined($start_exon) ) {
        throw("Translation does not define a start exon.");
      } elsif ( !defined($end_exon) ) {
        throw("Translation does not defined an end exon.");
      }

      if ( !defined( $start_exon->dbID() ) ) {
        my $key = $start_exon->hashkey();
        ($start_exon) = grep { $_->hashkey() eq $key } @{$exons};

        if ( defined($start_exon) ) {
          $alt_translation->start_Exon($start_exon);
        } else {
          throw(   "Translation's start_Exon does not appear "
                 . "to be one of the exon in"
                 . "its associated Transcript" );
        }
      }
      if ( !defined( $end_exon->dbID() ) ) {
        my $key = $end_exon->hashkey();
        ($end_exon) = grep { $_->hashkey() eq $key } @$exons;

        if ( defined($end_exon) ) {
          $alt_translation->end_Exon($end_exon);
        } else {
          throw(   "Translation's end_Exon does not appear "
                 . "to be one of the exons in "
                 . "its associated Transcript." );
        }
      }

      $db->get_TranslationAdaptor()
        ->store( $alt_translation, $transc_dbID );
    } ## end foreach my $alt_translation...
  } ## end if ( defined($alt_translations...))

  #
  # Store the xrefs/object xref mapping.
  #
  my $dbEntryAdaptor = $db->get_DBEntryAdaptor();

  foreach my $dbe ( @{ $transcript->get_all_DBEntries() } ) {
    $dbEntryAdaptor->store( $dbe, $transc_dbID, "Transcript", 1 );
  }

  #
  # Update transcript to point to display xref if it is set.
  #
  if ( my $dxref = $transcript->display_xref() ) {
    my $dxref_id;

    if ( $dxref->is_stored($db) ) {
      $dxref_id = $dxref->dbID();
    } else {
      $dxref_id = $dbEntryAdaptor->exists($dxref);
    }

    if ( defined($dxref_id) ) {
      my $sth =
        $self->prepare(   "UPDATE transcript "
                        . "SET display_xref_id = ? "
                        . "WHERE transcript_id = ?" );
      $sth->bind_param( 1, $dxref_id,    SQL_INTEGER );
      $sth->bind_param( 2, $transc_dbID, SQL_INTEGER );
      $sth->execute();
      $dxref->dbID($dxref_id);
      $dxref->adaptor($dbEntryAdaptor);
      $sth->finish();
    } else {
      warning(sprintf(
                     "Display_xref %s:%s is not stored in database.\n"
                       . "Not storing relationship to this transcript.",
                     $dxref->dbname(), $dxref->display_id() ) );
      $dxref->dbID(undef);
      $dxref->adaptor(undef);
    }
  } ## end if ( my $dxref = $transcript...)

  #
  # Link transcript to exons in exon_transcript table
  #
  my $etst = $self->prepare(
             "INSERT INTO exon_transcript (exon_id,transcript_id,rank) "
               . "VALUES (?,?,?)" );
  my $rank = 1;
  foreach my $exon ( @{ $transcript->get_all_Exons } ) {
    $etst->bind_param( 1, $exon->dbID,  SQL_INTEGER );
    $etst->bind_param( 2, $transc_dbID, SQL_INTEGER );
    $etst->bind_param( 3, $rank,        SQL_INTEGER );
    $etst->execute();
    $rank++;
  }

  $etst->finish();

  # Now the supporting evidence
  my $tsf_adaptor = $db->get_TranscriptSupportingFeatureAdaptor();
  $tsf_adaptor->store( $transc_dbID,
                       $transcript->get_all_supporting_features() );

  # store transcript attributes if there are any
  my $attr_adaptor = $db->get_AttributeAdaptor();

  $attr_adaptor->store_on_Transcript( $transc_dbID,
                                    $transcript->get_all_Attributes() );

  # store the IntronSupportingEvidence features
  my $ise_adaptor = $db->get_IntronSupportingEvidenceAdaptor();
  my $intron_supporting_evidence = $transcript->get_all_IntronSupportingEvidence();
  foreach my $ise (@{$intron_supporting_evidence}) {
    $ise_adaptor->store($ise);
    $ise_adaptor->store_transcript_linkage($ise, $transcript, $transc_dbID);
  }

  # Update the original transcript object - not the transfered copy that
  # we might have created.
  $original->dbID($transc_dbID);
  $original->adaptor($self);

  return $transc_dbID;
} ## end sub store


=head2 get_Interpro_by_transid

  Arg [1]    : String $trans_stable_id
               The stable if of the transcript to obtain
  Example    : @i = $tr_adaptor->get_Interpro_by_transid($trans->stable_id()); 
  Description: Gets interpro accession numbers by transcript stable id.
               A hack really - we should have a much more structured 
               system than this.
  Returntype : listref of strings (Interpro_acc:description)
  Exceptions : none 
  Caller     : domainview? , GeneView
  Status     : Stable

=cut

sub get_Interpro_by_transid {
   my ($self,$trans_stable_id) = @_;

   my $straight_join = $self->_can_straight_join ? 'STRAIGHT_JOIN' : '';
   my $sth = $self->prepare(qq(
      SELECT  ${straight_join} i.interpro_ac, x.description
      FROM    transcript t,
              translation tl,
              protein_feature pf,
	      interpro i,
              xref x
      WHERE   t.stable_id = ?
      AND     tl.transcript_id = t.transcript_id
      AND     tl.translation_id = pf.translation_id
      AND     i.id = pf.hit_name
      AND     i.interpro_ac = x.dbprimary_acc
      AND     t.is_current = 1
  ));

  $sth->bind_param(1, $trans_stable_id, SQL_VARCHAR);
  $sth->execute();

  my @out;
  my %h;
  while( (my $arr = $sth->fetchrow_arrayref()) ) {
     if( $h{$arr->[0]} ) { next; }
     $h{$arr->[0]}=1;
     my $string = $arr->[0] .":".$arr->[1];
     push(@out,$string);
  }

  return \@out;
}

=head2 is_Transcript_canonical()

  Arg [1]     : Bio::EnsEMBL::Transcript $transcript
                The transcript to query with
  Example     : $tr_adaptor->is_Transcript_canonical($transcript);
  Description : Returns a boolean if the given transcript is considered
                canonical with respect to a gene
  Returntype  : Boolean
  Exceptions  : None
  Caller      : Bio::EnsEMBL::Transcript
  Status      : Beta
  

=cut

sub is_Transcript_canonical {
  my ($self, $transcript) = @_;
  return $self->dbc()->sql_helper()->execute_single_result(
    -SQL => 'select count(*) from gene where canonical_transcript_id =?', 
    -PARAMS => [$transcript->dbID()]
  );
}


=head2 remove

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript
               The transcript to remove from the database
  Arg [2]    : Boolean, update Gene coordinates after removal. WARNING: this does not alter any other copies of the 
               gene currently in memory. Other copies will retain their original coordinates. Either refetch them
               or go directly through Gene->remove_Transcript first, then remove the Transcript here.
  Example    : $tr_adaptor->remove($transcript);
  Description: Removes a transcript completely from the database, and all
               associated information.
               This method is usually called by the GeneAdaptor::remove method
               because this method will not preform the removal of genes
               which are associated with this transcript. Do not call this
               method directly unless you know there are no genes associated
               with the transcript!
  Returntype : none
  Exceptions : throw on incorrect arguments
               warning if transcript is not in this database
  Caller     : GeneAdaptor::remove
  Status     : Stable

=cut

sub remove {
  my $self = shift;
  my $transcript = shift;
  my $update = shift;
  if(!ref($transcript) || !$transcript->isa('Bio::EnsEMBL::Transcript')) {
    throw("Bio::EnsEMBL::Transcript argument expected");
  }

  # sanity check: make sure nobody tries to slip past a prediction transcript
  # which inherits from transcript but actually uses different tables
  if($transcript->isa('Bio::EnsEMBL::PredictionTranscript')) {
    throw("TranscriptAdaptor can only remove Transcripts " .
          "not PredictionTranscripts");
  }

  if ( !$transcript->is_stored($self->db()) ) {
    warning("Cannot remove transcript ". $transcript->dbID .". Is not stored ".
            "in this database.");
    return;
  }

  # remove the supporting features of this transcript

  my $prot_adp = $self->db->get_ProteinAlignFeatureAdaptor;
  my $dna_adp = $self->db->get_DnaAlignFeatureAdaptor;

  my $sfsth = $self->prepare("SELECT feature_type, feature_id  " .
                             "FROM transcript_supporting_feature " .
                             "WHERE transcript_id = ?");

  $sfsth->bind_param(1, $transcript->dbID, SQL_INTEGER);
  $sfsth->execute();

  # statements to check for shared align_features
  my $sth1 = $self->prepare("SELECT count(*) FROM supporting_feature " .
			    "WHERE feature_type = ? AND feature_id = ?");
  my $sth2 = $self->prepare("SELECT count(*) " .
                            "FROM transcript_supporting_feature " .
			    "WHERE feature_type = ? AND feature_id = ?");

  SUPPORTING_FEATURE:
  while(my ($type, $feature_id) = $sfsth->fetchrow()){
    
    # only remove align_feature if this is the last reference to it
    $sth1->bind_param(1, $type, SQL_VARCHAR);
    $sth1->bind_param(2, $feature_id, SQL_INTEGER);
    $sth1->execute;
    $sth2->bind_param(1, $type, SQL_VARCHAR);
    $sth2->bind_param(2, $feature_id, SQL_INTEGER);
    $sth2->execute;
    my ($count1) = $sth1->fetchrow;
    my ($count2) = $sth2->fetchrow;
    if ($count1 + $count2 > 1) {
      #warn "transcript: shared feature, not removing $type|$feature_id\n";
      next SUPPORTING_FEATURE;
    }
    
    #warn "transcript: removing $type|$feature_id\n";
  
    if($type eq 'protein_align_feature'){
      my $f = $prot_adp->fetch_by_dbID($feature_id);
      $prot_adp->remove($f);
    }
    elsif($type eq 'dna_align_feature'){
      my $f = $dna_adp->fetch_by_dbID($feature_id);
      $dna_adp->remove($f);
    }
    else {
      warning("Unknown supporting feature type $type. Not removing feature.");
    }
  }
  $sfsth->finish();
  $sth1->finish();
  $sth2->finish();

  # delete the association to supporting features

  $sfsth = $self->prepare("DELETE FROM transcript_supporting_feature WHERE transcript_id = ?");
  $sfsth->bind_param(1, $transcript->dbID, SQL_INTEGER);
  $sfsth->execute();
  $sfsth->finish();
  
  # delete the associated IntronSupportingEvidence and if the ISE had no more
  # linked transcripts remove it
  my $ise_adaptor = $self->db->get_IntronSupportingEvidenceAdaptor();
  foreach my $ise (@{$transcript->get_all_IntronSupportingEvidence()}) {
    $ise_adaptor->remove_transcript_linkage($ise, $transcript);
    if(! $ise->has_linked_transcripts()) {
      $ise_adaptor->remove($ise);
    }
  }

  # remove all xref linkages to this transcript

  my $dbeAdaptor = $self->db->get_DBEntryAdaptor();
  foreach my $dbe (@{$transcript->get_all_DBEntries}) {
    $dbeAdaptor->remove_from_object($dbe, $transcript, 'Transcript');
  }

  # remove the attributes associated with this transcript
  my $attrib_adp = $self->db->get_AttributeAdaptor;  
  $attrib_adp->remove_from_Transcript($transcript);

  # remove the translation associated with this transcript

  my $translationAdaptor = $self->db->get_TranslationAdaptor();
  if( defined($transcript->translation()) ) {
    $translationAdaptor->remove( $transcript->translation );
  }

  # remove exon associations to this transcript

  my $exonAdaptor = $self->db->get_ExonAdaptor();
  foreach my $exon ( @{$transcript->get_all_Exons()} ) {
    # get the number of transcript references to this exon
    # only remove the exon if this is the last transcript to
    # reference it

    my $sth = $self->prepare( "SELECT count(*)
                               FROM   exon_transcript
                               WHERE  exon_id = ?" );
    $sth->bind_param(1, $exon->dbID, SQL_INTEGER);
    $sth->execute();
    my ($count) = $sth->fetchrow_array();
    $sth->finish();

    if($count == 1){
      $exonAdaptor->remove( $exon );
    }
  }

  my $sth = $self->prepare( "DELETE FROM exon_transcript
                             WHERE transcript_id = ?" );
  $sth->bind_param(1, $transcript->dbID, SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  my $gene = $transcript->get_Gene;

  $sth = $self->prepare( "DELETE FROM transcript
                          WHERE transcript_id = ?" );
  $sth->bind_param(1, $transcript->dbID, SQL_INTEGER);
  $sth->execute();
  $sth->finish();

  if ($update) {
    $gene->remove_Transcript($transcript);
  }

  $transcript->dbID(undef);
  $transcript->adaptor(undef);

  return;
}


=head2 update

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript
               The transcript to update
  Example    : $tr_adaptor->update($transcript);
  Description: Updates a transcript in the database.
  Returntype : None
  Exceptions : thrown if the $transcript is not a Bio::EnsEMBL::Transcript.
               warn if the method is called on a transcript that does not exist 
               in the database.
               Should warn if trying to update the number of attached exons, but
               this is a far more complex process and is not yet implemented.
  Caller     : general
  Status     : Stable

=cut

sub update {
  my ( $self, $transcript ) = @_;

  if (    !defined($transcript)
       || !ref($transcript)
       || !$transcript->isa('Bio::EnsEMBL::Transcript') )
  {
    throw("Must update a transcript object, not a $transcript");
  }

  my $update_transcript_sql = 
    sprintf "UPDATE transcript SET analysis_id = ?, display_xref_id = ?, description = ?,%s biotype = ?, is_current = ?, canonical_translation_id = ? WHERE transcript_id = ?", ($self->schema_version > 74)?" source = ?,":'';

  my $display_xref = $transcript->display_xref();
  my $display_xref_id;

  if ( defined($display_xref) && $display_xref->dbID() ) {
    $display_xref_id = $display_xref->dbID();
  } else {
    $display_xref_id = undef;
  }

  my $sth = $self->prepare($update_transcript_sql);
  my $i = 0;
  $sth->bind_param( ++$i, $transcript->analysis()->dbID(), SQL_INTEGER );
  $sth->bind_param( ++$i, $display_xref_id, SQL_INTEGER );
  $sth->bind_param( ++$i, $transcript->description(), SQL_LONGVARCHAR );

  $self->schema_version > 74 and 
    $sth->bind_param( ++$i,  $transcript->source(),      SQL_VARCHAR );

  $sth->bind_param( ++$i, $transcript->biotype(),     SQL_VARCHAR );
  $sth->bind_param( ++$i, $transcript->is_current(),  SQL_TINYINT );
  $sth->bind_param( ++$i, (
                      defined( $transcript->translation() )
                      ? $transcript->translation()->dbID()
                      : undef ),
                    SQL_INTEGER );
  $sth->bind_param( ++$i, $transcript->dbID(), SQL_INTEGER );

  $sth->execute();
} ## end sub update


=head2 list_dbIDs

  Example    : @transcript_ids = @{ $t_adaptor->list_dbIDs };
  Description: Gets a list of internal ids for all transcripts in the db.
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.  Returntype : Listref of Ints
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub list_dbIDs {
   my ($self, $ordered) = @_;

   return $self->_list_dbIDs("transcript",undef, $ordered);
}


=head2 list_stable_ids

  Example    : @stable_trans_ids = @{ $transcript_adaptor->list_stable_ids };
  Description: Gets a list of stable ids for all transcripts in the current
               database.
  Returntype : Listref of Strings
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub list_stable_ids {
   my ($self) = @_;

   return $self->_list_dbIDs("transcript", "stable_id");
}


#_objs_from_sth

#  Arg [1]    : StatementHandle $sth
#  Arg [2]    : Bio::EnsEMBL::AssemblyMapper $mapper
#  Arg [3]    : Bio::EnsEMBL::Slice $dest_slice
#  Description: PROTECTED implementation of abstract superclass method.
#               Responsible for the creation of Transcripts.
#  Returntype : Listref of Bio::EnsEMBL::Transcripts in target coord system
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

  my @transcripts;
  my %analysis_hash;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;

  my (
    $transcript_id,   $seq_region_id,      $seq_region_start,
    $seq_region_end,  $seq_region_strand,  $analysis_id,
    $gene_id,         $is_current,         $stable_id,
    $version,         $created_date,       $modified_date,
    $description,     $biotype,
    $external_db,     $external_status,    $external_db_name,
    $display_xref_id, $xref_display_label, $xref_primary_acc,
    $xref_version,    $xref_description,   $xref_info_type,
    $xref_info_text,  $external_release,   $source
  );

  if ($self->schema_version() > 74) {
    $sth->bind_columns(
           \(
       $transcript_id,   $seq_region_id,      $seq_region_start,
       $seq_region_end,  $seq_region_strand,  $analysis_id,
       $gene_id,         $is_current,         $stable_id,
       $version,         $created_date,       $modified_date,
       $description,     $biotype,
       $external_db,     $external_status,    $external_db_name,
       $display_xref_id, $xref_display_label, $xref_primary_acc,
       $xref_version,    $xref_description,   $xref_info_type,
       $xref_info_text,  $external_release,   $source
      ) );
  } else {
    $sth->bind_columns(
           \(
       $transcript_id,   $seq_region_id,      $seq_region_start,
       $seq_region_end,  $seq_region_strand,  $analysis_id,
       $gene_id,         $is_current,         $stable_id,
       $version,         $created_date,       $modified_date,
       $description,     $biotype,
       $external_db,     $external_status,    $external_db_name,
       $display_xref_id, $xref_display_label, $xref_primary_acc,
       $xref_version,    $xref_description,   $xref_info_type,
       $xref_info_text,  $external_release
      ) );    
  }

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
      $display_xref = Bio::EnsEMBL::DBEntry->new_fast( {
          'dbID'            => $display_xref_id,
          'adaptor'         => $dbEntryAdaptor,
          'display_id'      => $xref_display_label,
          'primary_id'      => $xref_primary_acc,
          'version'         => $xref_version,
          'description'     => $xref_description,
          'release'         => $external_release,
          'dbname'          => $external_db,
          'db_display_name' => $external_db_name,
          'info_type'       => $xref_info_type,
          'info_text'       => $xref_info_text
      });
      $display_xref->status($external_status);
    }

    # Finally, create the new Transcript.
    my $params = 
      {
       'analysis'              => $analysis,
       'biotype'               => $biotype,
       'start'                 => $seq_region_start,
       'end'                   => $seq_region_end,
       'strand'                => $seq_region_strand,
       'adaptor'               => $self,
       'slice'                 => $slice,
       'dbID'                  => $transcript_id,
       'stable_id'             => $stable_id,
       'version'               => $version,
       'created_date'          => $created_date || undef,
       'modified_date'         => $modified_date || undef,
       'description'           => $description,
       'external_name'         => $xref_display_label,

       'external_status'       => $external_status,
       'external_display_name' => $external_db_name,
       'external_db'           => $external_db,
       'display_xref'          => $display_xref,
       'is_current'            => $is_current,
       'edits_enabled'         => 1
      };
    
    $self->schema_version > 74 and $params->{'source'} = $source;
    push( @transcripts, 
    $self->_create_feature_fast(
          'Bio::EnsEMBL::Transcript',$params) );

  }

  return \@transcripts;
}


=head2 fetch_all_by_exon_supporting_evidence

  Arg [1]    : String $hit_name
               Name of supporting feature
  Arg [2]    : String $feature_type 
               one of "dna_align_feature" or "protein_align_feature"
  Arg [3]    : (optional) Bio::Ensembl::Analysis
  Example    : $tr = $tr_adaptor->fetch_all_by_exon_supporting_evidence
                  ('XYZ', 'dna_align_feature');
  Description: Gets all the transcripts with exons which have a specified hit
               on a particular type of feature. Optionally filter by analysis.
  Returntype : Listref of Bio::EnsEMBL::Transcript objects
  Exceptions : If feature_type is not of correct type.
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_exon_supporting_evidence {
  my ($self, $hit_name, $feature_type, $analysis) = @_;

  if($feature_type !~ /(dna)|(protein)_align_feature/) {
    throw("feature type must be dna_align_feature or protein_align_feature");
  }

  my $anal_from = "";
  $anal_from = ", analysis a " if ($analysis);
  my $anal_where = "";
  $anal_where = "AND a.analysis_id = f.analysis_id AND a.analysis_id=? "
    if ($analysis);

  my $sql = qq(
      SELECT DISTINCT(t.transcript_id)
        FROM transcript t,
             exon_transcript et,
             supporting_feature sf,
             $feature_type f
             $anal_from
       WHERE t.transcript_id = et.transcript_id
         AND t.is_current = 1
         AND et.exon_id = sf.exon_id
         AND sf.feature_id = f.${feature_type}_id
         AND sf.feature_type = ?
         AND f.hit_name=?
         $anal_where
  );

  my $sth = $self->prepare($sql);

  $sth->bind_param(1, $feature_type, SQL_VARCHAR);
  $sth->bind_param(2, $hit_name, SQL_VARCHAR);
  $sth->bind_param(3, $analysis->dbID(), SQL_INTEGER) if ($analysis);

  $sth->execute();

  my @transcripts;

  while( my $id = $sth->fetchrow_array ) {
    my $transcript = $self->fetch_by_dbID( $id  );
    push(@transcripts, $transcript) if $transcript;
  }

  return \@transcripts;
}


=head2 fetch_all_by_transcript_supporting_evidence

  Arg [1]    : String $hit_name
               Name of supporting feature
  Arg [2]    : String $feature_type 
               one of "dna_align_feature" or "protein_align_feature"
  Arg [3]    : (optional) Bio::Ensembl::Analysis
  Example    : $transcripts = $transcript_adaptor->fetch_all_by_transcript_supporting_evidence('XYZ', 'dna_align_feature');
  Description: Gets all the transcripts with evidence from a specified hit_name on a particular type of feature, stored in the  
               transcript_supporting_feature table. Optionally filter by analysis.  For hits stored in the supporting_feature 
               table (linked to exons) use fetch_all_by_exon_supporting_evidence instead.
  Returntype : Listref of Bio::EnsEMBL::Transcript objects
  Exceptions : If feature_type is not of correct type.
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_transcript_supporting_evidence {
  
  my ($self, $hit_name, $feature_type, $analysis) = @_;

  if($feature_type !~ /(dna)|(protein)_align_feature/) {
    throw("feature type must be dna_align_feature or protein_align_feature");
  }

  my $anal_from = "";
  $anal_from = ", analysis a " if ($analysis);
  my $anal_where = "";
  $anal_where = "AND a.analysis_id = f.analysis_id AND a.analysis_id=? "
    if ($analysis);

  my $sql = qq(
      SELECT DISTINCT(t.transcript_id)
        FROM transcript t,
             transcript_supporting_feature sf,
             $feature_type f
             $anal_from
       WHERE t.transcript_id = sf.transcript_id
         AND t.is_current = 1
         AND sf.feature_id = f.${feature_type}_id
         AND sf.feature_type = ?
         AND f.hit_name=?
         $anal_where
  );

  my $sth = $self->prepare($sql);

  $sth->bind_param(1, $feature_type, SQL_VARCHAR);
  $sth->bind_param(2, $hit_name, SQL_VARCHAR);
  $sth->bind_param(3, $analysis->dbID(), SQL_INTEGER) if ($analysis);

  $sth->execute();

  my @transcripts;

  while( my $id = $sth->fetchrow_array ) {
    my $transcript = $self->fetch_by_dbID( $id  );
    push(@transcripts, $transcript) if $transcript;
  }

  return \@transcripts;
}


1;
