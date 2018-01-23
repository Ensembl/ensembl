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

Bio::EnsEMBL::DBSQL::OperonAdaptor - Database adaptor for the retrieval and
storage of OperonTranscript objects

=head1 SYNOPSIS


my $operon_transcript_adaptor =  Bio::EnsEMBL::DBSQL::OperonTranscriptAdaptor->new($dba);
$operon_transcript_adaptor->store($operon_transcript);
my $operon_transcript2 = $operon_transcript_adaptor->fetch_by_dbID( $operon->dbID() );
my $operon_transcripts = $operon_transcript_adaptor->fetch_all_by_gene( $gene );

=head1 DESCRIPTION

This is a database aware adaptor for the retrieval and storage of operon
transcript objects.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::OperonTranscriptAdaptor;

use strict;

use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning );
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Operon;
use Bio::EnsEMBL::OperonTranscript;
use Bio::EnsEMBL::Utils::SqlHelper;

use vars '@ISA';
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);

# _tables
#  Arg [1]    : none
#  Description: PROTECTED implementation of superclass abstract method.
#               Returns the names, aliases of the tables to use for queries.
#  Returntype : list of listrefs of strings
#  Exceptions : none
#  Caller     : interna
#  Status     : Stable

sub _tables {
	return ( [ 'operon_transcript', 'o' ] );
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

	my $created_date =
	  $self->db()->dbc()->from_date_to_seconds("o.created_date");
	my $modified_date =
	  $self->db()->dbc()->from_date_to_seconds("o.modified_date");

	return ( 'o.operon_transcript_id', 'o.seq_region_id',
			 'o.seq_region_start',     'o.seq_region_end',
			 'o.seq_region_strand',    'o.display_label',
			 'o.analysis_id',          'o.stable_id',
			 'o.version',            $created_date,
			 $modified_date );
}

=head2 list_dbIDs

  Example    : @ot_ids = @{$ot_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all operon_transcripts in the current db
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.
  Returntype : Listref of Ints
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub list_dbIDs {
	my ( $self, $ordered ) = @_;

	return $self->_list_dbIDs( "operon_transcript", undef, $ordered );
}

=head2 list_stable_ids

  Example    : @stable_ot_ids = @{$ot_adaptor->list_stable_ids()};
  Description: Gets an listref of stable ids for all operon_transcripts in the current db
  Returntype : reference to a list of strings
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub list_stable_ids {
	my ($self) = @_;

	return $self->_list_dbIDs( "operon_transcript", "stable_id" );
}

sub list_seq_region_ids {
	my $self = shift;

	return $self->_list_seq_region_ids('operon');
}

=head2 fetch_by_stable_id

  Arg [1]    : String $id 
               The stable ID of the operon_transcript to retrieve
  Example    : $operon_transcript = $operon_transcript_adaptor->fetch_by_stable_id('T16152-16153-4840');
  Description: Retrieves a operon_transcript object from the database via its stable id.
               The operon_transcript will be retrieved in its native coordinate system (i.e.
               in the coordinate system it is stored in the database). It may
               be converted to a different coordinate system through a call to
               transform() or transfer(). If the operon_transcript is not found
               undef is returned instead.
  Returntype : Bio::EnsEMBL::OperonTranscript or undef
  Exceptions : if we cant get the operon_transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_stable_id {
	my ( $self, $stable_id ) = @_;

	my $constraint = "o.stable_id = ?";
	$self->bind_param_generic_fetch( $stable_id, SQL_VARCHAR );
	my ($operon_transcript) = @{ $self->generic_fetch($constraint) };

	# If we didn't get anything back, desperately try to see if there's
	# a version number in the stable_id
	if(!defined($operon_transcript) && (my $vindex = rindex($stable_id, '.'))) {
	    $operon_transcript = $self->fetch_by_stable_id_version(substr($stable_id,0,$vindex),
								   substr($stable_id,$vindex+1));
	}

	return $operon_transcript;
}

=head2 fetch_by_stable_id_version

  Arg [1]    : String $id 
               The stable ID of the operon_transcript to retrieve
  Arg [2]    : Integer $version
               The version of the stable_id to retrieve
  Example    : $operon_transcript = $operon_transcript_adaptor->fetch_by_stable_id('T16152-16153-4840', 2);
  Description: Retrieves an operon_transcript object from the database via its stable id and version.
               The operon_transcript will be retrieved in its native coordinate system (i.e.
               in the coordinate system it is stored in the database). It may
               be converted to a different coordinate system through a call to
               transform() or transfer(). If the operon_transcript is not found
               undef is returned instead.
  Returntype : Bio::EnsEMBL::OperonTranscript or undef
  Exceptions : if we cant get the operon_transcript in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_stable_id_version {
    my ($self, $stable_id, $version) = @_;

    # Enforce that version be numeric
    return unless($version =~ /^\d+$/);

    my $constraint = "o.stable_id = ? AND o.version = ?";
    $self->bind_param_generic_fetch($stable_id, SQL_VARCHAR);
    $self->bind_param_generic_fetch($version, SQL_INTEGER);
    my ($operon_transcript) = @{$self->generic_fetch($constraint)};

    return $operon_transcript;
}

=head2 fetch_by_name

  Arg [1]    : String $label - name of operon transcript to fetch
  Example    : my $operon_transcript = $operonAdaptor->fetch_by_name("ECK0012121342");
  Description: Returns the operon transcript which has the given display label or undef if
               there is none. If there are more than 1, only the first is
               reported.
  Returntype : Bio::EnsEMBL::OperonTranscript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_name {
	my $self  = shift;
	my $label = shift;

	my $constraint = "o.display_label = ?";
	$self->bind_param_generic_fetch( $label, SQL_VARCHAR );
	my ($operon) = @{ $self->generic_fetch($constraint) };

	return $operon;
}

=head2 fetch_all

  Example     : $operon_transcripts = $operon_adaptor->fetch_all();
  Description : Retrieves all operon transcripts stored in the database.
  Returntype  : listref of Bio::EnsEMBL::OperonTranscript
  Caller      : general
  Status      : At Risk

=cut

sub fetch_all {
	my ($self) = @_;

	my $constraint         = '';
	my @operon_transcripts = @{ $self->generic_fetch($constraint) };
	return \@operon_transcripts;
}

=head2 fetch_all_versions_by_stable_id 

  Arg [1]     : String $stable_id 
                The stable ID of the operon_transcript to retrieve
  Example     : $operon_transcript = $operon_transcript_adaptor->fetch_all_versions_by_stable_id
                  ('ENSG00000148944');
  Description : Similar to fetch_by_stable_id, but retrieves all versions of a
                operon_transcript stored in the database.
  Returntype  : listref of Bio::EnsEMBL::OperonTranscript
  Exceptions  : if we cant get the operon_transcript in given coord system
  Caller      : general
  Status      : At Risk

=cut

sub fetch_all_versions_by_stable_id {
	my ( $self, $stable_id ) = @_;

	my $constraint = "o.stable_id = ?";
	$self->bind_param_generic_fetch( $stable_id, SQL_VARCHAR );
	return $self->generic_fetch($constraint);
}

=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to fetch operon_transcripts on.
  Arg [2]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Arg [3]    : (optional) boolean $load_transcripts
               if true, transcripts will be loaded immediately rather than
               lazy loaded later.
  Arg [4]    : (optional) string $source
               the source name of the features to obtain.
  Arg [5]    : (optional) string biotype
                the biotype of the features to obtain.
  Example    : @operon_transcripts = @{$operon_transcript_adaptor->fetch_all_by_Slice()};
  Description: Overrides superclass method to optionally load transcripts
               immediately rather than lazy-loading them later.  This
               is more efficient when there are a lot of operon_transcripts whose
               transcripts are going to be used.
  Returntype : reference to list of operon_transcripts 
  Exceptions : thrown if exon cannot be placed on transcript slice
  Caller     : Slice::get_all_OperonTranscripts
  Status     : Stable

=cut

sub fetch_all_by_Slice {
  my ( $self, $slice, $logic_name, $load_transcripts ) = @_;

  my $constraint = '';

  my $operons =
    $self->SUPER::fetch_all_by_Slice_constraint( $slice, $constraint,
						 $logic_name );

  # If there are less than two operons, still do lazy-loading.
  if ( !$load_transcripts || @$operons < 2 ) {
    return $operons;
  }

  # Preload all of the transcripts now, instead of lazy loading later,
  # faster than one query per transcript.

  # First check if transcripts are already preloaded.
  # FIXME: Should check all transcripts.
  if ( exists( $operons->[0]->{'_operon_transcript_array'} ) ) {
    return $operons;
  }

  # Get extent of region spanned by transcripts.
  my ($min_start, $max_end);
  my $ext_slice;

  unless ($slice->is_circular()) {
    foreach my $o (@$operons) {
      if (!defined($min_start) || $o->seq_region_start() < $min_start) {
	$min_start = $o->seq_region_start();
      }
      if (!defined($max_end) || $o->seq_region_end() > $max_end) {
	$max_end = $o->seq_region_end();
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
    foreach my $o (@$operons) {
      if (!defined($min_start) || ($o->start() >= 0 && $o->start() < $min_start)) {
  	$min_start = $o->start();
  	$min_start_feature = $o;
      }
      if (!defined($max_end) || ($o->end() >= 0 && $o->end() > $max_end)) {
  	$max_end = $o->end();
  	$max_end_feature = $o;
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


  # Associate transcript identifiers with operon_transcripts.

  my %o_hash = map { $_->dbID => $_ } @{$operons};

  my $o_id_str = join( ',', keys(%o_hash) );

  my $sth =
    $self->prepare(   "SELECT operon_id, operon_transcript_id "
		      . "FROM   operon_transcript "
		      . "WHERE  operon_id IN ($o_id_str)" );

  $sth->execute();

  my ( $o_id, $tr_id );
  $sth->bind_columns( \( $o_id, $tr_id ) );

  my %tr_o_hash;

  while ( $sth->fetch() ) {
    $tr_o_hash{$tr_id} = $o_hash{$o_id};
  }

  my $ta = $self->db()->get_OperonTranscriptAdaptor();
  my $transcripts =
    $ta->fetch_all_by_Slice( $ext_slice,
			     1, undef,
			     sprintf( "ot.operon_transcript_id IN (%s)",
				      join( ',',
					    sort { $a <=> $b }
					    keys(%tr_o_hash) ) ) );

  # Move transcripts onto operon_transcript slice, and add them to operon_transcripts.
  foreach my $tr ( @{$transcripts} ) {
    if ( !exists( $tr_o_hash{ $tr->dbID() } ) ) {
      next;
    }

    my $new_tr;
    if ( $slice != $ext_slice ) {
      $new_tr = $tr->transfer($slice);
      if ( !defined($new_tr) ) {
	throw("Unexpected. "
	      . "Transcript could not be transfered onto OperonTranscript slice."
	     );
      }
    } else {
      $new_tr = $tr;
    }

    $tr_o_hash{ $tr->dbID() }->add_OperonTranscript($new_tr);
  }

  return $operons;
}				## end sub fetch_all_by_Slice

=head2 fetch_by_Operon

  Arg [1]    : Bio::EnsEMBL::Operon
  Example    : $ot = $ot_adaptor->fetch_by_Operon($operon);
  Description: Retrieves all operon transcripts belonging to an operon
  Returntype : arrayref of Bio::EnsEMBL::OperonTranscript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Operon {
	my ( $self, $operon ) = @_;
	return $self->fetch_by_operon_id( $operon->dbID() );
}

=head2 fetch_by_operon_id

  Arg [1]    : Int id
  Example    : $ot = $ot_adaptor->fetch_by_operon_transcript($operon);
  Description: Retrieves all operon transcripts belonging to an operon
  Returntype : arrayref of Bio::EnsEMBL::OperonTranscript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_operon_id {
	my ( $self, $operon_id ) = @_;

	my $constraint = "o.operon_id = ?";
	$self->bind_param_generic_fetch( $operon_id, SQL_INTEGER );
	return $self->generic_fetch($constraint);
}

=head2 fetch_genes_by_operon_transcript

  Arg [1]    : Bio::EnsEMBL::OperonTranscript
  Example    : $ot = $ot_adaptor->fetch_genes_by_operon_transcript($operon_transcript);
  Description: Retrieves all genes attached to an operon transcript
  Returntype : arrayref of Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_genes_by_operon_transcript {
	my ( $self, $operon_transcript ) = @_;
	assert_ref( $operon_transcript, 'Bio::EnsEMBL::OperonTranscript' );
	return $self->fetch_genes_by_operon_transcript_id(
												   $operon_transcript->dbID() );
}

=head2 fetch_genes_by_operon_transcript_id

  Arg [1]    : Int id
  Example    : $ot = $ot_adaptor->fetch_genes_by_operon_transcript($operon_transcript_id);
  Description: Retrieves all genes attached to an operon transcript
  Returntype : arrayref of Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_genes_by_operon_transcript_id {
	my ( $self, $operon_transcript_id ) = @_;
	my $helper =
	  Bio::EnsEMBL::Utils::SqlHelper->new( -DB_CONNECTION => $self->db->dbc() );

	my $gene_ids =
	  $helper->execute_simple(
		-SQL =>
'SELECT  gene_id FROM operon_transcript_gene tr WHERE  operon_transcript_id =?',
		-PARAMS => [$operon_transcript_id] );

	my $genes        = [];
	my $gene_adaptor = $self->db()->get_GeneAdaptor();
	for my $gene_id (@$gene_ids) {
		push @$genes, $gene_adaptor->fetch_by_dbID($gene_id);
	}
	return $genes;
}

=head2 fetch_all_by_gene

  Arg [1]    : Bio::EnsEMBL::Gene
  Example    : $ots = $ot_adaptor->fetch_all_by_gene($gene);
  Description: Retrieves all operon transcripts attached to a given gene
  Returntype : arrayref of Bio::EnsEMBL::OperonTranscript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_gene {
	my ( $self, $gene ) = @_;
	assert_ref( $gene, 'Bio::EnsEMBL::Gene' );
	return $self->fetch_all_by_gene_id( $gene->dbID() );
}

=head2 fetch_all_by_gene_id

  Arg [1]    : Int id of Bio::EnsEMBL::Gene
  Example    : $ots = $ot_adaptor->fetch_all_by_gene($gene);
  Description: Retrieves all operon transcripts attached to a given gene
  Returntype : arrayref of Bio::EnsEMBL::OperonTranscript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_gene_id {
	my ( $self, $gene_id ) = @_;
	my $helper =
	  Bio::EnsEMBL::Utils::SqlHelper->new( -DB_CONNECTION => $self->db->dbc() );

	my $ot_ids = $helper->execute_simple(
		-SQL =>
'SELECT operon_transcript_id FROM operon_transcript_gene tr WHERE gene_id =?',
		-PARAMS => [$gene_id] );

	my $ots = [];
	for my $ot_id (@$ot_ids) {
		push @$ots, $self->fetch_by_dbID($ot_id);
	}
	return $ots;
}

=head2 store

  Arg [1]    : Bio::EnsEMBL::OperonTranscript $gene
               The gene to store in the database
  Arg [2]    : ignore_release in xrefs [default 1] set to 0 to use release info 
               in external database references
  Example    : $gene_adaptor->store($gene);
  Description: Stores a gene in the database.
  Returntype : the database identifier (dbID) of the newly stored gene
  Exceptions : thrown if the $gene is not a Bio::EnsEMBL::OperonTranscript or if 
               $gene does not have an analysis object
  Caller     : general
  Status     : Stable

=cut

sub store {
	my ( $self, $operon_transcript, $operon_id ) = @_;

	assert_ref( $operon_transcript, 'Bio::EnsEMBL::OperonTranscript' );

	my $db = $self->db();

	if ( $operon_transcript->is_stored($db) ) {
		return $operon_transcript->dbID();
	}

	# ensure coords are correct before storing
	#$operon->recalculate_coordinates();

	my $seq_region_id;

	( $operon_transcript, $seq_region_id ) =
	  $self->_pre_store($operon_transcript);
	my $analysis = $operon_transcript->analysis();
	throw("OperonTranscripts must have an analysis object.")
	  if ( !defined($analysis) );
	my $analysis_id;
	if ( $analysis->is_stored($db) ) {
		$analysis_id = $analysis->dbID();
	} else {
		$analysis_id = $db->get_AnalysisAdaptor->store($analysis);
	}
	my @columns = qw(
               seq_region_id
               seq_region_start
               seq_region_end
               seq_region_strand
               display_label
               operon_id
               analysis_id
        );
  my @canned_columns;
  my @canned_values;
	if ( defined($operon_transcript->stable_id()) ) {
		push @columns, qw(
      stable_id
      version
    );
		my $created = $self->db->dbc->from_seconds_to_date($operon_transcript->created_date());
		my $modified = $self->db->dbc->from_seconds_to_date($operon_transcript->modified_date());

		if ($created) {
		  push @canned_columns, 'created_date';
		  push @canned_values,  $created;
		}
		if ($modified) {
		  push @canned_columns, 'modified_date';
		  push @canned_values,  $modified;
		}
	}
	my $i_columns = join(', ', @columns, @canned_columns);
  my $i_values  = join(', ', (('?') x scalar(@columns)), @canned_values);
  my $store_operon_transcript_sql = qq(
    INSERT INTO operon_transcript ( ${i_columns} ) VALUES ( $i_values )
  );

	my $sth = $self->prepare($store_operon_transcript_sql);
	$sth->bind_param( 1, $seq_region_id,                      SQL_INTEGER );
	$sth->bind_param( 2, $operon_transcript->start(),         SQL_INTEGER );
	$sth->bind_param( 3, $operon_transcript->end(),           SQL_INTEGER );
	$sth->bind_param( 4, $operon_transcript->strand(),        SQL_TINYINT );
	$sth->bind_param( 5, $operon_transcript->display_label(), SQL_VARCHAR );
	$sth->bind_param( 6, $operon_id,                          SQL_INTEGER );
	$sth->bind_param( 7, $analysis_id,                        SQL_INTEGER );

	if ( defined($operon_transcript->stable_id()) ) {
	    $sth->bind_param( 8, $operon_transcript->stable_id(), SQL_VARCHAR );
	    my $version = ($operon_transcript->version()) ? $operon_transcript->version() : 1;
	    $sth->bind_param( 9, $version,  SQL_INTEGER );
	}

	$sth->execute();
	$sth->finish();

	my $operon_transcript_dbID = $self->last_insert_id('operon_transcript_id', undef, 'operon_transcript');

	# store the dbentries associated with this gene
	my $dbEntryAdaptor = $db->get_DBEntryAdaptor();

	foreach my $dbe ( @{ $operon_transcript->get_all_DBEntries } ) {
		$dbEntryAdaptor->store( $dbe, $operon_transcript_dbID,
								"OperonTranscript" );
	}

	# store operon attributes if there are any
	my $attrs = $operon_transcript->get_all_Attributes();
	if ( $attrs && scalar @$attrs ) {
		my $attr_adaptor = $db->get_AttributeAdaptor();
		$attr_adaptor->store_on_OperonTranscript( $operon_transcript, $attrs );
	}

	# set the adaptor and dbID on the original passed in gene not the
	# transfered copy
	$operon_transcript->adaptor($self);
	$operon_transcript->dbID($operon_transcript_dbID);

	if ( defined $operon_transcript->{_gene_array} ) {
		$self->store_genes_on_OperonTranscript( $operon_transcript,
											$operon_transcript->{_gene_array} );
	}

	return $operon_transcript_dbID;
} ## end sub store

=head2 store_genes_on_OperonTranscript

  Arg [1]    : Bio::EnsEMBL::OperonTranscript $ot
               the operon_transcript to store genes on
  Arg [2]    : arrayref of Bio::EnsEMBL::Gene $gene
               the genes to store on operon transcript
  Example    : $ot_adaptor->store_genes_on_OperonTranscript(\@genes);
  Description: Associates genes with operon transcript
  Returntype : none
  Exceptions : throw on incorrect arguments 
               warning if operon_transcript is not stored in this database
  Caller     : general, store
  Status     : Stable

=cut

sub store_genes_on_OperonTranscript {
	my ( $self, $operon_transcript, $genes ) = @_;
	assert_ref( $operon_transcript, "Bio::EnsEMBL::OperonTranscript" );
	my $sth = $self->prepare(
'insert into operon_transcript_gene(operon_transcript_id,gene_id) values('
		  . $operon_transcript->dbID()
		  . ',?)' );
	for my $gene ( @{$genes} ) {
		assert_ref( $gene, "Bio::EnsEMBL::Gene" );
		$sth->bind_param( 1, $gene->dbID(), SQL_INTEGER );
		$sth->execute();
	}
	$sth->finish();
	return;
}

=head2 remove

  Arg [1]    : Bio::EnsEMBL::OperonTranscript $ot
               the operon_transcript to remove from the database
  Example    : $ot_adaptor->remove($ot);
  Description: Removes a operon transcript completely from the database.
  Returntype : none
  Exceptions : throw on incorrect arguments 
               warning if operon_transcript is not stored in this database
  Caller     : general
  Status     : Stable

=cut

sub remove {
	my $self              = shift;
	my $operon_transcript = shift;

	assert_ref( $operon_transcript, 'Bio::EnsEMBL::OperonTranscript' );

	if ( !$operon_transcript->is_stored( $self->db() ) ) {
		warning(   "Cannot remove operon transcript "
				 . $operon_transcript->dbID()
				 . ". Is not stored in "
				 . "this database." );
		return;
	}

	# remove all object xrefs associated with this gene

	my $dbe_adaptor = $self->db()->get_DBEntryAdaptor();
	foreach my $dbe ( @{ $operon_transcript->get_all_DBEntries() } ) {
		$dbe_adaptor->remove_from_object( $dbe, $operon_transcript,
										  'OperonTranscript' );
	}

	#	# remove the attributes associated with this transcript
	#	my $attrib_adaptor = $self->db->get_AttributeAdaptor;
	#	$attrib_adaptor->remove_from_OperonTranscript($operon_transcript);

	# remove from the database
	my $sth = $self->prepare(
			   "DELETE FROM operon_transcript WHERE operon_transcript_id = ? ");
	$sth->bind_param( 1, $operon_transcript->dbID, SQL_INTEGER );
	$sth->execute();
	$sth->finish();

	# unset the gene identifier and adaptor thereby flagging it as unstored

	$operon_transcript->dbID(undef);
	$operon_transcript->adaptor(undef);

	return;
} ## end sub remove

# _objs_from_sth

#  Arg [1]    : StatementHandle $sth
#  Arg [2]    : Bio::EnsEMBL::AssemblyMapper $mapper
#  Arg [3]    : Bio::EnsEMBL::Slice $dest_slice
#  Description: PROTECTED implementation of abstract superclass method.
#               responsible for the creation of OperonTranscripts
#  Returntype : listref of Bio::EnsEMBL::OperonTranscripts in target coordinate system
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

  my $sa = $self->db()->get_SliceAdaptor();
  my $aa = $self->db()->get_AnalysisAdaptor();

  my @operons;
  my %analysis_hash;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;

  my (
     $operon_transcript_id, $seq_region_id,     $seq_region_start,
     $seq_region_end,       $seq_region_strand, $display_label,
     $analysis_id,          $stable_id,         $version,
     $created_date,         $modified_date );

  $sth->bind_columns( \(
            $operon_transcript_id, $seq_region_id,     $seq_region_start,
            $seq_region_end,       $seq_region_strand, $display_label,
            $analysis_id,          $stable_id,         $version,
            $created_date,         $modified_date ));

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

  OPERON: while ( $sth->fetch() ) {

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
      next OPERON if ( !defined($seq_region_id) );

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
        next OPERON;
      }
      $slice = $dest_slice;
    }

    push( @operons,
        Bio::EnsEMBL::OperonTranscript->new(
                    -START         => $seq_region_start,
                    -END           => $seq_region_end,
                    -STRAND        => $seq_region_strand,
                    -SLICE         => $slice,
                    -DISPLAY_LABEL => $display_label,
                    -ADAPTOR       => $self,
                    -DBID          => $operon_transcript_id,
                    -STABLE_ID     => $stable_id,
                    -VERSION       => $version,
                    -CREATED_DATE  => $created_date || undef,
                    -MODIFIED_DATE => $modified_date || undef,
                    -ANALYSIS      => $analysis ) );

  } ## end while ( $sth->fetch() )

  return \@operons;
} ## end sub _objs_from_sth

1;

