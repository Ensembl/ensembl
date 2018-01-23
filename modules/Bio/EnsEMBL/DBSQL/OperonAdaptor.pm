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
storage of Operon objects

=head1 SYNOPSIS

my $operon_adaptor =  Bio::EnsEMBL::DBSQL::OperonAdaptor->new($dba);
$operon_adaptor->store($operon);
my $operon2 = $operon_adaptor->fetch_by_dbID( $operon->dbID() );

=head1 DESCRIPTION

This is a database aware adaptor for the retrieval and storage of operon
objects.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::OperonAdaptor;

use strict;

use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning );
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Operon;

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
	return ( [ 'operon', 'o' ] );
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

	return ( 'o.operon_id',      'o.seq_region_id',     'o.seq_region_start',
			 'o.seq_region_end', 'o.seq_region_strand', 'o.display_label',
			 'o.analysis_id',    'o.stable_id',       'o.version',
			 $created_date,      $modified_date );
}

=head2 list_dbIDs 

  Example    : @operon_ids = @{$operon_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all operons in the current db
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.
  Returntype : Listref of Ints
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub list_dbIDs {
	my ( $self, $ordered ) = @_;

	return $self->_list_dbIDs( "operon", undef, $ordered );
}

=head2 list_stable_ids

  Example    : @stable_operon_ids = @{$operon_adaptor->list_stable_ids()};
  Description: Gets an listref of stable ids for all operons in the current db
  Returntype : reference to a list of strings
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub list_stable_ids {
	my ($self) = @_;

	return $self->_list_dbIDs( "operon", "stable_id" );
}

sub list_seq_region_ids {
	my $self = shift;

	return $self->_list_seq_region_ids('operon');
}

=head2 fetch_by_name

  Arg [1]    : String $label - name of operon to fetch
  Example    : my $operon = $operonAdaptor->fetch_by_name("accBC");
  Description: Returns the operon which has the given display label or undef if
               there is none. If there are more than 1, only the first is
               reported.
  Returntype : Bio::EnsEMBL::Operon
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

=head2 fetch_by_stable_id

  Arg [1]    : String $id 
               The stable ID of the operon to retrieve
  Example    : $operon = $operon_adaptor->fetch_by_stable_id('ENSG00000148944');
  Description: Retrieves a operon object from the database via its stable id.
               The operon will be retrieved in its native coordinate system (i.e.
               in the coordinate system it is stored in the database). It may
               be converted to a different coordinate system through a call to
               transform() or transfer(). If the operon or exon is not found
               undef is returned instead.
  Returntype : Bio::EnsEMBL::Operon or undef
  Exceptions : if we cant get the operon in given coord system
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_stable_id {
	my ( $self, $stable_id ) = @_;

	my $constraint = "o.stable_id = ?";
	$self->bind_param_generic_fetch( $stable_id, SQL_VARCHAR );
	my ($operon) = @{ $self->generic_fetch($constraint) };

	# If we didn't get anything back, desperately try to see if there's
	# a version number in the stable_id
	if(!defined($operon) && (my $vindex = rindex($stable_id, '.'))) {
	    $operon = $self->fetch_by_stable_id_version(substr($stable_id,0,$vindex),
							substr($stable_id,$vindex+1));
	}

	return $operon;
}

=head2 fetch_all

  Example     : $operons = $operon_adaptor->fetch_all();
  Description : Similar to fetch_by_stable_id, but retrieves all
                operons stored in the database.
  Returntype  : listref of Bio::EnsEMBL::Operon
  Caller      : general
  Status      : At Risk

=cut

=head2 fetch_by_stable_id_version

  Arg [1]    : String $id 
               The stable ID of the operon to retrieve
  Arg [2]    : Integer $version
               The version of the stable_id to retrieve
  Example    : $operon = $operon_adaptor->fetch_by_stable_id('16152-16153-4840', 2);
  Description: Retrieves an operon object from the database via its stable id and version.
               The operon will be retrieved in its native coordinate system (i.e.
               in the coordinate system it is stored in the database). It may
               be converted to a different coordinate system through a call to
               transform() or transfer(). If the operon is not found
               undef is returned instead.
  Returntype : Bio::EnsEMBL::Operon or undef
  Exceptions : if we cant get the operon in given coord system
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
    my ($operon) = @{$self->generic_fetch($constraint)};

    return $operon;
}

sub fetch_all {
	my ($self) = @_;

	my $constraint = '';
	my @operons    = @{ $self->generic_fetch($constraint) };
	return \@operons;
}

=head2 fetch_all_versions_by_stable_id 

  Arg [1]     : String $stable_id 
                The stable ID of the operon to retrieve
  Example     : $operon = $operon_adaptor->fetch_all_versions_by_stable_id
                  ('ENSG00000148944');
  Description : Similar to fetch_by_stable_id, but retrieves all versions of a
                operon stored in the database.
  Returntype  : listref of Bio::EnsEMBL::Operon
  Exceptions  : if we cant get the operon in given coord system
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
               The slice to fetch operons on.
  Arg [2]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Arg [3]    : (optional) boolean $load_transcripts
               if true, transcripts will be loaded immediately rather than
               lazy loaded later.
  Arg [4]    : (optional) string $source
               the source name of the features to obtain.
  Arg [5]    : (optional) string biotype
                the biotype of the features to obtain.
  Example    : @operons = @{$operon_adaptor->fetch_all_by_Slice()};
  Description: Overrides superclass method to optionally load transcripts
               immediately rather than lazy-loading them later.  This
               is more efficient when there are a lot of operons whose
               transcripts are going to be used.
  Returntype : reference to list of operons 
  Exceptions : thrown if exon cannot be placed on transcript slice
  Caller     : Slice::get_all_operons
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


  # Associate transcript identifiers with operons.

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

  # Move transcripts onto operon slice, and add them to operons.
  foreach my $tr ( @{$transcripts} ) {
    if ( !exists( $tr_o_hash{ $tr->dbID() } ) ) {
      next;
    }

    my $new_tr;
    if ( $slice != $ext_slice ) {
      $new_tr = $tr->transfer($slice);
      if ( !defined($new_tr) ) {
	throw(   "Unexpected. "
		 . "Transcript could not be transfered onto operon slice."
	     );
      }
    } else {
      $new_tr = $tr;
    }

    $tr_o_hash{ $tr->dbID() }->add_OperonTranscript($new_tr);
  }

  return $operons;
}				## end sub fetch_all_by_Slice

=head2 fetch_by_transcript_id

  Arg [1]    : Int $trans_id
               Unique database identifier for the transcript whose operon should
               be retrieved. The operon is returned in its native coord
               system (i.e. the coord_system it is stored in). If the coord
               system needs to be changed, then tranform or transfer should
               be called on the returned object. undef is returned if the
               operon or transcript is not found in the database.
  Example    : $operon = $operon_adaptor->fetch_by_transcript_id(1241);
  Description: Retrieves a operon from the database via the database identifier
               of one of its transcripts.
  Returntype : Bio::EnsEMBL::Operon
  Exceptions : none
  Caller     : operonral
  Status     : Stable

=cut

sub fetch_by_operon_transcript_id {
	my ( $self, $trans_id ) = @_;

	# this is a cheap SQL call
	my $sth = $self->prepare(
		qq(
      SELECT tr.operon_id
      FROM operon_transcript tr
      WHERE tr.operon_transcript_id = ?
  ) );

	$sth->bind_param( 1, $trans_id, SQL_INTEGER );
	$sth->execute();

	my ($operonid) = $sth->fetchrow_array();

	$sth->finish();

	return undef if ( !defined $operonid );

	my $operon = $self->fetch_by_dbID($operonid);
	return $operon;
}

=head2 fetch_by_operon_transcript_stable_id

  Arg [1]    : string $trans_stable_id
               transcript stable ID whose operon should be retrieved
  Example    : my $operon = $operon_adaptor->fetch_by_operon_transcript_stable_id
                 ('ENST0000234');
  Description: Retrieves a operon from the database via the stable ID of one of
               its transcripts
  Returntype : Bio::EnsEMBL::Operon
  Exceptions : none
  Caller     : operonral
  Status     : Stable

=cut

sub fetch_by_operon_transcript_stable_id {
	my ( $self, $trans_stable_id ) = @_;

	my $sth = $self->prepare(
		qq(
        SELECT  operon_id
	FROM	operon_transcript
        WHERE   stable_id = ?
    ) );

	$sth->bind_param( 1, $trans_stable_id, SQL_VARCHAR );
	$sth->execute();

	my ($operonid) = $sth->fetchrow_array();
	$sth->finish;

	return undef if ( !defined $operonid );

	my $operon = $self->fetch_by_dbID($operonid);
	return $operon;
}

sub fetch_by_operon_transcript {
	my ( $self, $trans ) = @_;
	assert_ref( $trans, 'Bio::EnsEMBL::OperonTranscript' );
	$self->fetch_by_operon_transcript_id( $trans->dbID() );
}

=head2 store

  Arg [1]    : Bio::EnsEMBL::Operon $operon
               The operon to store in the database
  Arg [2]    : ignore_release in xrefs [default 1] set to 0 to use release info 
               in external database references
  Example    : $operon_adaptor->store($operon);
  Description: Stores a operon in the database.
  Returntype : the database identifier (dbID) of the newly stored operon
  Exceptions : thrown if the $operon is not a Bio::EnsEMBL::Operon or if 
               $operon does not have an analysis object
  Caller     : general
  Status     : Stable

=cut

sub store {
	my ( $self, $operon, $ignore_release ) = @_;

	if ( !ref $operon || !$operon->isa('Bio::EnsEMBL::Operon') ) {
		throw("Must store a operon object, not a $operon");
	}

	my $db = $self->db();

	if ( $operon->is_stored($db) ) {
		return $operon->dbID();
	}
		    my $analysis = $operon->analysis();
  throw("Operons must have an analysis object.") if(!defined($analysis));
	my $analysis_id;
	if ( $analysis->is_stored($db) ) {
		$analysis_id = $analysis->dbID();
	} else {
		$analysis_id = $db->get_AnalysisAdaptor->store( $analysis );
	}
	# ensure coords are correct before storing
	#$operon->recalculate_coordinates();

	my $seq_region_id;

	( $operon, $seq_region_id ) = $self->_pre_store($operon);

	my @columns = qw(
               seq_region_id
               seq_region_start
               seq_region_end
               seq_region_strand
               display_label
               analysis_id
        );
	my @canned_columns;
  my @canned_values;

	if ( defined($operon->stable_id()) ) {
		push @columns, qw(
		  stable_id
		  version
		);
		my $created = $self->db->dbc->from_seconds_to_date($operon->created_date());
		my $modified = $self->db->dbc->from_seconds_to_date($operon->modified_date());

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
  my $store_operon_sql = qq(
    INSERT INTO operon ( ${i_columns} ) VALUES ( $i_values )
  );

	my $sth = $self->prepare($store_operon_sql);
	$sth->bind_param( 1, $seq_region_id,           SQL_INTEGER );
	$sth->bind_param( 2, $operon->start(),         SQL_INTEGER );
	$sth->bind_param( 3, $operon->end(),           SQL_INTEGER );
	$sth->bind_param( 4, $operon->strand(),        SQL_TINYINT );
	$sth->bind_param( 5, $operon->display_label(), SQL_VARCHAR );
	$sth->bind_param( 6, $analysis_id,             SQL_INTEGER );

	if ( defined($operon->stable_id()) ) {
		$sth->bind_param( 7, $operon->stable_id(), SQL_VARCHAR );
		my $version = ($operon->version()) ? $operon->version() : 1;
		$sth->bind_param( 8, $version,             SQL_INTEGER );
	}
	$sth->execute();
	$sth->finish();

	my $operon_dbID = $self->last_insert_id('operon_id', undef, 'operon');

	my $transcripts = $operon->get_all_OperonTranscripts();

	if ( $transcripts && scalar @$transcripts ) {
		my $transcript_adaptor = $db->get_OperonTranscriptAdaptor();
		for my $transcript (@$transcripts) {
			$transcript_adaptor->store( $transcript, $operon_dbID );
		}
	}

	# store the dbentries associated with this operon
	my $dbEntryAdaptor = $db->get_DBEntryAdaptor();

	foreach my $dbe ( @{ $operon->get_all_DBEntries } ) {
		$dbEntryAdaptor->store( $dbe, $operon_dbID, "Operon", $ignore_release );
	}

	# store operon attributes if there are any
	my $attrs = $operon->get_all_Attributes();
	if ( $attrs && scalar @$attrs ) {
		my $attr_adaptor = $db->get_AttributeAdaptor();
		$attr_adaptor->store_on_Operon( $operon, $attrs );
	}

	# set the adaptor and dbID on the original passed in operon not the
	# transfered copy
	$operon->adaptor($self);
	$operon->dbID($operon_dbID);

	return $operon_dbID;
} ## end sub store

=head2 remove

  Arg [1]    : Bio::EnsEMBL::Operon $operon
               the operon to remove from the database
  Example    : $operon_adaptor->remove($operon);
  Description: Removes a operon completely from the database. All associated
               transcripts, exons, stable_identifiers, descriptions, etc.
               are removed as well. Use with caution!
  Returntype : none
  Exceptions : throw on incorrect arguments 
               warning if operon is not stored in this database
  Caller     : general
  Status     : Stable

=cut

sub remove {
	my $self   = shift;
	my $operon = shift;

	if ( !ref($operon) || !$operon->isa('Bio::EnsEMBL::Operon') ) {
		throw("Bio::EnsEMBL::Operon argument expected.");
	}

	if ( !$operon->is_stored( $self->db() ) ) {
		warning(   "Cannot remove operon "
				 . $operon->dbID()
				 . ". Is not stored in "
				 . "this database." );
		return;
	}

	# remove all object xrefs associated with this operon

	my $dbe_adaptor = $self->db()->get_DBEntryAdaptor();
	foreach my $dbe ( @{ $operon->get_all_DBEntries() } ) {
		$dbe_adaptor->remove_from_object( $dbe, $operon, 'Operon' );
	}

	# remove all of the transcripts associated with this operon
	my $transcriptAdaptor = $self->db->get_OperonTranscriptAdaptor();
	foreach my $trans ( @{ $operon->get_all_OperonTranscripts() } ) {
		$transcriptAdaptor->remove($trans);
	}

	# remove this operon from the database

	my $sth = $self->prepare("DELETE FROM operon WHERE operon_id = ? ");
	$sth->bind_param( 1, $operon->dbID, SQL_INTEGER );
	$sth->execute();
	$sth->finish();

	# unset the operon identifier and adaptor thereby flagging it as unstored

	$operon->dbID(undef);
	$operon->adaptor(undef);

	return;
} ## end sub remove

# _objs_from_sth

#  Arg [1]    : StatementHandle $sth
#  Arg [2]    : Bio::EnsEMBL::AssemblyMapper $mapper
#  Arg [3]    : Bio::EnsEMBL::Slice $dest_slice
#  Description: PROTECTED implementation of abstract superclass method.
#               responsible for the creation of Operons
#  Returntype : listref of Bio::EnsEMBL::Operon in target coordinate system
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
          $operon_id,      $seq_region_id,     $seq_region_start,
          $seq_region_end, $seq_region_strand, $display_label,
          $analysis_id,    $stable_id,         $version, 
          $created_date,   $modified_date
        );

  $sth->bind_columns( \(
                             $operon_id,         $seq_region_id,     $seq_region_start,
                             $seq_region_end,    $seq_region_strand, $display_label,
                             $analysis_id,       $stable_id,         $version,
                             $created_date,      $modified_date ));

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
      next OPERON if (!defined($seq_region_id));

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
    } ## end if ($dest_slice)

    push( @operons,
        Bio::EnsEMBL::Operon->new(
                    -START         => $seq_region_start,
                    -END           => $seq_region_end,
                    -STRAND        => $seq_region_strand,
                    -SLICE         => $slice,
                    -DISPLAY_LABEL => $display_label,
                    -ADAPTOR       => $self,
                    -DBID          => $operon_id,
                    -STABLE_ID     => $stable_id,
                    -VERSION       => $version,
                    -CREATED_DATE  => $created_date || undef,
                    -MODIFIED_DATE => $modified_date || undef,
                    -ANALYSIS      => $analysis ) );

  } ## end while ( $sth->fetch() )

  return \@operons;
} ## end sub _objs_from_sth

1;

