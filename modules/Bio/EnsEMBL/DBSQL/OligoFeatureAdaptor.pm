#
# Ensembl module for Bio::EnsEMBL::DBSQL::OligoFeatureAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::DBSQL::OligoFeatureAdaptor - A database adaptor for fetching and
storing OligoFeature objects.

=head1 SYNOPSIS

my $ofa = $db->get_OligoFeatureAdaptor();

my $features = $ofa->fetch_all_by_Probe($probe);
$features = $ofa->fetch_all_by_Slice_arrayname($slice, 'Array-1', 'Array-2');

=head1 DESCRIPTION

The OligoFeatureAdaptor is a database adaptor for storing and retrieving
OligoFeature objects.

=head1 AUTHOR

This module was created by Ian Sealy, but is almost entirely based on the
OligoFeatureAdaptor module written by Arne Stabenau.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::DBSQL::OligoFeatureAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::OligoFeature;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);


=head2 fetch_all_by_Probe

  Arg [1]    : Bio::EnsEMBL::OligoProbe
  Example    : my $features = $ofa->fetch_all_by_Probe($probe);
  Description: Fetchs all features that a given probe creates.
  Returntype : Listref of Bio::EnsEMBL::OligoFeature objects
  Exceptions : Throws if argument is not a stored OligoProbe object
  Caller     : OligoProbe->get_all_OligoFeatures()
  Status     : Medium Risk

=cut

sub fetch_all_by_Probe {
	my $self  = shift;
	my $probe = shift;
	
	if ( !ref($probe) && !$probe->isa('Bio::EnsEMBL::OligoProbe') ) {
		throw('fetch_all_by_Probe requires a Bio::EnsEMBL::OligoProbe object');
	}
	if ( !defined $probe->dbID() ) {
		throw('fetch_all_by_Probe requires a stored Bio::EnsEMBL::OligoProbe object');
	}
	
	return $self->generic_fetch( 'of.oligo_probe_id = ' . $probe->dbID() );
}

=head2 fetch_all_by_probeset

  Arg [1]    : string - probeset
  Example    : my $features = $ofa->fetch_all_by_probeset('Set-1');
  Description: Fetchs all features that a given probeset creates.
  Returntype : Listref of Bio::EnsEMBL::OligoFeature objects
  Exceptions : Throws if no probeset argument
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_all_by_probeset {
	my $self     = shift;
	my $probeset = shift;
	
	if (!$probeset) {
		throw('fetch_all_by_probeset requires a probeset argument');
	}
	
	return $self->generic_fetch( "op.probeset = '$probeset'" );
}

=head2 fetch_all_by_Slice_arrayname

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2...] : List of strings - array name(s)
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_by_Slice_arrayname($slice, '');
  Description: Retrieves a list of features on a given slice that are created
               by probes from the specified arrays.
  Returntype : Listref of Bio::EnsEMBL::OligoFeature objects
  Exceptions : Throws if no array name is provided
  Caller     : Slice->get_all_OligoFeatures()
  Status     : Medium Risk

=cut

sub fetch_all_by_Slice_arrayname {
	my ($self, $slice, @arraynames) = @_;
	
	throw('Need array name as parameter') if !@arraynames;
	
	my $constraint;
	if (scalar @arraynames == 1) {
		$constraint = qq( oa.name = '$arraynames[0]' );
	} else {
		$constraint = join q(','), @arraynames;
		$constraint = qq( oa.name IN ('$constraint') );
	}
	
	return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
}

=head2 fetch_all_by_Slice_type

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : string - type of array (e.g. AFFY or OLIGO)
  Arg [3]    : (optional) string - logic name
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_by_Slice_type($slice, 'OLIGO');
  Description: Retrieves a list of features on a given slice that are created
               by probes from the specified type of array.
  Returntype : Listref of Bio::EnsEMBL::OligoFeature objects
  Exceptions : Throws if no array type is provided
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_all_by_Slice_type {
	my ($self, $slice, $type, $logic_name) = @_;
	
	throw('Need type as parameter') if !$type;
	
	my $constraint = qq( oa.type = '$type' );
	
	return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint, $logic_name);
}

=head2 _tables

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns the names and aliases of the tables to use for queries.
  Returntype : List of listrefs of strings
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _tables {
	my $self = shift;
	
	return (
		[ 'oligo_feature', 'of' ], 
		[ 'oligo_probe',   'op' ], 
		[ 'oligo_array',   'oa' ]
	);
}

=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _columns {
	my $self = shift;
	
	return qw(
		of.oligo_feature_id  of.seq_region_id
		of.seq_region_start  of.seq_region_end
		of.seq_region_strand of.mismatches
		of.oligo_probe_id    of.analysis_id
		oa.name              op.probeset
		op.name
	);
}

=head2 _default_where_clause

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns an additional table joining constraint to use for
			   queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut
sub _default_where_clause {
	my $self = shift;
	
	return 'of.oligo_probe_id = op.oligo_probe_id AND op.oligo_array_id = oa.oligo_array_id';
}

=head2 _final_clause

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns an ORDER BY clause. Sorting by oligo_feature_id would be
			   enough to eliminate duplicates, but sorting by location might
			   make fetching features on a slice faster.
  Returntype : String
  Exceptions : None
  Caller     : generic_fetch
  Status     : Medium Risk

=cut


sub _final_clause {
	return ' ORDER BY of.seq_region_id, of.seq_region_start, of.oligo_feature_id';
}


=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates OligoFeature objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::OligoFeature objects
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _objs_from_sth {
	my ($self, $sth, $mapper, $dest_slice) = @_;
	
	# This code is ugly because caching is used to improve speed

	my $sa = $self->db->get_SliceAdaptor();
	my $aa = $self->db->get_AnalysisAdaptor();

	my @features;
	
	my (%analysis_hash, %slice_hash, %sr_name_hash, %sr_cs_hash);

	my (
		$oligo_feature_id,  $seq_region_id,
		$seq_region_start,  $seq_region_end,
		$seq_region_strand, $mismatches,
		$oligo_probe_id,    $analysis_id,
		$array_name,        $probeset,
		$oligo_probe_name,
	);
	$sth->bind_columns(
		\$oligo_feature_id,  \$seq_region_id,
		\$seq_region_start,  \$seq_region_end,
		\$seq_region_strand, \$mismatches,
		\$oligo_probe_id,    \$analysis_id,
		\$array_name,        \$probeset,
		\$oligo_probe_name,
	);

	my $asm_cs;
	my $cmp_cs;
	my $asm_cs_name;
	my $asm_cs_vers;
	my $cmp_cs_name;
	my $cmp_cs_vers;
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

	my $last_feature_id = -1;
	FEATURE: while ( $sth->fetch() ) {

		# This assumes that features come out sorted by ID
		next if ($last_feature_id == $oligo_feature_id);
		$last_feature_id = $oligo_feature_id;

		# Get the analysis object
		my $analysis = $analysis_hash{$analysis_id} ||= $aa->fetch_by_dbID($analysis_id);

		# Get the slice object
		my $slice = $slice_hash{'ID:'.$seq_region_id};

		if (!$slice) {
			$slice                            = $sa->fetch_by_seq_region_id($seq_region_id);
			$slice_hash{'ID:'.$seq_region_id} = $slice;
			$sr_name_hash{$seq_region_id}     = $slice->seq_region_name();
			$sr_cs_hash{$seq_region_id}       = $slice->coord_system();
		}

		my $sr_name = $sr_name_hash{$seq_region_id};
		my $sr_cs   = $sr_cs_hash{$seq_region_id};

		# Remap the feature coordinates to another coord system if a mapper was provided
		if ($mapper) {
			($seq_region_id, $seq_region_start, $seq_region_end, $seq_region_strand)
				= $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand, $sr_cs);

			# Skip features that map to gaps or coord system boundaries
			next FEATURE if !defined $seq_region_id;

			# Get a slice in the coord system we just mapped to
#			if ( $asm_cs == $sr_cs || ( $cmp_cs != $sr_cs && $asm_cs->equals($sr_cs) ) ) {
				$slice = $slice_hash{"ID:".$seq_region_id}
					||= $sa->fetch_by_seq_region_id($seq_region_id);
#			} else {
#				$slice = $slice_hash{"NAME:$sr_name:$asm_cs_name:$asm_cs_vers"}
#					||= $sa->fetch_by_region($asm_cs_name, $sr_name, undef, undef, undef, $asm_cs_vers);
#			}
		}

		# If a destination slice was provided convert the coords
		# If the destination slice starts at 1 and is forward strand, nothing needs doing
		if ($dest_slice) {
			unless ($dest_slice_start == 1 && $dest_slice_strand == 1) {
				if ($dest_slice_strand == 1) {
					$seq_region_start = $seq_region_start - $dest_slice_start + 1;
					$seq_region_end   = $seq_region_end   - $dest_slice_start + 1;
				} else {
					my $tmp_seq_region_start = $seq_region_start;
					$seq_region_start        = $dest_slice_end - $seq_region_end       + 1;
					$seq_region_end          = $dest_slice_end - $tmp_seq_region_start + 1;
					$seq_region_strand      *= -1;
				}
			}

			# Throw away features off the end of the requested slice
			next FEATURE if $seq_region_end < 1 || $seq_region_start > $dest_slice_length
				|| ( $dest_slice_sr_id ne $seq_region_id );

			$slice = $dest_slice;
		}

		push @features, $self->_new_fast( {
			'start'         => $seq_region_start,
			'end'           => $seq_region_end,
			'strand'        => $seq_region_strand,
			'slice'         => $slice,
			'analysis'      => $analysis,
			'adaptor'       => $self,
			'dbID'          => $oligo_feature_id,
			'mismatchcount' => $mismatches,
			'_probe_id'     => $oligo_probe_id,
			'probeset'      => $probeset,
			'_probe_name'   => $oligo_probe_name
		} );
	}

	return \@features;
}

=head2 _new_fast

  Args       : Hashref to be passed to OligoFeature->new_fast()
  Example    : None
  Description: Construct an OligoFeature object using quick and dirty new_fast.
  Returntype : Bio::EnsEMBL::OligoFeature
  Exceptions : None
  Caller     : _objs_from_sth
  Status     : Medium Risk

=cut

sub _new_fast {
	my $self = shift;
	
	my $hash_ref = shift;
	return Bio::EnsEMBL::OligoFeature->new_fast($hash_ref);
}

=head2 store

  Args       : List of Bio::EnsEMBL::OligoFeature objects
  Example    : $ofa->store(@features);
  Description: Stores given OligoFeature objects in the database. Should only
               be called once per feature because no checks are made for
			   duplicates. Sets dbID and adaptor on the objects that it stores.
  Returntype : None
  Exceptions : Throws if a list of OligoFeature objects is not provided or if
               an analysis is not attached to any of the objects
  Caller     : General
  Status     : Medium Risk

=cut

sub store{
	my ($self, @ofs) = @_;

	if (scalar(@ofs) == 0) {
		throw('Must call store with a list of OligoFeature objects');
	}

	my $sth = $self->prepare("
		INSERT INTO oligo_feature (
			seq_region_id,  seq_region_start,
			seq_region_end, seq_region_strand,
			oligo_probe_id,  analysis_id,
			mismatches
		) VALUES (?, ?, ?, ?, ?, ?, ?)
	");

	my $db = $self->db();
	my $analysis_adaptor = $db->get_AnalysisAdaptor();

	FEATURE: foreach my $of (@ofs) {

		if( !ref $of || !$of->isa('Bio::EnsEMBL::OligoFeature') ) {
			throw('Feature must be an OligoFeature object');
		}

		if ( $of->is_stored($db) ) {
			warning('OligoFeature [' . $of->dbID() . '] is already stored in the database');
			next FEATURE;
		}

		if ( !defined $of->analysis() ) {
			throw('An analysis must be attached to the OligoFeature objects to be stored.');
		}

		# Store the analysis if it has not been stored yet
		if ( !$of->analysis->is_stored($db) ) {
			$analysis_adaptor->store( $of->analysis() );
		}

		my $original = $of;
		my $seq_region_id;
		($of, $seq_region_id) = $self->_pre_store($of);

		$sth->bind_param(1, $seq_region_id,        SQL_INTEGER);
		$sth->bind_param(2, $of->start(),          SQL_INTEGER);
		$sth->bind_param(3, $of->end(),            SQL_INTEGER);
		$sth->bind_param(4, $of->strand(),         SQL_TINYINT);
		$sth->bind_param(5, $of->probe->dbID(),    SQL_INTEGER);
		$sth->bind_param(6, $of->analysis->dbID(), SQL_INTEGER);
		$sth->bind_param(7, $of->mismatchcount(),  SQL_TINYINT);

		$sth->execute();

		$original->dbID( $sth->{'mysql_insertid'} );
		$original->adaptor($self);
	}
}

=head2 list_dbIDs

  Args       : None
  Example    : my @feature_ids = @{$ofa->list_dbIDs()};
  Description: Gets an array of internal IDs for all OligoFeature objects in
               the current database.
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : Medium Risk

=cut

sub list_dbIDs {
  my ($self,$ordered) = shift;
	
  return $self->_list_dbIDs('oligo_feature',undef,$ordered);
}

1;

