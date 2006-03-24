#
# Ensembl module for Bio::EnsEMBL::DBSQL::OligoProbeAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::DBSQL::OligoProbeAdaptor - A database adaptor for fetching and
storing OligoProbe objects.

=head1 SYNOPSIS

my $opa = $db->get_OligoProbeAdaptor();

my $probe = $opa->fetch_by_array_probeset_probe('Array-1', undef, 'Probe-1');

=head1 DESCRIPTION

The OligoProbeAdaptor is a database adaptor for storing and retrieving
OligoProbe objects.

=head1 AUTHOR

This module was created by Ian Sealy, but is almost entirely based on the
OligoProbeAdaptor module written by Arne Stabenau.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::DBSQL::OligoProbeAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::OligoProbe;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 fetch_by_array_probeset_probe

  Arg [1]    : string - name of array
  Arg [2]    : (optional) string - name of probeset
  Arg [3]    : string - name of probe
  Example    : my $probe = $opa->fetch_by_array_probeset_probe('Array-1', 'Probeset-1', 'Probe-1');
  Description: Returns a probe given a combination of array name, probeset and
               probe name. This will uniquely define an Affy probe. Only one
			   probe is ever returned.
  Returntype : Bio::EnsEMBL::OligoProbe
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_by_array_probeset_probe {
	my $self          = shift;
	my $array_name    = shift;
	my $probeset_name = shift;
	my $probe_name    = shift;
	
	my $probeset_clause = 'AND op.probeset = ?';
	
	# Need to deal with non-Affy probes where probeset is NULL
	# (Also allow probeset to be empty string, just in case)
	if (!$probeset_name) {
		$probeset_name = '';
		$probeset_clause = 'AND (op.probeset IS NULL OR op.probeset = ?)';
	}
	
	my $sth = $self->prepare("
		SELECT oligo_probe_id
		FROM oligo_probe op, oligo_array oa
		WHERE op.oligo_array_id = oa.oligo_array_id
		AND oa.name = ?
		$probeset_clause
		AND op.name = ?
	");
	
	$sth->bind_param(1, $array_name,    SQL_VARCHAR);
	$sth->bind_param(2, $probeset_name, SQL_VARCHAR);
	$sth->bind_param(3, $probe_name,    SQL_VARCHAR);
	$sth->execute();
	
	my ($probe_id) = $sth->fetchrow();
	
	if ($probe_id) {
		return $self->fetch_by_dbID($probe_id);
	} else {
		return undef;
	}
}

=head2 fetch_all_by_probeset

  Arg [1]    : string - probeset name
  Example    : my @probes = @{$opa->fetch_all_by_probeset('Probeset-1')};
  Description: Fetch all probes in a particular probeset.
  Returntype : Listref of Bio::EnsEMBL::OligoProbe objects
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_all_by_probeset {
	my $self     = shift;
	my $probeset = shift;

	return $self->generic_fetch("op.probeset = '$probeset'");
}

=head2 fetch_all_by_Array

  Arg [1]    : Bio::EnsEMBL::OligoArray
  Example    : my @probes = @{$opa->fetch_all_by_Array($array)};
  Description: Fetch all probes on a particular array.
  Returntype : Listref of Bio::EnsEMBL::OligoProbe objects.
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_all_by_Array {
	my $self  = shift;
	my $array = shift;
	
	if ( !ref($array) || !$array->isa('Bio::EnsEMBL::OligoArray') ) {
		warning('fetch_all_by_Array requires a Bio::EnsEMBL::OligoArray object');
		return [];
	}
	
	my $array_id = $array->dbID();
	if (!defined $array_id) {
		warning('fetch_all_by_Array requires a stored Bio::EnsEMBL::OligoArray object');
		return [];
	}
	
	return $self->generic_fetch("op.oligo_array_id = $array_id");
}

=head2 fetch_by_OligoFeature

  Arg [1]    : Bio::EnsEMBL::OligoFeature
  Example    : my $probe = $opa->fetch_by_OligoFeature($feature);
  Description: Returns the probe that created a particular feature.
  Returntype : Bio::EnsEMBL::OligoProbe
  Exceptions : Throws if argument is not a Bio::EnsEMBL::OligoFeature object
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_by_OligoFeature {
	my $self    = shift;
	my $feature = shift;
	
	if (
		!ref($feature)
		|| !$feature->isa('Bio::EnsEMBL::OligoFeature')
		|| !$feature->{'_probe_id'}
	) {
		throw('fetch_by_OligoFeature requires a stored Bio::EnsEMBL::OligoFeature object');
	}
	
	return $self->fetch_by_dbID($feature->{'_probe_id'});
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

  	return [ 'oligo_probe', 'op' ];
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

  return qw( op.oligo_probe_id op.oligo_array_id op.probeset op.name op.description op.length );

}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates OligoProbe objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::OligoProbe objects
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _objs_from_sth {
	my ($self, $sth) = @_;
	
	my (@result, $current_dbid, $probe_id, $array_id, $probeset, $name, $description, $probelength);
	my ($array, %array_cache);
	
	$sth->bind_columns( \$probe_id, \$array_id, \$probeset, \$name, \$description, \$probelength );
	
	my $probe;
	while ( $sth->fetch() ) {
		$array = $array_cache{$array_id} || $self->db->get_OligoArrayAdaptor()->fetch_by_dbID($array_id);
		if (!$current_dbid || $current_dbid != $probe_id) {
			# New probe
			$probe = Bio::EnsEMBL::OligoProbe->new(
				-name        => $name,
				-array       => $array,
				-probeset    => $probeset,
				-description => $description,
				-probelength => $probelength,
				-dbID        => $probe_id,
				-adaptor     => $self,
			);
			push @result, $probe;
			$current_dbid = $probe_id;
		} else {
			# Extend existing probe
			$probe->add_Array_probename($array, $name);
		}
	}
	return \@result;
}

=head2 store

  Arg [1]    : List of Bio::EnsEMBL::OligoProbe objects
  Example    : $opa->store($probe1, $probe2, $probe3);
  Description: Stores given OligoProbe objects in the database. Should only be
               called once per probe because no checks are made for duplicates.
			   Sets dbID and adaptor on the objects that it stores.
  Returntype : None
  Exceptions : Throws if arguments aren't OligoProbe objects
  Caller     : General
  Status     : Medium Risk

=cut

sub store {
	my ($self, @probes) = @_;

	if (scalar @probes == 0) {
		throw('Must call store with a list of OligoProbe objects');
	}

	my $db = $self->db();

	PROBE: foreach my $probe (@probes) {

		if ( !ref $probe || !$probe->isa('Bio::EnsEMBL::OligoProbe') ) {
			throw('Probe must be an OligoProbe object');
		}

		if ( $probe->is_stored($db) ) {
			warning('OligoProbe [' . $probe->dbID() . '] is already stored in the database');
			next PROBE;
		}
		
		# Get all the arrays this probe is on and check they're all in the database
		my $arrays = $probe->get_all_Arrays();
		my @stored_arrays;
		for my $array (@$arrays) {
			if ( defined $array->dbID() ) {
				push @stored_arrays, $array;
			}
		}
		if ( !@stored_arrays ) {
			warning('Probes need attached arrays to be stored in the database');
			next PROBE;
		}

		# Insert separate entry (with same oligo_probe_id) in oligo_probe
		# for each array the probe is on
		my $dbID;
		for my $array (@stored_arrays) {
			if (defined $dbID) {
				# Probe we've seen already
				my $sth = $self->prepare("
					INSERT INTO oligo_probe
					(oligo_probe_id, oligo_array_id, name, probeset, description, length)
					VALUES (?, ?, ?, ?, ?, ?)
				");
				$sth->bind_param(1, $dbID,                                 SQL_INTEGER);
				$sth->bind_param(2, $array->dbID(),                        SQL_INTEGER);
				$sth->bind_param(3, $probe->get_probename($array->name()), SQL_VARCHAR);
				$sth->bind_param(4, $probe->probeset(),                    SQL_VARCHAR);
				$sth->bind_param(5, $probe->description(),                 SQL_VARCHAR);
				$sth->bind_param(6, $probe->probelength(),                 SQL_INTEGER);
				$sth->execute();
			} else {
				# New probe
				my $sth = $self->prepare("
					INSERT INTO oligo_probe
					(oligo_array_id, name, probeset, description, length)
					VALUES (?, ?, ?, ?, ?)
				");
				$sth->bind_param(1, $array->dbID,                          SQL_INTEGER);
				$sth->bind_param(2, $probe->get_probename($array->name()), SQL_VARCHAR);
				$sth->bind_param(3, $probe->probeset(),                    SQL_VARCHAR);
				$sth->bind_param(4, $probe->description(),                 SQL_VARCHAR);
				$sth->bind_param(5, $probe->probelength(),                 SQL_INTEGER);
				$sth->execute();
				$dbID = $sth->{'mysql_insertid'};
				$probe->dbID($dbID);
				$probe->adaptor($self);
			}
		}
	}
}

=head2 list_dbIDs

  Arg [1]    : none
  Example    : my @feature_ids = @{$opa->list_dbIDs()};
  Description: Gets an array of internal IDs for all OligoProbe objects in the
               current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : Medium Risk

=cut

sub list_dbIDs {
	my ($self) = @_;

	return $self->_list_dbIDs('oligo_probe');
}

1;

