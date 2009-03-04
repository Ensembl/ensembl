=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::AffyProbeAdaptor - A database adaptor for fetching and
storing AffyProbe objects.

=head1 SYNOPSIS

  my $apa = $db->get_AffyProbeAdaptor();

  my $probe =
    $opa->fetch_by_array_probeset_probe( 'Affy-1', 'Probeset-1',
    'Probe-1' );

=head1 DESCRIPTION

The AffyProbeAdaptor is a database adaptor for storing and retrieving
AffyProbe objects.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::DBSQL::AffyProbeAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning );
use Bio::EnsEMBL::AffyProbe;
use Bio::EnsEMBL::DBSQL::OligoProbeAdaptor;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::DBSQL::OligoProbeAdaptor);


=head2 fetch_all_by_AffyArray

  Arg [1]    : Bio::EnsEMBL::AffyArray
  Example    : my @probes = @{$apa->fetch_all_by_AffyArray($array)};
  Description: Fetch all probes on a particular array.
  Returntype : Listref of Bio::EnsEMBL::AffyProbe objects.
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_all_by_AffyArray {
	my $self  = shift;
	my $array = shift;

	return $self->fetch_all_by_Array($array);
}

=head2 fetch_by_AffyArray

  Arg [1]    : Bio::EnsEMBL::AffyArray
  Example    : my @probes = @{$apa->fetch_by_AffyArray($array)};
  Description: Fetch all probes on a particular array. This method was misnamed
               and is deprecated.
  Returntype : Listref of Bio::EnsEMBL::AffyProbe objects.
  Exceptions : None
  Caller     : General
  Status     : At Risk
             : Likely to be removed - use fetch_all_by_AffyArray instead

=cut

sub fetch_by_AffyArray {
  my $self  = shift;
  my $array = shift;

  deprecate('fetch_all_by_Array() should be used instead of the deprecated fetch_by_AffyArray().');

  return $self->fetch_all_by_Array($array);
}

=head2 fetch_by_AffyFeature

  Arg [1]    : Bio::EnsEMBL::AffyFeature
  Example    : my $probe = $apa->fetch_by_AffyFeature($feature);
  Description: Returns the probe that created a particular feature.
  Returntype : Bio::EnsEMBL::AffyProbe
  Exceptions : Throws if argument is not a Bio::EnsEMBL::AffyFeature object
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_by_AffyFeature {
  my $self    = shift;
  my $feature = shift;

  return $self->fetch_by_OligoFeature($feature);
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates OligoProbe objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::AffyProbe objects
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
		$array = $array_cache{$array_id} || $self->db->get_AffyArrayAdaptor()->fetch_by_dbID($array_id);
		if (!$current_dbid || $current_dbid != $probe_id) {
			# New probe
			$probe = Bio::EnsEMBL::AffyProbe->new(
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

=head2 list_dbIDs

  Arg [1]    : none
  Example    : my @probe_ids = @{$apa->list_dbIDs()};
  Description: Gets an array of internal IDs for all AffyProbe objects
               in the current database.  NOTE: In a multi-species
               database, this method will return the dbIDs of all
               AffyProbe objects, not just the ones associated with
               the current species.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : Medium Risk

=cut

sub list_dbIDs {
	my ($self) = @_;

	#return $self->_list_dbIDs('oligo_probe');
	# Can't use _list_dbIDs because only want OligoProbe objects on arrays of type AFFY
	
	my @out;

        # FIXME: This SQL will not work as expected on multi-species
        # databases.  It needs to be anchored in a coord_system entry
        # coord_system.species_id = $self->species_id(). /ak4@2008-07-15
	my $sql = "
		SELECT DISTINCT op.oligo_probe_id
		FROM   oligo_probe op, oligo_array oa
		WHERE  op.oligo_array_id=oa.oligo_array_id
		AND    oa.type='AFFY'
	";
	my $sth = $self->prepare($sql);
	$sth->execute;
	
	while (my ($id) = $sth->fetchrow() ) {
		push @out, $id;
	}
	
	$sth->finish;
	
	return \@out;
}

1;

