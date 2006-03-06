#
# Ensembl module for Bio::EnsEMBL::DBSQL::AffyFeatureAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::DBSQL::AffyFeatureAdaptor - A database adaptor for fetching and
storing AffyFeature objects.

=head1 SYNOPSIS

my $afa = $db->get_AffyFeatureAdaptor();

my $features = $afa->fetch_all_by_AffyProbe($probe);
$features = $afa->fetch_all_by_Slice_arrayname($slice, 'Affy-1', 'Affy-2');

=head1 DESCRIPTION

The AffyFeatureAdaptor is a database adaptor for storing and retrieving
AffyFeature objects.

=head1 AUTHOR

This module was originally written by Arne Stabenau, but was changed to be a
subclass of OligoFeatureAdaptor by Ian Sealy.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::DBSQL::AffyFeatureAdaptor;

use Bio::EnsEMBL::AffyFeature;
use Bio::EnsEMBL::DBSQL::OligoFeatureAdaptor;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::DBSQL::OligoFeatureAdaptor);


=head2 fetch_all_by_AffyProbe

  Arg [1]    : Bio::EnsEMBL::AffyProbe
  Example    : my $features = $afa->fetch_all_by_AffyProbe($probe);
  Description: Fetchs all features that a given probe creates.
  Returntype : Listref of Bio::EnsEMBL::AffyFeature objects
  Exceptions : None
  Caller     : AffyProbe->get_all_AffyFeatures()
  Status     : Medium Risk

=cut

sub fetch_all_by_AffyProbe {
	my $self  = shift;
	my $probe = shift;
	
	return $self->fetch_all_by_Probe($probe);
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
	
	return "
		    of.oligo_probe_id = op.oligo_probe_id
		AND op.oligo_array_id = oa.oligo_array_id
		AND oa.type='AFFY'
	";
}

=head2 _new_fast

  Args       : Hashref to be passed to AffyFeature->new_fast()
  Example    : None
  Description: Construct an AffyFeature object using quick and dirty new_fast.
  Returntype : Bio::EnsEMBL::AffyFeature
  Exceptions : None
  Caller     : _objs_from_sth
  Status     : Medium Risk

=cut

sub _new_fast {
	my $self = shift;
	
	my $hash_ref = shift;
	return Bio::EnsEMBL::AffyFeature->new_fast($hash_ref);
}

=head2 list_dbIDs

  Args       : None
  Example    : my @feature_ids = @{$afa->list_dbIDs()};
  Description: Gets an array of internal IDs for all AffyFeature objects in
               the current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : Medium Risk

=cut

sub list_dbIDs {
	my $self = shift;
	
	#return $self->_list_dbIDs('oligo_feature');
	# Can't use _list_dbIDs because only want OligoProbe objects on arrays of type AFFY
	
	my @out;
	my $sql = "
		SELECT DISTINCT of.oligo_feature_id
		FROM   oligo_feature of, oligo_probe op, oligo_array oa
		WHERE  of.oligo_probe_id=op.oligo_probe_id
		AND    op.oligo_array_id=oa.oligo_array_id
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

