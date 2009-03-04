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

Bio::EnsEMBL::AffyFeature - A module to represent an Affy probe's genomic
mapping.

=head1 SYNOPSIS

  use Bio::EnsEMBL::AffyFeature;

  my $feature = Bio::EnsEMBL::AffyFeature->new(
    -PROBE         => $probe,
    -MISMATCHCOUNT => 0,
    -SLICE         => $chr_1_slice,
    -START         => 1_000_000,
    -END           => 1_000_024,
    -STRAND        => -1,
  );

=head1 DESCRIPTION

An AffyFeature object represents the genomic placement of an AffyProbe
object. The data are stored in the oligo_feature table.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::AffyFeature;

use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::OligoFeature;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::OligoFeature);


=head2 new

  Arg [-PROBE]        : Bio::EnsEMBL::OligoProbe - probe
        An OligoFeature must have a probe. This probe must already be stored if
		you plan to store the feature.
  Arg [-MISMATCHCOUNT]: int
        Number of mismatches over the length of the probe. 
  Arg [-SLICE]        : Bio::EnsEMBL::Slice
        The slice on which this feature is.
  Arg [-START]        : int
        The start coordinate of this feature relative to the start of the slice
		it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-END]          : int
        The end coordinate of this feature relative to the start of the slice
		it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-STRAND]       : int
        The orientation of this feature. Valid values are 1, -1 and 0.
  Arg [-dbID]         : (optional) int
        Internal database ID.
  Arg [-ADAPTOR]      : (optional) Bio::EnsEMBL::DBSQL::BaseAdaptor
        Database adaptor.
  Example    : my $feature = Bio::EnsEMBL::OligoFeature->new(
				   -PROBE         => $probe,
				   -MISMATCHCOUNT => 0,
				   -SLICE         => $chr_1_slice,
				   -START         => 1_000_000,
				   -END           => 1_000_024,
				   -STRAND        => -1,
			   ); 
  Description: Constructor for OligoFeature objects.
  Returntype : Bio::EnsEMBL::OligoFeature
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub new {
	my $caller = shift;

	my $class = ref($caller) || $caller;

	my $self = $class->SUPER::new(@_);
	
	return $self;
}

=head2 probe

  Arg [1]    : Bio::EnsEMBL::AffyProbe - probe
  Example    : my $probe = $feature->probe();
  Description: Getter, setter and lazy loader of probe attribute for
               OligoFeature objects. Features are retrieved from the database
			   without attached probes, so retrieving probe information for a
			   feature will involve another query.
  Returntype : Bio::EnsEMBL::OligoProbe
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub probe {
    my $self = shift;
	my $probe = shift;
    if ($probe) {
		if ( !ref $probe || !$probe->isa('Bio::EnsEMBL::AffyProbe') ) {
			throw('Probe must be a Bio::EnsEMBL::AffyProbe object');
		}
		$self->{'probe'} = $probe;
    }
	if ( !defined $self->{'probe'} && $self->dbID() && $self->adaptor() ) {
	    $self->{'probe'} = $self->adaptor()->db()->get_AffyProbeAdaptor()->fetch_by_AffyFeature($self);
	}
    return $self->{'probe'};
}

1;

