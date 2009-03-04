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

Bio::EnsEMBL::OligoFeature - A module to represent an oligonucleotide probe's
genomic mapping.

=head1 SYNOPSIS

  use Bio::EnsEMBL::OligoFeature;

  my $feature = Bio::EnsEMBL::OligoFeature->new(
    -PROBE         => $probe,
    -MISMATCHCOUNT => 0,
    -SLICE         => $chr_1_slice,
    -START         => 1_000_000,
    -END           => 1_000_024,
    -STRAND        => -1,
  );

=head1 DESCRIPTION

An OligoFeature object represents the genomic placement of an OligoProbe
object. The data are stored in the oligo_feature table.

=head1 METHODS

=cut

package Bio::EnsEMBL::OligoFeature;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Feature;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Feature);


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
	
	my ($probe, $mismatchcount )
		= rearrange(['PROBE', 'MISMATCHCOUNT'], @_);
	
	$self->probe($probe);
	$self->mismatchcount($mismatchcount);
	
	return $self;
}

=head2 new_fast

  Args       : Hashref with all internal attributes set
  Example    : none
  Description: Quick and dirty version of new. Only works if the code is very
               disciplined.
  Returntype : Bio::EnsEMBL::OligoFeature
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub new_fast {
   my ($class, $hashref)  = @_;

   return bless ($hashref, $class);
}

=head2 probeset

  Arg [1]    : (optional) string - probeset
  Example    : my $probeset = $feature->probeset();
  Description: Getter and setter for the probeset for this feature. Shortcut
               for $feature->probe->probeset(), which should be used instead.
			   Probeset is not persisted if set with this method.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk
             : Use $feature->probe->probeset() because this may be removed

=cut

sub probeset {
    my $self = shift;
	
    $self->{'probeset'} = shift if @_;
	
    if ($self->{'probe'}) {
		$self->{'probeset'} = $self->probe()->probeset();
    }
	
    return $self->{'probeset'};
}

=head2 mismatchcount

  Arg [1]    : int - number of mismatches
  Example    : my $mismatches = $feature->mismatchcount();
  Description: Getter and setter for number of mismatches for this feature.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub mismatchcount {
    my $self = shift;
	
    $self->{'mismatchcount'} = shift if @_;
	
    return $self->{'mismatchcount'};
}

=head2 probelength

  Args       : None 
  Example    : my $probelength = $feature->probelength();
  Description: Getter for the length of the probe. Shortcut for
               $feature->probe->probelength(), which should be used instead.
			   Originally, this method returned the length of the feature,
			   which was often, but not always, the same as the length of the
			   probe.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : Medium Risk
             : Use $feature->probe->probelength() because this may be removed

=cut

sub probelength {
    my $self = shift;
	
    return $self->probe->probelength();
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
		if ( !ref $probe || !$probe->isa('Bio::EnsEMBL::OligoProbe') ) {
			throw('Probe must be a Bio::EnsEMBL::OligoProbe object');
		}
		$self->{'probe'} = $probe;
    }
	if ( !defined $self->{'probe'} && $self->dbID() && $self->adaptor() ) {
	    $self->{'probe'} = $self->adaptor()->db()->get_OligoProbeAdaptor()->fetch_by_OligoFeature($self);
	}
    return $self->{'probe'};
}

1;

