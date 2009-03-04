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

Bio::EnsEMBL::AffyProbe - A module to represent an Affymetrix probe.

=head1 SYNOPSIS

  use Bio::EnsEMBL::OligoProbe;

  my $probe = Bio::EnsEMBL::OligoProbe->new(
    -PROBENAME => 'Probe-1',
    -ARRAY     => $array,
    -PROBESET  => 'Probeset-1',
  );

=head1 DESCRIPTION

An AffyProbe object represents an Affy probe on a microarray.  The
data (currently the name, array, length, probeset and description) are
stored in the oligo_probe table.  The complete name of a probe is the
concatenation of the array name, the probeset and the probe name.

For Affy arrays, a probe can be part of more than one array, but
only part of one probeset.  On each Affy array the probe has a
slightly different name.  For example, two different complete names
for the same probe might be DrosGenome1:AFFX-LysX-5_at:535:35; and
Drosophila_2:AFFX-LysX-5_at:460:51;.  In the database, these two probes
will have the same oligo_probe_id.  Thus the same Affy probe can have a
number of different names and complete names depending on which array
it's on.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::AffyProbe;

use Bio::EnsEMBL::Utils::Exception qw( warning );
use Bio::EnsEMBL::OligoProbe;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::OligoProbe);


=head2 new

  Arg [-PROBENAME]  : string - probe name
        Used when the probe is on one array. Can also use -NAME.
  Arg [-PROBENAMES] : Listref of strings - probe names
        Used when the probe is on multiple arrays. Can also use -NAMES.
  Arg [-ARRAY]      : Bio::EnsEMBL::AffyArray
   or [-ARRAYNAME]  : string - array name
        Used when the probe is on one array. Either -ARRAY or -ARRAYNAME can be
		used. The latter is a convenience and -ARRAY must be used if the probe
		is to be stored.
  Arg [-ARRAYS]     : Listref of Bio::EnsEMBL::AffyArray objects
   or [-ARRAYNAMES] : Listref of strings - array names
        Used when the probe is on multiple arrays. Either -ARRAYS or
		-ARRAYNAMES can be used. The latter is a convenience and -ARRAYS must
		be used if the probes are to be stored.
  Arg [-PROBESET]   : string - probeset name
        Each probe is part of one (and only one) probeset.
  Arg [-PROBELENGTH]: int - probe length
        Like probesets, will obviously be the same for all probes if same probe
		is on multiple arrays.
  Arg [-DESCRIPTION]: string - probe description
        Like probesets, will be the same for all probes if same probe is on
		multiple arrays.
  Example    : my $probe = Bio::EnsEMBL::AffyProbe->new(
                   -PROBENAME   => 'Probe-1',
				   -ARRAY       => $array,
				   -PROBESET    => 'Probeset-1',
               );
  Description: Creates a new Bio::EnsEMBL::AffyProbe object.
  Returntype : Bio::EnsEMBL::AffyProbe
  Exceptions : Throws if not supplied with probe name(s) and array(s)
  Caller     : General
  Status     : Medium Risk

=cut

sub new {
	my $caller = shift;

	my $class = ref($caller) || $caller;

	my $self = $class->SUPER::new(@_);
	
	# Affy probes are 25mers (unless otherwise specified)
	if ( !$self->probelength() ) {
		$self->probelength(25);
	}
	
	return $self;
}

=head2 get_all_AffyFeatures

  Args       : None
  Example    : my $features = $probe->get_all_AffyFeatures();
  Description: Get all features produced by this probe. The probe needs to be
               database persistent.
  Returntype : Listref of Bio::EnsEMBL:AffyFeature objects
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_all_AffyFeatures {
	my $self = shift;
	if ( $self->adaptor() && $self->dbID() ) {
		return $self->adaptor()->db()->get_AffyFeatureAdaptor()->fetch_all_by_AffyProbe($self);
	} else {
		warning('Need database connection to retrieve Features');
		return [];
	}    
}

=head2 get_all_AffyArrays

  Args       : None
  Example    : my $arrays = $probe->get_all_AffyArrays();
  Description: Returns all arrays that this probe is part of. Only works if the
               probe was retrieved from the database or created using
			   add_Array_probename (rather than add_arrayname_probename).
  Returntype : Listref of Bio::EnsEMBL::AffyArray objects
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_all_AffyArrays {
	my $self = shift;
	
	return $self->get_all_Arrays();
}

1;

