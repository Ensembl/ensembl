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

Bio::EnsEMBL::OligoProbe - A module to represent an oligonucleotide probe.

=head1 SYNOPSIS

  use Bio::EnsEMBL::OligoProbe;

  my $probe = Bio::EnsEMBL::OligoProbe->new(
    -PROBENAME => 'Probe-1',
    -ARRAY     => $array,
  );

=head1 DESCRIPTION

An OligoProbe object represents an oligonucleotide probe on a
microarray. The data (currently the name, array, length, probeset and
description) are stored in the oligo_probe table. Probeset is only
really relevant for Affy probes.  The complete name of a probe is the
concatenation of the array name, the probeset (if relevant) and the
probe name.

For Affy arrays, a probe can be part of more than one array, but
only part of one probeset. On each Affy array the probe has a
slightly different name. For example, two different complete names
for the same probe might be DrosGenome1:AFFX-LysX-5_at:535:35; and
Drosophila_2:AFFX-LysX-5_at:460:51;. In the database, these two probes
will have the same oligo_probe_id. Thus the same Affy probe can have a
number of different names and complete names depending on which array
it's on.

=head1 METHODS

=cut

package Bio::EnsEMBL::OligoProbe;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw( rearrange ) ;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Storable;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

  Arg [-PROBENAME]  : string - probe name
        Used when the probe is on one array. Can also use -NAME.
  Arg [-PROBENAMES] : Listref of strings - probe names
        Used when the probe is on multiple arrays. Can also use -NAMES.
  Arg [-ARRAY]      : Bio::EnsEMBL::OligoArray
   or [-ARRAYNAME]  : string - array name
        Used when the probe is on one array. Either -ARRAY or -ARRAYNAME can be
		used. The latter is a convenience and -ARRAY must be used if the probe
		is to be stored.
  Arg [-ARRAYS]     : Listref of Bio::EnsEMBL::OligoArray objects
   or [-ARRAYNAMES] : Listref of strings - array names
        Used when the probe is on multiple arrays. Either -ARRAYS or
		-ARRAYNAMES can be used. The latter is a convenience and -ARRAYS must
		be used if the probes are to be stored.
  Arg [-PROBESET]   : string - probeset name
        Probably only needed for arrays of type AFFY. Each probe is part of one
		(and only one) probeset.
  Arg [-PROBELENGTH]: int - probe length
        Like probesets, will obviously be the same for all probes if same probe
		is on multiple arrays.
  Arg [-DESCRIPTION]: string - probe description
        Like probesets, will be the same for all probes if same probe is on
		multiple arrays.
  Example    : my $probe = Bio::EnsEMBL::OligoProbe->new(
                   -PROBENAME   => 'Probe-1',
				   -ARRAY       => $array,
				   -PROBELENGTH => 65,
               );
  Description: Creates a new Bio::EnsEMBL::OligoProbe object.
  Returntype : Bio::EnsEMBL::OligoProbe
  Exceptions : Throws if not supplied with probe name(s) and array(s)
  Caller     : General
  Status     : Medium Risk

=cut

sub new {
	my $caller = shift;
	
	my $class = ref($caller) || $caller;
	
	my $self = $class->SUPER::new(@_);
	
	my (
		$arrays,      $array,
		$arraynames,  $arrayname,
		$probenames,  $probename,
		$names,       $name,
		$probeset,    $description,
		$probelength,
	) = rearrange([
		'ARRAYS',     'ARRAY',
		'ARRAYNAMES', 'ARRAYNAME',
		'PROBENAMES', 'PROBENAME',
		'NAMES',      'NAME',
		'PROBESET',   'DESCRIPTION',
		'PROBELENGTH',
	], @_);
	
	# Originally arguments were inconsistent: NAME and PROBENAMES
	# Now can be NAME and NAMES or PROBENAME and PROBENAMES
	$probename  ||= $name;
	$probenames ||= $names;
	
	if ($probenames && ref($probenames) eq 'ARRAY') {
		# Multiple probes (all from the same probeset) have been specified
		my $probecount = scalar @$probenames;
		if ($arrays && ref($arrays) eq 'ARRAY') {
			if (scalar @$arrays != $probecount) {
				throw('Need a probename for each array');
			}
			# Array objects
			for (my $i = 0; $i < scalar @$probenames; $i++) {
				$self->add_Array_probename($arrays->[$i], $probenames->[$i]);
			}
		} elsif ($arraynames && ref($arraynames) eq 'ARRAY') {
			if (scalar @$arraynames != $probecount) {
				throw('Need a probename for each array');
			}
			# Named arrays
			for (my $i = 0; $i < scalar @$probenames; $i++) {
				$self->add_arrayname_probename($arraynames->[$i], $probenames->[$i]);
			}
		} else {
			throw('No arrays or arraynames have been supplied for this OligoProbe');
		}
	} elsif (defined $probename) {
		# Single probe specified
		if (defined $arrayname) {
			# Named array
			$self->add_arrayname_probename($arrayname, $probename);
		} elsif(defined $array) {
			# Array object
			$self->add_Array_probename($array, $probename);
		} else {
			throw('No array or arrayname has been supplied for this OligoProbe');
		}
	} else {
		throw('You need to provide a probe name (or names) to create an OligoProbe');
	}
		
	$self->probeset($probeset)       if defined $probeset;
	$self->description($description) if defined $description;
	$self->probelength($probelength) if defined $probelength;
	
	return $self;
}

=head2 add_Array_probename

  Arg [1]    : Bio::EnsEMBL::AffyArray - array
  Arg [2]    : string - probe name
  Example    : $probe->add_Array_probename($array, $probename);
  Description: Adds a probe name / array pair to a probe, allowing incremental
               generation of a probe.
  Returntype : None
  Exceptions : None
  Caller     : General,
               OligoProbe->new(),
               OligoProbeAdaptor->_obj_from_sth(),
			   AffyProbeAdaptor->_obj_from_sth()
  Status     : Medium Risk

=cut

sub add_Array_probename {
    my $self = shift;
    my ($array, $probename) = @_;
    $self->{ 'arrays'     } ||= {};
    $self->{ 'arrays'     }->{$array->name()} = $array;
    $self->{ 'probenames' }->{$array->name()} = $probename;
}

=head2 add_arrayname_probename

  Arg [1]    : string - array name
  Arg [2]    : string - probe name
  Example    : $probe->add_arrayname_probename($arrayname, $probename);
  Description: Adds a probe name / array pair to a probe. Such probes cannot be
               stored, so it's better to use add_Array_probename(). 
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub add_arrayname_probename {
    my $self = shift;
    my ($arrayname, $probename) = @_;
    $self->{'probenames'}->{$arrayname} = $probename;
}

=head2 get_all_OligoFeatures

  Args       : None
  Example    : my $features = $probe->get_all_OligoFeatures();
  Description: Get all features produced by this probe. The probe needs to be
               database persistent.
  Returntype : Listref of Bio::EnsEMBL:OligoFeature objects
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_all_OligoFeatures {
	my $self = shift;
	if ( $self->adaptor() && $self->dbID() ) {
		return $self->adaptor()->db()->get_OligoFeatureAdaptor()->fetch_all_by_Probe($self);
	} else {
		warning('Need database connection to retrieve Features');
		return [];
	}    
}

=head2 get_all_Arrays

  Args       : None
  Example    : my $arrays = $probe->get_all_Arrays();
  Description: Returns all arrays that this probe is part of. Only works if the
               probe was retrieved from the database or created using
			   add_Array_probename (rather than add_arrayname_probename).
  Returntype : Listref of Bio::EnsEMBL::OligoArray objects
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_all_Arrays {
    my $self = shift;
	# Do we have OligoArray objects for this probe?
    if (defined $self->{'arrays'}) {
		return [ values %{$self->{'arrays'}} ];
    } elsif ( $self->adaptor() && $self->dbID() ) { 
		# Only have names for arrays, so need to retrieve arrays from database
		warning('Not yet implemented');
		return [];
    } else {
		warning('Need database connection to get Arrays by name');
		return [];
    }
}

=head2 get_all_probenames

  Args       : None
  Example    : my @probenames = @{$probe->get_all_probenames()};
  Description: Retrieves all names for this probe. Only makes sense for probes
               that are part of a probeset (i.e. Affy probes), in which case
			   get_all_complete_names() would be more appropriate.
  Returntype : Listref of strings
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_all_probenames {
    my $self = shift;
    return [ values %{$self->{'probenames'}} ]
}

=head2 get_probename

  Arg [1]    : string - array name
  Example    : my $probename = $probe->get_probename('Array-1');
  Description: For a given array, retrieve the name for this probe.
  Returntype : string
  Exceptions : Throws if the array name is not known for this probe
  Caller     : General
  Status     : Medium Risk

=cut

sub get_probename {
    my $self = shift;
    my $arrayname = shift;
	
    my $probename = $self->{'probenames'}->{$arrayname};
    if (!defined $probename) {
		throw('Unknown array name');
    }
	
    return $probename;
}

=head2 get_all_complete_names

  Args       : None
  Example    : my @compnames = @{$probe->get_all_complete_names()};
  Description: Retrieves all complete names for this probe. The complete name
               is a concatenation of the array name, the probeset name and the
			   probe name.
  Returntype : Listref of strings
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_all_complete_names {
    my $self = shift;
	
    my @result = ();
	
    my $probeset = $self->probeset();
    $probeset .= ':' if $probeset;
    $probeset = "" if !$probeset;

    while ( my ($arrayname, $probename) = each %{$self->{'probenames'}} ) {
      push @result, "$arrayname:$probeset$probename";
    }
	
    return \@result;
}

=head2 get_complete_name

  Arg [1]    : string - array name
  Example    : my $compname = $probe->get_complete_name('Array-1');
  Description: For a given array, retrieve the complete name for this probe.
  Returntype : string
  Exceptions : Throws if the array name is not known for this probe
  Caller     : General
  Status     : Medium Risk

=cut

sub get_complete_name {
    my $self = shift;
    my $arrayname = shift;

    my $probename = $self->{'probenames'}->{$arrayname};
    if (!defined $probename) {
		throw('Unknown array name');
    }
	
	my $probeset = $self->probeset();
	$probeset .= ':' if $probeset;
	
    return "$arrayname:$probeset$probename";
}

=head2 probeset

  Arg [1]    : (optional) string - probeset
  Example    : my $probeset = $probe->probeset();
  Description: Getter and setter of probeset attribute for OligoProbe objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub probeset {
    my $self = shift;
    $self->{'probeset'} = shift if @_;
    return $self->{'probeset'};
}

=head2 description

  Arg [1]    : (optional) string - description
  Example    : my $description = $probe->description();
  Description: Getter and setter of description attribute for OligoProbe
               objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub description {
    my $self = shift;
    $self->{'description'} = shift if @_;
    return $self->{'description'};
}

=head2 probelength

  Arg [1]    : (optional) int - probe length
  Example    : my $probelength = $probe->probelength();
  Description: Getter and setter of probelength attribute for OligoProbe
               objects.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub probelength {
    my $self = shift;
    $self->{'probelength'} = shift if @_;
    return $self->{'probelength'};
}

1;

