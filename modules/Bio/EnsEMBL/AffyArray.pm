#
# Ensembl module for Bio::EnsEMBL::AffyArray
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::AffyArray - A module to represent an Affymetrix array.

=head1 SYNOPSIS

use Bio::EnsEMBL::AffyArray;

my $array = Bio::EnsEMBL::AffyArray->new(
	-NAME        => 'Affy-1',
	-INCLUDED_IN => $another_array,
	-SETSIZE     => 16,
);

=head1 DESCRIPTION

An AffyArray object represents an Affymetrix array. The data (currently the
name, probeset size and parent array) are stored in the oligo_array table.

Each array can have a parent array (another array that contains all the probes
of this array). This is rarely (if ever) used.

=head1 AUTHOR

This module was originally written by Arne Stabenau, but was changed to be a
subclass of OligoArray by Ian Sealy.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::AffyArray;

use Bio::EnsEMBL::OligoArray;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::OligoArray);


=head2 new

  Arg [-NAME]:
        string - the name of this array
  Arg [-INCLUDED_IN]:
        (optional) Bio::EnsEMBL::OligoArray - a possible superset array
  Arg [-SETSIZE]: 
        int - the number of probes in a probe set
  Example    : my $array = Bio::EnsEMBL::AffyArray->new(
                   -NAME        => 'Affy-1',
				   -INCLUDED_IN => $another_array,
				   -SETSIZE     => 16,
               );
  Description: Creates a new Bio::EnsEMBL::AffyArray object.
  Returntype : Bio::EnsEMBL::AffyArray
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub new {
	my $caller = shift;

	my $class = ref($caller) || $caller;

	my $self = $class->SUPER::new(@_);
	
	# All AffyArray objects are OligoArray objects of type AFFY
	$self->type('AFFY');
	
	return $self;
}


=head2 get_all_AffyProbes

  Args       : None
  Example    : my $probes = $array->get_all_AffyProbes();
  Description: Returns all probes on an Affy array. Needs a database
               connection.
  Returntype : Listref of Bio::EnsEMBL::AffyProbe objects
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_all_AffyProbes {
  my $self = shift;
  
  return $self->get_all_Probes();
}

1;

