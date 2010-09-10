package Bio::EnsEMBL::Utils::Scalar;

=pod

=head1 LICENSE

  Copyright (c) 1999-2010 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=pod

=head1 NAME

Bio::EnsEMBL::Utils::Scalar

=head1 SYNOPSIS

	use Bio::EnsEMBL::Utils::Scalar qw(check_ref assert_ref);
	
	check_ref([], 'ARRAY'); # Will return true
	check_ref({}, 'ARRAY'); # Will return false
	check_ref($dba, 'Bio::EnsEMBL::DBSQL::DBAdaptor'); #Returns true if $dba is a DBAdaptor
	
	assert_ref([], 'ARRAY'); #Returns true
	assert_ref({}, 'ARRAY'); #throws an exception
	assert_ref($dba, 'Bio::EnsEMBL::Gene'); #throws an exception if $dba is not a Gene  
	
=head1 DESCRIPTION

A collection of subroutines aimed to helping Scalar based operations

=head1 METHODS

See subroutines.

=head1 MAINTAINER

$Author$

=head1 VERSION

$Revision$

=cut

use strict;
use warnings;

use base qw(Exporter);

our @EXPORT_OK = qw(check_ref assert_ref);

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Scalar::Util qw(blessed);

=head2 check_ref()

  Arg [1]     : The reference to check
  Arg [2]     : The type we expect
  Description : A subroutine which checks to see if the given object/ref is 
                what you expect. If you give it a blessed reference then it 
                will perform an isa() call on the object after the defined 
                tests. If it is a plain reference then it will use ref().
                
                An undefined value will return a false.
  Returntype  : Boolean indicating if the reference was the type we 
                expect
  Example     : my $ok = check_ref([], 'ARRAY');
  Exceptions  : If the expected type was not set
  Status      : Stable

=cut

sub check_ref {
	my ($ref, $expected) = @_;
	throw('No expected type given') if ! defined $expected;
	if(defined $ref) {
		if(blessed($ref)) {
			return 1 if $ref->isa($expected);
		}
		else {
			my $ref_ref_type = ref($ref);
			return 1 if defined $ref_ref_type && $ref_ref_type eq $expected; 
		}
	}
	return 0;
}

=head2 assert_ref()

  Arg [1]     : The reference to check
  Arg [2]     : The type we expect
  Description : A subroutine which checks to see if the given object/ref is 
                what you expect. This behaves in an identical manner as
                C<check_ref()> does except this will raise exceptions when
                the values do not match rather than returning a boolean
                indicating the situation.
                
                Undefs cause exception circumstances.
  Returntype  : None
  Example     : assert_ref([], 'ARRAY');
  Exceptions  : If the expected type was not set and if the given reference
                was not assignable to the expected value
  Status      : Stable

=cut

sub assert_ref {
  my ($ref, $expected) = @_;
  throw('No expected type given') if ! defined $expected;
  my $class = ref($ref);
  throw('Given reference was undef') unless defined $ref;
  throw('Asking for the type of the reference produced no type; check your input is a reference') unless $class;
  if(blessed($ref)) {
    throw("Reference '${class}' is not an ISA of '${expected}'") if ! $ref->isa($expected);
  }
  else {    
    throw("'${expected}' expected class was not equal to actual class '${class}'") if $expected ne $class;
  }
  return 1;
}

1;
