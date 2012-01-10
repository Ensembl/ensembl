package Bio::EnsEMBL::Utils::Scalar;

=pod

=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
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

	use Bio::EnsEMBL::Utils::Scalar qw(check_ref assert_ref wrap_array check_ref_can assert_ref_can assert_numeric assert_integer);

	check_ref([], 'ARRAY'); # Will return true
	check_ref({}, 'ARRAY'); # Will return false
	check_ref($dba, 'Bio::EnsEMBL::DBSQL::DBAdaptor'); #Returns true if $dba is a DBAdaptor

	assert_ref([], 'ARRAY'); #Returns true
	assert_ref({}, 'ARRAY'); #throws an exception
	assert_ref($dba, 'Bio::EnsEMBL::Gene'); #throws an exception if $dba is not a Gene

	wrap_array([]); #Returns the same reference
	wrap_array($a); #Returns [$a] if $a was not an array
	wrap_array(undef); #Returns [] since incoming was undefined
	wrap_array(); #Returns [] since incoming was empty (therefore undefined)

	check_ref_can([], 'dbID'); #returns false as ArrayRef is not blessed
	check_ref_can($gene, 'dbID'); #returns true as Gene should implement dbID()
	check_ref_can(undef); #Throws an exception as we gave no method to test

	assert_ref_can([], 'dbID'); #throws an exception since ArrayRef is not blessed
	assert_ref_can($gene, 'dbID'); #returns true if gene implements dbID()
	assert_ref_can(undef); #Throws an exception as we gave no method to test

	asssert_integer(1, 'dbID'); #Passes
	asssert_integer(1.1, 'dbID'); #Fails
	asssert_numeric(1E-11, 'dbID'); #Passes
	asssert_numeric({}, 'dbID'); #Fails

	#Tags are also available for exporting
	use Bio::EnsEMBL::Utils::Scalar qw(:assert); # brings in all assert methods
	use Bio::EnsEMBL::Utils::Scalar qw(:check); #brings in all check methods
	use Bio::EnsEMBL::Utils::Scalar qw(:array); #brings in wrap_array
	use Bio::EnsEMBL::Utils::Scalar qw(:all); #import all methods

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

our %EXPORT_TAGS;
our @EXPORT_OK;

@EXPORT_OK = qw(
  check_ref check_ref_can
  assert_ref assert_ref_can assert_numeric assert_integer assert_boolean assert_strand assert_file_handle
  wrap_array
);
%EXPORT_TAGS = (
  assert  => [qw(assert_ref assert_ref_can assert_integer assert_numeric assert_boolean assert_strand assert_file_handle)],
  check   => [qw(check_ref check_ref_can)],
  array   => [qw/wrap_array/],
  all     => [@EXPORT_OK]
);

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Scalar::Util qw(blessed looks_like_number openhandle);

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
  Arg [3]     : The attribute name you are asserting; not required but allows
                for more useful error messages to be generated. Defaults to
                C<-Unknown->.
  Description : A subroutine which checks to see if the given object/ref is
                what you expect. This behaves in an identical manner as
                C<check_ref()> does except this will raise exceptions when
                the values do not match rather than returning a boolean
                indicating the situation.

                Undefs cause exception circumstances.
  Returntype  : Boolean; true if we managed to get to the return
  Example     : assert_ref([], 'ARRAY');
  Exceptions  : If the expected type was not set and if the given reference
                was not assignable to the expected value
  Status      : Stable

=cut

sub assert_ref {
  my ($ref, $expected, $attribute_name) = @_;
  $attribute_name ||= '-Unknown-';
  throw('No expected type given') if ! defined $expected;
  my $class = ref($ref);
  throw("The given reference for attribute $attribute_name was undef") unless defined $ref;
  throw("Asking for the type of the attribute $attribute_name produced no type; check it is a reference") unless $class;
  if(blessed($ref)) {
    throw("${attribute_name}'s type '${class}' is not an ISA of '${expected}'") if ! $ref->isa($expected);
  }
  else {
    throw("$attribute_name was expected to be '${expected}' but was '${class}'") if $expected ne $class;
  }
  return 1;
}

=head2 wrap_array()

  Arg         : The reference we want to wrap in an array
  Description : Takes in a reference and returns either the reference if it
                was already an array, the reference wrapped in an array or
                an empty array (if the given value was undefined).
  Returntype  : Array Reference
  Example     : my $a = wrap_array($input);
  Exceptions  : None
  Status      : Stable

=cut

sub wrap_array {
  my ($incoming_reference) = @_;
  if(defined $incoming_reference) {
    if(check_ref($incoming_reference, 'ARRAY')) {
      return $incoming_reference;
    }
    else {
      return [$incoming_reference];
    }
  }
  return [];
}

=head2 check_ref_can

  Arg [1]     : The reference to check
  Arg [2]     : The method we expect to run
  Description : A subroutine which checks to see if the given object/ref is
                implements the given method. This is very similar to the
                functionality given by C<UNIVERSAL::can()> but works
                by executing C<can()> on the object meaning we consult the
                object's potentially overriden version rather than Perl's
                default mechanism.
  Returntype  : CodeRef
  Example     : check_ref_can($gene, 'dbID');
  Exceptions  : If the expected type was not set.
  Status      : Stable

=cut

sub check_ref_can {
  my ($ref, $method) = @_;
  throw('No method given') if ! defined $method;
  return unless defined $ref && blessed($ref);
  return $ref->can($method);
}

=head2 assert_ref_can

  Arg [1]     : The reference to check
  Arg [2]     : The method we expect to run
  Arg [3]     : The attribute name you are asserting; not required but allows
                for more useful error messages to be generated. Defaults to
                C<-Unknown->.
  Description : A subroutine which checks to see if the given object/ref is
                implements the given method. Will throw exceptions.
  Returntype  : Boolean; true if we managed to get to the return
  Example     : assert_ref_can($gene, 'dbID');
  Exceptions  : If the reference is not defined, if the object does not
                implement the given method and if no method was given to check
  Status      : Stable

=cut

sub assert_ref_can {
  my ($ref, $method, $attribute_name) = @_;
  $attribute_name ||= '-Unknown-';
  throw('No method given') if ! defined $method;
  throw "The given reference $attribute_name is not defined" unless defined $ref;
  throw "The given reference $attribute_name is not blessed" unless blessed($ref);
  if(! $ref->can($method)) {
    my $str_ref = ref($ref);
    throw sprintf(q{The given blessed reference '%s' for attribute '%s' does not implement the method '%s'}, $str_ref, $attribute_name, $method);
  }
  return 1;
}

=head2 assert_numeric

  Arg [1]     : The Scalar to check
  Arg [2]     : The attribute name you are asserting; not required but allows
                for more useful error messages to be generated. Defaults to
                C<-Unknown->.
  Description : A subroutine which checks to see if the given scalar is
                number or not. If not then we raise exceptions detailing why
  Returntype  : Boolean; true if we had a numeric otherwise we signal failure
                via exceptions
  Example     : assert_numeric(1, 'dbID');
  Exceptions  : If the Scalar is not defined, if the Scalar was blessed and
                if the value was not a number
  Status      : Stable

=cut

sub assert_numeric {
  my ($integer, $attribute_name) = @_;
  $attribute_name ||= '-Unknown-';
  throw "$attribute_name attribute is undefined" if ! defined $integer;
  throw "The given attribute $attribute_name is blessed; cannot work with blessed values" if blessed($integer);
  if(! looks_like_number($integer)) {
    throw "Attribute $attribute_name was not a number";
  }
  return 1;
}

=head2 assert_integer

  Arg [1]     : The Scalar to check
  Arg [2]     : The attribute name you are asserting; not required but allows
                for more useful error messages to be generated. Defaults to
                C<-Unknown->.
  Description : A subroutine which checks to see if the given scalar is
                a whole integer; we delegate to L<assert_numeric> for number
                checking.
  Returntype  : Boolean; true if we had a numeric otherwise we signal failure
                via exceptions
  Example     : assert_integer(1, 'dbID');
  Exceptions  : See L<assert_numeric> and we raise exceptions if the value
                was not a whole integer
  Status      : Stable

=cut

sub assert_integer {
  my ($integer, $attribute_name) = @_;
  $attribute_name ||= '-Unknown-';
  assert_numeric($integer, $attribute_name);
  if($integer != int($integer)) {
    throw "Attribute $attribute_name was a number but not an Integer";
  }
  return 1;
}

=head2 assert_boolean

  Arg [1]     : The Scalar to check
  Arg [2]     : The attribute name you are asserting; not required but allows
                for more useful error messages to be generated. Defaults to
                C<-Unknown->.
  Description : A subroutine which checks to see if the given scalar is
                a boolean i.e. value is set to C<1> or C<0>
  Returntype  : Boolean; true if we were given a boolean otherwise we signal
                failure via exceptions
  Example     : assert_boolean(1, 'is_circular');
  Exceptions  : See L<assert_integer> and we raise exceptions if the value
                was not equal to the 2 valid states
  Status      : Stable

=cut

sub assert_boolean {
  my ($boolean, $attribute_name) = @_;
  $attribute_name ||= '-Unknown-';
  assert_numeric($boolean, $attribute_name);
  if($boolean != 0 && $boolean != 1) {
    throw "Attribute $attribute_name was an invalid boolean. Expected: 1 or 0. Got: $boolean";
  }
  return 1;
}

=head2 assert_strand

  Arg [1]     : The Scalar to check
  Arg [2]     : The attribute name you are asserting; not required but allows
                for more useful error messages to be generated. Defaults to
                C<-Unknown->.
  Description : A subroutine which checks to see if the given scalar is
                a whole integer and if the value is set to C<1>, C<0> or C<-1>
  Returntype  : Boolean; true if we had a strand integer otherwise we signal
                failure via exceptions
  Example     : assert_strand(1, 'strand');
  Exceptions  : See L<assert_integer> and we raise exceptions if the value
                was not equal to the 3 valid states
  Status      : Stable

=cut

sub assert_strand {
  my ($strand, $attribute_name) = @_;
  $attribute_name ||= '-Unknown-';
  assert_numeric($strand, $attribute_name);
  if($strand != -1 && $strand != 0 && $strand ne 1) {
    throw "Attribute $attribute_name was an invalid strand. Expected: 1, 0 or -1. Got: $strand";
  }
  return 1;
}


=head2 assert_file_handle

  Arg [1]     : The Scalar to check
  Arg [2]     : The attribute name you are asserting; not required but allows
                for more useful error messages to be generated. Defaults to
                C<-Unknown->.
  Description : A subroutine which checks to see if the given scalar is
                actually a file handle. This will handle those which are Glob
                references and those which inherit from C<IO::Handle>. It will
                also cope with a blessed Glob reference.
  Returntype  : Boolean;
  Example     : assert_file_handle($fh, '-FILE_HANDLE');
  Exceptions  : Raised if not defined, not a reference and was not a
                GLOB or did not inherit from IO::Handle   
  Status      : Stable

=cut

sub assert_file_handle {
  my ($file_handle, $attribute_name) = @_;
  $attribute_name ||= '-Unknown-';
  throw "Attribute $attribute_name was undefined" if ! defined $file_handle;
  my $ref = ref($file_handle);
  throw "Attribute $attribute_name was not a reference. Got: $file_handle" if ! $ref;
  if(!openhandle($file_handle)) {
    if(blessed($file_handle)) {
      if(! $file_handle->isa('IO::Handle')) {
        throw "Attribute $attribute_name was blessed but did not inherit from IO::Handle. Ref was: $ref";
      }
    }
    else {
      throw "Attribute $attribute_name was not a file handle. Ref was: $ref";
    }
  }
  return 1;
}

1;
