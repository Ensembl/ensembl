=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package Bio::EnsEMBL::Utils::Scalar;

=pod


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=pod

=head1 NAME

Bio::EnsEMBL::Utils::Scalar

=head1 SYNOPSIS

	use Bio::EnsEMBL::Utils::Scalar qw(check_ref assert_ref wrap_array check_ref_can assert_ref_can assert_numeric assert_integer scope_guard);

	check_ref([], 'ARRAY'); # Will return true
	check_ref({}, 'ARRAY'); # Will return false
	check_ref($dba, 'Bio::EnsEMBL::DBSQL::DBAdaptor'); #Returns true if $dba is a DBAdaptor

  # Returns true if all array contents are of the given type
  check_array_contents([$dba], 'Bio::EnsEMBL::DBSQL::DBAdaptor'); 

	assert_ref([], 'ARRAY'); #Returns true
	assert_ref({}, 'ARRAY'); #throws an exception
	assert_ref($dba, 'Bio::EnsEMBL::Gene'); #throws an exception if $dba is not a Gene

  # Throws an exception if all array contents are not of the given type
  assert_array_contents([$dba], 'Bio::EnsEMBL::Gene'); #throws an exception if $dba is not a Gene

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
	
	#Scope guards
	my $v = 'wibble'; 
  {
    #Build a guard to reset $v to wibble
    my $guard = scope_guard(sub { $v = 'wibble'});
    $v = 'wobble';
    warn $v; # prints wobble
  }
  # $guard is out of scope; sub is triggered and $v is reset
  warn $v; # prints wibble

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

#
# Interface with some of the module function XS reimplementation
# 
# If Bio::EnsEMBL::XS is installed, assign the function glob to
# the XS counterpart, otherwise assign to the original function
#
BEGIN {

  if (eval { require Bio::EnsEMBL::XS; 1 }) {
    *check_ref = \&Bio::EnsEMBL::XS::Utils::Scalar::check_ref;
    *assert_ref = \&Bio::EnsEMBL::XS::Utils::Scalar::assert_ref;
    # *assert_numeric = \&Bio::EnsEMBL::XS::Utils::Scalar::assert_numeric;
    # *assert_integer = \&Bio::EnsEMBL::XS::Utils::Scalar::assert_integer;
  } else {
    *check_ref = \&check_ref_pp;
    *assert_ref = \&assert_ref_pp;
    # *assert_numeric = \&assert_numeric_pp;
    # *assert_integer = \&assert_integer_pp;

  } 

  *assert_numeric = \&assert_numeric_pp;
  *assert_integer = \&assert_integer_pp;
}


use base qw(Exporter);

our %EXPORT_TAGS;
our @EXPORT_OK;

@EXPORT_OK = qw(
  check_ref check_ref_can check_array_contents check_hash_contents
  assert_ref assert_ref_can assert_numeric assert_integer assert_boolean assert_strand assert_file_handle assert_array_contents assert_hash_contents
  wrap_array
  scope_guard
  split_array
);
%EXPORT_TAGS = (
  assert  => [qw(assert_ref assert_ref_can assert_integer assert_numeric assert_boolean assert_strand assert_file_handle assert_array_contents assert_hash_contents)],
  check   => [qw(check_ref check_ref_can check_array_contents check_hash_contents)],
  array   => [qw/wrap_array split_array/],
  all     => [@EXPORT_OK]
);

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Scalar::Util qw(blessed looks_like_number openhandle);

our $ASSERTIONS = 1;

=head2 check_ref_pp()

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

sub check_ref_pp {
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

=head2 assert_ref_pp()

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
                
                You can turn assertions off by using the global variable
                $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS = 0
  Returntype  : Boolean; true if we managed to get to the return
  Example     : assert_ref([], 'ARRAY');
  Exceptions  : If the expected type was not set and if the given reference
                was not assignable to the expected value
  Status      : Stable

=cut

sub assert_ref_pp {
  my ($ref, $expected, $attribute_name) = @_;
  return 1 unless $ASSERTIONS;
  $attribute_name ||= '-Unknown-';
  throw('No expected type given') if ! defined $expected;
  my $class = ref($ref);
  throw("The given reference for attribute $attribute_name was undef. Expected '$expected'") unless defined $ref;
  throw("Asking for the type of the attribute $attribute_name produced no type; check it is a reference. Expected '$expected'") unless $class;
  if(blessed($ref)) {
    throw("${attribute_name}'s type '${class}' is not an ISA of '${expected}'") if ! $ref->isa($expected);
  }
  else {
    throw("$attribute_name was expected to be '${expected}' but was '${class}'") if $expected ne $class;
  }
  return 1;
}

=head2 assert_array_contents

 Arg [1]     : ArrayRef references to check
 Arg [2]     : The type we expect
 Arg [3]     : The attribute name you are asserting; not required but allows
               for more useful error messages to be generated. Defaults to
               C<-Unknown->.
 Description : A subroutine which checks to see if the given objects/refs are
               what you expect. This behaves in an identical manner as
               C<assert_ref> does works on an array ref of references

               You can turn assertions off by using the global variable
               $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS = 0
 Returntype  : Boolean; true if we managed to get to the return
 Example     : assert_array_contents([[],[],[]], 'ARRAY');
 Exceptions  : Throws is references argument is not an ArrayRef, also
               if the expected type was not set and if the given reference
               was not assignable to the expected value.
 Status      : Stable

=cut

sub assert_array_contents {
  my ($array, $expected, $attribute_name) = @_;
  return 1 unless $ASSERTIONS;
  throw('No expected type given') if ! defined $expected;
  $attribute_name ||= '-Unknown-';
  assert_ref($array, 'ARRAY', $attribute_name);
  my $count = scalar(@{$array});
  for(my $i = 0; $i<$count; $i++) {
    my $ref = $array->[$i];
    my $class = ref($ref);
    throw("The given reference for attribute $attribute_name was undef (at position ${i}). Expected '$expected'") unless defined $ref;
    throw("Asking for the type of the attribute $attribute_name produced no type; check it is a reference (at position ${i}). Expected '$expected'") unless $class;
    if(blessed($ref)) {
      throw("${attribute_name}'s type '${class}' is not an ISA of '${expected}' (at position ${i})") if ! $ref->isa($expected);
    }
    else {
      throw("$attribute_name was expected to be '${expected}' but was '${class}' (at position ${i})") if $expected ne $class;
    }
  }
  return 1;
}

=head2 check_array_contents

 Arg [1]     : ArrayRef references to check
 Arg [2]     : The type we expect
 Arg [3]     : The attribute name you are asserting; not required but allows
               for more useful error messages to be generated. Defaults to
               C<-Unknown->.
 Description : A subroutine which checks to see if the given objects/refs are
               what you expect. 
 Returntype  : Boolean; true if all contents were as expected
 Example     : check_array_contents([[],[],[]], 'ARRAY');
 Exceptions  : Thrown if no type was given
 Status      : Stable

=cut

sub check_array_contents {
  my ($array, $expected, $attribute_name) = @_;
  return 0 if ! check_ref($array, 'ARRAY');
  throw('No expected type given') if ! defined $expected;
  my $contents_ok = 1;
  my $count = scalar(@{$array});
  for(my $i = 0; $i<$count; $i++) {
    my $ref = $array->[$i];
    if(!$ref) {
      $contents_ok = 0;
      last;
    }
    my $class = ref($ref);
    if(!$class) {
      $contents_ok = 0;
      last;
    }
    if(blessed($ref)) {
      if(! $ref->isa($expected)) {
        $contents_ok = 0;
        last;
      }
    }
    elsif($expected ne $class) {
      $contents_ok = 0;
      last;
    }
  }
  return $contents_ok;
}

=head2 assert_hash_contents

 Arg [1]     : HashRef references to check
 Arg [2]     : The type we expect
 Arg [3]     : The attribute name you are asserting; not required but allows
               for more useful error messages to be generated. Defaults to
               C<-Unknown->.
 Description : A subroutine which checks to see if the given objects/refs are
               what you expect. This behaves in an identical manner as
               C<assert_ref> does works on a HashRef of references. Hash keys 
               are always Strings so do not need asserting.

               You can turn assertions off by using the global variable
               $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS = 0
 Returntype  : Boolean; true if we managed to get to the return
 Example     : assert_hash_contents({a => [], b => []}, 'ARRAY');
 Exceptions  : Throws is references argument is not an ArrayRef, also
               if the expected type was not set and if the given reference
               was not assignable to the expected value.
 Status      : Stable

=cut

sub assert_hash_contents {
  my ($hash, $expected, $attribute_name) = @_;
  return 1 unless $ASSERTIONS;
  throw('No expected type given') if ! defined $expected;
  $attribute_name ||= '-Unknown-';
  assert_ref($hash, 'HASH', $attribute_name);
  my @keys = keys %{$hash};
  while(my $key = shift @keys) {
    my $ref = $hash->{$key};
    my $class = ref($ref);
    throw("The given reference for attribute $attribute_name was undef (with key ${key}). Expected '$expected'") unless defined $ref;
    throw("Asking for the type of the attribute $attribute_name produced no type; check it is a reference (with key ${key}). Expected '$expected'") unless $class;
    if(blessed($ref)) {
      throw("${attribute_name}'s type '${class}' is not an ISA of '${expected}' (with key ${key})") if ! $ref->isa($expected);
    }
    else {
      throw("$attribute_name was expected to be '${expected}' but was '${class}' (with key ${key})") if $expected ne $class;
    }
  }
  return 1;
}

=head2 check_hash_contents

 Arg [1]     : HashRef references to check
 Arg [2]     : The type we expect
 Arg [3]     : The attribute name you are asserting; not required but allows
               for more useful error messages to be generated. Defaults to
               C<-Unknown->.
 Description : A subroutine which checks to see if the given objects/refs are
               what you expect. 
 Returntype  : Boolean; true if all contents were as expected
 Example     : check_hash_contents({a => [], b => []}, 'ARRAY');
 Exceptions  : Thrown if no type was given
 Status      : Stable

=cut

sub check_hash_contents {
  my ($hash, $expected, $attribute_name) = @_;
  throw('No expected type given') if ! defined $expected;
  return 0 if ! check_ref($hash, 'HASH');
  my $contents_ok = 1;
  my @keys = keys %{$hash};
  while(my $key = shift @keys) {
    my $ref = $hash->{$key};
    if(!$ref) {
      $contents_ok = 0;
      last;
    }
    my $class = ref($ref);
    if(!$class) {
      $contents_ok = 0;
      last;
    }
    if(blessed($ref)) {
      if(! $ref->isa($expected)) {
        $contents_ok = 0;
        last;
      }
    }
    elsif($expected ne $class) {
      $contents_ok = 0;
      last;
    }
  }
  return $contents_ok;
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
                
                You can turn assertions off by using the global variable
                $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS = 0
  Returntype  : Boolean; true if we managed to get to the return
  Example     : assert_ref_can($gene, 'dbID');
  Exceptions  : If the reference is not defined, if the object does not
                implement the given method and if no method was given to check
  Status      : Stable

=cut

sub assert_ref_can {
  my ($ref, $method, $attribute_name) = @_;
  return 1 unless $ASSERTIONS;
  $attribute_name ||= '-Unknown-';
  throw('No method given') if ! defined $method;
  throw "The given reference $attribute_name is not defined. Expected method '$method'" unless defined $ref;
  throw "The given reference $attribute_name is not blessed. Expected method '$method'" unless blessed($ref);
  if(! $ref->can($method)) {
    my $str_ref = ref($ref);
    throw sprintf(q{The given blessed reference '%s' for attribute '%s' does not implement the method '%s'}, $str_ref, $attribute_name, $method);
  }
  return 1;
}

=head2 assert_numeric_pp

  Arg [1]     : The Scalar to check
  Arg [2]     : The attribute name you are asserting; not required but allows
                for more useful error messages to be generated. Defaults to
                C<-Unknown->.
  Description : A subroutine which checks to see if the given scalar is
                number or not. If not then we raise exceptions detailing why
                
                You can turn assertions off by using the global variable
                $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS = 0
  Returntype  : Boolean; true if we had a numeric otherwise we signal failure
                via exceptions
  Example     : assert_numeric(1, 'dbID');
  Exceptions  : If the Scalar is not defined, if the Scalar was blessed and
                if the value was not a number
  Status      : Stable

=cut

sub assert_numeric_pp {
  my ($integer, $attribute_name) = @_;
  return 1 unless $ASSERTIONS;
  $attribute_name ||= '-Unknown-';
  throw "$attribute_name attribute is undefined. Expected a number" if ! defined $integer;
  throw "The given attribute $attribute_name is blessed; cannot work with blessed values. Expected a number" if blessed($integer);
  if(! looks_like_number($integer)) {
    throw "Attribute $attribute_name was not a number";
  }
  return 1;
}

=head2 assert_integer_pp

  Arg [1]     : The Scalar to check
  Arg [2]     : The attribute name you are asserting; not required but allows
                for more useful error messages to be generated. Defaults to
                C<-Unknown->.
  Description : A subroutine which checks to see if the given scalar is
                a whole integer; we delegate to L<assert_numeric> for number
                checking.
                
                You can turn assertions off by using the global variable
                $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS = 0
  Returntype  : Boolean; true if we had a numeric otherwise we signal failure
                via exceptions
  Example     : assert_integer(1, 'dbID');
  Exceptions  : See L<assert_numeric> and we raise exceptions if the value
                was not a whole integer
  Status      : Stable

=cut

sub assert_integer_pp {
  my ($integer, $attribute_name) = @_;
  return 1 unless $ASSERTIONS;
  $attribute_name ||= '-Unknown-';
  throw "$attribute_name attribute is undefined. Expected an Integer" if ! defined $integer;
  throw "The given attribute $attribute_name is blessed; cannot work with blessed values. Expected an Integer" if blessed($integer);
  if(! looks_like_number($integer)) {
    throw "Attribute $attribute_name was not a number. Expected an Integer";
  }
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
                
                You can turn assertions off by using the global variable
                $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS = 0
  Returntype  : Boolean; true if we were given a boolean otherwise we signal
                failure via exceptions
  Example     : assert_boolean(1, 'is_circular');
  Exceptions  : See L<assert_integer> and we raise exceptions if the value
                was not equal to the 2 valid states
  Status      : Stable

=cut

sub assert_boolean {
  my ($boolean, $attribute_name) = @_;
  return 1 unless $ASSERTIONS;
  $attribute_name ||= '-Unknown-';
  throw "$attribute_name attribute is undefined. Expected a boolean" if ! defined $boolean;
  throw "The given attribute $attribute_name is blessed; cannot work with blessed values. Expected an Integer" if blessed($boolean);
  if(! looks_like_number($boolean)) {
    throw "Attribute $attribute_name was not a number. Expected a boolean";
  }
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
                
                You can turn assertions off by using the global variable
                $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS = 0
  Returntype  : Boolean; true if we had a strand integer otherwise we signal
                failure via exceptions
  Example     : assert_strand(1, 'strand');
  Exceptions  : See L<assert_integer> and we raise exceptions if the value
                was not equal to the 3 valid states
  Status      : Stable

=cut

sub assert_strand {
  my ($strand, $attribute_name) = @_;
  return 1 unless $ASSERTIONS;
  $attribute_name ||= '-Unknown-';
  throw "$attribute_name attribute is undefined. Expected: 1, 0 or -1" if ! defined $strand;
  throw "The given attribute $attribute_name is blessed; cannot work with blessed values. Expected: 1, 0 or -1" if blessed($strand);
  if(! looks_like_number($strand)) {
    throw "Attribute $attribute_name was not a number. Expected: 1, 0 or -1";
  }
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
                
                You can turn assertions off by using the global variable
                $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS = 0
  Returntype  : Boolean;
  Example     : assert_file_handle($fh, '-FILE_HANDLE');
  Exceptions  : Raised if not defined, not a reference and was not a
                GLOB or did not inherit from IO::Handle   
  Status      : Stable

=cut

sub assert_file_handle {
  my ($file_handle, $attribute_name) = @_;
  return 1 unless $ASSERTIONS;
  $attribute_name ||= '-Unknown-';
  throw "Attribute $attribute_name was undefined. Expected a FileHandle" if ! defined $file_handle;
  my $ref = ref($file_handle);
  throw "Attribute $attribute_name was not a reference. Got: $file_handle. Expected a FileHandle" if ! $ref;
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

=head2 split_array

  Arg [1]     : Integer Maximum size of an array produced
  Arg [2]     : ArrayRef The array to split
  Description : Takes an array of values and splits the array into multiple 
                arrays where the maximum size of each array is as specified
  Example     : my $split_arrays = split_array($large_array, 10);
  Returntype  : ArrayRef of ArrayRefs where each element is a split list
=cut

sub split_array {
  my ($amount, $array) = @_;
  assert_ref($array, 'ARRAY', 'array');
  my @split;
  my $counter = 0;
  my $index = 0;
  foreach my $e (@$array) {
    if($counter == $amount) {
      $index++;
      $counter = 0;
    }
    push(@{$split[$index]}, $e);
    $counter++;
  }
  return \@split;
}

=head2 scope_guard

  Arg [1]     : CodeRef The block of code to exit once it escapes out of scope
  Description : Simple subroutine which blesses your given code reference into
                a L<Bio::EnsEMBL::Utils::Scalar::ScopeGuard> object. This has
                a DESTROY implemented which will cause the code reference
                to execute once the object goes out of scope and its reference
                count hits 0.
  Returntype  : Bio::EnsEMBL::Utils::Scalar::ScopeGuard
  Example     : my $v = 'wibble'; 
                {
                  #Build a guard to reset $v to wibble
                  my $guard = scope_guard(sub { $v = 'wibble'});
                  $v = 'wobble';
                  warn $v;
                }
                # $guard is out of scope; sub is triggered and $v is reset
                warn $v;
  Exceptions  : Raised if argument was not a CodeRef   
  Status      : Stable

=cut

sub scope_guard {
  my ($callback) = @_;
  assert_ref($callback, 'CODE', 'callback');
  return bless($callback, 'Bio::EnsEMBL::Utils::Scalar::ScopeGuard');
}

1;

#### SUPER SECRET PACKAGE. IGNORE ME
package Bio::EnsEMBL::Utils::Scalar::ScopeGuard;
sub DESTROY {
  my ($self) = @_;
  $self->();
  return;
}

1;
