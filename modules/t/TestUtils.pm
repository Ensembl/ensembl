use strict;

package TestUtils;

require Exporter;
use vars qw( @ISA @EXPORT_OK );
@ISA=('Exporter');
@EXPORT_OK=qw(&debug &test_getter_setter);

=head2 test_getter_setter

  Arg [1]    : Object $object
               The object to test the getter setter on
  Arg [2]    : string $method
               The name of the getter setter method to test
  Arg [3]    : $test_val
               The value to use to test the set behavior of the method.
  Example    : ok(&TestUtils::test_getter_setter($object, 'type', 'value'));
  Description: Tests a getter setter method by attempting to set a value
               and verifying that the newly set value can be retrieved.  The
               old value of the the attribute is restored after the test 
               (providing the method functions correctly).
  Returntype : boolean - true value on success, false on failure
  Exceptions : none
  Caller     : test scripts

=cut

sub test_getter_setter {
    my ($object, $method, $test_val) = @_;
    
    my $ret_val = 0;
    
    #save the old value
    my $old_val = $object->$method;
    
    $object->$method($test_val);
    
    #verify value was set
    $ret_val = (!defined($test_val) && !defined($object->$method)) || 
	       ($object->$method eq $test_val);

    #restore the old value
    $object->$method($old_val);
    
    return $ret_val;
}

sub debug {
  my $txt = shift;
  if( $::verbose ) {
    print STDERR $txt,"\n";
  }
}




1;

    
