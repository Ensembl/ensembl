use strict;
use warnings;


BEGIN { $| = 1;  
	use Test;
	plan tests => 9;
}


use Bio::EnsEMBL::Container;

#
#1 TEST - Container Compiles
#
ok(1); 

#
#2-5 TEST new and isa
#
my $test_obj = new TestObj;

ok(!$test_obj->deleteObj_called);

my $container = new Bio::EnsEMBL::Container($test_obj);
ok($container->isa('TestObj'));
ok(!$container->isa('Cruft'));
ok($container->isa('Bio::EnsEMBL::Container'));

#
# 6 TEST _obj method
#
ok($container->_obj == $test_obj);

#
# 7-8 TEST AUTOLOAD (and symbol table caching mechanism)
#
ok($container->test_method(5) == 5);
ok($container->test_method(6) == 6);

#
# 9 test destroy
#
$container = undef;
#sleep(1);
ok($test_obj->deleteObj_called);



package TestObj;

sub new {
  my $class = shift;

  my $self = bless {}, $class;

  $self->{'deleteObjCalled'} = 0;

  return $self;

}

sub test_method {
  my ($self, $val) = @_;

  return $val;
}

sub deleteObj_called {
  my $self = shift;

  return $self->{'deleteObjCalled'};
}


sub deleteObj {
  my $self = shift;

  $self->{'deleteObjCalled'} = 1;
}

1;


