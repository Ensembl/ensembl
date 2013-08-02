package Bio::EnsEMBL::DBSQL::Driver;

use warnings;
use strict;

use Scalar::Util qw(weaken);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

sub new {
    my ($class, $parent) = @_;

    my $self = bless {}, $class;
    $self->parent($parent);

    return $self;
}

sub parent {
    my ($self, @args) = @_;
    if (@args) {
        ( $self->{'_parent'} ) = @args;
        weaken $self->{'_parent'};
    }
    return $self->{'_parent'};
}

sub last_insert_id_args {
    return;
}

sub can_straight_join {
    return;
}

sub set_wait_timeout {
    my $self = shift;
    my $class = ref $self;
    warning("'set_wait_timeout()' is not implemented in ${class}");
    return;
}

sub AUTOLOAD {
    my ($self, @args) = @_;
    my $class = ref $self;
    my $method = $Bio::EnsEMBL::DBSQL::Driver::AUTOLOAD;
    $method =~ s/^${class}:://;
    throw("'${method}()' has not been implemented in ${class}");
    return;
}

sub DESTROY { }                 # prevent DESTROY being handled by AUTOLOAD

1;
