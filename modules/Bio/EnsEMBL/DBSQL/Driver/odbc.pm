package Bio::EnsEMBL::DBSQL::Driver::odbc;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Exception qw(warning);

use base 'Bio::EnsEMBL::DBSQL::Driver';

sub new {
    my ($self, @args) = @_;
    warning(__PACKAGE__ . ' is untested and is being used only to collect odbc-specific code');
    return $self->SUPER::new(@args);
}

sub from_date_to_seconds {
    my ($self, $column) = @_;
    return "DATEDIFF(second,'JAN 1 1970',$column)";
}

sub from_seconds_to_date {
    my ($self, $seconds) = @_;
    return "DATEDIFF(date,'JAN 1 1970',$seconds)";
}

1;
