package Bio::EnsEMBL::DBSQL::Driver::Oracle;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Exception qw(warning);

use base 'Bio::EnsEMBL::DBSQL::Driver';

sub new {
    my ($self, @args) = @_;
    warning(__PACKAGE__ . ' is untested and is being used only to collect Postgres-specific code');
    return $self->SUPER::new(@args);
}

1;
