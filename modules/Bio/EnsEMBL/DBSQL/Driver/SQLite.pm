package Bio::EnsEMBL::DBSQL::Driver::SQLite;

use warnings;
use strict;

use base 'Bio::EnsEMBL::DBSQL::Driver';

sub from_date_to_seconds {
    my ($self, $column) = @_;
    return "STRFTIME('%s', $column)";
}

sub from_seconds_to_date {
    my ($self, $seconds) = @_;
    return "DATETIME($seconds)";
}

1;
