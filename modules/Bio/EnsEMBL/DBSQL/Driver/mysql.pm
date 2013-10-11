package Bio::EnsEMBL::DBSQL::Driver::mysql;

use warnings;
use strict;

use base 'Bio::EnsEMBL::DBSQL::Driver';

sub from_date_to_seconds {
    my ($self, $column) = @_;
    return "UNIX_TIMESTAMP($column)";
}

sub from_seconds_to_date {
    my ($self, $seconds) = @_;
    return "from_unixtime($seconds)";
}

1;
