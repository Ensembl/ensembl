package Bio::EnsEMBL::DBSQL::Driver::SQLite;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Exception qw/throw/;

use base 'Bio::EnsEMBL::DBSQL::Driver';

sub connect_params {
    my ($self, $conn) = @_;

    my $dbname = $conn->dbname();
    throw "We require a dbname to connect to a SQLite database"
      if !$dbname;

    return {
        dsn        => sprintf( "DBI:SQLite:%s", $dbname ),
        username   => '',
        password   => '',
        attributes => { 'RaiseError' => 1 },
    };
}

sub from_date_to_seconds {
    my ($self, $column) = @_;
    return "STRFTIME('%s', $column)";
}

sub from_seconds_to_date {
    my ($self, $seconds) = @_;
    return "DATETIME($seconds)";
}

sub insert_ignore_clause {
    return 'INSERT OR IGNORE';
}

1;
