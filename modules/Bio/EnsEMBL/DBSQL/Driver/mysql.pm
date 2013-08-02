package Bio::EnsEMBL::DBSQL::Driver::mysql;

use warnings;
use strict;

use base 'Bio::EnsEMBL::DBSQL::Driver';

sub connect_params {
    my ($self, $conn) = @_;

    my $dbname = $conn->dbname();
    my $dbparam = ($dbname) ? "database=${dbname};" : q{};

    my $dsn = sprintf( "DBI:%s:%shost=%s;port=%s",
                       $conn->driver(), $dbparam,
                       $conn->host(),   $conn->port() );

    if ( $conn->{'disconnect_when_inactive'} ) {
      $conn->{'count'}++;
      if ( $conn->{'count'} > 1000 ) {
        sleep 1;
        $conn->{'count'} = 0;
      }
    }

    return {
        dsn        => $dsn,
        username   => $conn->username(),
        password   => $conn->password(),
        attributes => { 'RaiseError' => 1 },
    };
}

sub from_date_to_seconds {
    my ($self, $column) = @_;
    return "UNIX_TIMESTAMP($column)";
}

sub from_seconds_to_date {
    my ($self, $seconds) = @_;
    return "from_unixtime($seconds)";
}

sub last_insert_id_args {
    return (undef, undef, undef, undef);
}

sub insert_ignore_clause {
    return 'INSERT IGNORE';
}

sub can_straight_join {
    return 1;
}

sub set_wait_timeout {
    my ($self, $dbh, $timeout) = @_;
    $dbh->do( "SET SESSION wait_timeout=${timeout}" );
    return;
}

1;
