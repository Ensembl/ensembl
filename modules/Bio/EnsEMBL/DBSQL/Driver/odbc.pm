# FIXME: should this be ::ODBC rather than ::odbc ??

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

sub connect_params {
    my ($self, $conn) = @_;

    my $dsn = sprintf( "DBI:ODBC:%s", $conn->dbname() );

    return {
        dsn        => $dsn,
        username   => $conn->username(),
        password   => $conn->password(),
        attributes => {
            'LongTruncOk'     => 1,
            'LongReadLen'     => 2**16 - 8,
            'RaiseError'      => 1,
            'PrintError'      => 0,
            'odbc_cursortype' => 2,
        },
    };
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
