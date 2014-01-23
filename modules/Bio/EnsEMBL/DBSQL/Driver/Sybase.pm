package Bio::EnsEMBL::DBSQL::Driver::Sybase;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Exception qw(warning);

use base 'Bio::EnsEMBL::DBSQL::Driver';

sub new {
    my ($self, @args) = @_;
    warning(__PACKAGE__ . ' is untested and is being used only to collect Sybase-specific code');
    return $self->SUPER::new(@args);
}

sub connect_params {
    my ($self, $conn) = @_;

    my $dbname = $conn->dbname();

    my $dbparam = ($dbname) ? ";database=${dbname}" : q{};
    my $dsn = sprintf( "DBI:Sybase:server=%s%s;tdsLevel=CS_TDS_495",
                       $conn->host(), $dbparam );

    return {
        dsn        => $dsn,
        username   => $conn->username(),
        password   => $conn->password(),
        attributes => {
            'LongTruncOk' => 1,
            'RaiseError'  => 1,
            'PrintError'  => 0,
        },
    };
}

1;
