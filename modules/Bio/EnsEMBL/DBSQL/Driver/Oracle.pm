package Bio::EnsEMBL::DBSQL::Driver::Oracle;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Exception qw(warning);

use base 'Bio::EnsEMBL::DBSQL::Driver';

sub new {
    my ($self, @args) = @_;
    warning(__PACKAGE__ . ' is untested and is being used only to collect Oracle-specific code');
    return $self->SUPER::new(@args);
}

sub connect_params {
    my ($self, $conn) = @_;

    my $dbname = $conn->dbname();

    my $dsn      = "DBI:Oracle:";
    my $username = sprintf("%s@%s", $conn->username(), $dbname );

    return {
        dsn        => $dsn,
        username   => $username,
        password   => $conn->password(),
        attributes => { 'RaiseError' => 1, 'PrintError' => 0 },
    };
}

1;
