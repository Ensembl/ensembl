=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

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
