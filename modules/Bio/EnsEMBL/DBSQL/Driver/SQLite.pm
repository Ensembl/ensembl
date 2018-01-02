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
    return "DATETIME($seconds, 'unixepoch')";
}

sub insert_ignore_clause {
    return 'INSERT OR IGNORE';
}

1;
