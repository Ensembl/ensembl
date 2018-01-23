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

package Bio::EnsEMBL::DBSQL::Driver::mysql;

use warnings;
use strict;

use base 'Bio::EnsEMBL::DBSQL::Driver';

#
# override parent's method to enable MySQL local load data in case DBD::mysql 
# has been compiled against a C client library which has been built with
# no support for this feature
#
sub connect_params {
  my ($self, $conn) = @_;

  my $params = $self->SUPER::connect_params($conn);
  $params->{attributes}{mysql_local_infile} = 1;

  return $params;
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
