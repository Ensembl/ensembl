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

package Bio::EnsEMBL::DBSQL::Driver;

use warnings;
use strict;

use Scalar::Util qw(weaken);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

sub new {
    my ($class, $parent) = @_;

    my $self = bless {}, $class;
    $self->parent($parent);

    return $self;
}

sub parent {
    my ($self, @args) = @_;
    if (@args) {
        ( $self->{'_parent'} ) = @args;
        weaken $self->{'_parent'};
    }
    return $self->{'_parent'};
}

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

sub last_insert_id_args {
    return;
}

sub can_straight_join {
    return;
}

sub set_wait_timeout {
    my $self = shift;
    my $class = ref $self;
    warning("'set_wait_timeout()' is not implemented in ${class}");
    return;
}

sub AUTOLOAD {
    my ($self, @args) = @_;
    my $class = ref $self;
    my $method = $Bio::EnsEMBL::DBSQL::Driver::AUTOLOAD;
    $method =~ s/^${class}:://;
    throw("'${method}()' has not been implemented in ${class}");
    return;
}

sub DESTROY { }                 # prevent DESTROY being handled by AUTOLOAD

1;
