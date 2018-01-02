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

package XrefParser::DatabaseParser;

use strict;

use DBI;

use base qw( XrefParser::BaseParser );

# Base class for parsers that parse from databases rather than files

# Format of $dsn (in populate_metadata.sql)
# mysql:HOST:PORT:DATABASE:USERNAME:PASSWORD
# e.g.
# mysql:ecs4:3350:glenn_elegans_xrefs:ensro::

my $db;

# Parse the DSN.
# Return hash with keys: host, port, database, user, pass

sub parse_dsn {

  my ($self,$dsn) = @_;

  my %hash;

  my @bits = split /:/, $dsn;

  $hash{'host'} = $bits[1];
  $hash{'port'} = $bits[2];
  $hash{'database'} = $bits[3];
  $hash{'user'} = $bits[4];
  $hash{'pass'} = $bits[5];

  return %hash;

}

# Connect to a database. Return DB connection.
sub connect {

  my ($self, $dsn) = @_;

  my %dsn = $self->parse_dsn($dsn);

  $db = DBI->connect( "DBI:mysql:host=" . $dsn{'host'} . ":port=" . $dsn{'port'} . ";database=" . $dsn{'database'}, $dsn{'user'}, $dsn{'pass'},
		      {'RaiseError' => 1}) || die "Can't connect to database\n";

  return $db;

}

sub db {

  return $db;

}

1;

