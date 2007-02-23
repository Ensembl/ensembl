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

