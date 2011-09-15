use strict;
use warnings;

use Test::More;
use Test::Exception;

use Bio::EnsEMBL::Utils::URI qw/parse_uri/;

sub assert_parsing {
  my ($msg, $url, $expected) = @_;
  my $actual = parse_uri($url);
  ok(! $actual->is_db_scheme(), 'Checking we do not have a DB scheme');
  is_deeply($actual, $expected, $msg) or diag explain $actual;
  return;
}

sub assert_db_parsing {
  my ($msg, $url, $expected, $expected_params) = @_;
  my $actual = parse_uri($url);
  ok($actual->is_db_scheme(), 'Checking we have a DB scheme');
  my $actual_dbsql_params = {$actual->generate_dbsql_params()};
  if(! is_deeply($actual_dbsql_params, $expected, $msg)) {
    diag explain $actual;
    diag explain $actual_dbsql_params;
  }
  if($expected_params) {
    is_deeply($actual->{params}, $expected_params, "Checking paramaters for $msg") or diag explain $actual->{params};
  }

  my $roundtrip_uri = $actual->generate_uri();
  is($roundtrip_uri, $url, "Checking that '$msg' will roundtrip its URI");

  return;
}

assert_db_parsing('Basic full URL','mysql://user:pass@host:51/dbname',{
  -DRIVER => 'mysql',
  -USER => 'user',
  -PASS => 'pass',
  -HOST => 'host',
  -PORT => 51,
  -DBNAME => 'dbname'
});

assert_db_parsing('Basic full URL with table','mysql://user:pass@host:51/dbname/table',{
  -DRIVER => 'mysql',
  -USER => 'user',
  -PASS => 'pass',
  -HOST => 'host',
  -PORT => 51,
  -DBNAME => 'dbname',
  -TABLE => 'table'
});

assert_db_parsing('URL with no pass','mysql://user@host:51/dbname',{
  -DRIVER => 'mysql',
  -USER => 'user',
  -HOST => 'host',
  -PORT => 51,
  -DBNAME => 'dbname'
});

assert_db_parsing('URL with no user but a password','mysql://:pass@host:51/dbname',{
  -DRIVER => 'mysql',
  -PASS => 'pass',
  -HOST => 'host',
  -PORT => 51,
  -DBNAME => 'dbname'
});


assert_db_parsing('URL no user credentials with table','mysql://host:51/dbname/table',{
  -DRIVER => 'mysql',
  -HOST => 'host',
  -PORT => 51,
  -DBNAME => 'dbname',
  -TABLE => 'table'
});

assert_db_parsing('URL no port','mysql://host/dbname',{
  -DRIVER => 'mysql',
  -HOST => 'host',
  -DBNAME => 'dbname'
});

assert_db_parsing('URL no host','mysql:///dbname',{
  -DRIVER => 'mysql',
  -DBNAME => 'dbname'
});

assert_db_parsing('URL table only','mysql:////table',{
  -DRIVER => 'mysql',
  -TABLE => 'table',
});

assert_db_parsing('URL table with paramaters','mysql:////table?insert_scheme=INSERT_IGNORE',{
  -DRIVER => 'mysql',
  -TABLE => 'table',
}, {
  insert_scheme => ['INSERT_IGNORE']
});

assert_db_parsing('URL params','mysql://host/db?param1=val1;param2',{
  -DRIVER => 'mysql',
  -HOST => 'host',
  -DBNAME => 'db',
}, {
  param1 => ['val1'],
  param2 => [undef]
});

assert_parsing('http with params','http://host/path?param1=val1',{
  scheme => 'http',
  param_keys => [qw/param1/],
  params => {param1 => ['val1']},
  path => '/path',
  host => 'host',
  db_params => {}
});

assert_parsing('http real life params','http://www.google.co.uk:80/search?hl=en&source=hp&q=testing+%3D&hl=en',{
  scheme => 'http',
  param_keys => [qw/hl source q/],
  params => {
    hl => [qw/en en/],
    source => ['hp'],
    q => [qw/testing+%3D/]
  },
  path => '/search',
  host => 'www.google.co.uk',
  port => 80,
  db_params => {}
});

assert_parsing('File path','file:///my/path',{
  scheme => 'file',
  params => {},
  param_keys => [],
  path => '/my/path',
  db_params => {},
});

{
  my $uri = Bio::EnsEMBL::Utils::URI->new('http');
  $uri->host('www.google.co.uk');
  $uri->port(80);
  $uri->path('/search');
  $uri->add_param('q', 'testing');
  is($uri->generate_uri(), 'http://www.google.co.uk:80/search?q=testing', 'Checking URI generation for example usage');
}

{
  my $u = Bio::EnsEMBL::Utils::URI->new('http');
  dies_ok { $u->port(-1) } 'Checking we die if port is given an integer less than 1';
  dies_ok { $u->port(80.265) } 'Checking we die if port is given a floating point number';
  dies_ok { $u->port('my port') } 'Checking we die if port is given a string';
  lives_ok { $u->port('80') } 'Checking we live if port is given a string which decodes to an integer';
  lives_ok { $u->port(80) } 'Checking we live if port is given an integer';
}

done_testing();
