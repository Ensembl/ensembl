# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

use Test::More;
use Test::Warnings;
use Test::Exception;
use File::Spec;
use Bio::EnsEMBL::Utils::IO qw/work_with_file/;
use Bio::EnsEMBL::Utils::URI qw/parse_uri is_uri/;

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

assert_db_parsing('URL with no user but a password','mysql://:pass@host:51',{
  -DRIVER => 'mysql',
  -PASS => 'pass',
  -HOST => 'host',
  -PORT => 51
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

assert_db_parsing('SQLite simple', 'sqlite://db', {
  -DRIVER => 'sqlite',
  -DBNAME => 'db'
});

assert_db_parsing('SQLite simple', 'sqlite:///db', {
  -DRIVER => 'sqlite',
  -DBNAME => '/db'
});

{
  my $path = File::Spec->catfile(File::Spec->tmpdir(), 'ensembl_utilsuri.sqlite.tmp');
  my $db = $path.q{/tab};
  assert_db_parsing('SQLite with file and table but no file on disk', 'sqlite://'.$db, {
    -DRIVER => 'sqlite',
    -DBNAME => $db
  });
  work_with_file($path, 'w', sub { my $fh = shift @_; print $fh '1'; });
  assert_db_parsing('SQLite with file and with table', 'sqlite://'.$db, {
    -DRIVER => 'sqlite',
    -DBNAME => $path,
    -TABLE  => 'tab'
  });
  unlink $path;
}

assert_db_parsing('URL params','mysql://host/db?param1=val1;param2',{
  -DRIVER => 'mysql',
  -HOST => 'host',
  -DBNAME => 'db',
}, {
  param1 => ['val1'],
  param2 => [undef]
});

#URI without a scheme
assert_parsing('URL no scheme',':////table',{
  scheme => '',
  db_params => { table => 'table', dbname => undef},
  params => {},
  param_keys => []
});

assert_parsing('http with params','http://host/path?param1=val1',{
  scheme => 'http',
  param_keys => [qw/param1/],
  params => {param1 => ['val1']},
  path => '/path',
  host => 'host',
  db_params => {}
});

SKIP : {
  skip 'Skipping real-life tests as we expect URI encoding/decoding', 2 if ! $Bio::EnsEMBL::Utils::URI::URI_ESCAPE;
  assert_parsing('http real life params','http://www.google.co.uk:80/search?hl=en&source=hp&q=testing+%3D&hl=en',{
    scheme => 'http',
    param_keys => [qw/hl source q/],
    params => {
      hl => [qw/en en/],
      source => ['hp'],
      q => [qw/testing+=/]
    },
    path => '/search',
    host => 'www.google.co.uk',
    port => 80,
    db_params => {}
  });
}

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

is(is_uri('mysql://user@host'), 1, 'Checking mysql scheme is detected as a URI');
is(is_uri('file:///path'), 1, 'Checking file scheme is ok detected as a URI');
is(is_uri('/path'), 0, 'Checking plain path is not detected as a URI');
is(is_uri(''), 0, 'Checking empty path is not detected as a URI');
is(is_uri(undef), 0, 'Checking undef is not detected as a URI');

done_testing();
