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
use Test::MockObject::Extends;
use Time::HiRes qw/usleep/;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::ProxyDBConnection;

#
# 1 DBConnection compiles
#
ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;
my $db    = $multi->get_DBAdaptor('core');


#
# 2 new
#
my $dbc;
my %dbc_args;
{
  my $db_name = $db->dbc->dbname;
  my $port    = $db->dbc->port;
  my $user    = $db->dbc->username;
  my $pass    = $db->dbc->password;
  my $host    = $db->dbc->host;
  my $driver  = $db->dbc->driver;
  %dbc_args = (
    -DBNAME => $db_name,
    -USER   => $user,
    -PASS   => $pass,
    -PORT   => $port,
    -HOST   => $host,
    -DRIVER => $driver
  );
  $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(%dbc_args);
}

ok($dbc->isa('Bio::EnsEMBL::DBSQL::DBConnection'));

#
# 3 driver
#
ok(test_getter_setter($dbc, 'driver'  , 'oracle'));

#
# 4 port
#
ok(test_getter_setter($dbc, 'port'    , 6666));

#
# 5 dbname
#
ok(test_getter_setter($dbc, 'dbname'  , 'ensembl_db_name'));

#
# 6 username
#
ok(test_getter_setter($dbc, 'username', 'ensembl_user'));
is($dbc->user(), $dbc->username(), 'Checking user() returns same as username()');

#
# 7 password
#
ok(test_getter_setter($dbc, 'password', 'ensembl_password'));
is($dbc->pass(), $dbc->password(), 'Checking pass() returns same as password()');

#
# 10 dbhandle
#
ok(test_getter_setter($dbc, 'db_handle', $dbc->db_handle));

#
# dbname
#
ok(test_getter_setter($dbc, 'host', $dbc->host));
is($dbc->hostname(), $dbc->host(), 'Checking hostname() returns same as host()');

# Check the to_hash() works
is_deeply($dbc->to_hash(), \%dbc_args, 'Checking to_hash() can roundtrip a DBConnection');

{
  #
  # 11 prepare
  #
  my $sth = $dbc->prepare('SELECT * from gene limit 1');
  $sth->execute;
  my @row = $sth->fetchrow_array;
  ok($sth->rows);
  $sth->finish;
}

#
# try the database with the disconnect_when_inactive flag set.
# this should automatically disconnect from the db
#
my $dbh = $dbc->db_handle();

$dbc->disconnect_when_inactive(1);

ok(!$dbh->ping());

{
  # reconnect should happen now
  my $sth2 = $dbc->prepare('SELECT * from gene limit 1');
  $sth2->execute;
  $sth2->fetchrow_array;
  ok($sth2->rows);
  $sth2->finish;

  # disconnect should occur now
}

ok(!$dbh->ping());

$dbh = $dbc->db_handle();

#
# try the same thing but with 2 connections at a time
#

{
  my $sth1 = $dbc->prepare('SELECT * from gene limit 1');
  $sth1->execute;

  {
    my $sth2 = $dbc->prepare('SELECT * from gene limit 1');
    $sth2->execute();
    $sth2->finish();
  }

  my @gene = $sth1->fetchrow_array();
  ok($sth1->rows);
  ok(@gene);

  $sth1->finish;
}


ok(!$dbh->ping());

$dbc->disconnect_when_inactive(0);

ok($dbc->db_handle->ping());

{
  my $sth1 = $dbc->prepare('SELECT * from gene limit 1');
  $sth1->execute();
  $sth1->fetchrow_array();
  ok($sth1->rows());
  $sth1->finish();
}

#should not have disconnected this time
ok($dbc->db_handle->ping());
$dbc->disconnect_when_inactive(1);

# construct a second database handle using a first one:
my $dbc2 = Bio::EnsEMBL::DBSQL::DBConnection->new(-DBCONN => $dbc);

#equality checks
ok($dbc->equals($dbc2));
ok(! $dbc->equals());
is($dbc2->dbname(),   $dbc->dbname());
is($dbc2->host(),     $dbc->host());
is($dbc2->username(), $dbc->username());
is($dbc2->password(), $dbc->password());
is($dbc2->port(),     $dbc->port());
is($dbc2->driver(),   $dbc->driver());

ok(! $dbc->equals(Bio::EnsEMBL::DBSQL::DBConnection->new(-DRIVER => 'TestDummy', -DBNAME => 'bill')));

#Test spanning db_handle() methods
my $db_handle_ref;
is($dbc->disconnect_when_inactive(), 1, 'Ensuring disconnect_when_inactive() is on');
$dbc->prevent_disconnect(sub {
  is($dbc->disconnect_when_inactive(), 0, 'Ensuring disconnect_when_inactive() has been turned off');
  my $sql = 'select 1';
  $db_handle_ref = $dbc->db_handle();
  my $do_result_1 = $dbc->do($sql);
  ok($do_result_1, "Asserting do returns true [$do_result_1]");
  is($dbc->db_handle(), $db_handle_ref, 'Checking DBH is the same as it was at the beginning');
  my $do_result_2 = $dbc->do($sql);
  ok($do_result_2, "Asserting do returns true [$do_result_2]");
  is($dbc->db_handle(), $db_handle_ref, 'Checking DBH is the same as it was at the beginning');
  my $do_result_3 = $dbc->do($sql, undef, 1);
  ok($do_result_3, "Asserting do accepts bind_values arguments [$do_result_3]");
  is($dbc->db_handle(), $db_handle_ref, 'Checking DBH is the same as it was at the beginning');
  my $sth1 = $dbc->prepare($sql);
  $sth1->execute();
  my $expected = $sth1->fetchall_arrayref()->[0]->[0];
  $sth1->finish();
  is($expected, 1, 'Asserting prepare returns 1');
  is($dbc->db_handle(), $db_handle_ref, 'Checking DBH is the same as it was at the beginning');
});
is($dbc->disconnect_when_inactive(), 1, 'Ensuring disconnect_when_inactive() is on');
isnt($dbc->db_handle(), $db_handle_ref, 'Checking if a new DBH has been handed out');
$dbc->disconnect_when_inactive();

#Testing quote identifiers
my $quote_result = $dbc->quote_identifier(qw/a b c/, [undef, qw/db table/], [1]);
my $expected_quote_result;
if ( $dbc->driver() eq 'SQLite' ) {
  $expected_quote_result = [qw/"a" "b" "c"/, '"db"."table"', '"1"'];
} else {
  $expected_quote_result = [qw/`a` `b` `c`/, '`db`.`table`', '`1`'];
}
is_deeply($quote_result, $expected_quote_result, 'Checking quote identifier will quote everything') or diag explain $quote_result;

#Testing reconnection via proxy
$db->dbc()->disconnect_if_idle();
#my $dbc_copy = bless({%{$db->dbc()}}, ref($db->dbc()));
#$dbc_copy = Test::MockObject::Extends->new($dbc_copy);
my $dbc_copy = mock_object($dbc);
{
  my $pdbc = Bio::EnsEMBL::DBSQL::ProxyDBConnection->new(-DBC => $dbc_copy, -DBNAME => $dbc->dbname(), -RECONNECT_INTERVAL => 0);
  is($pdbc->sql_helper()->execute_single_result(-SQL => 'select 1'), 1, 'Checking we get a 1 back from the DB');
  usleep(1000); #sleep 1ms
  is($pdbc->sql_helper()->execute_single_result(-SQL => 'select 1'), 1, 'Checking we get a 1 back from the DB');
  $dbc_copy->__is_called('reconnect', 0, "No need to reconnect as we have no interval set");
  
  #Disconnect, set the interval to 1ms and sleep for 10ms. Skips reconnection as the connection is not active
  $dbc_copy->disconnect_if_idle();
  $pdbc->reconnect_interval(1);
  usleep((10*1000));
  $pdbc->check_reconnection();
  $dbc_copy->__is_called('connected', 1, "connected() was called since we had to check if the DBConnection was active");
  $dbc_copy->__is_called('reconnect', 0, "No need to reconnect as we have no interval set");
  
  #Connect, set the interval to 1ms and sleep for 10ms. Reconnect
  $dbc_copy->connect();
  $dbc_copy->__clear();
  usleep((10*1000));
  $pdbc->check_reconnection();
  $dbc_copy->__is_called('connected', 1, "connected() was called since we had to check if the DBConnection was active");
  $dbc_copy->__is_called('reconnect', 1, "reconnect() as we had gone beyond our normal timeout interval");
}

# Test error reporting when we connect using a bogus driver
{
  throws_ok {
    my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(
      -HOST => '127.0.0.1', -PORT => 3306, -USER => 'usr', -DRIVER => 'bogusdriver'
    );
    local $SIG{__WARN__} = sub {
      #swallow the warning. we get it raised via an error anyway
    };
    $dbc->db_handle;
  } qr/Cannot load.+::bogusdriver/, 'Checking for an error message from DBI detailing missing driver problems';
  
  if($db->dbc->driver() eq 'mysql') {
    throws_ok {
      my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(
        -HOST => $db->dbc->host, -PORT => $db->dbc->port, -USER => 'arandomuser', -PASS => 'sobogusthisis', -DRIVER => 'mysql'
      );
      local $SIG{__WARN__} = sub {
        #swallow the warning. we get it raised via an error anyway
      };
      $dbc->db_handle;
    } qr/Access denied for user/, 'Checking we raise an error about a bogus connection details';
  }
}

done_testing();

1;
