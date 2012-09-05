use strict;
use warnings;

use Test::More;
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
{
  my $db_name = $db->dbc->dbname;
  my $port    = $db->dbc->port;
  my $user    = $db->dbc->username;
  my $pass    = $db->dbc->password;
  my $host    = $db->dbc->host;
  my $driver  = $db->dbc->driver;

  $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(-dbname => $db_name,
						-user   => $user,
						-pass   => $pass,
						-port   => $port,
						-host   => $host,
						-driver => $driver);
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

#
# 7 password
#
ok(test_getter_setter($dbc, 'password', 'ensembl_password'));

#
# 8-9 _get_adaptor NO LONGER ALLOWED MUST GO VIA DBAdaptor
#
#my $adaptor_name = 'Bio::EnsEMBL::DBSQL::GeneAdaptor';
#my $adaptor = $dbc->_get_adaptor($adaptor_name);
#ok($adaptor->isa($adaptor_name));
#ok($adaptor == $dbc->_get_adaptor($adaptor_name)); #verify cache is used

#
# 10 dbhandle
#
ok(test_getter_setter($dbc, 'db_handle', $dbc->db_handle));


{
  #
  # 11 prepare
  #
  my $sth = $dbc->prepare('SELECT * from gene limit 1');
  $sth->execute;
  ok($sth->rows);
  $sth->finish;
}

# AGAIN NOW DONE VIA DBADAPTOR
#
# 12 add_db_adaptor
#
#$dbc->add_db_adaptor('core', $db);
#
#my $db1 = $dbc->get_all_db_adaptors->{'core'};
#my $db2 = $db->_obj;
#debug("\n\ndb1=[$db1] db2=[$db2]\n\n"); 
#ok($db1 == $db2);
#
##
## 13 get_db_adaptor
##
#ok($dbc->get_db_adaptor('core')->isa('Bio::EnsEMBL::DBSQL::DBConnection'));
#
##
## 14-15 remove_db_adaptor
##
#$dbc->remove_db_adaptor('core');
#ok(!defined $dbc->get_db_adaptor('core'));
#ok(!defined $dbc->get_all_db_adaptors->{'core'});


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

  ok($sth1->rows);
  my @gene = $sth1->fetchrow_array();
  ok(@gene);

  $sth1->finish;
}


ok(!$dbh->ping());

$dbc->disconnect_when_inactive(0);

ok($dbc->db_handle->ping());

{
  my $sth1 = $dbc->prepare('SELECT * from gene limit 1');
  $sth1->execute();
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
ok($dbc2->dbname eq $dbc->dbname());
ok($dbc2->host   eq $dbc->host());
ok($dbc2->username()   eq $dbc->username());
ok($dbc2->password eq $dbc->password());
ok($dbc2->port  == $dbc->port());
ok($dbc2->driver eq $dbc->driver());

ok(! $dbc->equals(Bio::EnsEMBL::DBSQL::DBConnection->new(-DRIVER => 'sqlite', -DBNAME => 'bill')));

#Test spanning db_handle() methods
my $db_handle_ref;
is($dbc->disconnect_when_inactive(), 1, 'Ensuring disconnect_when_inactive() is on');
$dbc->prevent_disconnect(sub {
  is($dbc->disconnect_when_inactive(), 0, 'Ensuring disconnect_when_inactive() has been turned off');
  my $sql = 'select 1';
  $db_handle_ref = $dbc->db_handle();
  is($dbc->do($sql), 1, 'Asserting do returns 1');
  is($dbc->db_handle(), $db_handle_ref, 'Checking DBH is the same as it was at the beginning');
  is($dbc->do($sql), 1, 'Asserting do returns 1');
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
is_deeply($quote_result, [qw/`a` `b` `c`/, '`db`.`table`', '`1`'], 'Checking quote identifier will quote everything') or diag explain $quote_result;

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

done_testing();

1;
