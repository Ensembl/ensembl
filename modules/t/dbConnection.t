use lib 't';
use strict;
use warnings;


BEGIN { $| = 1;
	use Test;
	plan tests => 24;
}

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::DBSQL::DBConnection;


our $verbose = 0;

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

# construct a second database handle using a first one:
my $dbc2 = Bio::EnsEMBL::DBSQL::DBConnection->new(-DBCONN => $dbc);

ok($dbc2->dbname eq $dbc->dbname());
ok($dbc2->host   eq $dbc->host());
ok($dbc2->username()   eq $dbc->username());
ok($dbc2->password eq $dbc->password());
ok($dbc2->port  == $dbc->port());
ok($dbc2->driver eq $dbc->driver());

