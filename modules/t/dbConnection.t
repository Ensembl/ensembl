use lib 't';
use strict;
use warnings;


BEGIN { $| = 1;
	use Test;
	plan tests => 15;
}

use MultiTestDB;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use TestUtils qw(test_getter_setter debug);
use Bio::EnsEMBL::DBSQL::DBConnection;


our $verbose = 0;

#
# 1 DBConnection compiles
#
ok(1);

my $multi = MultiTestDB->new;
my $db    = $multi->get_DBAdaptor('core');


#
# 2 new
#
my $dbc;
{
  my $db_name = $db->dbname;
  my $port    = $db->port;
  my $user    = $db->username;
  my $pass    = $db->password;
  my $host    = $db->host;
  my $driver  = $db->driver;

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
# 8-9 _get_adaptor
#
my $adaptor_name = 'Bio::EnsEMBL::DBSQL::SliceAdaptor';
my $adaptor = $dbc->_get_adaptor($adaptor_name);
ok($adaptor->isa($adaptor_name));
ok($adaptor == $dbc->_get_adaptor($adaptor_name)); #verify cache is used

#
# 10 dbhandle
#
ok(test_getter_setter($dbc, 'db_handle', $db->db_handle));

#
# 11 prepare
#
my $sth = $dbc->prepare('SELECT * from gene limit 1');
$sth->execute;
ok($sth->rows);
$sth->finish();


#
# 12 add_db_adaptor
#
$dbc->add_db_adaptor('core', $db);

my $db1 = $dbc->get_all_db_adaptors->{'core'};
my $db2 = $db->_obj;
debug("\n\ndb1=[$db1] db2=[$db2]\n\n"); 
ok($db1 == $db2);

#
# 13 get_db_adaptor
#
ok($dbc->get_db_adaptor('core')->isa('Bio::EnsEMBL::DBSQL::DBConnection'));

#
# 14-15 remove_db_adaptor
#
$dbc->remove_db_adaptor('core');
ok(!defined $dbc->get_db_adaptor('core'));
ok(!defined $dbc->get_all_db_adaptors->{'core'});







