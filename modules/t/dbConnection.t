use lib 't';
use strict;
use warnings;


BEGIN { $| = 1;
	use Test;
	plan tests => 14;
}

use MultiTestDB;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use TestUtils qw(test_getter_setter debug);
use Bio::EnsEMBL::DBSQL::DBConnection;


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

  print STDERR "Using port [$port]\n";

  $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(-dbname => $db_name,
						-user   => $user,
						-pass   => $pass,
						-port   => $port,
						-host   => $host,
						-driver => $driver);
}

ok($dbc->isa('Bio::EnsEMBL::DBSQL::DBConnection'));

#
# 2 driver
#
ok(test_getter_setter($dbc, 'driver'  , 'oracle'));

#
# 3 port
#
ok(test_getter_setter($dbc, 'port'    , 6666));

#
# 4 dbname
#
ok(test_getter_setter($dbc, 'dbname'  , 'ensembl_db_name'));

#
# 5 username
#
ok(test_getter_setter($dbc, 'username', 'ensembl_user'));

#
# 6 password
#
ok(test_getter_setter($dbc, 'password', 'ensembl_password'));

#
# 7-8 _get_adaptor
#
my $adaptor_name = 'Bio::EnsEMBL::DBSQL::SliceAdaptor';
my $adaptor = $dbc->_get_adaptor($adaptor_name);
ok($adaptor->isa($adaptor_name));
ok($adaptor == $dbc->_get_adaptor($adaptor_name)); #verify cache is used

#
# 9 dbhandle
#
ok(test_getter_setter($dbc, 'db_handle', $db->db_handle));

#
# 10 prepare
#
my $sth = $dbc->prepare('SELECT * from gene limit 1');
$sth->execute;
ok(scalar($sth->fetchrow_array) == 1);

#
#  add_db_adaptor
#
$dbc->add_db_adaptor('core', $db);


ok($dbc->get_all_db_adaptors->{'core'} == $db->_obj);

#
# 12 get_db_adaptor
#
ok($dbc->get_db_adaptor('core') == $db->_obj);

#
# 13-14 remove_db_adaptor
#
$dbc->remove_db_adaptor('core');
ok(!defined $dbc->get_db_adaptor('core'));
ok(!defined $dbc->get_all_db_adaptors->{'core'});







