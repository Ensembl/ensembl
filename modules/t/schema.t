use strict;
use warnings;

use FindBin qw/$Bin/;
use Test::More;
use Bio::EnsEMBL::Test::MultiTestDB;
use File::Spec::Functions qw/updir catfile catdir/;

my $db = Bio::EnsEMBL::Test::MultiTestDB->new();
my $dba = $db->get_DBAdaptor('core');
my $dbc = $dba->dbc();

# Get file location
my $sql_dir = catdir($Bin, updir(), updir(), 'sql');
my $sql_file = catfile($sql_dir, 'table.sql');

my $result = ok(-f $sql_file, 'Checking SQL schema file exists');

SKIP: {
  skip 'Skipping DB creation tests as schema file cannot be found at '.$sql_file, 1 unless $result;
  
  #Create DB & load schema
  my $new_db_name = $db->create_db_name('schematemp');
  note 'Creating database '.$new_db_name;
  $dba->dbc()->do("create database $new_db_name");
  my %args = ( host => $dbc->host(), port => $dbc->port(), user => $dbc->username(), password => $dbc->password());
  my $cmd_args = join(q{ }, map { "--${_}=$args{$_}" } keys %args);
  my $cmd = "mysql $cmd_args $new_db_name < $sql_file 2>&1";
  my $output = `$cmd`;
  my $ec = ($? >> 8);
  if($ec != 0) {
    note($output);
    fail("MySQL command failed with error code '$ec'");
  }
  else {
    pass("MySQL was able to load the Ensembl core schema");
  }
  
  note 'Dropping database '.$new_db_name;
  $dba->dbc()->do("drop database $new_db_name");
}

done_testing();