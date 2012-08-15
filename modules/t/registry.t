use strict;
use warnings;

use Config;
use Test::More;
use File::Temp qw/tempfile/;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Test::MultiTestDB;

my $threads;
if($Config{useithreads}) {
  note 'Using threaded tests';
  require threads;
  $threads = 1;
}
else {
  note 'Using non-threaded tests';
}

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db = $multi_db->get_DBAdaptor('core');
my $dbc = $db->dbc();

my $reg = 'Bio::EnsEMBL::Registry';

my $registry_template = <<'TMPL';
{
  package Reg;
  use Bio::EnsEMBL::DBSQL::DBAdaptor;
  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -HOST => '%s',
    -PORT => %d,
    -USER => '%s',
    -PASSWORD => '%s',
    -DBNAME => '%s',
    -DRIVER => 'mysql',
    -SPECIES => 'new'
  );
}
1;
TMPL

{
  my ($fh, $filename) = tempfile();
  my $final = sprintf($registry_template, $dbc->host(), $dbc->port(), $dbc->username(), $dbc->password(), $dbc->dbname());
  print $fh $final;
  close $fh;
  
  my $call = sub {
    my @results;
    foreach my $inc (0..9) {
      push(@results, $reg->load_all($filename));
    }
    return \@results;
  };
  
  if($threads) {
    my @thrds;
    foreach my $thr (0..9) {
      push(@thrds, threads->create($call));
    }
    foreach my $thr (@thrds) {
      my $results = $thr->join();
      my $msg = sprintf('THREAD %s: Checking first call loaded 1 DBAdaptor and the remaining 9 did nothing', $thr->tid());
      is_deeply($results, [1, ((0)x9)], $msg);
    }
  }
  else {
    foreach my $itr (0..9) {
      $call->();
    }
  }
}

done_testing();