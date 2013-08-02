use strict;
use warnings;

use Config;
use Test::More;
use File::Temp qw/tempfile/;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Test::MultiTestDB;

my $threads;
if($Config{useithreads} && ! $ENV{ENS_FORCE_NOTHREADS}) {
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
    -DRIVER => '%s',
    -SPECIES => 'new'
  );
}
1;
TMPL

{
  my ($fh, $filename) = tempfile();
  my $final = sprintf($registry_template,
                      $dbc->host() || '',
                      $dbc->port() || 0,
                      $dbc->username() || '',
                      $dbc->password() || '',
                      $dbc->dbname(),
                      $dbc->driver(),
      );
  print $fh $final;
  close $fh;
  
  my $call = sub {
    my @results;
    foreach my $inc (0..9) {
      push(@results, $reg->load_all($filename));
    }
    return \@results;
  };
  
  my @results;
  my @expected;
  if($threads) {
    $ENV{RUNTESTS_HARNESS_NORESTORE} = 1;
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
    foreach my $itr (1..10) {
      my $res = $call->();
      is_deeply($res, [1, (0)x9], "Testing iteration 1 where we can load an adaptor") if $itr == 1;
      is_deeply($res, [(0)x10], "Testing iteration $itr where we can no longer load adaptors") if $itr > 1;
    }
    ok("Calling of single-threaded load went off without any problems");
  }
}

done_testing();
