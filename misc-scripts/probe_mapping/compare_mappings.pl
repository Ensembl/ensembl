use strict;

use Getopt::Long;

use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Mapper::RangeRegistry;

my ($s_host, $s_port, $s_user, $s_pass, $s_dbname, $t_host, $t_port, $t_user, $t_pass, $t_dbname, $array);

GetOptions('source_host=s'      => \$s_host,
           'source_user=s'      => \$s_user,
           'source_port=i'      => \$s_port,
           'source_pass=s'      => \$s_pass,
           'source_dbname=s'    => \$s_dbname,
	   'target_host=s'      => \$t_host,
           'target_user=s'      => \$t_user,
           'target_port=i'      => \$t_port,
           'target_pass=s'      => \$t_pass,
           'target_dbname=s'    => \$t_dbname,
	   'array=s'            => \$array,
           'help'               => sub { usage(); exit(0); });

$s_port ||= 3306; $t_port ||= 3306;

usage() if(!$s_user || !$s_dbname || !$s_host || !$t_dbname);

$t_host = $s_host if (!$t_host);
$t_port = $s_port if (!$t_port);
$t_user = $s_user if (!$t_user);

my $s_db = new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'   => $s_host,
					      '-port'   => $s_port,
					      '-user'   => $s_user,
					      '-pass'   => $s_pass,
					      '-dbname' => $s_dbname,
					      '-species'=> $s_dbname);

my $t_db = new Bio::EnsEMBL::DBSQL::DBAdaptor('-host'   => $t_host,
					      '-port'   => $t_port,
					      '-user'   => $t_user,
					      '-pass'   => $t_pass,
					      '-dbname' => $t_dbname,
					      '-species'=> $t_dbname);

# TODO fix this when external db is sorted out
my $restrict_sql = "x.external_db_id > 3000 AND x.external_db_id < 3210";
$restrict_sql = " e.db_name='$array'" if ($array);

run();

sub run {

  print "Using only array $array\n" if ($array);

  # compare total counts first
  my $count_sql = "SELECT COUNT(*) FROM xref x, object_xref ox, external_db e WHERE e.external_db_id=x.external_db_id AND x.xref_id=ox.xref_id AND ox.ensembl_object_type='Transcript' AND $restrict_sql";

  my $s_sth = $s_db->dbc()->prepare($count_sql);
  $s_sth->execute();
  my $count = ($s_sth->fetchrow_array())[0];
  print "Total oligo xrefs in $s_dbname: $count\n";

  my $t_sth = $t_db->dbc()->prepare($count_sql);
  $t_sth->execute();
  $count = ($t_sth->fetchrow_array())[0];
  print "Total oligo xrefs in $t_dbname: $count\n";

  # cache all mappings from each database
  # key format: transcript_stable_id:array:probe
  my %source_mappings = cache_mappings($s_db);
  my %target_mappings = cache_mappings($t_db);

  # compare mappings
  my $source_only = 0;
  my $target_only = 0;
  my $both = 0;

  open(SOURCE_ONLY, ">source_only.txt");
  open(TARGET_ONLY, ">target_only.txt");
  open(BOTH,        ">both.txt");

  foreach my $key (keys %source_mappings) {

    if ($target_mappings{$key}) {
      $both++;
      print BOTH "$key\n";
    } else {
      $source_only++;
      print SOURCE_ONLY "$key\n";
    }

  }

  print "$both unique mappings in both databases\n";
  print "$source_only mappings in source only\n";

  foreach my $key (keys %target_mappings) {

    if (!$source_mappings{$key}) {
      $target_only++;
      print TARGET_ONLY "$key\n";
    }

  }

  print "$target_only mappings in target only\n";

}

# ----------------------------------------------------------------------

sub cache_mappings {

  my ($db, $print) = @_;

  my %mappings;

  # TODO - fix SQL when external_db changes
  my $sth = $db->dbc()->prepare("SELECT tsi.stable_id, e.db_name, x.dbprimary_acc, t.biotype FROM xref x, external_db e, object_xref ox, transcript_stable_id tsi, transcript t WHERE x.xref_id=ox.xref_id AND x.external_db_id=e.external_db_id AND ox.ensembl_object_type='Transcript' AND ox.ensembl_id=t.transcript_id AND t.transcript_id=tsi.transcript_id AND $restrict_sql");
  $sth->execute();
  while (my @row = $sth->fetchrow_array()) {

    my $key = $row[0] . "\t" . $row[3] . "\t". $row[1]. "\t" . $row[2];
    $mappings{$key} = $key;

    print "$key\n" if ($print);

  }

  $sth->finish();

  return %mappings;

}

# ----------------------------------------------------------------------

sub usage {

  print << "EOF";

  Compare oligo xrefs between 2 databases.

  perl $0 {options}

  Options ([..] indicates optional):

  ORIGINAL XREFS:

   --source_host          The database server to read the first set of xrefs from.
		
   [--source_port]        The port to use. Defaults to 3306.
		
   --source_user          Database username.
		
   --source_pass          Password for source_user, if required.
		
   --source_dbname        Database name.

  NEW XREFS:

   --target_host          The database server to read the second set of xrefs from.
                          Defaults to source_host.
		
   [--target_port]        The port to use. Defaults to 3306.
		
   --target_user          Database username. Defaults to source_user.
		
   --target_pass          Password for target_user, if required.
		
   --target_dbname        Database name.

  MISCELLANEOUS:

  [--array]           Just compare results from this array.

  [--help]            This text.


EOF

  exit(0);

}
