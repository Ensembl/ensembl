use strict;

use Getopt::Long;

use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Mapper::RangeRegistry;

my ($s_host, $s_port, $s_user, $s_pass, $s_dbname, $t_host, $t_port, $t_user, $t_pass, $t_dbname, );

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
           'help'               => sub { usage(); exit(0); });

$s_port ||= 3306; $t_port ||= 3306;

usage() if(!$s_user || !$s_dbname || !$s_host || !$t_user || !$t_dbname || !$t_host);

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

run();

sub run {

  # compare total counts first
  my $count_sql = "SELECT COUNT(*) FROM xref x, object_xref ox WHERE x.xref_id=ox.xref_id AND $restrict_sql";
  my $s_sth = $s_db->dbc()->prepare($count_sql);
  $s_sth->execute();
  my $count = ($s_sth->fetchrow_array())[0];
  print "Total oligo xrefs in $s_dbname: $count\n";

  my $t_sth = $t_db->dbc()->prepare($count_sql);
  $t_sth->execute();
  $count = ($t_sth->fetchrow_array())[0];
  print "Total oligo_xrefs in $t_dbname: $count\n";

  # cache all mappings from each database
  # key format: transcript_stable_id:array:probe
  my %source_mappings = cache_mappings($s_db);
  my %target_mappings = cache_mappings($t_db);

  # compare mappings
  my $source_only = 0;
  my $target_only = 0;
  my $both = 0;

  foreach my $key (keys %source_mappings) {

    if ($target_mappings{$key}) {
      $both++;
    } else {
      $source_only++;
    }

  }

  print "$both mappings in both databases\n";
  print "$source_only mappings in source only\n";

  foreach my $key (keys %target_mappings) {

    $target_only++ if (!$source_mappings{$key});

  }

  print "$target_only mappings in target only\n";

}

# ----------------------------------------------------------------------

sub cache_mappings {

  my ($db) = @_;

  my %mappings;

  # TODO - fix SQL when external_db changes
  my $sth = $db->dbc()->prepare("SELECT tsi.stable_id, e.db_name, x.dbprimary_acc FROM xref x, external_db e, object_xref ox, transcript_stable_id tsi WHERE x.xref_id=ox.xref_id AND x.external_db_id=e.external_db_id AND ox.ensembl_object_type='Transcript' AND ox.ensembl_id=tsi.transcript_id AND $restrict_sql");
  $sth->execute();
  while (my @row = $sth->fetchrow_array()) {

    #my $key = $row[0] . ":" . $row[1]. ":" . $row[2]; # TODO - add array back in
    my $key = $row[0] . ":" . $row[2];
    $mappings{$key} = $key;

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
		
   [--target_port]        The port to use.. Defaults to 3306.
		
   --target_user          Database username.
		
   --target_pass          Password for target_user, if required.
		
   --target_dbname        Database name.

  MISCELLANEOUS:

  [--help]            This text.


EOF

  exit(0);

}
