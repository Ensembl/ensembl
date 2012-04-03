#!/usr/bin/env perl

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBConnection;
use Getopt::Long;

my $help = 0;
my ( $host, $port, $user, $pass, $dbpattern );

$port = '3306';

my @table_names = qw(
  assembly_exception
  density_feature
  ditag
  ditag_feature
  dna_align_feature
  exon
  gene
  karyotype
  marker_feature
  misc_feature
  prediction_exon
  prediction_transcript
  protein_align_feature
  qtl_feature
  repeat_feature
  simple_feature
  splicing_event
  transcript
);

sub usage {
  print <<USAGE_END;
USAGE:

  $0 --dbhost=ens-staging1 [--dbport=3306] \\
  \t--dbuser=ensadmin --dbpass=XXX \\
  \t--dbpattern=core

  $0 --help

  --dbpattern   Specifies a regular expression for (possibly) matching
                multiple databases.

  --help        Displays this help text.

This script will dump the current meta_coord table to a backup file in
the current directory.  Then it will update the meta_coord table for the
data in the following tables:

USAGE_END

  print( "\t", join( "\n\t", @table_names ), "\n" );

}

if ( scalar(@ARGV) == 0 ) {
  usage();
  exit 0;
}

if ( !GetOptions( 'help!'                  => \$help,
                  'dbhost|host=s'          => \$host,
                  'dbport|port=i'          => \$port,
                  'dbuser|user=s'          => \$user,
                  'dbpass|password|pass=s' => \$pass,
                  'dbpattern=s'            => \$dbpattern
     ) ||
     $help )
{
  usage();
  exit;
}

my $dsn = "DBI:mysql:host=$host;port=$port";

my $db = DBI->connect( $dsn, $user, $pass );

my @dbnames =
  map { $_->[0] } @{ $db->selectall_arrayref("SHOW DATABASES") };

for my $dbname (@dbnames) {

  if ( $dbname !~ /$dbpattern/ ) { next }

  print("==> Looking at $dbname...\n");

  my $dbc =
    new Bio::EnsEMBL::DBSQL::DBConnection( -host   => $host,
                                           -port   => $port,
                                           -user   => $user,
                                           -pass   => $pass,
                                           -dbname => $dbname );

  if ( system( "mysql --host=$host --port=$port " .
                 "--user=$user --password='$pass' " .
                 "--database=$dbname --skip-column-names " .
                 "--execute='SELECT * FROM meta_coord'" .
                 ">$dbname.meta_coord.backup" ) )
  {
    warn( "Can't dump the original meta_coord for back up " .
          "(skipping this database)\n" );
    next;
  }
  else {
    print STDERR "Original meta_coord table backed up in " .
      "$dbname.meta_coord.backup\n";
  }

  foreach my $table_name (@table_names) {
    print("Updating $table_name table entries... ");

    my $sql = "DELETE FROM meta_coord WHERE table_name = '$table_name'";
    $dbc->do($sql);

    $sql =
      "INSERT INTO meta_coord " .
      "SELECT '$table_name', s.coord_system_id, " .
      "MAX( t.seq_region_end - t.seq_region_start + 1 ) " .
      "FROM $table_name t JOIN seq_region s USING (seq_region_id) " .
      "GROUP BY s.coord_system_id";
    $dbc->do($sql);

    print("done\n");
  }

  print("==> Done with $dbname\n");
} ## end for my $dbname (@dbnames)

print("==> All done.\n");
