#!/software/bin/perl -w

# Populate meta table with (e.g.) genebuild.level = toplevel if all genes are
# top level. Using v41 API code this can speed fetching & dumping greatly.
#

use strict;
use DBI;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my ( $host, $user, $pass, $port, $dbpattern, $print, $sing_db_name );

$port = 3306;

GetOptions( "dbhost|host=s",       \$host,
            "dbuser|user=s",       \$user,
            "dbpass|pass=s",       \$pass,
            "dbport|port=i",       \$port,
            "dbpattern|pattern=s", \$dbpattern,
            "dbname=s",            \$sing_db_name,
            "print",               \$print,
            "help",                \&usage );

if ( !$host || ( !$dbpattern && !$sing_db_name ) ) {
  usage();
}

$dbpattern = $sing_db_name unless $dbpattern;

my @feature_types = qw[gene transcript exon repeat_feature dna_align_feature protein_align_feature simple_feature prediction_transcript prediction_exon];

# loop over databases

my $dsn = "DBI:mysql:host=$host";
$dsn .= ";port=$port" if ($port);

my $db = DBI->connect( $dsn, $user, $pass );

my @dbnames =
  map { $_->[0] } @{ $db->selectall_arrayref("show databases") };

for my $dbname (@dbnames) {

  next if ( $dbname !~ /$dbpattern/ );

  my $dba =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( '-host'    => $host,
                                        '-port'    => $port,
                                        '-user'    => $user,
                                        '-pass'    => $pass,
                                        '-dbname'  => $dbname,
                                        '-species' => $dbname );

  my $ma = $dba->get_MetaContainer();

  my @inserted;
  my @not_inserted;

  foreach my $type (@feature_types) {

    delete_existing( $ma, $type ) if ( !$print );

    if ( can_use_key( $dba, $type ) ) {

      insert_key( $ma, $type ) if ( !$print );
      push @inserted, $type;

    } else {

      push @not_inserted, $type;

    }

  }

  print "$dbname inserted keys for " . join( ", ", @inserted ) . ".\n"
    if (@inserted);
  print "$dbname did not insert keys for "
    . join( ", ", @not_inserted ) . ".\n"
    if (@not_inserted);

} ## end for my $dbname (@dbnames)

#------------------------------------------------------------------------------

sub delete_existing {

  my ( $ma, $type ) = @_;

  $ma->delete_key( $type . "build.level" );

}

#------------------------------------------------------------------------------

sub can_use_key {

  my ( $dba, $type ) = @_;

# compare total count of typewith number of toplevel type, if they're the same,
# then we can use the key

  my $sth = $dba->dbc()->prepare("SELECT COUNT(*) FROM $type");
  $sth->execute();
  my $total = ( $sth->fetchrow_array() )[0];

  $sth =
    $dba->dbc()
    ->prepare(   "SELECT COUNT(*) "
               . "FROM $type t, seq_region_attrib sra, attrib_type at "
               . "WHERE t.seq_region_id=sra.seq_region_id "
               . "AND sra.attrib_type_id=at.attrib_type_id "
               . "AND at.code='toplevel'" );
  $sth->execute();
  my $toplevel = ( $sth->fetchrow_array() )[0];

  return $total == $toplevel;
}

#------------------------------------------------------------------------------

sub insert_key {
  my ( $ma, $type ) = @_;
  $ma->store_key_value( $type . "build.level", "toplevel" );
}

#------------------------------------------------------------------------------

sub usage {
  print <<EOF; exit(0);

Populate meta table with (e.g.) genebuild.level = toplevel if all genes
are top level. Using v41 API code this can speed fetching and dumping
greatly.

Usage: perl $0 <options>

  --host|dbhost     Database host to connect to.

  --port|dbport     Database port to connect to (default is 3306).

  --dbpattern       Database name regexp

  --user|dbuser     Database username.

  --pass|dbpass     Password for user.

  --print           Just print, don't insert or delete keys.

  --help            This message.

  --verbose         Prints out the name of the database
                    which is going to be patched.



   If you like to patch just a single database you can also use this
   command line :

   perl meta_levels.pl --host ... --user ... --pass ... --port ... --name ...
EOF
} ## end sub usage
