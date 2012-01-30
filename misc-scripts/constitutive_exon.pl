#!/usr/local/ensembl/bin/perl -w

### This script is used to mark exons which are used in more than one
### transcript linked to a single gene. The script does not limit
### by differing biotypes. Exons are identified as being the same if they
### occupy identical genomic locations if the -uselocations flag is specified

use strict;
use warnings;

use Getopt::Long;
use DBI;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

sub usage {
  print("Command line switches:\n");
  print("\t--dbhost=<host>      \tDatabase server host\n");
  print("\t--dbport=<port>      \tDatabase server port (optional)\n");
  print("\t--dbuser=<user>      \tDatabase user\n");
  print("\t--dbpass=<password>  \tUser password (optional)\n");
  print("\t--dbpattern=<regex>  "
      . "\tRegular expresssion to match database names\n" );
  print("\t--uselocations       "
      . "\tUse locations to find indentical exons rather than DB identifiers(optional)\n" );
  print("\t--help               "
      . "\tDisplays this info and exits (optional)\n" );
  exit;
}

my ( $dbhost, $dbport, $dbuser, $dbpass, $dbpattern, $uselocations );

GetOptions(
  'dbhost|host=s'       => \$dbhost,
  'dbport|port=i'       => \$dbport,
  'dbuser|user=s'       => \$dbuser,
  'dbpass|pass=s'       => \$dbpass,
  'dbpattern|pattern=s' => \$dbpattern,
  'uselocations'        => \$uselocations,
  'help|h'              => \&usage
);

$dbport ||= 3306;

if ( !( defined($dbhost) && defined($dbuser) && defined($dbpattern) ) )
{
  usage();
}

my $dsn = sprintf( 'dbi:mysql:host=%s;port=%d', $dbhost, $dbport );
my $dbh = DBI->connect( $dsn, $dbuser, $dbpass );

my @dbnames =
  map { $_->[0] } @{ $dbh->selectall_arrayref('SHOW DATABASES') };

foreach my $dbname (@dbnames) {
  if ( $dbname !~ /$dbpattern/ ) { next }

  printf( "Connecting to '%s'\n", $dbname );

  my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    '-host'    => $dbhost,
    '-port'    => $dbport,
    '-user'    => $dbuser,
    '-pass'    => $dbpass,
    '-dbname'  => $dbname,
    '-species' => $dbname
  );

  my $gene_adaptor = $dba->get_GeneAdaptor();

  # Retrieve all Gene dbIDs
  my @gene_dbIDs = sort { $a <=> $b } @{ $gene_adaptor->list_dbIDs() };
  my $gene_count = scalar(@gene_dbIDs);
  printf( "There are %d genes to go through...\n", $gene_count );

  my ( $update_0,     $update_1 )     = ( 0, 0 );
  my ( $gene_counter, $exon_counter ) = ( 0, 0 );

  while (
    my @gene_list = @{
      $gene_adaptor->fetch_all_by_dbID_list(
        [ splice( @gene_dbIDs, 0, 1000 ) ] ) } )
  {
    while ( my $gene = shift(@gene_list) ) {
      my @transcripts      = @{ $gene->get_all_Transcripts() };
      my $transcript_count = scalar(@transcripts);

      my %exon_count;
      my %exon_object;

      foreach my $transcript (@transcripts) {
        foreach my $exon ( @{ $transcript->get_all_Exons() } ) {
          my $key = exon_key($exon);
          ++$exon_count{$key};
          $exon_object{$key} = $exon;
        }
      }

      foreach my $exon_key ( keys(%exon_count) ) {
        my $exon = $exon_object{$exon_key};
        my $is_constitutive = $exon->is_constitutive(); 
        if ( $exon_count{$exon_key} == $transcript_count ) {
          if ( !$exon->is_constitutive() ) {
            $is_constitutive = 1;
            ++$update_1;
          }
        } else {
          if ( $exon_object{$exon_key}->is_constitutive() ) {
            $is_constitutive = 0;
            ++$update_0;
          }
        }
        $gene_adaptor->dbc()->sql_helper()->excute_update(
          -SQL => 'update exon set is_constitutive = ? where exon_id =?',
          -PARAMS => [$is_constitutive, $exon->dbID()]
        );

        if ( ( ++$exon_counter % 1000 ) == 0 ) {
          printf(
            "After %d exons (%d genes, %.3g%%):\n"
              . "\tupdated %d to constitutive, "
              . "%d to non-constitutive\n",
            $exon_counter, $gene_counter, 100*$gene_counter/$gene_count,
            $update_1, $update_0
          );
        }
      } ## end foreach my $exon_dbID ( keys...)

      ++$gene_counter;

    } ## end while ( my $gene = shift(...))
  } ## end while ( my @gene_list = @...)

  print( '-' x 4, '{ ', $dbname, ' }',
    '-' x ( 80 - ( length($dbname) + 4 ) ), "\n" );
  print("Summary:\n");
  printf( "\t%d exons in total, %d genes\n",
    $exon_counter, $gene_counter );
  printf( "\tUpdated %d exons to constitutive\n",     $update_1 );
  printf( "\tUpdated %d exons to non-constitutive\n", $update_0 );
  print( '-' x 80, "\n" );
  print("\n");
} ## end foreach my $dbname (@dbnames)

sub exon_key {
  my ($e) = @_;
  if($uselocations) {
    return join(q{:}, 
      ($e->slice()->name() ? $e->slice()->name() : 'undef'),
      $e->seq_region_start(),
      $e->seq_region_end(),
      $e->seq_region_strand()
    );
  }
  return $e->dbID();
}
