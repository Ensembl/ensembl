# Script that selects the best candidate for canonical transcripts on each
# gene.
# For usage instructions, run ./select_canonical_transcripts.pl -help

#!/usr/bin/env perl

use strict;
use warnings;

use Bio::EnsEMBL::Utils::CliHelper;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::TranscriptSelector;

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();

# get the basic options for connecting to a database server
my $opt_definitions = [ @{ $cli_helper->get_dba_opts() },
              @{ $cli_helper->get_dba_opts('dna') },
              @{ $cli_helper->get_dba_opts('ccds') } ];
push( @{$opt_definitions}, "coord_system_name:s" );
push( @{$opt_definitions}, "logic_name:s" );
push( @{$opt_definitions}, "write!" );
push( @{$opt_definitions}, "include_non_ref!" );
push( @{$opt_definitions}, "include_duplicates" );
push( @{$opt_definitions}, "verbose!" );
push( @{$opt_definitions}, "seq_region_name:s");
# process the command line with the supplied options plus a help subroutine
my $opts = $cli_helper->process_args( $opt_definitions, \&usage );

$opts->{'write'}           ||= 0;
$opts->{'include_non_ref'} ||= 1;
$opts->{'verbose'}         ||= undef;

unless ( $opts->{'write'} ) {
  print "You have not used the -write option "
    . "so results will not be written into the database\n";
}

my @db_args = @{ $cli_helper->get_dba_args_for_opts($opts) };
my @dnadb_args = @{ $cli_helper->get_dba_args_for_opts( $opts, 'dna' ) };
my @ccdsdb_args = @{ $cli_helper->get_dba_args_for_opts( $opts, 'ccds' ) };

if ( defined $dnadb_args[0] && scalar(@dnadb_args) != scalar(@db_args) ) {
  throw "Different number of DBAs found for DB and DNADB";
}
if ( defined $ccdsdb_args[0] && scalar(@ccdsdb_args) != scalar(@db_args) ) {
  throw "Different number of DBAs found for DB and CCDSDB";
}

while (my $db_args = shift (@db_args)) {
    # synchronise removal of dnadb info from array
    my $dna_args = shift (@dnadb_args);
    my $ccds_args = $ccdsdb_args[0]; # Always only have one CCDS source
    
    my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%$db_args);
    
    if (!check_if_DB_contains_DNA($dba)) {
        if ($dna_args) {
            my $dna_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new( %{ $dna_args } );
            $dba->dnadb($dna_dba);
        } else {
            throw("Your gene DB contains no DNA. You must provide a DNA_DB to connect to");
        }
    }
    my $ccds_dba;
    if ($ccds_args) {
        $ccds_dba =  Bio::EnsEMBL::DBSQL::DBAdaptor->new( %{$ccds_args} );
    } else {
        $ccds_dba = undef;
    }
    
    my $transcript_selector = Bio::EnsEMBL::Utils::TranscriptSelector->new($ccds_dba);
    
    my $slice_adaptor = $dba->get_SliceAdaptor;
    my $slices;
    if ($opts->{'seq_region_name'}) {
        $slices = [$slice_adaptor->fetch_by_region($opts->{'coord_system_name'},
                                  $opts->{'seq_region_name'},
                                  $opts->{'include_non_ref'},
                                  $opts->{'include_duplicates'}) ];
    } else {
        if (!$opts->{'coord_system_name'}) {throw 'Requires a coordinate system name to function in this mode';}
        $slices = $slice_adaptor->fetch_all($opts->{'coord_system_name'},
                                            '',
                                            $opts->{'include_non_ref'},
                                            $opts->{'include_duplicates'});
    }
    my $canonical_changes = 0;
    my $total_genes = 0;
    
    foreach my $slice (@$slices) {
        my $genes = $slice->get_all_Genes($opts->{logic_name}, undef, 1);
        
        foreach my $gene (@$genes) {
            $total_genes++;
            my $canonical = $transcript_selector->select_canonical_transcript_for_Gene($gene);
            
            my $old_canonical = $gene->canonical_transcript;
            
            if ($canonical->stable_id ne $old_canonical->stable_id) {
                printf "%s changed transcript from %s to %s\n",
                    $gene->stable_id,$canonical->stable_id,$old_canonical->stable_id;
                $canonical_changes++;
                
                if ($opts->{'verbose'}) {
                    printf "Old transcript: [%s,%s,%s,%s,%s,%s]\n",
                        @{ $transcript_selector->encode_transcript($old_canonical) };
                    printf "New transcript: [%s,%s,%s,%s,%s,%s]\n",
                        @{ $transcript_selector->encode_transcript($canonical) };
                }
            }
        }
    }

    print "Canonical transcript alterations: ".$canonical_changes." from ".$total_genes." genes\n";
    
    
    
}

sub check_if_DB_contains_DNA {
  my ($dba)        = @_;
  my $sql_command = "select count(*) from dna";
  my $sth         = $dba->dbc->prepare($sql_command);
  $sth->execute();
  my @dna_array = $sth->fetchrow_array;
  if ( $dna_array[0] > 0 ) {
    print "Your DB "
      . $dba->dbc->dbname
      . " contains DNA sequences. No need to attach a "
      . "DNA_DB to it.\n"
      if ( $opts->{verbose} );
    return 1;
  } else {
    print "Your DB " . $dba->dbc->dbname . " does not contain DNA sequences.\n"
      if ( $opts->{verbose} );
    return 0;
  }
}

sub usage {
print "
Example usage: perl set_canonical_transcripts.pl -dbhost host -dbuser user 
     -dbpass *** -dbname dbname -dbport 3306 -coord_system toplevel -write

Script options:

    -dbname       Database name

    -dbhost       Database host

    -dbport       Database port

    -dbuser       Database user

    -dbpass       Database password

Optional DB connection arguments:

    -dnadbname    DNA Database name

    -dnadbhost    DNA Database host

    -dnadbuser    DNA Database user

    -ccdsdbname  CCDS database name

    -ccdshost    CCDS database host

    -ccdsuser    CCDS database user

Other optional arguments:

    -coord_system_name    Coordinate system to use

    -include_non_ref      Specify if the non_reference regions should 
                          be _excluded_. (default: include) 

    -include_duplicates   Specify if the duplicate regions should be 
                          _included_. eg. Human PAR on Y (default: exclude) 

    -seq_region_name      Chromosome name if running a single seq_region

    -write                Specify if results should be written to the database

    -verbose              Increase verbosity of output messages


To check the script has run correctly you can run the
CoreForeignKeys healthcheck:

./run-healthcheck.sh -d dbname -output problem CoreForeignKeys
";
    
}