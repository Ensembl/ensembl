#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Script that selects the best candidate for canonical transcripts on
# each gene.

# For usage instructions, run ./select_canonical_transcripts.pl -help


use strict;
use warnings;

use Getopt::Long;
use IO::File;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::TranscriptSelector;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;

my ($host, $port, $dbname, $user, $pass);
my ($dnahost, $dnaport, $dnadbname, $dnauser, $dnapass);
my ($ccds_host, $ccds_dbname, $ccds_user, $ccds_port, $ccds_pass);
my ($log_path, $help);

Bio::EnsEMBL::Registry->no_cache_warnings(1);

my $coord_system_name;
my $seq_region_name;

# keep as undefined unless you only want to run on a specific analysis
my $logic_name;

my $write = 0;
my $include_non_ref = 1;
my $include_duplicates;
my $verbose = 0;

GetOptions( 'dbhost:s'            => \$host,
            'dbport:n'            => \$port,
            'dbname:s'            => \$dbname,
            'dbuser:s'            => \$user,
            'dbpass:s'            => \$pass,
            'dnadbhost:s'         => \$dnahost,
            'dnadbname:s'         => \$dnadbname,
            'dnadbuser:s'         => \$dnauser,
            'dnadbpass:s'         => \$dnapass,
            'dnadbport:s'         => \$dnaport,
            'ccdshost:s'          => \$ccds_host,
            'ccdsdbname:s'        => \$ccds_dbname,
            'ccdsuser:s'          => \$ccds_user,
            'ccdsport:s'          => \$ccds_port,
            'ccdspass:s'          => \$ccds_pass,

            'coord_system_name:s' => \$coord_system_name,
            'seq_region_name:s'   => \$seq_region_name,

            'logic_name:s'        => \$logic_name,
            'write!'              => \$write,
            'include_non_ref!'    => \$include_non_ref,
            'include_duplicates'  => \$include_duplicates,
            'verbose!'            => \$verbose, 

            # log file used for analysing choices in bulk
            'log:s'               => \$log_path,
            'help!'               => \$help,
            'h!'                  => \$help
    ) or die "check options\n";

if ($help) { &usage; exit; }

unless ($write) {
  print "You have not used the -write option "
      . "so results will not be written into the database\n";
}

my $dba =
  new Bio::EnsEMBL::DBSQL::DBAdaptor(
      -host   => $host,
      -user   => $user,
      -port   => $port,
      -dbname => $dbname,
      -pass   => $pass,
      -species => 'default',
      -no_cache => 1,
    );

if($dnadbname) {
  if(!$dnauser || !$dnahost) {
    throw ("You must provide user, host and dbname details ".
           "to connect to DNA DB!");
  }
  my $dna_db =
      new Bio::EnsEMBL::DBSQL::DBAdaptor(
          -host   => $dnahost,
          -user       => $dnauser,
          -port       => $dnaport,
          -dbname     => $dnadbname,
          -pass       => $dnapass,
          -species    => 'dna_'.$dba->species(),
      );
  $dba->dnadb($dna_db);
}
else {
  my $dna = check_if_DB_contains_DNA($dba);
  if(!$dna) {
    throw ("Your gene DB contains no DNA. ".
           "You must provide a DNA_DB to connect to");
  }
}

my $ccds_dba;

if ($ccds_dbname) {
  if (!$ccds_user || !$ccds_host) {
    throw ("You must provide user, host and dbname details ".
           "to connect to CCDS DB!");
  }
  $ccds_dba =
  new Bio::EnsEMBL::DBSQL::DBAdaptor(
      -host   => $ccds_host,
      -user   => $ccds_user,
      -port   => $ccds_port,
      -pass   => $ccds_pass,
      -dbname => $ccds_dbname,
      -species => 'ccds_'.$dba->species(),
      -no_cache => 1,
  );
}

my $log_fh;
if ($log_path) {
    $log_fh = IO::File->new($log_path,"w")
        or throw ("Could not open logging file.");
}

my $transcript_selector =
    Bio::EnsEMBL::Utils::TranscriptSelector->new($ccds_dba);

my $slice_adaptor = $dba->get_SliceAdaptor;
my $slices;

if ($seq_region_name) {
    $slices = [
        $slice_adaptor->
          fetch_by_region(
              $coord_system_name,
              $seq_region_name,
              $include_non_ref,
              $include_duplicates
          )
        ];
}
else {
    if (!$coord_system_name) {
        throw 'Requires a coordinate system name to function in this mode';
    }
    $slices = $slice_adaptor->
        fetch_all($coord_system_name,
                  '',
                  $include_non_ref,
                  $include_duplicates
        );
}

my $canonical_changes = 0;
my $total_genes = 0;

my @change_list;

foreach my $slice (@$slices) {
    my $genes = $slice->
        get_all_Genes($logic_name, undef, 1);

    while (my $gene = shift @$genes) {
        $total_genes++;
        my $new_canonical = $transcript_selector->
            select_canonical_transcript_for_Gene($gene);

        my $old_canonical = $gene->canonical_transcript;

        if ( !defined $old_canonical ) {
            # Original canonical transcript is now absent, or never set.
            if ($log_fh) {
                print $log_fh "//\n";
                print $log_fh "Old=[undef,undef,undef,undef,undef,undef,undef]\n";
                printf $log_fh "New=[%s,%s,%s,%s,%s,%s,'%s']\n",
                  @{ $transcript_selector->encode_transcript($new_canonical) };
            }
            push @change_list, [ $gene->dbID, $new_canonical->dbID ];
            $canonical_changes++;
        }

        elsif ($new_canonical->dbID != $old_canonical->dbID) {
            no warnings 'uninitialized';
            printf
                "%s (%s) changed transcript from %s (%s) to %s (%s)\n",
                $gene->stable_id,
                $gene->dbID,
                $old_canonical->stable_id,
                $old_canonical->dbID,
                $new_canonical->stable_id,
                $new_canonical->dbID;

            push @change_list, [ $gene->dbID, $new_canonical->dbID ];
            $canonical_changes++;

            if ($verbose) {
                printf "Old transcript: [%s,%s,%s,%s,%s,%s,%s]\n",
                    @{ $transcript_selector->encode_transcript($old_canonical) };
                printf "New transcript: [%s,%s,%s,%s,%s,%s,%s]\n",
                    @{ $transcript_selector->encode_transcript($new_canonical) };
            }

            if ($log_fh) {
                print $log_fh "//\n";
                printf $log_fh "Old=[%s,%s,%s,%s,%s,%s,'%s']\n",
                  @{ $transcript_selector->encode_transcript($old_canonical) };
                printf $log_fh "New=[%s,%s,%s,%s,%s,%s,'%s']\n",
                  @{ $transcript_selector->encode_transcript($new_canonical) };
            }
        }
    }
}

print
    "Canonical transcript alterations: $canonical_changes ".
    "from $total_genes genes\n";
if ($log_fh) {$log_fh->close}



## Change database entries.
if ($write) {
    my $gene_update_sql = "UPDATE gene SET canonical_transcript_id = ? where gene_id = ?";
    my $gene_sth = $dba->dbc->prepare($gene_update_sql);

    print "Updating database with new canonical transcripts...\n";
    foreach my $change (@change_list) {
        print "Changin' ". $change->[1]. " on ". $change->[0]."\n" if $verbose;
        $gene_sth->execute( $change->[1], $change->[0]);
    }
    print "Done\n";
}

print "Done\n";



## Subroutines

sub check_if_DB_contains_DNA {
  my $dba = shift;
  my $sql_command = "SELECT COUNT(*) FROM dna";

  my $sth = $dba->dbc->prepare($sql_command);
  $sth->execute();

  my @dna_array = $sth->fetchrow_array;

  if ( $dna_array[0] > 0 ) {
    print
      "Your DB ". $dba->dbc->dbname. " contains DNA sequences. ".
      "No need to attach a DNA_DB to it.\n"
        if $verbose;
    return 1;
  }
  else {
    print
      "Your DB " . $dba->dbc->dbname . " does not contain DNA sequences.\n"
        if $verbose;
    return 0;
  }
}


## POD anyone?

sub usage {
print <<EOS

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
    
    -dnadbport    DNA Database port
    
    -dnadbpass    DNA Database pass

    -ccdsdbname  CCDS database name

    -ccdshost    CCDS database host

    -ccdsuser    CCDS database user
    
    -ccdspass    CCDS database pass
    
    -ccdsport    CCDS database port

Other optional arguments:

    -coord_system_name    Coordinate system to use

    -include_non_ref      Specify if the non_reference regions should 
                          be _excluded_. (default: include) 

    -include_duplicates   Specify if the duplicate regions should be 
                          _included_. eg. Human PAR on Y (default: exclude) 

    -seq_region_name      Chromosome name if running a single seq_region

    -write                Specify if results should be written to the database

    -verbose              Increase verbosity of output messages

    -log                  Dump decision matrices into a log file for analysis

To check the script has run correctly you can run the
CoreForeignKeys healthcheck:

./run-healthcheck.sh -d dbname -output problem CoreForeignKeys

A warning about not using CCDS is perfectly acceptible when not
running on Human, Mouse and Zebrafish.

EOS
;

}
