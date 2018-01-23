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


=head1 NAME



=head1 SYNOPSIS

run_all.pl [arguments]

Required arguments:

  --dbname, db_name=NAME              database name NAME
  --host, --dbhost, --db_host=HOST    database host HOST
  --port, --dbport, --db_port=PORT    database port PORT
  --user, --dbuser, --db_user=USER    database username USER
  --pass, --dbpass, --db_pass=PASS    database passwort PASS

Optional arguments:

  --conffile, --conf=FILE             read parameters from FILE
                                      (default: conf/Conversion.ini)

  --logfile, --log=FILE               log to FILE (default: *STDOUT)
  --logpath=PATH                      write logfile to PATH (default: .)
  --logappend, --log_append           append to logfile (default: truncate)
  --loglevel=LEVEL                    define log level (default: INFO)

  -i, --interactive=0|1               run script interactively (default: true)
  -n, --dry_run, --dry=0|1            don't write results to database
  -h, --help, -?                      print help (this message)

=head1 DESCRIPTION



=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<http://lists.ensembl.org/mailman/listinfo/dev>

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin);
use Bio::EnsEMBL::Utils::ConfParser;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::Utils::ScriptUtils qw(path_append);
use Bio::EnsEMBL::DBSQL::DBAdaptor;

# parse configuration and commandline arguments
my $conf = new Bio::EnsEMBL::Utils::ConfParser(
  -SERVERROOT => "$Bin/../../../..",
  -DEFAULT_CONF => "$Bin/../default.conf"
);

$conf->parse_options(
  'sourcehost|source_host=s' => 1,
  'sourceport|source_port=n' => 1,
  'sourceuser|source_user=s' => 1,
  'sourcepass|source_pass=s' => 0,
  'sourcedbname|source_dbname=s' => 1,
  'targethost|target_host=s' => 1,
  'targetport|target_port=n' => 1,
  'targetuser|target_user=s' => 1,
  'targetpass|target_pass=s' => 0,
  'targetdbname|target_dbname=s' => 1,
  'althost|alt_host=s' => 0,
  'altport|alt_port=n' => 0,
  'altuser|alt_user=s' => 0,
  'altpass|alt_pass=s' => 0,
  'altdbname|alt_dbname=s' => 1,
  'basedir|basedir=s' => 1,
  'lsf!' => 0,
  'lsf_opt|lsfopt=s' => 0,
  'suffix|sfx=s' => 0,
  'debug1|d1=s' => 1,
  'debug2|d2=s' => 1,
  'type|t=s' => 0,
);

# set default logpath
unless ($conf->param('logpath')) {
  $conf->param('logpath', path_append($conf->param('basedir'), 'log'));
}

# assume both dbs are on same host unless specified otherwise
foreach my $p (qw(host port user pass)) {
  unless ($conf->param("alt$p")) {
    $conf->param("alt$p", $conf->param("target$p"));
  }
}

# get log filehandle and print header and parameters to logfile
my $logger = new Bio::EnsEMBL::Utils::Logger(
  -LOGFILE      => $conf->param('logfile'),
  -LOGAUTO      => $conf->param('logauto'),
  -LOGAUTOBASE  => 'compare_results',
  -LOGPATH      => $conf->param('logpath'),
  -LOGAPPEND    => $conf->param('logappend'),
  -LOGLEVEL     => $conf->param('loglevel'),
);

# if user wants to run via lsf, submit script with bsub (this will exit this
# instance of the script)
&bsubmit if ($conf->param('lsf'));

# initialise log
$logger->init_log($conf->list_param_values);

# connect to dbs
my $dba_s = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -host   => $conf->param("sourcehost"),
        -port   => $conf->param("sourceport"),
        -user   => $conf->param("sourceuser"),
        -pass   => $conf->param("sourcepass"),
        -dbname => $conf->param("sourcedbname"),
        -group  => 'source',
);
$dba_s->dnadb($dba_s);

my $dba1 = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -host   => $conf->param("targethost"),
        -port   => $conf->param("targetport"),
        -user   => $conf->param("targetuser"),
        -pass   => $conf->param("targetpass"),
        -dbname => $conf->param("targetdbname"),
        -group  => 'target',
);
$dba1->dnadb($dba1);
my $dbh1 = $dba1->dbc->db_handle;

my $dba2 = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -host   => $conf->param("althost"),
        -port   => $conf->param("altport"),
        -user   => $conf->param("altuser"),
        -pass   => $conf->param("altpass"),
        -dbname => $conf->param("altdbname"),
        -group  => 'alt',
);
$dba2->dnadb($dba2);
my $dbh2 = $dba2->dbc->db_handle;

# compare mapping results
my $type = $conf->param('type') || 'gene';
my $function = "compare_${type}s";
no strict 'refs';
&$function;

# finish logfile
$logger->finish_log;

### END main ###


sub compare_genes {
  $logger->info("Comparing genes...\n\n", 0, 'stamped');
  &compare_features('gene');
  $logger->info("Done\n\n");
}


sub compare_transcripts {
  $logger->info("Comparing transcripts...\n\n", 0, 'stamped');
  &compare_features('transcript');
  $logger->info("Done\n\n");
}


sub compare_translations {
  $logger->info("Comparing translations...\n\n", 0, 'stamped');
  &compare_features('translation');
  $logger->info("Done\n\n");
}


sub compare_exons {
  $logger->info("Comparing exons...\n\n", 0, 'stamped');
  &compare_features('exon');
  $logger->info("Done\n\n");
}


sub compare_features {
  my $ftype = shift;

  # get a filehandle to write results for debugging
  my $path = path_append($conf->param('basedir'), 'debug');
  my $file = "$path/${ftype}_diff.txt";
  open(my $fh, '>', $file) or die "Can't open $file for writing: $!\n";

  # read scores from files
  my $scores = {};
  unless ($ftype eq 'translation') {
    foreach my $path1 (qw(debug1 debug2)) {
      my $p1 = $conf->param($path1);
      my $file1 = "$p1/${ftype}_scores.txt";
      open(my $fh1, '<', $file1) or die "Can't open $file1 for reading: $!\n";

      while (my $line = <$fh1>) {
        chomp $line;
        my ($old_id, $new_id, $score) = split(/\s+/, $line);
        $score = sprintf("%.6f", $score);

        # remember the highest score for each new_id
        if ($score > $scores->{$path1}->{$new_id}) {
          $scores->{$path1}->{$new_id} = $score;
        }
      }

      close($fh1);
    }
  }
  
  #
  # fetch all features from both runs and create lookup hash by stable_id
  #
  $logger->info("Fetching ${ftype} data from dbs...\n", 0, 'stamped');

  # db 2
  my $sql1 = qq(SELECT ${ftype}_id, stable_id FROM ${ftype}_stable_id);
  my $sth1 = $dbh1->prepare($sql1);
  $sth1->execute;
  
  my %fsi1 = ();
  my %fii1 = ();
  
  while (my $r = $sth1->fetchrow_arrayref) {
    # create lookup hashes of dbID to stable ID and vice versa
    $fsi1{$r->[1]} = $r->[0];
    $fii1{$r->[0]} = $r->[1];
  }

  $sth1->finish;

  # db 2
  my $suffix = $conf->param('suffix');
  my $sql2 = qq(SELECT ${ftype}_id, stable_id FROM ${ftype}_stable_id${suffix});
  my $sth2 = $dbh2->prepare($sql2);
  $sth2->execute;
  
  my %fsi2 = ();
  my %fii2 = ();
  
  while (my $r = $sth2->fetchrow_arrayref) {
    # create lookup hashes of dbID to stable ID and vice versa
    $fsi2{$r->[1]} = $r->[0];
    $fii2{$r->[0]} = $r->[1];
  }

  $sth2->finish;
  
  $logger->info("Done.\n\n", 0, 'stamped');

  #
  # get max(gene_stable_id) from source db
  #
  my $dbh = $dba_s->dbc->db_handle;
  my $sql = qq(SELECT max(stable_id) FROM ${ftype}_stable_id);
  my $sth = $dbh->prepare($sql);
  $sth->execute;
  my ($max_stable_id) = $sth->fetchrow_array;
  $sth->finish;

  #
  # now loop over dbIDs in db 1 and compare results with db2
  #
  $logger->info("Comparing results...\n", 0, 'stamped');

  my @stat_keys = qw(TOT NN II NE EN EE);
  my %stats = map { $_ => 0 } @stat_keys;
  my $fmt = "%-3s%6d %-20s %-20s %-10s %-10s\n";

  foreach my $dbID1 (sort { $a <=> $b } keys %fii1) {
    $stats{TOT}++;
    my $status;
    
    my $sid1 = $fii1{$dbID1};
    my $sid2 = $fii2{$dbID1};

    # db 1 has new stable ID
    if (($max_stable_id cmp $sid1) == -1) {
      
      # db 2 has new stable ID too
      if (($max_stable_id cmp $sid2) == -1) {
        $status = 'NN';

      # db 2 reuses an existing stable ID
      } else {
        $status = 'NE';
      }
      
    # else db 1 reused an existing stable ID
    } else {
      
      # db 2 uses the same stable ID
      if ($sid1 eq $sid2) {
        $status = 'II';

      # db 2 has a new stable ID
      } elsif (($max_stable_id cmp $sid2) == -1) {
        $status = 'EN';
      
      # db 2 reuses an existing (but different from db 1) stable ID
      } else {
        $status = 'EE';
      }
    }

    # stats
    $stats{$status}++;

    # print result line (status dbID sid1 sid2)
    print $fh sprintf($fmt, $status, $dbID1, $sid1, $sid2,
      $scores->{'debug1'}->{$dbID1}, $scores->{'debug2'}->{$dbID1});
  }

  close($fh);
  
  $logger->info("Done.\n\n", 0, 'stamped');

  # print stats
  $logger->info("Stats:\n");
  foreach my $key (@stat_keys) {
    $logger->info(sprintf("  %-5s%8d (%6s)\n", $key, $stats{$key}, sprintf("%3.1f%%", 100*$stats{$key}/$stats{'TOT'})));
  }

}


sub bsubmit {
  #
  # build bsub commandline
  #

  # automatically create a filename for lsf output
  my $cmd = 'bsub -o '.$conf->param('logpath');
  $cmd .= "/lsf_compare_".$logger->log_auto_id.'.out';

  # add extra lsf options as configured by the user
  $cmd .= ' '.$conf->param('lsf_opt');

  # this script's name
  $cmd .= " $0";

  # options for this script
  my $options = $conf->create_commandline_options(
    logautoid => $logger->log_auto_id,
    interactive   => 0,
    lsf       => 0,
  );
  $cmd .= " $options";

  #
  # execute bsub
  #
  print "\nRe-executing via lsf:\n";
  print "$cmd\n\n";

  exec($cmd) or die "Could not exec $0 via lsf: $!\n";
  #exit;
}

