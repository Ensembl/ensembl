#!/software/bin/perl

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


=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

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
  'dumppath|dump_path=s' => 1,
  'lsf!' => 0,
  'lsf_opt|lsfopt=s' => 0,
);

# set default logpath
unless ($conf->param('logpath')) {
  $conf->param('logpath', path_append($conf->param('dumppath'), 'log'));
}

# assume both dbs are on same host unless specified otherwise
foreach my $p (qw(host port user pass)) {
  unless ($conf->param("alt$p")) {
    $conf->param("alt$p", $conf->param("target$p"));
  }
}

# get log filehandle and print heading and parameters to logfile
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

my $dba2 = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -host   => $conf->param("althost"),
        -port   => $conf->param("altport"),
        -user   => $conf->param("altuser"),
        -pass   => $conf->param("altpass"),
        -dbname => $conf->param("altdbname"),
        -group  => 'alt',
);
$dba2->dnadb($dba2);

# compare gene mapping results
&compare_genes;

# finish logfile
$logger->finish_log;

### END main ###


sub compare_genes {

  $logger->info("Comparing genes...\n\n", 0, 'stamped');

  # get a filehandle to write results for debugging
  my $path = path_append($conf->param('dumppath'), 'debug');
  my $file = "$path/genes_diff.txt";
  open(my $fh, '>', $file) or die "Can't open $file for writing: $!\n";
  
  #
  # fetch all genes from both dbs and create lookup hash by stable_id
  #
  $logger->info("Fetching genes...\n", 0, 'stamped');
  my $ga1 = $dba1->get_GeneAdaptor;
  my $ga2 = $dba2->get_GeneAdaptor;

  my @genes1 = @{ $ga1->fetch_all(undef, undef, 0) };
  my @genes2 = @{ $ga2->fetch_all(undef, undef, 0) };
  
  $logger->info("Done.\n\n", 0, 'stamped');

  my %gsi1 = map { $_->stable_id => $_ } @genes1;
  my %gi1 = map { $_->dbID => $_ } @genes1;
  my %gsi2 = map { $_->stable_id => $_ } @genes2;
  my %gi2 = map { $_->dbID => $_ } @genes2;

  # get max(gene_stable_id) from source db
  my $dbh = $dba_s->dbc->db_handle;
  my $sql = qq(SELECT max(stable_id) FROM gene_stable_id);
  my $sth = $dbh->prepare($sql);
  $sth->execute;
  my ($max_stable_id) = $sth->fetchrow_array;
  $sth->finish;

  my $fmt1 = "%-20s%-8s%-40s%-1s\n";

  my @stat_keys = qw(TOT OK S I D N);
  my %stats = map { $_ => 0 } @stat_keys;

  #
  # now loop over genes in db 1 and print information about genes not found in 
  # db 2
  #
  foreach my $gsid1 (sort keys %gsi1) {

    my $gene1 = $gsi1{$gsid1};
    my $gene2 = $gsi2{$gsid1};

    $stats{TOT}++;

    # next if gene in db 2 is the same (same stable and internal ID)
    if ($gene2 and ($gene1->dbID == $gene2->dbID)) {
      $stats{OK}++;
      next;
    }

    my $flag;
    my $gene1a;
    my $gene2a;

    if ($gene2) {
      # we found the stable_id, but it was assigned to a different gene in db 2
      $gene1a = $gi1{$gene2->dbID};
      $gene2a = $gi2{$gene1->dbID};
      $flag = 'S';
    } else {
      # the gene was assigned a different stable_id
      $gene2 = $gi2{$gene1->dbID};

      if ($gsi1{$gene2->stable_id}) {
        # stable ID used in db 1, but for different gene
        $flag = 'I';
      } elsif (($max_stable_id cmp $gene2->stable_id) == 1) {
        # stable ID not used in db 1, gene deleted from db 1
        $flag = 'D';
      } else {
        # new stable ID used in db 2
        $flag = 'N';
      }
    }

    my $slice1 = $gene1->feature_Slice;
    my $ss1 = join(':', map { $slice1->$_ }
      qw(seq_region_name start end strand length));
      
    my $slice2 = $gene2->feature_Slice;
    my $ss2 = join(':', map { $slice2->$_ }
      qw(seq_region_name start end strand length));

    # print debug statement
    my $txt1 = sprintf($fmt1,
                       $gene1->stable_id,
                       $gene1->dbID,
                       $ss1,
                       $flag
                      );

    my $txt1a;
    
    if ($gene1a) {
      my $slice1a = $gene1a->feature_Slice;
      my $ss1a = join(':', map { $slice1a->$_ }
        qw(seq_region_name start end strand length));
        
      my $txt1a = sprintf($fmt1,
                         $gene1a->stable_id,
                         $gene1a->dbID,
                         $ss1a,
                         undef
      );
    } else {
      $txt1a = "none\n";
    }

    my $txt2 = sprintf($fmt1,
                       $gene2->stable_id,
                       $gene2->dbID,
                       $ss2,
                       $flag
                      );
    
    my $txt2a;

    if ($gene2a) {
      my $slice2a = $gene2a->feature_Slice;
      my $ss2a = join(':', map { $slice2a->$_ }
        qw(seq_region_name start end strand length));
        
      my $txt2a = sprintf($fmt1,
                         $gene2a->stable_id,
                         $gene2a->dbID,
                         $ss2a,
                         undef
      );
    } else {
      $txt2a = "none\n";
    }

    $logger->info($txt1, 1);
    $logger->info($txt1a, 1) if ($flag eq 'S');
    $logger->info($txt2, 1);
    $logger->info($txt2a, 1) if ($flag eq 'S');
    $logger->info("\n");

    $stats{$flag}++;
  }

  close($fh);
  
  $logger->info("\nDone.\n\n", 0, 'stamped');

  # print stats
  $logger->info("Stats:\n");
  foreach my $key (@stat_keys) {
    $logger->info(sprintf("  %-5s%8.0f\n", $key, $stats{$key}));
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

