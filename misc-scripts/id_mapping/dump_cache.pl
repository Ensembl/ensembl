#!/software/bin/perl
#!/usr/local/ensembl/bin/perl

=head1 NAME


=head1 SYNOPSIS

.pl [arguments]

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

Use --sourceschema and --targetschema to specify a schema version (default:
latest). This will be used to determine the subroutine to build the cache. By
default, &build_cache_latest() is run which uses Bio::EnsEMBL::IdMapping::Cache
to read from the database and write the cache.  An alternative subroutine can
use a different module for that, which will usually inherit from the former and
overwrite Cache->build_cache(). This is useful for backwards compatibility with
older schema versions. Once the cache is built, no API access is needed,
therefore the ID mapping application is independent of the underlying database
schema.


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
use Bio::EnsEMBL::Utils::ScriptUtils qw(inject path_append);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

# parse configuration and commandline arguments
my $conf = new Bio::EnsEMBL::Utils::ConfParser(
  -SERVERROOT => "$Bin/../../..",
  -DEFAULT_CONF => "$Bin/default.conf"
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
  'basedir|basedir=s' => 1,
  'chromosomes|chr=s@' => 0,
  'region=s' => 0,
  'biotypes=s@' => 0,
  'lsf_opt_dump_cache|lsfoptdumpcache=s' => 0,
  'cache_method=s' => 0,
  'build_cache_auto_threshold=n' => 0,
  'build_cache_concurrent_jobs=n' => 0,
);

# set default logpath
unless ($conf->param('logpath')) {
  $conf->param('logpath', path_append($conf->param('basedir'), 'log'));
}

# get log filehandle and print heading and parameters to logfile
my $logger = new Bio::EnsEMBL::Utils::Logger(
  -LOGFILE      => $conf->param('logfile'),
  -LOGAUTO      => $conf->param('logauto'),
  -LOGAUTOBASE  => 'dump_cache',
  -LOGAUTOID    => $conf->param('logautoid'),
  -LOGPATH      => $conf->param('logpath'),
  -LOGAPPEND    => $conf->param('logappend'),
  -LOGLEVEL     => $conf->param('loglevel'),
  -IS_COMPONENT => $conf->param('is_component'),
);

# initialise log
$logger->init_log($conf->list_param_values);

# determin cache method to use.
# this can be used to support different caching strategies or access to old
# database schemas.
my $cache_method = $conf->param('cache_method') || 'build_cache_auto';
no strict 'refs';
my $retval = &$cache_method;

# finish logfile
$logger->finish_log;

exit($retval);


### END main ###


sub build_cache_auto {
  # load the cache implementation
  my $cache_impl = 'Bio::EnsEMBL::IdMapping::Cache';
  inject($cache_impl);

  # EG - populate cache for each species in turn
  $logger->debug("\nChecking number of toplevel seq_regions...\n");
  my $max         = 0;
  my @species_ids = @{ get_species_ids("target") };
  foreach my $dbtype (qw(source target)) {

    # populate the cache for each species in turn
    for my $species (@species_ids) {
      $conf->param( 'species_id',   $$species[1] );
      $conf->param( 'species_name', $$species[0] );

      my $cache = $cache_impl->new( -LOGGER => $logger,
                                    -CONF   => $conf, );

      my $num =
        scalar( @{ $cache->slice_names( $dbtype, @$species ) } );

      $max = $num if ( $num > $max );
      $logger->debug( "$dbtype: $num.\n", 1 );
    }
  }

  my $threshold = $conf->param('build_cache_auto_threshold') || 100;
  my $retval;

  if ($max > $threshold) {
    $logger->debug("\nWill use build_cache_all.\n");
    $retval = &build_cache_all;
  } else {
    $logger->debug("\nWill use build_cache_by_seq_region.\n");
    $retval = &build_cache_by_seq_region;
  }

  return $retval;
}


sub build_cache_by_seq_region {

  my %jobs = ();

  # create empty directory for logs
  my $logpath = path_append($conf->param('logpath'), 'dump_by_seq_region');
  system("rm -rf $logpath") == 0 or
    $logger->error("Unable to delete lsf log dir $logpath: $!\n");
  system("mkdir -p $logpath") == 0 or
    $logger->error("Can't create lsf log dir $logpath: $!\n");


  # EG get the list of species IDs for sources and targets
  my @source_species_ids = @{ get_species_ids("source") };
  my @species_ids        = @{ get_species_ids("target") };

  # load the cache implementation
  my $cache_impl = 'Bio::EnsEMBL::IdMapping::Cache';
  inject($cache_impl);

  # EG store the base directory onto which the species ID will be added
  my $basedir = $conf->param('basedir');

  # submit jobs to lsf
  foreach my $dbtype (qw(source target)) {

    # EG iterate over individual species for source and target

    # determine which slices need to be done
    my $filename = "$dbtype.dump_cache.slices.txt";
    open(my $fh, '>', "$logpath/$filename") or
      throw("Unable to open $logpath/$filename for writing: $!");
    
    my $num_jobs = 0;
    for my $species (@species_ids) {
      # EG set config based on species ID in turn
      $conf->param( 'basedir', path_append( $basedir, $$species[1] ) );
      $conf->param( 'species_id',   $$species[1] );
      $conf->param( 'species_name', $$species[0] );
      # EG load cache for current species ID
      my $cache = $cache_impl->new( -LOGGER => $logger,
                                    -CONF   => $conf, );
      foreach my $slice_name (
                        @{ $cache->slice_names( $dbtype, @$species ) } )
      {
        my $type = "$dbtype.$slice_name";
        my $src_species_id;
        for my $src_id (@source_species_ids) {
          if ( $$species[1] == $$src_id[1] ) {
            $src_species_id = $$src_id[1];
            last;
          }
        }
        $logger->info( "\n" . ucfirst($dbtype) . " db...\n",
                       0, 'stamped' );

        foreach my $slice_name ( @{ $cache->slice_names($dbtype) } ) {
          my $type = "$dbtype.$slice_name";
          unless ( $cache->cache_file_exists($type) ) {
            print $fh "$slice_name\n";
            print $fh "$slice_name,$$species[0],$$species[1],"
              . $src_species_id . "\n";
            $num_jobs++;
          }
        }
      }

      unless ($num_jobs) {
        $logger->info("All cache files for $dbtype exist.\n");
        next;
      }
    } ## end for my $species (@species_ids)
    close($fh);
    # EG reset original basedir
    $conf->param( 'basedir', $basedir );

    # build lsf command
    my $lsf_name = 'dump_by_seq_region_'.time;
    my $concurrent = $conf->param('build_cache_concurrent_jobs') || 200;

    my $options = $conf->create_commandline_options(
        logauto       => 1,
        logautobase   => "dump_by_seq_region",
        interactive   => 0,
        is_component  => 1,
        dbtype        => $dbtype,
        cache_impl    => $cache_impl,
    );

    # EG invoke perl with correct path rather than relying on shebang
    my $cmd =
        qq{perl -I ./modules }
      . qq{./misc-scripts/id_mapping/dump_by_seq_region.pl }
      . qq{$options --index \$LSB_JOBINDEX};


    my $pipe =
        '|bsub '
      . $conf->param('lsf_opt_run')
      . qq{ -J '$lsf_name\[1-$num_jobs\]\%$concurrent' }
      . qq{-o $logpath/dump_by_seq_region.$dbtype.\%I.out }
      . qq{-e $logpath/dump_by_seq_region.$dbtype.\%I.err }
      . $conf->param('lsf_opt_dump_cache');

    # run lsf job array
    $logger->info("\nSubmitting $num_jobs jobs to lsf.\n");
    $logger->debug("$cmd\n\n");
    $logger->debug("$pipe\n\n");

    local *BSUB;
    open BSUB, $pipe or
      $logger->error("Could not open open pipe to bsub: $!\n");

    print BSUB $cmd;
    $logger->error("Error submitting jobs: $!\n")
      unless ($? == 0); 
    close BSUB;

    # submit dependent job to monitor finishing of jobs
    $logger->info("Waiting for jobs to finish...\n", 0, 'stamped');

    my $dependent_job =
        qq{bsub -K -w "ended($lsf_name)" }
      . $conf->param('lsf_opt_run_small')
      . qq{ -o $logpath/dump_cache.$dbtype.depend.out /bin/true};

    system($dependent_job) == 0 or
      $logger->error("Error submitting dependent job: $!\n");

    $logger->info("All jobs finished.\n", 0, 'stamped');

    # check for lsf errors
    sleep(5);
    my $err;
    foreach my $i (1..$num_jobs) {
      $err++ unless (-e "$logpath/dump_by_seq_region.$dbtype.$i.success");
    }

    if ($err) {
      $logger->error("At least one of your jobs failed.\nPlease check the logfiles at $logpath for errors.\n");
      return 1;
    }

  }

  return 0;
}


sub build_cache_all {

  # load the cache implementation
  my $cache_impl = 'Bio::EnsEMBL::IdMapping::Cache';
  inject($cache_impl);

  my $cache = $cache_impl->new(
    -LOGGER       => $logger,
    -CONF         => $conf,
  );

  # submit jobs to lsf
  foreach my $dbtype (qw(source target)) {

    $logger->info("\n".ucfirst($dbtype)." db...\n", 0, 'stamped');
    $logger->info("Building cache for whole genome...\n");

    my $i = 0;
    my $size = 0;
    ($i, $size) = $cache->build_cache_all($dbtype);

    $logger->info("Done with $dbtype (genes: $i, filesize: $size).\n", 0,
      'stamped');
  }

  return 0;
}


# EG new method for getting species IDs
sub get_species_ids {

  my ($prefix) = @_;
  my @speciesIds;
  my $dsn =
      "DBI:mysql:database="
    . $conf->param("${prefix}dbname")
    . ";host="
    . $conf->param("${prefix}host")
    . ";port="
    . $conf->param("${prefix}port");

  my $ensemblCoreDbh = DBI->connect( $dsn,
                                     $conf->param("${prefix}user"),
                                     $conf->param("${prefix}pass") )
    || die "Cannot connect to server: $DBI::errstr\n";

  my $query = "SELECT DISTINCT meta_value, species_id
               FROM meta WHERE meta_key = 'species.production_name'";

  my $psmt = $ensemblCoreDbh->prepare($query);
  $psmt->execute();

  while ( my (@results) = $psmt->fetchrow() ) {
    push @speciesIds, [ $results[0], $results[1] ];
  }
  return \@speciesIds;
} ## end sub get_species_ids
