use strict;

use DBI;
use File::Basename;
use IPC::Open3;

# Use exonerate (or other program) to find xref-ensembl obejct mappings


# XXX
my $queryfile = "xref_dna.fasta";
my $targetfile = "ensembl_transcripts.fasta";

run_mapping($queryfile, $targetfile, ".");
store($targetfile);

sub run_mapping {

  my ($queryfile, $targetfile, $root_dir) = @_;

  # get list of methods
  my @methods = ("ExonerateBasic", "ExonerateBest1"); # TODO get from Ian, maybe files as well

  # foreach method, submit the appropriate job & keep track of the job name
  my @job_names;

  foreach my $method (@methods) {

    my $obj_name = "XrefMapper::Methods::$method";
    # check that the appropriate object exists
    eval "require $obj_name";
    if($@) {

      warn("Could not find object $obj_name corresponding to mapping method $method, skipping\n$@");

    } else {

      my $obj = $obj_name->new();
      my $job_name = $obj->run($queryfile, $targetfile);
      push @job_names, $job_name;
      print "Submitted LSF job $job_name to list\n";
      sleep 1; # make sure unique names really are unique

    }

  } # foreach method

  # submit depend job to wait for all mapping jobs
  submit_depend_job($root_dir, @job_names);


} # run_exonerate


sub submit_depend_job {

  my ($root_dir, @job_names) = @_;

  # Submit a job that does nothing but wait on the main jobs to
  # finish. This job is submitted interactively so the exec does not
  # return until everything is finished.

  # build up the bsub command; first part
  my @depend_bsub = ('bsub', '-K');

  # one -wended clause for each main job
  foreach my $job (@job_names) {
    push @depend_bsub, "-wended($job)";
  }

  # rest of command
  push @depend_bsub, ('-q', 'small', '-o', "$root_dir/depend.out", '-e', "$root_dir/depend.err", '/bin/true');

  #print "depend bsub:\n" . join (" ", @depend_bsub) . "\n";

  my ($depend_wtr, $depend_rtr, $depend_etr, $depend_pid);
  $depend_pid = open3($depend_wtr, $depend_rtr, $depend_etr, @depend_bsub);
  my $depend_jobid;
  while (<$depend_rtr>) {
    if (/Job <([0-9]+)> is/) {
      $depend_jobid = $1;
      print "LSF job ID for depend job: $depend_jobid \n" ;
    }
  }
  if (!defined($depend_jobid)) {
    print STDERR "Error: could not get depend job ID\n";
  }



}

=head2 store

  Arg[1]     : The target file used in the exonerate run. Used to work out the Ensembl object type.
  Arg[2]     : 
  Example    : none
  Description: Parse exonerate output files and build files for loading into target db tables.
  Returntype : List of strings
  Exceptions : none
  Caller     : general

=cut

sub store {

  my ($target_file_name) = @_;

  my $type = get_ensembl_object_type($target_file_name);

  # get or create the appropriate analysis ID
  my $analysis_id = get_analysis_id($type);

  # TODO - get this from config
  my $dbi = DBI->connect("dbi:mysql:host=ecs1g;port=3306;database=arne_core_20_34",
			 "ensadmin",
			 "ensembl",
			 {'RaiseError' => 1}) || die "Can't connect to database";

  # get current highest internal ID from xref and object_xref
  my $max_xref_id = 0;
  my $sth = $dbi->prepare("SELECT MAX(xref_id) FROM xref");
  $sth->execute();
  my $max_xref_id = ($sth->fetchrow_array())[0];
  if (!defined $max_xref_id) {
    print "Can't get highest existing xref_id, using 0\n)";
  } else {
    print "Maximum existing xref_id = $max_xref_id\n";
  }
  my $max_object_xref_id = 0;
  $sth = $dbi->prepare("SELECT MAX(object_xref_id) FROM object_xref");
  $sth->execute();
  my $max_object_xref_id = ($sth->fetchrow_array())[0];
  if (!defined $max_object_xref_id) {
    print "Can't get highest existing object_xref_id, using 1\n)";
  } else {
    print "Maximum existing object_xref_id = $max_object_xref_id\n";
  }


  #my $ox_sth = $dbi->prepare("INSERT INTO object_xref(ensembl_id, ensembl_object_type, xref_id) VALUES(?,?,?)");

  #my $ix_sth = $dbi->prepare("INSERT INTO identity_xref VALUES(?,?,?,?,?,?,?,?,?,?,?)");

  # files to write table data to
  open (OBJECT_XREF, ">object_xref.txt");
  open (IDENTITY_XREF, ">identity_xref.txt");

  my $total_lines = 0;
  my $total_files = 0;

  my $xref_id = $max_xref_id + 1;
  my $object_xref_id = $max_object_xref_id + 1;

  # TODO - store xrefs in a file as well

  foreach my $file (glob("*.map")) {

    print "Parsing results from $file \n";
    open(FILE, $file);
    $total_files++;

    while (<FILE>) {

      $total_lines++;
      chomp();
      my ($label, $query_id, $target_id, $query_start, $query_end, $target_start, $target_end, $cigar_line, $score) = split(/:/, $_);
      $cigar_line =~ s/ //;

      # TODO make sure query & target are the right way around

      print OBJECT_XREF "$object_xref_id\t$target_id\t$type\t$query_id\n";
      print IDENTITY_XREF "$object_xref_id\t$query_id\t$target_id\t$query_start\t$query_end\t$target_start\t$target_end\t$cigar_line\t$score\t\\N\t$analysis_id\n";
      # TODO - evalue?
      $object_xref_id++;

      # Store in database
      # create entry in object_xref and get its object_xref_id
      #$ox_sth->execute($target_id, $type, $query_id) || warn "Error writing to object_xref table";
      #my $object_xref_id = $ox_sth->{'mysql_insertid'};

      # create entry in identity_xref
      #$ix_sth->execute($object_xref_id, $query_id, $target_id, $query_start, $query_end, $target_start, $target_end, $cigar_line, $score, undef, $analysis_id) || warn "Error writing to identity_xref table";

    }

    close(FILE);

  }

  close(IDENTITY_XREF);
  close(OBJECT_XREF);

  print "Read $total_lines lines from $total_files exonerate output files\n";

}


sub get_ensembl_object_type {

  my $filename = shift;
  my $type;

  if ($filename =~ /gene/i) {

    $type = "Gene";

  } elsif ($filename =~ /transcript/i) {

    $type = "Transcript";

  } elsif ($filename =~ /translation/i) {

    $type = "Translation";

  } else {

    print STDERR "Cannot deduce Ensembl object type from filename $filename";
  }

  return $type;

}


sub get_analysis_id {

  my $ensembl_type = shift;

  my %typeToLogicName = ( 'transcript' => 'XrefExonerateDNA',
			  'translation' => 'XrefExonerateProtein' );

  my $logic_name = $typeToLogicName{lc($ensembl_type)};

  # TODO - get these details from Config
  my $host = "ecs1g";
  my $port = 3306;
  my $database = "arne_core_20_34";
  my $user = "ensadmin";
  my $password = "ensembl";

  my $dbi = DBI->connect("dbi:mysql:host=$host;port=$port;database=$database",
			 "$user",
			 "$password",
			 {'RaiseError' => 1}) || die "Can't connect to database";


  my $sth = $dbi->prepare("SELECT analysis_id FROM analysis WHERE logic_name='" . $logic_name ."'");
  $sth->execute();

  my $analysis_id;

  if (my @row = $sth->fetchrow_array()) {

    $analysis_id = $row[0];
    print "Found exising analysis ID ($analysis_id) for $logic_name\n";

  } else {

    print "No analysis with logic_name $logic_name found, creating ...\n";
    $sth = $dbi->prepare("INSERT INTO analysis (logic_name, created) VALUES ('" . $logic_name. "', NOW())");
    # TODO - other fields in analysis table
    $sth->execute();
    $analysis_id = $sth->{'mysql_insertid'};
    print "Done (analysis ID=" . $analysis_id. ")\n";

  }

  return $analysis_id;

}
