use strict;

use DBI;
use File::Basename;
use IPC::Open3;

# Use exonerate (or other program) to find xref-ensembl obejct mappings

# XXX
my $queryfile = "xref_dna.fasta";
my $targetfile = "ensembl_transcripts.fasta";

#run_exonerate($queryfile, $targetfile);
parse_and_store($targetfile);


=head2 run_exonerate

  Arg[1]     : Query filename to pass to exonerate; this should be the XREF file.
  Arg[2]     : Target filename to pass to exonerate; this should be the ENSEMBL file.
  Example    : none
  Description: submit mapping jobs to LSF via bsub.
  Returntype : List of strings
  Exceptions : none
  Caller     : general

=cut

sub run_exonerate {

  my ($query, $target) = @_;

  my $num_jobs = calculate_num_jobs($query);

  # TODO - get root_dir from config
  my $root_dir = "/nfs/acari/gp1/work/ensembl/misc-scripts/xref_mapping";

  # TODO - get exonerate executable from config
  my $exonerate_path = "/usr/local/ensembl/bin/exonerate-0.8.3";

  # DNA xref - transcript mapping

  my $dna_exonerate_options = join(" ", get_dna_exonerate_options());

  my $unique_name = "mapXrefDNA" . time();

  my $prefix = $root_dir . "/" . basename($query);
  $prefix =~ s/\.\w+$//;

  my $output = "\$LSB_JOBINDEX.map";

  my @dna_bsub = ( 'bsub', '-J', $unique_name . "[1-$num_jobs]", '-o', "$prefix.%I.out", '-e', "$prefix.%I.err");

  #print "bsub command: " . join(" ", @dna_bsub) . "\n\n";

  # Open a pipe to the stdin of a bsub command with the appropriate options
  #open( BSUB, '|-' ) or exec (@dna_bsub);

  # Use IPC::Open3 to open the process so that we can read and write from/to its stdout/stderr
  my ($wtr, $rtr, $etr, $pid);
  $pid = open3($wtr, $rtr, $etr, @dna_bsub);

  # Create actual execute script to be executed with LSF, and write to pipe
  my $dna_job = <<EOF;
. /usr/local/lsf/conf/profile.lsf

cd /tmp

rm -f /tmp/\$LSB_JOBINDEX.$query /tmp/\$LSB_JOBINDEX.$target /tmp/$output

lsrcp ecs1a:$root_dir/$target /tmp/\$LSB_JOBINDEX.$target
lsrcp ecs1a:$root_dir/$query  /tmp/\$LSB_JOBINDEX.$query

$exonerate_path /tmp/\$LSB_JOBINDEX.$query /tmp/\$LSB_JOBINDEX.$target --querychunkid \$LSB_JOBINDEX --querychunktotal $num_jobs --showvulgar false --showalignment FALSE --ryo "xrefdna %qi %ti %qab %qae %tab %tae %C %s\n" $dna_exonerate_options | grep '^xrefdna ' > /tmp/$output

lsrcp /tmp/$output ecs1a:$root_dir/$output

rm -f /tmp/\$LSB_JOBINDEX.$query /tmp/\$LSB_JOBINDEX.$target /tmp/$output
EOF

# TODO make sure query/target are the right way round
# TODO analysis ID

  print $wtr $dna_job;

  #print "Job:\n" . $dna_job . "\n\n";

  close($wtr);

  # Wait until bsub has actually run - will print to its stdout ($rtr) and then close it
  my $jobid;
  while (<$rtr>) {
    if (/Job <([0-9]+)> is/) {
      $jobid = $1;
      print "LSF job ID for main mapping job: $jobid \n"
    }
  }
  if (!defined($jobid)) {
    print STDERR "Error: could not get LSF job ID for mapping job\n";
  }

  # TODO same for protein

  # The above calls to bsub return immediately. We now submit a job that does
  # nothing but wait on them to finish. This job is submitted interactively
  # so the exec does not return until everything is finished.
  my @depend_bsub = ('bsub', '-K', "-wended($unique_name)", '-q', 'small', '-o', "$root_dir/depend.out", '-e', "$root_dir/depend.err", '/bin/true');

  #print "depend bsub:\n" . join (" ", @depend_bsub) . "\n";

  my ($depend_wtr, $depend_rtr, $depend_etr, $depend_pid);
  $depend_pid = open3($depend_wtr, $depend_rtr, $depend_etr, @depend_bsub);
  my $depend_jobid;
  while (<$depend_rtr>) {
    if (/Job <([0-9]+)> is/) {
      $jobid = $1;
      print "LSF job ID for depend job: $jobid \n" ;
    }
  }
  if (!defined($depend_jobid)) {
    print STDERR "Error: could not get depend job ID\n";
  }


} # run_exonerate


sub parse_and_store {

  my $target_file_name = shift;
  my $type = get_ensembl_object_type($target_file_name);

  # files to write table data to
  open (OBJECT_XREF, ">object_xref.txt");
  open (IDENTITY_XREF, ">identity_xref.txt");

  my $total_lines = 0;
  my $total_files = 0;

  foreach my $file (glob("*.map")) {

    print "Parsing results from $file \n";
    open(FILE, $file);
    $total_files++;

    while (<FILE>) {

      $total_lines++;
      my ($label, $query_id, $target_id, $query_start, $query_end, $target_start, $target_end, $cigar_line, $score) = split();
      # TODO make sure query & target are the right way around

      print OBJECT_XREF "$target_id $type $query_id\n"; # XXX object_xref_id
      print IDENTITY_XREF "$query_id $target_id $query_start $query_end $target_start $target_end $cigar_line $score\n"; # XXX object_xref_id

    }

    close(FILE);

  }

  close(IDENTITY_XREF);
  close(OBJECT_XREF);

  print "Read $total_lines lines from $total_files exonerate output files\n";

}

=head2 get_protein_exonerate_options

  Args       : none
  Example    : none
  Description: Return additional options to pass to exonerate when running protein mapping
  Returntype : List of strings
  Exceptions : none
  Caller     : general

=cut

sub get_protein_exonerate_options {

   return();

}

=head2 get_dna_exonerate_options

  Args       : none
  Example    : none
  Description: Return additional options to pass to exonerate when running DNA mapping
  Returntype : List of strings
  Exceptions : none
  Caller     : general

=cut

sub get_dna_exonerate_options {

   return("--bestn", "1");

}


sub calculate_num_jobs {

  my $query = shift;

  my $bytes_per_job = 250000;

  my $size = (stat $query)[7];

  return int($size/$bytes_per_job);

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
