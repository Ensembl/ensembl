use strict;

use DBI;
use File::Basename;

# Use exonerate (or other program) to find xref-ensembl obejct mappings

# XXX
run_matching();

sub run_matching {

   run_exonerate("queryfile.fasta", "targetfile.fasta");

}

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

  # TODO - calculate number of jobs from file sizes
  my $num_jobs = 20;

  # TODO - get root_dir from config
  my $root_dir = ".";

  # TODO - get exonerate executable from config
  my $exonerate_path = "/usr/local/ensembl/bin/exonerate-0.8.3";

  # DNA xref - transcript mapping

  my $dna_exonerate_options = join(" ", get_dna_exonerate_options());

  my $unique_name = "mapXrefDNA" . time();

  # TODO - prefix = working dir + basename of query file
  my $prefix = basename($query);
  $prefix =~ s/\.\w+$//;

  my $output = "\$LSB_JOBINDEX.map";

  my @dna_bsub = ( 'bsub', '-J', $unique_name . "[1-$num_jobs]", '-o', "$prefix.%I.out", '-e', "$prefix.%I.err");

  print "bsub command: " . join(" ", @dna_bsub) . "\n\n";

  # Open a pipe to the stdin of a bsub command with the appropriate options
  #open( BSUB, '|-' ) or exec (@dna_bsub);

  # Create actual execute script to be executed with LSF, and write to pipe
  my $dna_job = <<EOF;
. /usr/local/lsf/conf/profile.lsf

cd /tmp

rm -f /tmp/\$LSB_JOBINDEX.$query /tmp/\$LSB_JOBINDEX.$target /tmp/$output

lsrcp ecs1a:$root_dir/$target /tmp/\$LSB_JOBINDEX.$target
lsrcp ecs1a:$root_dir/$query  /tmp/\$LSB_JOBINDEX.$query

$exonerate_path /tmp/\$LSB_JOBINDEX.$query /tmp/\$LSB_JOBINDEX.$target --querychunkid \$LSB_JOBINDEX --querychunktotal $num_jobs --model affine:local -M 900 --showalignment FALSE --subopt no --ryo "xrefdna %qi %ti %qab %qae %tab %tae %C %s\n" $dna_exonerate_options | grep '^xrefdna ' > /tmp/$output

lsrcp /tmp/$output ecs1a:$root_dir/$output

rm -f /tmp/\$LSB_JOBINDEX.$query /tmp/\$LSB_JOBINDEX.$target /tmp/$output
EOF

# TODO make sure query/target are the right way round
# TODO analysis ID

  #print BSUB $dna_job;

  print "job:\n" . $dna_job . "\n\n";

  #close(BSUB);

  # TODO same for protein

  # The above calls to bsub return immediately. We now submit a job that does
  # nothing but wait on them to finish. This job is submitted interactively
  # so the exec does not return until everything is finished.
  my @depend_bsub = ('bsub', '-K', "-wended($unique_name)", '-q', 'small', '-o', "$root_dir/depend.out", '-e', "$root_dir/depend.err", '/bin/true');

  print "depend bsub:\n" . join (" ", @depend_bsub) . "\n";

  # TODO - system or exec this

} # run_exonerate


sub parse_and_store {

  # files to write table data to
  open (OBJECT_XREF, ">object_xref.txt");
  open (IDENTITY_XREF, ">identity_xref.txt");

  my $total_lines = 0;
  foreach my $file (glob <*.map>) {

    # XXX figure out type of ensembl objects we're dealing with
    $type = 'Transcript';
    print "Parsing results from $file \n";
    open(FILE, file);

    while (<FILE>) {

      $total_lines++;
      my ($label, $query_id, $target_id, $query_start, $query_end, $target_start, $target_end, $cigar_line, $score) = split();
      # TODO make sure query & target are the right way around

      print OBJECT_XREF "$target_id $type $query_id\n"; # XXX object_xref_id
      print IDENTITY_XREF "$query_id $target_id $query_start $query_end $target_start $target_end $cigar_line $score"; # XXX object_xref_id

    }

    close(FILE);

  }

  close(IDENTITY_XREF);
  close(OBJECT_XREF);

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

   return('-q', 'normal');

}
