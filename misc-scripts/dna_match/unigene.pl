#!/usr/local/ensembl/bin/perl -w

use strict;
use constant QUEUENAME => "acari";	# queue name for bsub
use constant CHUNK_SIZE => 40;		# transcripts per farm job
use Bio::SeqIO;
use Bio::EnsEMBL::Pipeline::GeneComparison::CdnaComparison;

# Perform ENSEMBL transcript stable ID to Unigene cluster ID mapping.
# This script should be launched on a head node. It will distribute
# the tasks over the farm. It does not, however, wait for them to
# complete: the user must do that manually by watching the LSF queues.

# Usage: mapper.pl <ensembl_transcript_file_name.fa>

# Implementation details. This same script is used to control the
# analysis on the head node and perform sub-analyses on the farm. It
# decides which to do according to the first command-line argument.
# The master image creates a temporary, fairly small transcript file
# for each slave, and launches the slaves through LSF. Each slave
# runs a mapping for the transcripts in the file it is given.

if ($ARGV[0] eq "SLAVE") {	# slave: do some work

  die "File name argument required" if @ARGV < 2;
  
  # read the transcripts we are going to analyse into memory
  my $transcript_fnam = $ARGV[1];
  print STDERR "$0 slave: analysing $transcript_fnam\n";
  my $in  = Bio::SeqIO->new(-file     => $transcript_fnam,
                            '-format' => 'Fasta');
  my @transcript_seqs;
  while (my $seq = $in->next_seq) {
    push @transcript_seqs, $seq;
  }

  # run comparison
  my $cdna_comp = Bio::EnsEMBL::Pipeline::GeneComparison::CdnaComparison->new(
                                 -TRANSCRIPT_ARRAYREF => \@transcript_seqs,
		                 -CDNA_DATABASE       => 'unigene.seq');
  $cdna_comp->run_transcript_mapping;
  my %name_map = $cdna_comp->get_transcript_mapping;

  # output to stdout
  foreach my $transcript_id (keys %name_map) {
    print "$transcript_id";
    foreach my $mapped_name (@{$name_map{$transcript_id}}) {
      print "\t$mapped_name";
    }
    print "\n";
  }

} else {			# master: prepare for others to do the work

  die "File name argument required" if @ARGV < 1;

  # read all transcripts into memory
  my $full_transcript_file_name = $ARGV[0];
  my $in = Bio::SeqIO->new(-file     => $full_transcript_file_name,
                           '-format' => 'Fasta');
  my @transcript_seqs;
  while (my $seq = $in->next_seq) {
    push @transcript_seqs, $seq;
  }

  # prepare and launch small-scale analyses

  my @bjob_nos;
  my $seq_no = 0;
  OUTER:
  while ($seq_no < @transcript_seqs) {

    # create a file containing a small number of transcripts
    my $chunk_fnam = "mapper_pl_tmp_" . $$ . "_" . ($seq_no + 1);
    open (MAPPER_PL_OUT_FH, ">$chunk_fnam") or die "file open error";
    my $out = new Bio::SeqIO(-fh => \*MAPPER_PL_OUT_FH, '-format' => 'Fasta');
    my $lower = $seq_no;		# first transcript to output
    my $upper = $seq_no + CHUNK_SIZE;	# 1 beyond last transcript to output
    $upper = @transcript_seqs if $upper > @transcript_seqs;
    for (my $i = $lower; $i < $upper; $i++) { 
      $out->write_seq($transcript_seqs[$i]);
    }
    $seq_no = $upper;			# write this one on next iteration
    close MAPPER_PL_OUT_FH or die 'file close error';

    # LSF launch: first 500 all at once, then with 1 sec delay
    my $command = "bsub -q " . QUEUENAME . " -C0 "
      . "\"-f $chunk_fnam > /tmp/$chunk_fnam\" "
      . "-o ./$chunk_fnam" . "_mapping -e /dev/null "
      . "$0 SLAVE /tmp/$chunk_fnam";
    print STDERR "$0 master: about to launch: $command\n";
    my $bsub_output=`$command`;
    sleep 1
      unless $seq_no < 500;
    
  }

}
