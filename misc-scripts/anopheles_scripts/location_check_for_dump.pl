#!/usr/local/ensembl/bin/perl -w

## check of old and new protein positions, to see where old & new versions of the same transcript are in different (non-overlapping) positions - identifies 50 on same chr of which just 16 are same scaffold (identified manually - should incorporate this in script!)  Currently set so STDOUT prints just those that have jumped location within a chromosome (where chromosome includes UNKN)


use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Bio::EnsEMBL::Exon;



my $db1 = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    -host   => 'ecs1a',
					    -user   => 'ensro',
					    -dbname => 'anopheles_gambiae_core_19_2b',
					   );

print STDERR "Connecting to ecs1a, anopheles_gambiae_core_19_2b\n";

my $slice_adapt1 = $db1->get_SliceAdaptor();
my $trans_adapt1 = $db1->get_TranscriptAdaptor();


my $db2 = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    -host   => 'ecs3a',
					    -user   => 'ensro',
					    -port => '3304',
					    -dbname => 'anopheles_gambiae_core_16_2',
					   );

print STDERR "Connecting to ecs3a, anopheles_gambiae_core_16_2\n";

my $slice_adapt2 = $db2->get_SliceAdaptor();
my $trans_adapt2 = $db2->get_TranscriptAdaptor();


# for testing purposes, get all the new transcript internal ids

my $fetch_count=0;
my ($chr_mismatch,$loc_mismatch,$hurrah,$not_found);

my $sth_transcripts = $db1->prepare("select transcript_id from transcript");
$sth_transcripts->execute();
while (my $transcript_id=$sth_transcripts->fetchrow) {
  $fetch_count++;
  print STDERR "$fetch_count so far\n" if ($fetch_count % 100 == 0);
  my $transcript1=$trans_adapt1->fetch_by_dbID($transcript_id);
  my $tsi =$transcript1->stable_id;
#  print "$tsi\n";
  my $slice1=$slice_adapt1->fetch_by_transcript_id($transcript1->dbID);
  my $chr1 = $slice1->chr_name;
  my $start1= $slice1->chr_start;
  my $end1= $slice1->chr_end;
#  print "$chr1\n";

# Now get old transcript with same stable id

  if (! defined($trans_adapt2->fetch_by_stable_id($tsi))) {
    print "$tsi\tNot found in old\n";
   $not_found++;
    next;
  }
  my $transcript2=$trans_adapt2->fetch_by_stable_id($tsi);
  my $slice2=$slice_adapt2->fetch_by_transcript_stable_id($transcript2->stable_id);
  my $chr2 = $slice2->chr_name;
  my $start2= $slice2->chr_start;
  my $end2= $slice2->chr_end;

  if (($chr1) ne ($chr2)) {
    print "$tsi\tDifferent chromosomes\t$chr1  $chr2\n";
    $chr_mismatch++;
    next;
  }

  unless ( (($start1>=$start2)&&($start1<=$end2))||(($end1>=$start2)&&($end1<=$end2))||(($start1<=$start2)&&($end1>=$end2)) ) {
    print "$tsi\tSame chr but different place\t$start1..$end1  $start2..$end2\n";
    $loc_mismatch++;
    next;
  }
#  print "$tsi\t Overlap hurrah  $start1...$end1  $start2...$end2\n";
  $hurrah++;
}
print "fetched: $fetch_count, of which not found in old: $not_found\n";
print "ok: $hurrah\tloc bad: $loc_mismatch\tchr bad: $chr_mismatch\n";

