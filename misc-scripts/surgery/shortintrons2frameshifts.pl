use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Getopt::Long;


my ($host, $port, $user, $pass, $dbname);

GetOptions('host=s' => \$host,
           'user=s' => \$user,
           'port=i' => \$port,
           'pass=s' => \$pass,
           'dbname=s' => \$dbname);

$port ||= 3306;

usage("-user, -host, -dbname args are required") 
  if(!$user || !$dbname || !$host);

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new
  (-host => $host,
   -user => $user,
   -pass => $pass,
   -dbname => $dbname,
   -port => $port);


print STDERR "Finding short introns\n";

my $find_introns_sth = $db->dbc()->prepare
  (qq{SELECT e1.exon_id as exon1_id, e2.exon_id as exon2_id,
             IF(e1.seq_region_strand = 1,
                e2.seq_region_start - e1.seq_region_end - 1,
                e1.seq_region_start - e2.seq_region_end - 1) AS intron_len
      FROM exon e1, exon e2, exon_transcript et1, exon_transcript et2
      WHERE et1.exon_id = e1.exon_id
      AND et2.exon_id = e2.exon_id
      AND et1.transcript_id = et2.transcript_id
      AND et1.rank = et2.rank - 1
      GROUP BY exon1_id, exon2_id
      HAVING intron_len < 3});

$find_introns_sth->execute();

my %DELETED_EXONS;
my %MULT_FRAMESHIFTS;

my $BAD_COUNT = 0;

print STDERR "Processing frameshifts";
process_frameshifts($find_introns_sth);

$find_introns_sth->finish();


# keep reprocessing the transcripts that had multiple frameshifts
# until there are no transcripts with frameshifts left

while(keys %MULT_FRAMESHIFTS) {
  my $id_list = join(',',keys %MULT_FRAMESHIFTS);

  print STDERR "\nReprocessing transcripts which had multiple frameshifts\n";

  %DELETED_EXONS    = ();
  %MULT_FRAMESHIFTS = ();

  $find_introns_sth = $db->dbc()->prepare
    (qq{SELECT e1.exon_id as exon1_id, e2.exon_id as exon2_id,
             IF(e1.seq_region_strand = 1,
                e2.seq_region_start - e1.seq_region_end - 1,
                e1.seq_region_start - e2.seq_region_end - 1) AS intron_len
      FROM exon e1, exon e2, exon_transcript et1, exon_transcript et2
      WHERE et1.exon_id = e1.exon_id
      AND et2.exon_id = e2.exon_id
      AND et1.transcript_id = et2.transcript_id
      AND et1.rank = et2.rank - 1
      AND et2.exon_id in ($id_list)
      GROUP BY exon1_id, exon2_id
      HAVING intron_len < 3});

  $find_introns_sth->execute();

  process_frameshifts($find_introns_sth);

  $find_introns_sth->finish();
}


print STDERR "\nFound $BAD_COUNT bad introns\n";
print STDERR "\ndone\n";



#
# takes an executed statement handle and processes the frameshifts
# (short introns) which are the results
#
sub process_frameshifts {
  my $sth = shift;

  my ($ex1_id, $ex2_id, $intron_len);
  $find_introns_sth->bind_columns(\$ex1_id, \$ex2_id, \$intron_len);

  while($find_introns_sth->fetch()) {
    if($intron_len < 1) {
      die("Unexpected: intron_len less than 1 between exons " .
          "$ex1_id and $ex2_id");
    }

    print STDERR ".";

    if($DELETED_EXONS{$ex1_id}) {
      # multiple frameshift transcript, first exon already merged
      # with another...
      $MULT_FRAMESHIFTS{$ex2_id} = 1;
      next;
    }

    next if(!check_introns($db, $ex1_id, $ex2_id));
    add_rna_edits($db, $ex1_id, $ex2_id, $intron_len);
    merge_exon($db, $ex1_id, $ex2_id);
  }
}


#
# stores frameshifts as TranscriptAttribs in the database
#
sub add_rna_edits {
  my $db = shift;
  my $ex1_id = shift;
  my $ex2_id = shift;
  my $intron_len = shift;

  my $sth = $db->dbc->prepare(qq{SELECT et1.transcript_id
                                 FROM  exon_transcript et1, exon_transcript et2
                                 WHERE et1.transcript_id = et2.transcript_id
                                 AND   et1.exon_id = ?
                                 AND   et2.exon_id = ?});

  $sth->execute($ex1_id, $ex2_id);

  my @tr_ids = map {$_->[0]} @{$sth->fetchall_arrayref()};

  my $tra = $db->get_TranscriptAdaptor();
  my $aa = $db->get_AttributeAdaptor();

  foreach my $tr_id (@tr_ids) {
    my $tr = $tra->fetch_by_dbID($tr_id,1);

    my $exons = $tr->get_all_Exons();

    my $cdna_start = 1;

    my $prev_was_ex1 = 0;
    foreach my $ex (@$exons) {

      if($ex->dbID() == $ex2_id) {
        if(!$prev_was_ex1) {
          die("Unexpected: exon1 $ex1_id was not exon before exon2 $ex2_id\n");
        }

        my $seqed = Bio::EnsEMBL::SeqEdit->new
          (-CODE    => '_rna_edit',
           -NAME    => 'RNA Edit',
           -DESC    => 'Post transcriptional RNA edit',
           -START   => $cdna_start,
           -END     => $cdna_start + $intron_len - 1,
           -ALT_SEQ => '');

        $aa->store_on_Transcript($tr, [$seqed->get_Attribute]);

        last;
      }

      $cdna_start += $ex->length();
      $prev_was_ex1 = ($ex->dbID() == $ex1_id) ? 1 : 0;
    }
  }

  return;
}


#
# makes sure that all transcripts that have one of the exons have
# both of the exons that are on either side of the short intron
#
sub check_introns {
  my $db = shift;
  my $ex1_id = shift;
  my $ex2_id = shift;

  my $sth = $db->dbc()->prepare
    (qq{SELECT count(*)
        FROM   exon_transcript et
        WHERE  et.exon_id = ?});

  $sth->execute($ex1_id);
  my $ex1_count = $sth->fetchall_arrayref()->[0]->[0];
  $sth->execute($ex2_id);
  my $ex2_count = $sth->fetchall_arrayref()->[0]->[0];

  $sth->finish();

  $sth = $db->dbc()->prepare
    (qq{SELECT count(*)
        FROM exon_transcript et1, exon_transcript et2
        WHERE et1.transcript_id = et2.transcript_id
        AND   et1.exon_id = ?
        AND   et2.exon_id = ?});

  $sth->execute($ex1_id, $ex2_id);

  my $both_count = $sth->fetchall_arrayref()->[0]->[0];

  if($ex1_count != $ex2_count || $ex1_count != $both_count) {
    warn("Exons $ex1_id and $ex2_id define a small intron but are not " .
        "both shared by all transcripts. Skipping transcript.");
    $BAD_COUNT++;
    return 0;
  }


  return 1;
}


#
# Removes the short intron from the database by enlarging the first exon
# and deleting the second exon.
#
sub merge_exon {
  my $db = shift;
  my $ex1_id = shift;
  my $ex2_id = shift;

  my $exa = $db->get_ExonAdaptor();

  my $ex1 = $exa->fetch_by_dbID($ex1_id);
  my $ex2 = $exa->fetch_by_dbID($ex2_id);


  # update the size of the first exon

  my $new_start = ($ex1->strand() == 1) ? $ex1->start() : $ex2->start();
  my $new_end   = ($ex1->strand() == 1) ? $ex2->end()   : $ex2->end();

  my $new_end_phase = $ex2->end_phase();

  my $sth = $db->dbc->prepare(qq{UPDATE exon
                                 SET seq_region_start = ?,
                                     seq_region_end = ?,
                                     end_phase = ?
                                 WHERE exon_id = ?});

  $sth->execute($new_start, $new_end, $new_end_phase, $ex1->dbID());
  $sth->finish();

  # delete the second exon

  $DELETED_EXONS{$ex2->dbID()} = 1;

  $sth = $db->dbc->prepare(qq{DELETE FROM exon WHERE exon_id = ?});
  $sth->execute($ex2->dbID());
  $sth->finish();

  ### should supporting evidence be deleted or merged???
  $sth = $db->dbc->prepare(qq{DELETE FROM supporting_feature
                              WHERE exon_id = ?});
  $sth->execute($ex2->dbID());
  $sth->finish();

  $sth = $db->dbc->prepare(qq{DELETE FROM exon_stable_id
                              WHERE exon_id = ?});
  $sth->execute($ex2->dbID());
  $sth->finish();

  # update the rank of other exons in these transcripts

  my $update_rank_sth = $db->dbc->prepare(qq{UPDATE exon_transcript
                                             SET rank = rank-1
                                             WHERE exon_id = ?
                                             AND   transcript_id = ?});


  $sth = $db->dbc->prepare(qq{SELECT et1.exon_id, et1.transcript_id
                              FROM exon_transcript et1, exon_transcript et2
                              WHERE et1.transcript_id = et2.transcript_id
                              AND et2.exon_id = ?
                              AND et1.rank > et2.rank});

  $sth->execute($ex2->dbID());

  while(my @array = $sth->fetchrow_array()) {
    $update_rank_sth->execute(@array);
  }

  $update_rank_sth->finish();
  $sth->finish();


  $sth = $db->dbc->prepare(qq{DELETE FROM exon_transcript
                              WHERE exon_id = ?});

  $sth->execute($ex2->dbID());

  $sth->finish();

  return;
}




sub usage {
  print STDERR qq {
This program merges exons which are seperated by introns of length 1 or 2
and replaces them with a frameshift mrna edit.

Usage:
perl shortintrons2frameshift.pl -user <user> -host <host> [-pass <password>] \
                                -dbname <dbname> [-port <port>]
};

  exit;
}
