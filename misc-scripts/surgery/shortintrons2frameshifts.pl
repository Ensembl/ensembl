use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Getopt::Long;
use POSIX qw(ceil);

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


# algorithm:
# - find genes with transcripts that have short (1-2bp) introns)
# - remember the exon positions of the exons of each transcript
# - loop through each transcript:
#      * merge exons seperated by short introns
#      * store frameshifts as attributes for transcripts
# - compare set of exons now used by updated transcripts
#      * exons which have same positions as before as left alone
#      * old exons which are no longer used are deleted
#      * new exons are given new stable ids and stored
# - update the exon_transcript table with new transcript exon composition


print STDERR "Finding short introns\n";


my $del_et_sth = $db->dbc->prepare(qq{DELETE FROM exon_transcript
                                      WHERE transcript_id = ?});

my $ins_et_sth = $db->dbc->prepare(qq{INSERT INTO exon_transcript
                                      SET exon_id = ?,
                                      transcript_id = ?,
                                      rank = ?});

my $sth = $db->dbc->prepare(qq{SELECT max(stable_id) FROM exon_stable_id});
$sth->execute();
my $ex_stable_id  = $sth->fetchall_arrayref->[0]->[0];
$ex_stable_id++;
$sth->finish();

$sth = $db->dbc()->prepare
  (qq{SELECT t.gene_id,
             MIN(IF(e1.seq_region_strand = 1,
                e2.seq_region_start - e1.seq_region_end - 1,
                e1.seq_region_start - e2.seq_region_end - 1)) AS intron_len
      FROM exon e1, exon e2, exon_transcript et1, exon_transcript et2,
           transcript t
      WHERE et1.exon_id = e1.exon_id
      AND et2.exon_id = e2.exon_id
      AND et1.transcript_id = et2.transcript_id
      AND et1.rank = et2.rank - 1
      AND et1.transcript_id = t.transcript_id
      GROUP BY t.gene_id
      HAVING intron_len < 3});

$sth->execute();

my $total_rows = $sth->rows();
my $cur_row = 0;
my $last_percent = undef;

my $ga = $db->get_GeneAdaptor();
my $ea = $db->get_ExonAdaptor();
my $aa = $db->get_AttributeAdaptor();

print STDERR "Merging exons and storing frameshifts\n";
while(my $array = $sth->fetchrow_arrayref()) {

  my $percent = ceil(($cur_row++ / $total_rows) * 100);
  if(($percent % 5) == 0 && 
     (!defined($last_percent) || $percent != $last_percent)) {
    $last_percent = $percent;
    print STDERR "$percent% complete\n";
  }

  my $g = $ga->fetch_by_dbID($array->[0]);

  my %old_exons = map {$_->hashkey() => $_} @{$g->get_all_Exons()};
  my %new_exons = ();

  foreach my $tr (@{$g->get_all_Transcripts}) {
    my @frameshifts;

    foreach my $intron (@{$tr->get_all_Introns()}) {
      push(@frameshifts, $intron) if($intron->length() < 3);
    }

    my $fs = shift(@frameshifts);

    my $merging_exon;
    my $cdna_start = 1;
    my @exons;
    my %seen_evidence;

    my @old_exons = @{$tr->get_all_Exons()};

    while(@old_exons) {
      my $ex = shift(@old_exons);
      $cdna_start += $ex->length();

      if($fs && $fs->prev_Exon->stable_id() eq $ex->stable_id()) {
        $merging_exon = undef;

        while($fs && $fs->prev_Exon->stable_id() eq $ex->stable_id()) {
          # this exon has a frameshift, so merge it with next exon
          # and any subsequent exons serperated only by frameshifts

          if(!$merging_exon) {
            $merging_exon = {};
            %{$merging_exon} = %$ex;
            bless $merging_exon, ref($ex);
            push @exons, $merging_exon;

            %seen_evidence = ();
          }

          if($merging_exon->strand() == 1) {
            $merging_exon->end($fs->next_Exon()->end());
          } else {
            $merging_exon->start($fs->next_Exon()->start());
          }

          # merge supporting evidence
          foreach my $ev (@{$ex->get_all_supporting_features()}) {
            if(!$seen_evidence{$ev->dbID()}) {
              $merging_exon->add_supporting_features($ev);
              $seen_evidence{$ev->dbID()} = 1;
            }
          }

          # store frameshift as transcript attrib
          my $seqed = Bio::EnsEMBL::SeqEdit->new
            (-CODE    => '_rna_edit',
             -NAME    => 'RNA Edit',
             -DESC    => 'Post transcriptional RNA edit',
             -START   => $cdna_start,
             -END     => $cdna_start + $fs->length() - 1,
             -ALT_SEQ => '');
          $aa->store_on_Transcript($tr, [$seqed->get_Attribute]);

          # take this frameshift and the next exon
          # (which was merged into the current exon) off the list
          $fs = shift(@frameshifts);
          $ex = shift(@old_exons);
          $cdna_start += $ex->length();
        }
      } else {
        # this exon does not have a frameshift so just add it to the list

        push @exons, $ex;
      }
    }

    $tr->flush_Exons();

    foreach my $ex (@exons) {
      $new_exons{$ex->hashkey} = $ex;
      $tr->add_Exon($ex);
    }
  }

  # determine which exons can be deleted
  foreach my $old_key (keys %old_exons) {
    if(!$new_exons{$old_key}) {
      $ea->remove($old_exons{$old_key});
      delete $old_exons{$old_key};
    }
  }

  # determine which exons are brand new and need storing
  foreach my $new_key (keys %new_exons) {
    if(!$old_exons{$new_key}) {
      my $new_ex = $new_exons{$new_key};
      # "unstore" the merged exon by unsetting adaptor and dbID
      $new_ex->dbID(undef);
      $new_ex->adaptor(undef);

      # assign a new stable id
      $new_ex->stable_id($ex_stable_id++); 
      $ea->store($new_exons{$new_key});
    }
  }

  foreach my $tr (@{$g->get_all_Transcripts()}) {
    $del_et_sth->execute($tr->dbID());

    my $rank = 1;
    foreach my $ex (@{$tr->get_all_Exons()}) {
      my $ex_id;

      if($old_exons{$ex->hashkey()}) {
        $ex_id = $old_exons{$ex->hashkey()}->dbID();
      } else {
        $ex_id = $new_exons{$ex->hashkey()}->dbID();
      }

      $ins_et_sth->execute($ex_id, $tr->dbID(), $rank);

      $rank++;
    }
  }
}

$del_et_sth->finish();
$ins_et_sth->finish();
$sth->finish();

print STDERR "\nAll Done.\n";


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
