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

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Getopt::Long;
use POSIX qw(ceil);

my ($host, $port, $user, $pass, $dbname, $testid);

GetOptions('host=s' => \$host,
           'user=s' => \$user,
           'port=i' => \$port,
           'pass=s' => \$pass,
           'dbname=s' => \$dbname,
           'testid=i' => \$testid);

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

my $upd_tl_sth = $db->dbc->prepare(qq{UPDATE translation
                                      SET start_exon_id = ?,
                                          end_exon_id   = ?,
                                          seq_start    = ?,
                                          seq_end      = ?
                                      WHERE translation_id = ?});


my $sth = $db->dbc->prepare(qq{SELECT max(stable_id) FROM exon});
$sth->execute();
my $ex_stable_id  = $sth->fetchall_arrayref->[0]->[0];
$ex_stable_id++;
$sth->finish();

if(!$testid) {

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

}

my @gene_ids =
  ($testid) ? ($testid) : map{$_->[0]} @{$sth->fetchall_arrayref()};

my $total_rows = ($testid) ?  1 : $sth->rows();
my $cur_row = 0;
my $last_percent = undef;

my $ga = $db->get_GeneAdaptor();
my $ea = $db->get_ExonAdaptor();
my $aa = $db->get_AttributeAdaptor();

print STDERR "Merging exons and storing frameshifts\n";
foreach my $gene_id (@gene_ids) {

  my $percent = ceil(($cur_row++ / $total_rows) * 100);
  if(($percent % 5) == 0 && 
     (!defined($last_percent) || $percent != $last_percent)) {
    $last_percent = $percent;
    print STDERR "$percent% complete\n";
  }

  my $g = $ga->fetch_by_dbID($gene_id);

  my %old_exons = map {$_->hashkey() => $_} @{$g->get_all_Exons()};
  my %new_exons = ();

  my %tl_cdna_starts;
  my %tl_cdna_ends;

  foreach my $tr (@{$g->get_all_Transcripts}) {
    # keep track of translation cdna coordinates before messing with transcript
    my ($cdna_coding_start, $cdna_coding_end);
    if($tr->translation()) {
      $cdna_coding_start = $tr->cdna_coding_start();
      $cdna_coding_end   = $tr->cdna_coding_end();
    }

    my @frameshifts;

    # collect list of frameshifts
    foreach my $intron (@{$tr->get_all_Introns()}) {
      push(@frameshifts, $intron) if($intron->length() < 3);
    }

    my %seen_evidence;
    my @exons;
    my $fs = shift(@frameshifts);
    my $cdna_start = 1;

    # merge together exons which are seperated by frameshift introns
    foreach my $ex (@{$tr->get_all_Exons()}) {

      # was the previous exon seperated from this one by a frameshift intron?
      if($fs && $fs->next_Exon()->stable_id() eq $ex->stable_id()) {

        # merge this exon with the previous one
        if($ex->strand() == 1) {
          $exons[$#exons]->end($ex->end());
        } else {
          $exons[$#exons]->start($ex->start());
        }

        # merge supporting evidence
        foreach my $sf (@{$ex->get_all_supporting_features()}) {
          if(!$seen_evidence{ref($sf) . ':' . $sf->dbID()}) {
            $ex->add_supporting_feature($sf);
            $seen_evidence{ref($sf . ':' . $sf->dbID())} = 1;
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
        $aa->store_on_Transcript($tr->dbID, [$seqed->get_Attribute]);

        # adjust cdna coordinates for frameshifted basepairs
        if($tr->translation) {
          if($cdna_coding_start >= $cdna_start) {
            $cdna_coding_start += $fs->length();
          }
          if($cdna_coding_end >= $cdna_start) {
            $cdna_coding_end += $fs->length();
          }
        }
        $cdna_start += $fs->length();

        # look at the next frameshift
        $fs = shift(@frameshifts);
      } else {
        push @exons, $ex;

        %seen_evidence = ();
        foreach my $sf (@{$ex->get_all_supporting_features()}) {
          $seen_evidence{ref($_) . ':' . $_->dbID()} = 1;
        }
      }

      $cdna_start += $ex->length();
    }

    # rebuild transcript from new set of exons
    $tr->flush_Exons();

    foreach my $ex (@exons) {
      $new_exons{$ex->hashkey} = $ex;
      $tr->add_Exon($ex);
    }

    $tl_cdna_starts{$tr->dbID()} = $cdna_coding_start;
    $tl_cdna_ends{$tr->dbID()} = $cdna_coding_end;
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

  # update the transcript composition by updating the exon transcript table
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

  # update the translations to use the new exons
  foreach my $tr (@{$g->get_all_Transcripts()}) {
    my $tl = $tr->translation();

    next if(!$tl);

    my $tl_start = $tl_cdna_starts{$tr->dbID()};
    my $tl_end   = $tl_cdna_ends{$tr->dbID()};

    foreach my $ex (@{$tr->get_all_Exons()}) {
      if($tl_start > 0 && $tl_start <= $ex->length()) {
        $tl->start_Exon($ex);
        $tl->start($tl_start);
      }

      if($tl_end > 0 && $tl_end <= $ex->length()) {
        $tl->end_Exon($ex);
        $tl->end($tl_end);
      }

      $tl_start -= $ex->length();
      $tl_end   -= $ex->length();
    }

    my ($start_ex_id, $end_ex_id);

    # use consolidated exon ids. Exons may have been duplicated across
    # multiple transcripts but only one has been stored and had id set.
    if($old_exons{$tl->start_Exon->hashkey()}) {
      $start_ex_id = $old_exons{$tl->start_Exon->hashkey()}->dbID();
    } else {
      $start_ex_id = $new_exons{$tl->start_Exon->hashkey()}->dbID();
    }

    if($old_exons{$tl->end_Exon->hashkey()}) {
      $end_ex_id = $old_exons{$tl->end_Exon->hashkey()}->dbID();
    } else {
      $end_ex_id = $new_exons{$tl->end_Exon->hashkey()}->dbID();
    }


    $upd_tl_sth->execute($start_ex_id, $end_ex_id,
                         $tl->start(), $tl->end(), $tl->dbID());
  }
}

$del_et_sth->finish();
$ins_et_sth->finish();
$sth->finish();
$upd_tl_sth->finish();

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
