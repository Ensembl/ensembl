=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::VegaCuration::Transcript;

use strict;
use warnings;
no warnings 'uninitialized';
use vars qw(@ISA);

use Bio::EnsEMBL::Utils::VegaCuration::Gene;
use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Utils::VegaCuration::Gene);


=head2 find_non_overlaps

   Args       : arrayref of B::E::Transcripts
   Example    : find_non_overlaps($all_transcripts)
   Description: identifies any non-overlapping transcripts
   Returntype : array refs of stable IDs
   Exceptions : none

=cut

sub find_non_overlaps {
  my $self = shift;
  my ($all_transcripts) = @_;
  my $non_overlaps = [];
  foreach my $transcript1 (@{$all_transcripts}) {
    foreach my $transcript2 (@{$all_transcripts}) {
      if ($transcript1->end < $transcript2->start) {
	push @{$non_overlaps}, $transcript1->stable_id;
	push @{$non_overlaps}, $transcript2->stable_id;
      }
    }
  }
  return $non_overlaps;
}

=head2 check_remarks_and_update_names

   Arg[1]     : B::E::Gene (with potentially duplicated transcript names)
   Arg[2]     : counter 1 (no. of patched genes)
   Arg[3]     : counter 2 (no. of patched transcripts)
   Example    : $support->update_names($gene,\$c1,\$c2)
   Description: - checks remarks and patches transcripts with identical names according to
                CDS and length
   Returntype : true | false (depending on whether patched or not), counter1, counter2

=cut

sub check_remarks_and_update_names {
  my $self = shift;
  my ($gene,$gene_c,$trans_c) = @_;
  my $action = ($self->param('dry_run')) ? 'Would add' : 'Added';
  my $aa  = $gene->adaptor->db->get_AttributeAdaptor;
  my $dbh = $gene->adaptor->db->dbc->db_handle;

  #get list of IDs that have previously been sent to annotators
  my $seen_genes = $self->get_havana_fragmented_loci_comments;

  my $gsi    = $gene->stable_id;
  my $gid    = $gene->dbID;
  my $g_name;
  my $study_more = 1;
  eval {
    $g_name = $gene->display_xref->display_id;
  };	
  if ($@) {
    $g_name = $gene->get_all_Attributes('name')->[0]->value;
  }

  #get existing gene remarks
  my $remarks = [ map {$_->value} @{$gene->get_all_Attributes('remark')} ];

  #shout if there is no remark to identify this as being fragmented
  if ( grep {$_ eq 'fragmented locus' } @$remarks) {
    $study_more = 0;
  }
  else {
    $self->log_warning("Gene $gsi should have a fragmented locus remark\n");
  }

  ##patch transcript names according to length and CDS
  $gene_c++;

  #separate coding and non_coding transcripts
  my $coding_trans = [];
  my $noncoding_trans = [];
  foreach my $trans ( @{$gene->get_all_Transcripts()} ) {
    if ($trans->translate) {
      push @$coding_trans, $trans;
    }
    else {
      push @$noncoding_trans, $trans;
    }
  }

  #sort transcripts coding > non-coding, then on length
  my $c = 0;
  $self->log("\nPatching names according to CDS and length:\n",1);
  foreach my $array_ref ($coding_trans,$noncoding_trans) {
    foreach my $trans ( sort { $b->length <=> $a->length } @$array_ref ) {
      $trans_c++;
      my $tsi = $trans->stable_id;
      my $t_name;
      eval {
	$t_name = $trans->display_xref->display_id;
      };	
      if ($@) {
	$t_name = $trans->get_all_Attributes('name')->[0]->value;
      }
      $c++;
      my $ext = sprintf("%03d", $c);
      my $new_name = $g_name.'-'.$ext;
      $self->log(sprintf("%-20s%-3s%-20s", "$t_name ", "-->", "$new_name")."\n",1);
      if (! $self->param('dry_run')) {
	
	# update transcript display xref
	$dbh->do(qq(UPDATE xref x, external_db edb
                       SET x.display_label  = "$new_name"
                     WHERE x.external_db_id = edb.external_db_id
                       AND x.dbprimary_acc  = "$tsi"
                       AND edb.db_name      = "Vega_transcript"));
      }
    }
  }
  return ($study_more,$gene_c,$trans_c);
}

=head2 check_names_and_overlap

   Arg[1]     : arayref of arrayrefs of duplicated names
   Arg[2]     : B::E::Gene (with potentially duplicated transcript names)
   Arg[3]     : FH (to log new duplicates)
   Example    : $support->check_names_and_overlap($transcripts,$gene,$fh)
   Description: checks pairs of transcripts identified as having duplicate Vega names:
                - to see if they have identical names in loutre (shouldn't have)
                - distinguish between overlapping and non overlapping transcripts
   Returntype : none

=cut

sub check_names_and_overlap {
  my $self = shift;
  my ($transcript_info,$gene,$n_flist_fh) = @_;
  my $ta  = $gene->adaptor->db->get_TranscriptAdaptor;
  my $gsi = $gene->stable_id;
  my $g_name = $gene->get_all_Attributes('name')->[0]->value;
  foreach my $set (values %{$transcript_info} ) {
    next if (scalar @{$set} == 1);
    my $transcripts = [];
    my $all_t_names;
    my %ids_to_names;
    foreach my $id1 (@{$set}) {
      my ($name1,$tsi1) = split /\|/, $id1;
      $ids_to_names{$tsi1} = $name1;
      $all_t_names .= "$tsi1 [$name1] ";
      my $t = $ta->fetch_by_stable_id($tsi1);
      push @{$transcripts}, $t;
    }

    my $non_overlaps;
    eval {
      $non_overlaps = $self->find_non_overlaps($transcripts);
    };
    if ($@) {
      $self->log_warning("Problem looking for overlapping transcripts for gene $gsi (is_current = 0 ?). Skipping this bit\n");
    }

    #if the transcripts don't overlap
    elsif (@{$non_overlaps}) {
      my $tsi_string;
      foreach my $id (@{$non_overlaps}) {
	my $string = " $id [ $ids_to_names{$id} ] ";
	$tsi_string .= $string;
      }

      $self->log_warning("NEW: Non-overlapping: $gsi ($g_name) has non-overlapping transcripts ($tsi_string) with duplicated Vega names, and it has no \'fragmented locus\' gene remark. Neither has it been OKeyed by Havana before. Transcript names are being patched but this needs checking by Havana.\n");
      #log gsi (to be sent to Havana)
      print $n_flist_fh "$gsi\n";
    }
    #...otherwise if the transcripts do overlap
    else {
      $self->log_warning("NEW: Overlapping: $gsi ($g_name) has overlapping transcripts ($all_t_names) with duplicated Vega names and it has no \'fragmented locus\' gene_remark. Neither has it been OKeyed by Havana before. Transcript names are being patched but this could be checked by Havana if they were feeling keen.\n");
      print $n_flist_fh "$gsi\n";
    }
  }
}		

=head2 get_havana_fragmented_loci_comments

   Args       : none
   Example    : my $results = $support->get_havana_fragmented_loci_comments
   Description: parses the HEREDOC containing Havana comments in this module
   Returntype : hashref

=cut

sub get_havana_fragmented_loci_comments {
  my $seen_genes;
  while (<DATA>) {
    next if /^\s+$/ or /#+/;
    my ($obj,$comment) = split /=/;
    $obj =~ s/^\s+|\s+$//g;
    $comment =~ s/^\s+|\s+$//g;
    $seen_genes->{$obj} = $comment;
  }
  return $seen_genes;
}



#details of genes with duplicated transcript names that have already been reported to Havana
#identified as either fragmented or as being OK to patch
__DATA__

OTTMUSG00000005478 = fragmented
OTTMUSG00000001936 = fragmented
OTTMUSG00000017081 = fragmented
OTTMUSG00000011441 = fragmented
OTTMUSG00000013335 = fragmented
OTTMUSG00000011654 = fragmented
OTTMUSG00000001835 = fragmented
OTTHUMG00000035221 = fragmented
OTTHUMG00000037378 = fragmented
OTTHUMG00000060732 = fragmented
OTTHUMG00000132441 = fragmented
OTTHUMG00000031383 = fragmented
OTTHUMG00000012716 = fragmented
OTTHUMG00000031102 = fragmented
OTTHUMG00000148816 = fragmented
OTTHUMG00000149059 = fragmented
OTTHUMG00000149221 = fragmented
OTTHUMG00000149326 = fragmented
OTTHUMG00000149644 = fragmented
OTTHUMG00000149574 = fragmented
OTTHUMG00000058101 = fragmented

OTTHUMG00000150119 = OK
OTTHUMG00000149850 = OK
OTTHUMG00000058101 = OK
OTTHUMG00000058907 = OK

OTTMUSG00000011654 = fragmented
OTTMUSG00000019369 = fragmented
OTTMUSG00000017081 = fragmented
OTTMUSG00000001835 = fragmented
OTTMUSG00000011499 = fragmented
OTTMUSG00000013335 = fragmented
OTTMUSG00000008023 = fragmented
OTTMUSG00000019369 = fragmented


OTTMUSG00000022266
OTTMUSG00000006697





OTTMUSG00000012302 =
OTTMUSG00000013368 =
OTTMUSG00000015766 =
OTTMUSG00000016025 =
OTTMUSG00000001066 =
OTTMUSG00000016331 =
OTTMUSG00000006935 =
OTTMUSG00000007263 =
OTTMUSG00000000304 =
OTTMUSG00000009150 =
OTTMUSG00000008023 =
OTTMUSG00000017077 =
OTTMUSG00000003440 =
OTTMUSG00000016310 =
OTTMUSG00000026199 =
OTTMUSG00000028423 =
OTTMUSG00000007427 =
