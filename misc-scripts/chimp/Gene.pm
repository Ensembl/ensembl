use strict;
use warnings;



#
# Utility methods for the storing of chimp genes
#

package Gene;


use Bio::EnsEMBL::Utils::Exception qw(info throw warning);
use Bio::EnsEMBL::DBEntry;

use constant NEAR => 1.5e6;

our %KEEP_XREF = ('SWISSPROT' => 1,
                  'SPTREMBL'  => 1,
                  'HUGO'      => 1);

###############################################################################
# store gene
#
# Builds Ensembl genes from the generated chimp transcripts and stores them
# in the database.
#
###############################################################################

sub store_gene {
  my $db = shift;
  my $hum_gene = shift; # human gene
  my $ctranscripts = shift; # chimp transcripts

  my $MIN_AA_LEN = 15;
  my $MIN_NT_LEN = 600;

  my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name('ensembl');

  # Look at the translations and convert any transcripts with stop codons
  # into pseudogenes
  foreach my $ct (@$ctranscripts) {
    if($ct->translation && $ct->translate->seq() =~ /\*/) {
      $ct->translation(undef);
    }
  }

  # create xrefs to reference the human transcripts and translations
  create_ensembl_xrefs($ctranscripts);

  # transfer xrefs from human transcripts/translations
  transfer_xrefs($hum_gene, $ctranscripts);

  # compact duplicate transcripts into the same transcripts
  $ctranscripts = compact_transcripts($ctranscripts);

  # group all close together transcripts on the same strand into clusters
  my $clusters = cluster_transcripts($ctranscripts);

  my $gene_adaptor = $db->get_GeneAdaptor();

  foreach my $cluster (@$clusters) {

    # keep genes only if there is a minimum amount of nucleotide
    # OR amino acid sequence in transcripts in the same region
    if($cluster->{'nt_len'} < $MIN_NT_LEN &&
       $cluster->{'aa_len'} < $MIN_AA_LEN) {
      next;
    }

    # one gene for each cluster
    my $cgene = Bio::EnsEMBL::Gene->new();

    # add reference to the original human gene
    $cgene->add_DBEntry(Bio::EnsEMBL::DBEntry->new
                        (-primary_id => $hum_gene->stable_id(),
                         -version    => $hum_gene->version(),
                         -dbname     => 'Ens_Hs_gene',
                         -release    => 1,
                         -display_id => $hum_gene->stable_id()));

    generate_stable_id($cgene);

    # rename transcripts and add to gene
    foreach my $ctrans (@{$cluster->{'transcripts'}}) {
      generate_stable_id($ctrans);

      # rename translation
      if($ctrans->translation) {
        generate_stable_id($ctrans->translation);
      }

      $cgene->add_Transcript($ctrans);
    }

    # rename all of the exons
    # but watch out because duplicate exons will be merged and we do not
    # want to generate multiple names
    my %ex_stable_ids;
    foreach my $ex (@{$cgene->get_all_Exons()}) {
      if($ex_stable_ids{$ex->hashkey()}) {
        $ex->stable_id($ex_stable_ids{$ex->hashkey()});
      } else {
        generate_stable_id($ex);
        $ex_stable_ids{$ex->hashkey()} = $ex->stable_id();
      }
    }


    foreach my $gx (@{$hum_gene->get_all_DBEntries}) {
      $cgene->add_DBEntry($gx) if($KEEP_XREF{uc($gx->dbname())});
    }

    if($hum_gene->display_xref &&
       $KEEP_XREF{uc($hum_gene->display_xref->dbname)}){
      $cgene->display_xref($hum_gene->display_xref);
    }

    # set the analysis on the gene object
    $cgene->analysis($analysis);

    my $name = $cgene->stable_id();

    $name .= '/'.$cgene->display_xref->display_id() if($cgene->display_xref());

    $cgene->type('ensembl');

    # store the bloody thing
    print STDERR "Storing gene: $name\n";
    $gene_adaptor->store($cgene);
  }

  return;
}



###############################################################################
# generate_stable_id
#
# Generates a stable_id for a gene, transcript, translation or exon and sets
# it on the object.
#
###############################################################################


my ($TRANSCRIPT_NUM, $GENE_NUM, $EXON_NUM, $TRANSLATION_NUM);


sub generate_stable_id {
  my $object = shift;

  my $SPECIES_PREFIX = 'PTR';
  my $PAD            = 18;

  my $type_prefix;
  my $num;

  if($object->isa('Bio::EnsEMBL::Exon')) {
    $type_prefix = 'E';
    $EXON_NUM       ||= 0;
    $num = ++$EXON_NUM;
  } elsif($object->isa('Bio::EnsEMBL::Transcript')) {
    $type_prefix = 'T';
    $TRANSCRIPT_NUM ||= 0;
    $num = ++$TRANSCRIPT_NUM;
  } elsif($object->isa('Bio::EnsEMBL::Gene')) {
    $type_prefix = 'G';
    $GENE_NUM       ||= 0;
    $num = ++$GENE_NUM;
  } elsif($object->isa('Bio::EnsEMBL::Translation')) {
    $type_prefix = 'P';
    $TRANSLATION_NUM ||= 0;
    $num = ++$TRANSLATION_NUM;
  } else {
    throw('Unknown object type '.ref($object).'. Cannot create stable_id.');
  }

  my $prefix = "ENS${SPECIES_PREFIX}${type_prefix}";

  my $pad = $PAD - length($prefix) - length($num);

  $object->version(1);
  $object->stable_id($prefix . ('0'x$pad) . $num);
}





###############################################################################
# compact_transcripts
#
# Given a set of transcripts this function removes duplicate transcripts and
# returns the unique set
#
###############################################################################

sub compact_transcripts {
  my $transcripts = shift;

  my %unique_hash;
  my @unique_list;

  info("Compacting transcripts");

  foreach my $transcript (@$transcripts) {
    my $hashkey = 'exons:';

    foreach my $exon (@{$transcript->get_all_Exons}) {
      $hashkey .= '('.$exon->hashkey.')';
    }

    # include the translation in the hashkey because sometimes the
    # same transcripts may end up with different translations
    # (this is especially possible when split transcripts become partially the
    # same)
    if($transcript->translation) {
      $hashkey .= 'translation:' .
        $transcript->translation->start() . '-' .
        $transcript->translation->end() . '(' .
        $transcript->translation->start_Exon->hashkey() . ')('.
        $transcript->translation->end_Exon->hashkey() . ')';
    }


    $unique_hash{$hashkey} ||= [];
    push @{$unique_hash{$hashkey}}, $transcript;
  }

  # merge xrefs from duplicated transcripts

  foreach my $key (keys %unique_hash) {
    # choose one of the duplicates arbitrarily
    my $duplicates = $unique_hash{$key};
    my $transcript = pop(@$duplicates);

    # merge all of the xrefs and add them to the one transcript that will
    # be kept
    merge_xrefs($transcript, $duplicates);
    push @unique_list, $transcript;
  }

  return \@unique_list;
}



sub merge_xrefs {
  my $kept_transcript = shift;
  my $duplicates      = shift;

  return if(!@$duplicates);

  info('Merging xrefs from duplicate transcripts');

  # construct a set of xrefs the transcript to keep already has

  my %existing_tl_xrefs;
  my %existing_tr_xrefs;
  foreach my $xref (@{$kept_transcript->get_all_DBEntries()}) {
    $existing_tr_xrefs{$xref->dbname().':'.$xref->primary_id()} = 1;
  }

  if($kept_transcript->translation()) {
    foreach my $xref (@{$kept_transcript->translation->get_all_DBEntries()}) {
      $existing_tl_xrefs{$xref->dbname().':'.$xref->primary_id()} = 1;
    }
  }

  # collect a list of additional xrefs which should be kept

  my %tl_xrefs;
  my %tr_xrefs;
  foreach my $dup (@$duplicates) {
    foreach my $xref (@{$dup->get_all_DBEntries()}) {
      my $key = $xref->dbname().':'.$xref->primary_id();
      $tr_xrefs{$key} = $xref if(!$existing_tr_xrefs{$key});
    }

    if($dup->translation()) {
      foreach my $xref (@{$dup->translation->get_all_DBEntries()}) {
        my $key = $xref->dbname().':'.$xref->primary_id();
        $tl_xrefs{$key} = $xref if(!$existing_tl_xrefs{$key});
      }
    }
  }

  # add any new xrefs which were found to the kept translation

  foreach my $xref (values %tr_xrefs) {
    $kept_transcript->add_DBEntry($xref);
  }

  my @new_tl_xrefs = values(%tl_xrefs);
  my $tl = $kept_transcript->translation();

  if(@new_tl_xrefs && !$tl) {
    throw("Some duplicate transcripts have translations, and others do not?");
    return;
  }

  foreach my $xref (@new_tl_xrefs) {
    $tl->add_DBEntry($xref);
  }

  return;
}


###############################################################################
# cluster_transcripts
#
# Given a set of transcripts this function will cluster them into groups
# which are near to each other. Note this may be more efficient to use
# a sorted list, but for the small size of the clusters it should not matter.
#
##############################################################################

sub cluster_transcripts {
  my $transcripts = shift;

  my @clusters;

  info("Clustering transcripts");

  foreach my $tr (@$transcripts) {
    my $cl = undef;

    foreach my $c (@clusters) {
      if(is_near($tr->start(), $tr->end(), $tr->strand(), $tr->slice(),
                 $c->{'start'}, $c->{'end'}, $c->{'strand'}, $c->{'slice'})) {

        $cl = $c;
      }
    }

    if($cl) {
      push @{$cl->{'transcripts'}}, $tr;

      $cl->{'end'} =
        ( $tr->end > $cl->{'end'} )     ? $tr->{'end'}   : $cl->{'end'};
      $cl->{'start'} =
        ( $tr->start < $cl->{'start'} ) ? $tr->{'start'} : $cl->{'start'};

    } else {
      $cl = {'start'  => $tr->start(),
             'end'    => $tr->end(),
             'strand' => $tr->strand(),
             'slice'  => $tr->slice(),
             'transcripts' => [$tr]};

      push @clusters, $cl;
    }
  }

  # now cluster clusters

  for(my $i = 0; $i < @clusters; $i++) {
    for(my $j = $i+1; $j < @clusters; $j++) {
      my $c1 = $clusters[$i];
      my $c2 = $clusters[$j];
      if(is_near($c1->{'start'}, $c1->{'end'}, $c1->{'strand'}, $c1->{'slice'},
              $c2->{'start'}, $c2->{'end'}, $c2->{'strand'}, $c2->{'slice'})) {
        # merge clusters and take one of them out of the list
        splice(@clusters, $j, 1);

        $c1->{'start'} =
          ($c1->{'start'} < $c2->{'start'}) ? $c1->{'start'} : $c2->{'start'};

        $c1->{'end'} =
          ($c1->{'end'} > $c2->{'end'}) ? $c1->{'end'} : $c2->{'end'};

        push @{$c1->{'transcripts'}}, @{$c2->{'transcripts'}};

        # start over again
        $i = -1;
        $j = -1;
      }
    }
  }

  # sum the nucleotide length and amino acid lengths of each of the
  # transcripts in a given cluster
  foreach my $cl (@clusters) {
    $cl->{'nt_len'} = 0;
    $cl->{'aa_len'} = 0;

    foreach my $tr (@{$cl->{'transcripts'}}) {
      $cl->{'nt_len'} += length($tr->spliced_seq());
      if($tr->translation) {
        $cl->{'aa_len'} += length($tr->translate->seq());
      }
    }
  }


  return \@clusters;
}



sub is_near {
  my ($start1,$end1,$strand1, $slice1, $start2, $end2, $strand2, $slice2) = @_;

  # cannot be clustered if not on same strand
  if($strand1 != $strand2) {
    return 0;
  }

  # cannot be clustered if not on same slice
  if($slice1->name() ne $slice2->name()) {
    return 0;
  }

  # check if the regions overlap
  if($end1 >= $start2 && $start1 <= $end2) {
    return 1;
  }

  if($start1 > $end2) {
    return (($start1 - $end2) < NEAR) ? 1 : 0;
  }

  return (($start2 - $end1) < NEAR) ? 1 : 0;
}






sub transfer_xrefs {
  my $hum_gene = shift;
  my $chimp_transcripts = shift;

  my %chimp_transcripts;
  my %chimp_translations;

  foreach my $tr (@$chimp_transcripts) {
    $chimp_transcripts{$tr->stable_id()} ||= [];
    push @{$chimp_transcripts{$tr->stable_id()}}, $tr;

    my $tl = $tr->translation();

    if($tl) {
      $chimp_translations{$tl->stable_id()} ||= [];
      push @{$chimp_translations{$tl->stable_id()}}, $tl;
    }
  }

  foreach my $tr (@{$hum_gene->get_all_Transcripts()}) {
    foreach my $chimp_tr (@{$chimp_transcripts{$tr->stable_id}}) {
      foreach my $xref (@{$tr->get_all_DBEntries}) {
        $chimp_tr->add_DBEntry($xref) if($KEEP_XREF{uc($xref->dbname())});
      }

      if($tr->display_xref() && $KEEP_XREF{uc($tr->display_xref->dbname)}) {
        $chimp_tr->display_xref($tr->display_xref);
      }
    }

    my $tl = $tr->translation();
    if($tl) {
      foreach my $xref (@{$tl->get_all_DBEntries}) {
        foreach my $chimp_tl (@{$chimp_translations{$tl->stable_id()}}) {
          $chimp_tl->add_DBEntry($xref) if($KEEP_XREF{uc($xref->dbname())});
        }
      }
    }
  }

  return;
}



sub create_ensembl_xrefs {
  my $chimp_transcripts = shift;

  foreach my $transcript (@$chimp_transcripts) {
    my $dbe = Bio::EnsEMBL::DBEntry->new
      (-primary_id => $transcript->stable_id(),
       -version    => $transcript->version(),
       -dbname     => 'Ens_Hs_transcript',
       -release    => 1,
       -display_id => $transcript->stable_id());
    $transcript->add_DBEntry($dbe);

    if($transcript->translation()) {
      $dbe = Bio::EnsEMBL::DBEntry->new
        (-primary_id => $transcript->translation->stable_id(),
         -version    => $transcript->translation->version(),
         -dbname     => 'Ens_Hs_translation',
         -release    => 1,
         -display_id => $transcript->translation->stable_id());
      $transcript->translation->add_DBEntry($dbe);
    }
  }
}

1;
