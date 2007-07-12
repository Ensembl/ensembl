package Bio::EnsEMBL::IdMapping::InternalIdMapper;

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS


=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http:#www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut


use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::BaseObject;
our @ISA = qw(Bio::EnsEMBL::IdMapping::BaseObject);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


# scores are considered the same if (2.0 * (s1-s2))/(s1 + s2) < this
use constant SIMILAR_SCORE_RATIO => 0.01;

    
sub map_genes {
  my $self = shift;
  my $gene_scores = shift;
  my $transcript_scores = shift;
  my $gsb = shift;

  # argument checks
  unless ($gene_scores and
          $gene_scores->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a gene Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }
  
  unless ($transcript_scores and
          $transcript_scores->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a transcript Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  unless ($gsb and
          $gsb->isa('Bio::EnsEMBL::IdMapping::GeneScoreBuilder')) {
    throw('Need a Bio::EnsEMBL::IdMapping::GeneScoreBuilder.');
  }
  
  $self->logger->info("Starting gene mapping...\n\n", 0, 'stamped');

  my $mappings = Bio::EnsEMBL::IdMapping::MappingList->new(
    -DUMP_PATH   => $self->conf->param('dumppath'),
    -CACHE_FILE  => 'gene_mappings.ser',
  );

  my $mapping_cache = $mappings->cache_file;

  if (-s $mapping_cache) {
    
    # read from file
    $self->logger->info("Reading gene mappings from file...\n", 0, 'stamped');
    $self->logger->debug("Cache file $mapping_cache.\n", 1);
    $matrix->read_from_file;
    $self->logger->info("Done.\n\n", 0, 'stamped');
    
  } else {
    
    # create gene mappings
    $self->logger->info("No gene mappings found. Will calculate then now.\n");

    $self->logger->info("Mapping genes...\n", 0, 'stamped');

    #
    # basic mapping
    #
    $self->basic_mapping($gene_scores, $mappings);
    my $synteny_gene_scores = $gsb->create_shrinked_matrix($gene_scores,
      $mappings, 'synteny_gene_matrix.ser');

    #
    # build the synteny from unambiguous mappings
    #
    $self->logger->info("Synteny Framework building\n");
    my $sf = Bio::EnsEMBL::IdMapping::SyntenyFramework->new(
      -LOGGER       => $self->logger,
      -CONF         => $self->conf,
      -CACHE        => $self->cache
    );
    $sf->build_synteny($mappings);

    # use it to rescore the genes
    $self->logger->info("Synteny assisted mapping\n");
    $sf->rescore_gene_matrix($synteny_gene_scores);

    my $synteny_mappings = Bio::EnsEMBL::IdMapping::MappingList->new(
      -DUMP_PATH   => $self->conf->param('dumppath'),
      -CACHE_FILE  => 'gene_synteny_mappings.ser',
    );
    
    $self->basic_mapping($synteny_gene_scores, $synteny_mappings);
    
    $self->logger->info("Found ".$synteny_mappings->get_entry_count.
      " additional mappings.\n");


    #
    # rescore with simple scoring function and try again
    #
    $self->logger->info("Retry with simple best transcript score\n");
    
    my $simple_gene_scores = $gsb->create_shrinked_matrix($synteny_gene_scores,
      $synteny_mappings, 'simple_gene_matrix.ser');
    
    $gsb->simple_gene_rescore($simple_gene_scores, $transcript_scores);
    
    my $simple_mappings = Bio::EnsEMBL::IdMapping::MappingList->new(
      -DUMP_PATH   => $self->conf->param('dumppath'),
      -CACHE_FILE  => 'gene_simple_mappings.ser',
    );
    
    $self->basic_mapping($simple_gene_scores, $simple_mappings);
    
    $self->logger->info("Found ".$simple_mappings->get_entry_count.
      " additional mappings.\n");


    #
    # rescore by penalising scores between genes with different biotypes  
    #
    $self->logger->info("Retry with biotype disambiguation\n");
    
    my $biotype_gene_scores = $gsb->create_shrinked_matrix($simple_gene_scores, 
      $simple_mappings, 'biotype_gene_matrix.ser');

    $gsb->biotype_gene_rescore($biotype_gene_scores);

    my $biotype_mappings = Bio::EnsEMBL::IdMapping::MappingList->new(
      -DUMP_PATH   => $self->conf->param('dumppath'),
      -CACHE_FILE  => 'gene_biotype_mappings.ser',
    );
    
    $self->basic_mapping($biotype_gene_scores, $biotype_mappings);
    
    $self->logger->info("Found ".$biotype_mappings->get_entry_count.
      " additional mappings.\n");
    

    #
    # selectively rescore by penalising scores between genes with different
    # internalIDs  
    #
    $self->logger->info("Retry with internalID disambiguation\n");
    
    my $internal_id_gene_scores = $gsb->create_shrinked_matrix(
      $biotype_gene_scores, $biotype_mappings, 'internal_id_gene_matrix.ser');

    $gsb->internal_id_gene_rescore($internal_id_gene_scores);

    my $internal_id_mappings = Bio::EnsEMBL::IdMapping::MappingList->new(
      -DUMP_PATH   => $self->conf->param('dumppath'),
      -CACHE_FILE  => 'gene_internal_id_mappings.ser',
    );
    
    $self->basic_mapping($internal_id_gene_scores, $internal_id_mappings);
    
    $self->logger->info("Found ".$internal_id_mappings->get_entry_count.
      " additional mappings.\n");


    #
    # report remaining ambiguities
    #
    my $remaining_gene_scores = $gsb->create_shrinked_matrix(
      $internal_id_gene_scores, $internal_id_mappings,
      'remaining_gene_matrix.ser');

    $self->logger->info("\n".$remaining_gene_scores->get_source_count.
      " source genes are ambiguous with ".
      $remaining_gene_scores->get_target_count." target genes.\n\n");

    $self->log_ambiguous($remaining_gene_scores, 'gene');

    
    #
    # merge mappings and write to file
    #
    $mappings->add_all($synteny_mappings, $simple_mappings, $biotype_mappings,
                       $internal_id_mappings);

    $mappings->write_to_file;

    if ($self->logger->loglevel eq 'debug') {
      $mappings->log('gene');
    }

    $self->logger->info("Done.\n\n", 0, 'stamped');

  }

  return $mappings;

}

#
# find the highest unambiguous score for all sources and targets in a scoring
# matrix
#
sub basic_mapping {
  my $self = shift;
  my $matrix = shift;
  my $mappings = shift;

  # argument checks
  unless ($matrix and
          $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }
  
  unless ($mappings and
          $mappings->isa('Bio::EnsEMBL::IdMapping::MappingList')) {
    throw('Need a gene Bio::EnsEMBL::IdMapping::MappingList.');
  }

  my $sources_done = {};
  my $targets_done = {};

  # sort scoring matrix entries by descending score
  my @sorted_entries = sort { $b->score <=> $a->score }
    @{ $matrix->get_all_Entries };

  while (my $entry = shift(@sorted_entries)) {
    
    # we already found a mapping for either source or target yet
    next if ($sources_done->{$entry->source} or
             $targets_done->{$entry->target});
    
    # there's a better mapping for either source or target
    next if ($self->higher_score_exists($entry, $matrix, $sources_done,
      $targets_done));

    # check for ambiguous mappings; they are dealt with later
    my $other_sources = [];
    my $other_targets = [];

    if ($self->ambiguous_mapping($entry, $matrix, $other_sources, $other_targets)) {
      $self->filter_sources($other_sources, $sources_done);
      $self->filter_targets($other_targets, $targets_done);

      next if (scalar(@$other_sources) or scalar(@$other_targets));
    }

    # this is the best mapping, add it
    $mappings->add_Entry($entry);

    $sources_done->{$entry->source} = 1;
    $targets_done->{$entry->target} = 1;
  }
  
}


sub higher_score_exists {
  my ($self, $entry, $matrix, $sources_done, $targets_done) = @_;

  my $source = $entry->source;
  my $target = $entry->target;
  my $score = $entry->score;

  foreach my $other_source @{ $matrix->get_sources_for_target($target) }) {
    if ($other_source != $source and !$sources_done->{$other_source} and
        $score < $matrix->get_score($other_source, $target)) {
          return 1;
    }
  }

  foreach my $other_target @{ $matrix->get_targets_for_source($source) }) {
    if ($other_target != $target and !$targets_done->{$other_target} and
        $score < $matrix->get_score($source, $other_target)) {
          return 1;
    }
  }

  return 0;
}


#
# find ambiguous mappings (see scores_similar() for definition)
#
sub ambiguous_mapping {
  my ($self, $entry, $matrix, $other_sources, $other_targets) = @_;

  my $source = $entry->source;
  my $target = $entry->target;
  my $score = $entry->score;

  my $retval = 0;

  foreach my $other_source @{ $matrix->get_sources_for_target($target) }) {
    my $other_score = $matrix->get_score($other_source, $target);
    
    if ($other_source != $source and
      ($self->scores_similar($score, $other_score) or $score < $other_score)) {
        retval = 1;
        push @{ $other_sources }, $other_source;
    }
  }

  foreach my $other_target @{ $matrix->get_targets_for_source($source) }) {
    my $other_score = $matrix->get_score($source, $other_target);
    
    if ($other_target != $target and
      ($self->scores_similar($score, $other_score) or $score < $other_score)) {
        retval = 1;
        push @{ $other_targets }, $other_target;
    }
  }

  return $retval;
}


# 
# rule for similarity taken from java code...
#
sub scores_similar {
  my ($self, $s1, $s2) = @_;

  # always give priority to exact matches over very similar ones
  return 0 if ($s1 == 1 and $s2 < 1);

  my $diff = $s1 -$s2;
  $diff = -$diff if ($diff < 0);
  
  my $pc = 2 * $diff / ($s1 + $s2);
  
  return ($pc < SIMILAR_SCORE_RATIO);
}


sub filter_sources {
  my ($self, $other_sources, $sources_done) = @_;

  return 0 unless (scalar(@$other_sources) and scalar(@$sources_done));

  my @tmp = ();

  foreach my $e (@{ $other_sources }) {
    push @tmp, $e unless ($sources_done->{$e}); 
  }

  $other_sources = \@tmp;
}


sub filter_targets {
  my ($self, $other_targets, $targets_done) = @_;

  return 0 unless (scalar(@{ $other_targets) and scalar(@$targets_done));

  my @tmp = ();

  foreach my $e (@{ $other_targets }) {
    push @tmp, $e unless ($targets_done->{$e}); 
  }

  $other_targets = \@tmp;
}


1;

