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


package Bio::EnsEMBL::IdMapping::InternalIdMapper;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::BaseObject;
our @ISA = qw(Bio::EnsEMBL::IdMapping::BaseObject);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(inject path_append);
use Bio::EnsEMBL::IdMapping::Entry;
use Bio::EnsEMBL::IdMapping::MappingList;
use Bio::EnsEMBL::IdMapping::SyntenyFramework;


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
  
  $self->logger->info("== Internal ID mapping for genes...\n\n", 0, 'stamped');

  my $dump_path = path_append($self->conf->param('basedir'), 'mapping');

  my $mappings = Bio::EnsEMBL::IdMapping::MappingList->new(
    -DUMP_PATH   => $dump_path,
    -CACHE_FILE  => 'gene_mappings.ser',
  );

  my $mapping_cache = $mappings->cache_file;

  if (-s $mapping_cache) {
    
    # read from file
    $self->logger->info("Reading gene mappings from file...\n", 0, 'stamped');
    $self->logger->debug("Cache file $mapping_cache.\n", 1);
    $mappings->read_from_file;
    $self->logger->info("Done.\n\n", 0, 'stamped');
    
  } else {
    
    # create gene mappings
    $self->logger->info("No gene mappings found. Will calculate them now.\n");

    # determine which plugin methods to run
    my @default_plugins = (qw(
      Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblGeneGeneric::init_basic
      Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblGeneGeneric::location
      Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblGeneGeneric::best_transcript
      Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblGeneGeneric::biotype
      Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblGeneGeneric::synteny
      Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblGeneGeneric::internal_id
    ));

    my @plugins = $self->conf->param('plugin_internal_id_mappers_gene');
    @plugins = @default_plugins unless (defined($plugins[0]));

    my $new_mappings = Bio::EnsEMBL::IdMapping::MappingList->new(
      -DUMP_PATH   => $dump_path,
      -CACHE_FILE  => 'gene_mappings0.ser',
    );
    my @mappings = ();
    my $i = 0;

    #
    # run the scoring chain
    #
    foreach my $plugin (@plugins) {
      ($gene_scores, $new_mappings) = $self->delegate_to_plugin($plugin, $i++,
        $gsb, $new_mappings, $gene_scores, $transcript_scores);

      push(@mappings, $new_mappings);
    }

    # report remaining ambiguities
    $self->logger->info($gene_scores->get_source_count.
      " source genes are ambiguous with ".
      $gene_scores->get_target_count." target genes.\n\n");

    $self->log_ambiguous($gene_scores, 'gene');
    
    # merge mappings and write to file
    $mappings->add_all(@mappings);
    $mappings->write_to_file;

    if ($self->logger->loglevel eq 'debug') {
      $mappings->log('gene', $self->conf->param('basedir'));
    }

    $self->logger->info("Done.\n\n", 0, 'stamped');

  }

  return $mappings;
}


sub map_transcripts {
  my $self = shift;
  my $transcript_scores = shift;
  my $gene_mappings = shift;
  my $tsb = shift;

  # argument checks
  unless ($transcript_scores and
      $transcript_scores->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a transcript Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  unless ($gene_mappings and
          $gene_mappings->isa('Bio::EnsEMBL::IdMapping::MappingList')) {
    throw('Need a gene Bio::EnsEMBL::IdMapping::MappingList.');
  }
  
  unless ($tsb and
          $tsb->isa('Bio::EnsEMBL::IdMapping::TranscriptScoreBuilder')) {
    throw('Need a Bio::EnsEMBL::IdMapping::TranscriptScoreBuilder.');
  }
  
  $self->logger->info("== Internal ID mapping for transcripts...\n\n", 0, 'stamped');

  my $dump_path = path_append($self->conf->param('basedir'), 'mapping');

  my $mappings = Bio::EnsEMBL::IdMapping::MappingList->new(
    -DUMP_PATH   => $dump_path,
    -CACHE_FILE  => 'transcript_mappings.ser',
  );

  my $mapping_cache = $mappings->cache_file;

  if (-s $mapping_cache) {
    
    # read from file
    $self->logger->info("Reading transcript mappings from file...\n", 0,
      'stamped');
    $self->logger->debug("Cache file $mapping_cache.\n", 1);
    $mappings->read_from_file;
    $self->logger->info("Done.\n\n", 0, 'stamped');
    
  } else {
    
    # create transcript mappings
    $self->logger->info("No transcript mappings found. Will calculate them now.\n");

    # determine which plugin methods to run
    my @default_plugins = (qw(
      Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblTranscriptGeneric::init_basic
      Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblTranscriptGeneric::non_exact_translation
      Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblTranscriptGeneric::biotype
      Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblTranscriptGeneric::mapped_gene
      Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblTranscriptGeneric::single_gene
      Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblTranscriptGeneric::internal_id
    ));

    my @plugins = $self->conf->param('plugin_internal_id_mappers_transcript');
    @plugins = @default_plugins unless (defined($plugins[0]));

    my $new_mappings = Bio::EnsEMBL::IdMapping::MappingList->new(
      -DUMP_PATH   => $dump_path,
      -CACHE_FILE  => 'transcript_mappings0.ser',
    );
    my @mappings = ();
    my $i = 0;

    #
    # run the scoring chain
    #
    foreach my $plugin (@plugins) {
      ($transcript_scores, $new_mappings) = $self->delegate_to_plugin($plugin,
        $i++, $tsb, $new_mappings, $transcript_scores, $gene_mappings);

      push(@mappings, $new_mappings);
    }

    # report remaining ambiguities
    $self->logger->info($transcript_scores->get_source_count.
      " source transcripts are ambiguous with ".
      $transcript_scores->get_target_count." target transcripts.\n\n");

    $self->log_ambiguous($transcript_scores, 'transcript');

    # merge mappings and write to file
    $mappings->add_all(@mappings);
    $mappings->write_to_file;

    if ($self->logger->loglevel eq 'debug') {
      $mappings->log('transcript', $self->conf->param('basedir'));
    }

    $self->logger->info("Done.\n\n", 0, 'stamped');

  }

  return $mappings;

}


sub map_exons {
  my $self = shift;
  my $exon_scores = shift;
  my $transcript_mappings = shift;
  my $esb = shift;

  # argument checks
  unless ($exon_scores and
      $exon_scores->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix of exons.');
  }

  unless ($transcript_mappings and
          $transcript_mappings->isa('Bio::EnsEMBL::IdMapping::MappingList')) {
    throw('Need a Bio::EnsEMBL::IdMapping::MappingList of transcripts.');
  }
  
  unless ($esb and
          $esb->isa('Bio::EnsEMBL::IdMapping::ExonScoreBuilder')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ExonScoreBuilder.');
  }
  
  $self->logger->info("== Internal ID mapping for exons...\n\n", 0, 'stamped');

  my $dump_path = path_append($self->conf->param('basedir'), 'mapping');

  my $mappings = Bio::EnsEMBL::IdMapping::MappingList->new(
    -DUMP_PATH   => $dump_path,
    -CACHE_FILE  => 'exon_mappings.ser',
  );

  my $mapping_cache = $mappings->cache_file;

  if (-s $mapping_cache) {
    
    # read from file
    $self->logger->info("Reading exon mappings from file...\n", 0,
      'stamped');
    $self->logger->debug("Cache file $mapping_cache.\n", 1);
    $mappings->read_from_file;
    $self->logger->info("Done.\n\n", 0, 'stamped');
    
  } else {
    
    # create exon mappings
    $self->logger->info("No exon mappings found. Will calculate them now.\n");

    # determine which plugin methods to run
    my @default_plugins = (qw(
      Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblExonGeneric::init_basic
      Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblExonGeneric::mapped_transcript
      Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblExonGeneric::single_transcript
      Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblExonGeneric::internal_id
    ));

    my @plugins = $self->conf->param('plugin_internal_id_mappers_exon');
    @plugins = @default_plugins unless (defined($plugins[0]));

    my $new_mappings = Bio::EnsEMBL::IdMapping::MappingList->new(
      -DUMP_PATH   => $dump_path,
      -CACHE_FILE  => 'exon_mappings0.ser',
    );
    my @mappings = ();
    my $i = 0;

    #
    # run the scoring chain
    #
    foreach my $plugin (@plugins) {
      ($exon_scores, $new_mappings) = $self->delegate_to_plugin($plugin, $i++,
        $esb, $new_mappings, $exon_scores, $transcript_mappings);

      push(@mappings, $new_mappings);
    }

    # report remaining ambiguities
    $self->logger->info($exon_scores->get_source_count.
      " source exons are ambiguous with ".
      $exon_scores->get_target_count." target exons.\n\n");

    $self->log_ambiguous($exon_scores, 'exon');

    # merge mappings and write to file
    $mappings->add_all(@mappings);
    $mappings->write_to_file;

    if ($self->logger->loglevel eq 'debug') {
      $mappings->log('exon', $self->conf->param('basedir'));
    }

    $self->logger->info("Done.\n\n", 0, 'stamped');

  }

  return $mappings;

}


#
# this is not implemented as a plugin, since a) it's too simple and b) it's
# tied to transcripts so there are no translation scores or score builder.
#
sub map_translations {
  my $self = shift;
  my $transcript_mappings = shift;

  # argument checks
  unless ($transcript_mappings and
          $transcript_mappings->isa('Bio::EnsEMBL::IdMapping::MappingList')) {
    throw('Need a Bio::EnsEMBL::IdMapping::MappingList of transcripts.');
  }
  
  $self->logger->info("== Internal ID mapping for translations...\n\n", 0, 'stamped');

  my $dump_path = path_append($self->conf->param('basedir'), 'mapping');

  my $mappings = Bio::EnsEMBL::IdMapping::MappingList->new(
    -DUMP_PATH   => $dump_path,
    -CACHE_FILE  => 'translation_mappings.ser',
  );

  my $mapping_cache = $mappings->cache_file;

  if (-s $mapping_cache) {
    
    # read from file
    $self->logger->info("Reading translation mappings from file...\n", 0,
      'stamped');
    $self->logger->debug("Cache file $mapping_cache.\n", 1);
    $mappings->read_from_file;
    $self->logger->info("Done.\n\n", 0, 'stamped');
    
  } else {
    
    # create translation mappings
    $self->logger->info("No translation mappings found. Will calculate them now.\n");

    $self->logger->info("Translation mapping...\n", 0, 'stamped');

    #
    # map translations for mapped transcripts
    #
    my $i = 0;

    foreach my $entry (@{ $transcript_mappings->get_all_Entries }) {

      my $source_tl = $self->cache->get_by_key('transcripts_by_id',
        'source', $entry->source)->translation;
      my $target_tl = $self->cache->get_by_key('transcripts_by_id',
        'target', $entry->target)->translation;

      if ($source_tl and $target_tl) {
      
        # add mapping for the translations; note that the score is taken from
        # the transcript mapping
        my $tl_entry = Bio::EnsEMBL::IdMapping::Entry->new_fast([
          $source_tl->id, $target_tl->id, $entry->score
        ]);
        $mappings->add_Entry($tl_entry);
      
      } else {
        $i++;
      }

    }

    $self->logger->debug("Skipped transcripts without translation: $i\n", 1);
    $self->logger->info("New mappings: ".$mappings->get_entry_count."\n\n");

    $mappings->write_to_file;

    if ($self->logger->loglevel eq 'debug') {
      $mappings->log('translation', $self->conf->param('basedir'));
    }

    $self->logger->info("Done.\n\n", 0, 'stamped');

  }

  return $mappings;

}


sub delegate_to_plugin {
  my $self = shift;
  my $plugin = shift;
  my $num = shift;
  my $score_builder = shift;
  my $mappings = shift;
  my $scores = shift;

  # argument checks
  unless ($score_builder and
          $score_builder->isa('Bio::EnsEMBL::IdMapping::ScoreBuilder')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoreBuilder.');
  }

  unless ($mappings and
          $mappings->isa('Bio::EnsEMBL::IdMapping::MappingList')) {
    throw('Need a Bio::EnsEMBL::IdMapping::MappingList.');
  }
  
  unless ($scores and
          $scores->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  # split plugin name into module and method
  $plugin =~ /(.*)::(\w+)$/;
  my $module = $1;
  my $method = $2;

  unless ($module and $method) {
    throw("Unable to determine module and method name from $plugin.\n");
  }

  # instantiate the plugin unless we already have an instance
  my $plugin_instance;
  if ($self->has_plugin($module)) {
    
    # re-use an existing plugin instance
    $plugin_instance = $self->get_plugin($module);
  
  } else {
    
    # inject and instantiate the plugin module
    inject($module);
    $plugin_instance = $module->new(
        -LOGGER       => $self->logger,
        -CONF         => $self->conf,
        -CACHE        => $self->cache
    );
    $self->add_plugin($plugin_instance);

  }

  # run the method on the plugin
  #
  # pass in a sequence number (number of method run, used for generating
  # checkpoint files), the scores used for determining the mapping, and all
  # other arguments passed to this method (these will vary for different object
  # types)
  #
  # return the scores and mappings to feed into the next plugin in the chain
  return $plugin_instance->$method($num, $score_builder, $mappings, $scores, @_);
}


sub has_plugin {
  my $self = shift;
  my $module = shift;

  defined($self->{'_plugins'}->{$module}) ? (return 1) : (return 0);
}


sub get_plugin {
  my $self = shift;
  my $module = shift;

  return $self->{'_plugins'}->{$module};
}


sub add_plugin {
  my $self = shift;
  my $plugin_instance = shift;

  $self->{'_plugins'}->{ref($plugin_instance)} = $plugin_instance;
}


sub log_ambiguous {
  my $self = shift;
  my $matrix = shift;
  my $type = shift;

  unless ($matrix and
          $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }
  
  # create dump directory if it doesn't exist
  my $debug_path = $self->conf->param('basedir').'/debug';
  unless (-d $debug_path) {
    system("mkdir -p $debug_path") == 0 or
      throw("Unable to create directory $debug_path.\n");
  }
  
  my $logfile = "$debug_path/ambiguous_${type}.txt";
  
  open(my $fh, '>', $logfile) or
    throw("Unable to open $logfile for writing: $!");

  my @low_scoring = ();
  my @high_scoring = ();
  my $last_id;

  # log by source
  foreach my $entry (sort { $a->source <=> $b->source }
                        @{ $matrix->get_all_Entries }) {
    
    $last_id ||= $entry->target;

    if ($last_id != $entry->source) {
      $self->write_ambiguous($type, 'source', $fh, \@low_scoring,
        \@high_scoring);
      $last_id = $entry->source;
    }
    
    if ($entry->score < 0.5) {
      push @low_scoring, $entry;
    } else {
      push @high_scoring, $entry;
    }
  }

  # write last source
  $self->write_ambiguous($type, 'source', $fh, \@low_scoring, \@high_scoring);

  # now do the same by target
  $last_id = undef;
  foreach my $entry (sort { $a->target <=> $b->target }
                        @{ $matrix->get_all_Entries }) {

    $last_id ||= $entry->target;

    if ($last_id != $entry->target) {
      $self->write_ambiguous($type, 'target', $fh, \@low_scoring,
        \@high_scoring);
      $last_id = $entry->target;
    }
    
    if ($entry->score < 0.5) {
      push @low_scoring, $entry;
    } else {
      push @high_scoring, $entry;
    }
  }

  # write last target
  $self->write_ambiguous($type, 'target', $fh, \@low_scoring, \@high_scoring);

  close($fh);
}


sub write_ambiguous {
  my ($self, $type, $db_type, $fh, $low, $high) = @_;
  
  # if only source or target are ambiguous (i.e. you have only one mapping from
  # this perspective) then log from the other perspective
  if (scalar(@$low) + scalar(@$high) <= 1) {
    @$low = ();
    @$high = ();
    return;
  }

  my $first_id;
  if (@$low) {
    $first_id = $low->[0]->$db_type;
  } else {
    $first_id = $high->[0]->$db_type;
  }

  my $other_db_type;
  if ($db_type eq 'source') {
    $other_db_type = 'target';
  } else {
    $other_db_type = 'source';
  }

  print $fh "$db_type $type $first_id scores ambiguously:\n";

  # high scorers
  if (@$high) {
    print $fh "  high scoring ${other_db_type}s\n";

    while (my $e = shift(@$high)) {
      print $fh "    ", $e->$other_db_type, " ", $e->score, "\n";
    }
  }

  # low scorers
  if (@$low) {
    print $fh "  low scoring ${other_db_type}s\n    ";

    my $i = 1;

    while (my $e = shift(@$low)) {
      print $fh "\n    " unless (($i++)%10);
      print $fh $e->$other_db_type, ", ";
    }
    print $fh "\n";
  }

  print $fh "\n";
}


1;

