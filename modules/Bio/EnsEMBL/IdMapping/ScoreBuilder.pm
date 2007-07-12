package Bio::EnsEMBL::IdMapping::ScoreBuilder;

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS


=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

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
use Bio::EnsEMBL::IdMapping::ScoredMappingMatrix;

#
# create a shrinked matrix which doesn't contain entries which were already 
# mapped
#
sub create_shrinked_matrix {
  my $self = shift;
  my $matrix = shift;
  my $mappings = shift;
  my $cache_file = shift;

  # argument checks
  unless ($matrix and
          $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }
  
  unless ($mappings and
          $mappings->isa('Bio::EnsEMBL::IdMapping::MappingList')) {
    throw('Need a gene Bio::EnsEMBL::IdMapping::MappingList.');
  }

  throw('Need a cache file name.') unless ($cache_file);

  my $shrinked_matrix = Bio::EnsEMBL::IdMapping::ScoredMappingMatrix->new(
    -DUMP_PATH   => $self->conf->param('dumppath'),
    -CACHE_FILE  => $cache_file,
  );

  # create lookup hashes for sources and targets in the MappingList
  my %sources = ();
  my %targets = ();

  foreach my $entry (@{ $mappings->get_all_Entries }) {
    $sources{$entry->source} = 1;
    $targets{$entry->target} = 1;
  }

  # add all entries to shrinked matrix which are not in the MappingList
  foreach my $entry (@{ $matrix->get_all_Entries }) {
    unless ($sources{$entry->source} or $targets{$entry->target}) {
      $shrinked_matrix->add_Entry($entry);
    }
  }

  # log shrinking stats
  $self->logger->info('Sources '.$matrix->get_source_count.' --> '.
    $shrinked_matrix->get_source_count."\n");
  $self->logger->info('Targets '.$matrix->get_target_count.' --> '.
    $shrinked_matrix->get_target_count."\n");
  $self->logger->info('Entries '.$matrix->get_entry_count.' --> '.
    $shrinked_matrix->get_entry_count."\n");

  return $shrinked_matrix;
}


sub log_matrix_stats {
  my $self = shift;
  my $matrix = shift;

  unless ($matrix and
          $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('You must provide a ScoredMappingMatrix.');
  }

  my $fmt1 = "%-40s%10.0f\n";
  my $fmt2 = "%-40s%10.2f\n";
  
  $self->logger->info(sprintf($fmt1, "Scoring matrix entries:",
    $matrix->get_entry_count), 1);
  
  $self->logger->info(sprintf($fmt1, "Scoring matrix sources:",
    $matrix->get_source_count), 1);
  
  $self->logger->info(sprintf($fmt1, "Scoring matrix targets:",
    $matrix->get_target_count), 1);
  
  $self->logger->info(sprintf($fmt2, "Average score:",
    $matrix->get_average_score), 1);
  
  my ($min, $max) = @{ $matrix->get_min_max_scores };
  $self->logger->info(sprintf($fmt2, "Min. score:", $min), 1);
  $self->logger->info(sprintf($fmt2, "Max. score:", $max), 1);
}


1;

