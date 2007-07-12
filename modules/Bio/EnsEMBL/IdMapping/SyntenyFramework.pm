package Bio::EnsEMBL::IdMapping::SyntenyFramework;

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS


=head1 REALTED MODULES


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


sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  # initialise internal datastructure
  $self->{'merged_syntenies'} = [];

  return $self;
}


sub build_synteny {
  my $self = shift;
  my $mappings = shift;
  
  unless ($mappings and
          $mappings->isa('Bio::EnsEMBL::IdMapping::MappingList')) {
    throw('Need a gene Bio::EnsEMBL::IdMapping::MappingList.');
  }

  # create a synteny region for each mapping
  my @synteny_regions = ();

  foreach my $entry (@{ $mappings->get_all_Entries }) {
    
    my $source_gene = $self->cache->get_by_key('genes_by_id', 'source',
      $entry->source);
    my $target_gene = $self->cache->get_by_key('genes_by_id', 'target',
      $entry->target);

    my $sr = Bio::EnsEMBL::IdMapping::SyntenyRegion->new_fast([
      $source_gene->start,
      $source_gene->end,
      $source_gene->strand,
      $source_gene->seq_region_name,
      $target_gene->start,
      $target_gene->end,
      $target_gene->strand,
      $target_gene->seq_region_name,
      $entry->score,
    ]);

    push @synteny_regions, $sr;
  }

  # sort synteny regions
  my @sorted = sort { $a->seq_region_name cmp $b->seq_region_name or
                      $a->start <=> $b->start } @synteny_regions;
  
  # now create merged regions from overlapping syntenies, but only merge a
  # maximum of 2 regions (otherwise you end up with large synteny blocks which
  # won't contain much information in this context)
  my $last_merged = 0;
  my $last_sr = shift(@sorted);

  while (my $sr = shift(@sorted)) {
  
    my $merged_sr = $last_sr->merge($sr);

    if (! $merged_sr) {
      unless ($last_merged) {
        $self->add_SyntenyRegion($last_sr->stretch(2));
        $last_merged = 0;
      }
    } else {
      $self->add_SyntenyRegion($merged_sr->stretch(2));
      $last_merged = 1;
    }

    $last_sr = $sr;
  }

  # deal with last synteny region in @sorted
  unless ($last_merged) {
    $self->add_SyntenyRegion($last_sr->stretch(2));
    $last_merged = 0;
  }

}


sub add_SyntenyRegion {
  push @{ $_[0]->{'merged_syntenies'} }, $_[1];
}


sub get_all_SyntenyRegions {
  return $_[0]->{'merged_syntenies'};
}


# 
# retain 70% of old score and build other 30% from synteny match
#
sub rescore_gene_matrix {
  my $self = shift;
  my $matrix = shift;

  unless ($matrix and
          $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  my $retain_factor = 0.7;

  foreach my $entry (@{ $matrix->get_all_Entries }) {

    my $source_gene = $self->cache->get_by_key('genes_by_id', 'source',
      $entry->source);

    my $target_gene = $self->cache->get_by_key('genes_by_id', 'target',
      $entry->target);

    my $highest_score = 0;

    foreach my $sr (@{ $self->get_all_SyntenyRegions }) {
      my $score = $sr->score_location_relationship($source_gene, $target_gene);
      $highest_score = $score if ($score > $highest_score);
    }

    $entry->score($entry_score * 0.7 + $highest_score * 0.3);
  }

}


1;

