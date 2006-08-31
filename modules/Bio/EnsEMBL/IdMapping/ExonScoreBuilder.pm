package Bio::EnsEMBL::IdMapping::ExonScoreBuilder;

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

use Bio::EnsEMBL::IdMapping::ScoreBuilder;
our @ISA = qw(Bio::EnsEMBL::IdMapping::ScoreBuilder);

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ScriptUtils qw(parse_bytes);
use Bio::EnsEMBL::IdMapping::ScoredMappingMatrix;


sub score_exons {
  my $self = shift;
  
  my $matrix = Bio::EnsEMBL::IdMapping::ScoredMappingMatrix->new(
    -DUMP_PATH   => $self->conf('dumppath')
  );

  my $cache_file = $matrix->cache_file;

  if (-s $cache_file) {
    # read from file
    $self->logger->log_stamped("Reading exon scoring matrix from file...\n");
    $self->logger->log("Cache file $cache_file.\n", 1);
    $matrix->read_from_file;
    $self->logger->log_stamped("Done.\n");
  } else {
    #
    # build scoring matrix
    #

    # direct mapping (by overlap, if common coord_system exists)
    if ($self->cache->highest_common_cs) {
      $matrix = $self->build_overlap_scores($matrix);
    }


    # map the remaining exons using exonerate
    
    # dump exons to fasta files
    my $dump_size = $self->cache->dump_filtered_exons;

    if ($dump_size) {
      # run exonerate
      my $exonerate_matrix;

      # merge matrices
      $matrix->merge($exonerate_matrix);
    }


    #
    # write scoring matrix to file
    #
    $matrix->write_to_file;

    return $matrix;
  }

  return $matrix;
}


sub build_overlap_scores {
  my $self = shift;
  my $matrix = shift;

  unless ($matrix and
          $matrix->isa('Bio::EnsEMBL::IdMapping::ScoredMappingMatrix')) {
    throw('Need a Bio::EnsEMBL::IdMapping::ScoredMappingMatrix.');
  }

  my @s_exons = $self->sort_exons(
    [values %{ $self->cache->get_by_name('exons_by_id', 'source') }
  );
  my @t_exons = $self->sort_exons(
    [values %{ $self->cache->get_by_name('exons_by_id', 'target') }
  );

  
}


#
# Return a list of exons, sorted by seq_region_name, then location (where
# location is either start-1 or end, so each exon is in the list twice)
#
# TODO: this implementation isn't good enough, since you'll need location later
# to compare source and target exons again. best add this as an instance
# variable to the TinyExon returned
#
sub sort_exons {
  my $self = shift;
  my $exons = shift;

  return
    map { $_->[2] }
      sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] }
        ( 
          map { [$_->common_name, $_->common_start - 1, $_] } @$exons,
          map { [$_->common_name, $_->common_end, $_] } @$exons
        );

}


sub compare_exons {
  my $self = shift;
  my $e1 = shift;
  my $e2 = shift;

  return ( $e1->
}

