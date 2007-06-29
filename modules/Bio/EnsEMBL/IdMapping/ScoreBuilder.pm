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

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);


=head2 new

  Arg[1]      : 
  Example     : 
  Description : constructor
  Return type : 
  Exceptions  : 
  Caller      : general

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($logger, $conf, $cache) = rearrange(['LOGGER', 'CONF', 'CACHE'], @_);

  unless ($logger->isa('Bio::EnsEMBL::Utils::Logger')) {
    throw("You must provide a Bio::EnsEMBL::Utils::Logger for logging.");
  }
  
  unless ($conf->isa('Bio::EnsEMBL::Utils::ConfParser')) {
    throw("You must provide configuration as a Bio::EnsEMBL::Utils::ConfParser object.");
  }
  
  unless ($cache->isa('Bio::EnsEMBL::IdMapping::Cache')) {
    throw("You must provide configuration as a Bio::EnsEMBL::IdMapping::Cache object.");
  }
  
  my $self = {};
  bless ($self, $class);

  # initialise
  $self->logger($logger);
  $self->conf($conf);
  $self->cache($cache);
  
  return $self;
}


sub logger {
  my $self = shift;
  $self->{'_logger'} = shift if (@_);
  return $self->{'_logger'};
}


sub conf {
  my $self = shift;
  $self->{'_conf'} = shift if (@_);
  return $self->{'_conf'};
}


sub cache {
  my $self = shift;
  $self->{'_cache'} = shift if (@_);
  return $self->{'_cache'};
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
  
  $self->logger->info(sprintf($fmt2, "Average score:",
    $matrix->get_average_score), 1);
  
  my ($min, $max) = @{ $matrix->get_min_max_scores };
  $self->logger->info(sprintf($fmt2, "Min. score:", $min), 1);
  $self->logger->info(sprintf($fmt2, "Max. score:", $max), 1);
}


1;

