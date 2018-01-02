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

Bio::EnsEMBL::SubSlicedFeature

=head1 SYNOPSIS

my $truncated_gene = Bio::EnsEMBL::Utils::SubSlicedFeature->new(
    -start => 300, 
    -end => 10000, 
    -feature => $gene);
my $transcripts = $truncated_gene->get_all_Transcripts();
# list of transcripts is limited to those within the coordinates, rather
# than the original feature Slice. 

=head1 DESCRIPTION

  Alters the behaviour of a normal Feature object to act within a user-specified 
  sub-slice of of its boundaries. As it stands, this only affects get_all_*
  methods, meaning that seq_region_start() and project() will work on original
  coordinates.

=cut

package Bio::EnsEMBL::SubSlicedFeature;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw/rearrange/;
use base qw/Bio::EnsEMBL::Utils::Proxy/;

sub new {
  my ($class, @args) = @_;
  my ($start,$end,$original_feature) = rearrange([qw/start end feature/], @args);
  my $self = $class->SUPER::new($original_feature);
  $self->{'start'} = $start;
  $self->{'end'} = $end;
  return $self;
}

# Required by Proxy to control scope of the autoloaded methods.
# Also intercepts calls to get methods that would access Slice

sub __resolver {
  my ($self, $package_name, $method) = @_;
  
  if ($method =~ /^get_all_/) {
      # Call original method and filter results to Proxy coordinates
      return sub {
          my ($local_self, @args) = @_;
          my $feature_list = $local_self->__proxy()->$method(@args);
          my @short_list;
          foreach my $feature (@$feature_list) {
              if ($feature->start > $local_self->{'start'}
                && $feature->end < $local_self->{'end'}) {
                    push @short_list,$feature;
                }
          }
          return \@short_list;
      }
  } else {
      # No intervention required, call original object method
      return sub {
        my ($local_self, @args) = @_;
        return $local_self->__proxy()->$method(@args);
      };
  }
  
}

1;
