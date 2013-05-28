=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

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