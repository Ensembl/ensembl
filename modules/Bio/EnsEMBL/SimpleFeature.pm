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

Bio::EnsEMBL::SimpleFeature - A simple feature with a location and label

=head1 SYNOPSIS

  use Bio::EnsEMBL::SimpleFeature;

  $feature = Bio::EnsEMBL::SimpleFeature->new(
    -start         => 100,
    -end           => 220,
    -strand        => -1,
    -slice         => $slice,
    -analysis      => $analysis,
    -score         => 58,
    -display_label => 'EponineTSS',
    -dbID          => 1230,
    -adaptor       => $adaptor
  );

=head1 DESCRIPTION

This is a simple feature which extends the Feature class to add
display_label and score attributes.

=head1 METHODS

=cut

use strict;

package Bio::EnsEMBL::SimpleFeature;

use vars qw(@ISA);

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Scalar::Util qw(weaken isweak);

@ISA = qw(Bio::EnsEMBL::Feature);


=head2 new

  Arg [DISPLAY_LABEL]: The label assigned to this simple feature
  Arg [...]  : Named arguments passed to superclass
  Example    : $feature = Bio::EnsEMBL::SimpleFeature->new
                        (-start    => 1,
                         -end      => 100,
                         -strand   => 1,
                         -slice    => $slice,
                         -analysis => $analysis,
                         -adaptor => $adaptor,
                         -dbID    => 10,
                         -display_label => 'EponineTSS',
                         -score => 100);
  Description: Constructs a new Bio::EnsEMBL::Feature.  Generally subclasses
               of this method are instantiated, rather than this class itself.
  Returntype : Bio::EnsEMBL::Feature
  Exceptions : Thrown on invalid -SLICE, -ANALYSIS, -STRAND arguments
  Caller     : general, subclass constructors
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  #allow this to be called as class or object method
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($display_label, $score) = rearrange(['DISPLAY_LABEL','SCORE'],@_);

  $self->{'display_label'} = $display_label;
  $self->{'score'} = $score;

  return $self;
}


=head2 display_label

  Arg [1]    : (optional) string $value
  Example    : $label = $simple_feature->display_label();
  Description: Getter/Setter for the display label associated with this
               feature.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub display_label{
  my $self = shift;

  $self->{'display_label'} = shift if(@_);

  return $self->{'display_label'};
}


=head2 display_id

  Arg [1]    : none
  Example    : print $rf->display_id();
  Description: This method returns a string that is considered to be
               the 'display' identifier.  For simple features this is the 
               display_label if it is available otherwise it is an empty 
               string.
  Returntype : string
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub display_id {
  my $self = shift;
  return $self->{'display_label'} || '';
}



=head2 score

  Arg [1]    : (optional) string $value
  Example    : $score = $simple_feature->score();
  Description: Getter/Setter for the score associated with this
               feature.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub score {
  my $self = shift;
  $self->{'score'} = shift if(@_);
  return $self->{'score'};
}

=head2 summary_as_hash

  Example    : my $hash = $simple_feature->summary_as_hash();
  Description: Generates a HashRef compatible with GFFSerializer. Adds
               score, external_name and logic_name to the basic Feature
               hash
  Returntype : Hash
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub summary_as_hash {
  my ($self) = @_;
  my $hash = $self->SUPER::summary_as_hash();
  $hash->{score} = $self->score() if $self->score();
  $hash->{'external_name'} = $self->display_label() if $self->display_label();
  $hash->{'logic_name'} = $self->analysis->logic_name();
  delete $hash->{id};
  return $hash;
}


1;
