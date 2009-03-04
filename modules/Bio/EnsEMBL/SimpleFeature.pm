=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

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


=head2 new_fast

  Arg [1]    : hashref to be blessed
  Description: Construct a new Bio::EnsEMBL::Feature using the hashref.
  Exceptions : none
  Returntype : Bio::EnsEMBL::Feature
  Caller     : general, subclass constructors
  Status     : Stable

=cut


sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
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


1;
