#
# Ensembl module for Bio::EnsEMBL::SimpleFeature
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::SimpleFeature - A simple feature with a location and label

=head1 SYNOPSIS

use Bio::EnsEMBL::SimpleFeature;

$feature = Bio::EnsEMBL::SimpleFeature->new(-start    => 100,
                                            -end      => 220,
                                            -strand   => -1,
                                            -slice    => $slice,
                                            -analysis => $analysis,
                                            -display_label => 'EponineTSS',
                                            -dbID     => 1230,
                                            -adaptor  => $adaptor);

=head1 DESCRIPTION

This is a simple feature which extends the Feature class to add a display_label
attribute

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a _

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

=cut

sub score {
  my $self = shift;
  $self->{'score'} = shift if(@_);
  return $self->{'score'};
}


1;
