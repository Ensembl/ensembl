#
# Ensembl module for Bio::EnsEMBL::RegulatorySearchRegion
#
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::RegulatorySearchRegion - A feature representing a part of the genome that was searched for regulatory features.

=head1 SYNOPSIS

use Bio::EnsEMBL::RegulatorySearchRegion;

$region  = Bio::EnsEMBL::RegulatorySearchRegion->new(-start            	  => 100,
                                                     -end              	  => 220,
                                                     -strand           	  => -1,
                                                     -slice            	  => $slice,
                                                     -analysis         	  => $analysis,
                                                     -ensembl_object_type => "Gene",
                                                     -ensembl_object_id   => 12340,
                                                     -dbID                => 1230,
                                                     -adaptor             => $adaptor);

=head1 DESCRIPTION

A feature representing a part of the genome that was searched for regulatory features.

=head1 AUTHOR - Glenn Proctor

This module is part of the Ensembl project http://www.ensembl.org

=head1 CONTACT

Post comments/questions to the ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;

package Bio::EnsEMBL::RegulatorySearchRegion;

use vars qw(@ISA);

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw(Bio::EnsEMBL::Feature);


=head2 new

  Arg [...]  : Named arguments passed to superclass
  Example    : $feature = Bio::EnsEMBL::RegulatoryFeatureSearchRegion->new
                                               (-name                => 'CisRed_search_1',
						-start               => 100,
                                                -end                 => 220,
                                                -strand              => -1,
                                                -slice               => $slice,
                                                -analysis            => $analysis,
                                                -ensembl_object_type => "Gene",
                                                -ensembl_object_id   => 12340,
                                                -dbID                => 1230,
                                                -adaptor             => $adaptor);
  Description: Constructs a new Bio::EnsEMBL::RegulatorySearchRegion.
  Exceptions : Thrown on invalid -SLICE, -ANALYSIS, -STRAND arguments
  Caller     : general, subclass constructors
  Status     : At Risk
             : under development

=cut

sub new {
  my $caller = shift;

  #allow this to be called as class or object method
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  my ($name, $factor, $ensembl_object_type, $ensembl_object_id) = rearrange(['NAME', 'ENSEMBL_OBJECT_TYPE', 'ENSEMBL_OBJECT_ID'],@_);

  $self->{'name'} = $name;
  $self->{'ensembl_object_type'} = $ensembl_object_type;
  $self->{'ensembl_object_id'} = $ensembl_object_id;

  return $self;
}


sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}

=head2 name

  Arg [1]    : none
  Example    : print $rsr->name();
  Description: Getter/setter for this search region's name.
  Returntype : string
  Exceptions : none
  Caller     : web drawing code
  Status     : At Risk
             : under development

=cut

sub name {
  my $self = shift;
  $self->{'name'} = shift if(@_);
  return $self->{'name'};
}

=head2 ensembl_object_type

  Arg [1]    : none
  Example    : print $rsr->ensembl_object_type();
  Description: Getter/setter for this search region's ensembl_object_type.
  Returntype : string
  Exceptions : none
  Caller     : web drawing code
  Status     : At Risk
             : under development

=cut

sub ensembl_object_type {
  my $self = shift;
  $self->{'ensembl_object_type'} = shift if(@_);
  return $self->{'ensembl_object_type'};
}

=head2 ensembl_object_id

  Arg [1]    : none
  Example    : print $rsr->ensembl_object_id();
  Description: Getter/setter for this search region's ensembl_object_id.
  Returntype : string
  Exceptions : none
  Caller     : web drawing code
  Status     : At Risk
             : under development

=cut

sub ensembl_object_id {
  my $self = shift;
  $self->{'ensembl_object_id'} = shift if(@_);
  return $self->{'ensembl_object_id'};
}

1;
