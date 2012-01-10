=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
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

Bio::EnsEMBL::SeqFeatureI

=head1 DESCRIPTION

Do not use this class. It is deprecated and has been replaced by 
Bio::EnsEMBL::Feature.

=head1 METHODS

=cut


# Let the code begin...


package Bio::EnsEMBL::SeqFeatureI;

use vars qw(@ISA);
use strict;
use Carp;

# Object preamble - inherits from Bio::Root::Object

use Bio::SeqFeatureI;
use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Root Bio::SeqFeatureI);


=head1 Abstract methods

These methods must be implemented in all subclasses.


=head2 analysis

 Title   : analysis
 Usage   : $sf->analysis();
 Function: Store details of the program/database
           and versions used to create this feature.
           
 Example :
 Returns : 
 Args    :


=cut

sub analysis {
   my ($self,$value) = @_;

   $self->throw("Have not implemeneted analysis");

}

=head2 validate

 Title   : validate
 Usage   : $sf->validate;
 Function: Checks whether all the data is present in the
           object.
 Example :
 Returns : 
 Args    :


=cut

sub validate {
   my ($self,$value) = @_;

   $self->throw("Have not implemeneted validate");


}


sub id {
    my ($self,$value) = @_;

    $self->throw("Have not implemented id");
}

=head2 percent_id

 Title   : percent_id
 Usage   : $pid = $feat->percent_id()
           $feat->percent_id($pid)
 Function: get/set on percentage identity information
 Returns : float
 Args    : none if get, the new value if set

=cut

sub percent_id {
    my ($self) = @_;
    $self->throw("percent_id() not yet implemented");
}

=head2 e_value

 Title   : p_value
 Usage   : $p_val = $feat->p_value()
           $feat->p_value($p_val)
 Function: get/set on p value information
 Returns : float
 Args    : none if get, the new value if set

=cut

sub e_value {
    my ($self) = @_;
    $self->throw("e value() not yet implemented");
}

=head2 phase

 Title   : phase
 Usage   : $phase = $feat->phase()
           $feat->phase($phase)
 Function: get/set on start phase of predicted exon feature
 Returns : [0,1,2]
 Args    : none if get, 0,1 or 2 if set. 

=cut

sub phase {
    my ($self) = @_;
    $self->throw("phase() not yet implemented");
}

=head2 end_phase

 Title   : end_phase
 Usage   : $end_phase = $feat->end_phase()
           $feat->end_phase($end_phase)
 Function: get/set on end phase of predicted exon feature
 Returns : [0,1,2]
 Args    : none if get, 0,1 or 2 if set. 

=cut

sub end_phase {
    my ($self) = @_;
    $self->throw("end_phase() not yet implemented");
}


# this is a bit too sneaky. 
sub location {
    my ($self)= @_;
    return $self;
}


1;
