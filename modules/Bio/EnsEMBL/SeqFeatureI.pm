=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
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

   deprecate('SeqFeatureI->analysis is deprecated and will be removed in e84');

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

   deprecate('SeqFeatureI->validate is deprecated and will be removed in e84');

}


sub id {
    my ($self,$value) = @_;

    deprecate('SeqFeatureI->id is deprecated and will be removed in e84');
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

    deprecate('SeqFeatureI->percent_id is deprecated and will be removed in e84');
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

    deprecate('SeqFeatureI->percent_id is deprecated and will be removed in e84');
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
    deprecate('SeqFeatureI->phase is deprecated and will be removed in e84');
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
    deprecate('SeqFeatureI->end_phase is deprecated and will be removed in e84');
}


# this is a bit too sneaky. 
sub location {
    my ($self)= @_;
    deprecate('SeqFeatureI->location is deprecated and will be removed in e84');
    return $self;
}


1;
