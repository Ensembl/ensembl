
#
# BioPerl module for Bio::EnsEMBL::SeqFeature
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::FeaturePairI - Ensembl specific feature pair interface.

=head1 SYNOPSIS

    # Bio::SeqFeatureI methods can be used
    my $start = $feat->start;
    my $end   = $feat->end;


    # Bio::SeqFeature::FeaturePair methods can be used
    my $hid = $feat->hid;

    # New Bio::EnsEMBL::SeqFeatureI specific methods are :
    my $analysis = $feat->analysis;

    # Validate all the data in the object
    $feat->validate  || $feat->throw("Invalid data in $feat");

=head1 DESCRIPTION

This is an extension of the bioperl Bio::SeqFeatureI interface.  Extra
methods are to store details of the analysis program/database/version used
to create this data and also a method to validate all data in the object is
present and of the right type.  This is useful before writing into
a relational database for example.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::FeaturePairI;

use vars qw(@ISA);
use strict;
use Carp;


use Bio::SeqFeatureI;
use Bio::EnsEMBL::SeqFeatureI;
use Bio::Root::RootI;

@ISA = qw(Bio::EnsEMBL::SeqFeatureI);


=head

Waiting for some definitions here...

=cut

=head2 percent_id

 Title   : percent_id
 Usage   : $percent_id = $featpair->percent_id
           $featpair->percent_id($pid)
 Function: Get/set on the percent_id of feature1
 Returns : integer
 Args    : none

=cut

sub percent_id {
    my ($self) = @_;
    $self->throw("Method not impelemented here");
}

=head2 hpercent_id

 Title   : hpercent_id
 Usage   : $percent_id = $featpair->hpercent_id
           $featpair->hpercent_id($pid)
 Function: Get/set on the percent_id of feature2
 Returns : integer
 Args    : none

=cut

sub hpercent_id {
    my ($self) = @_;
    $self->throw("Method not impelemented here");
}

=head2 p_value

 Title   : p_value
 Usage   : $p_value = $featpair->p_value
           $featpair->p_value($p_value)
 Function: Get/set on the p_value of feature1
 Returns : integer
 Args    : none

=cut

sub p_value {
    my ($self) = @_;
    $self->throw("Method not impelemented here");
}
=head2 hp_value

 Title   : hp_value
 Usage   : $p_value = $featpair->hp_value
           $featpair->hp_value($p_value)
 Function: Get/set on the p_value of feature2
 Returns : integer
 Args    : none

=cut

sub hp_value {
    my ($self) = @_;
    $self->throw("Method not impelemented here");
}

=head2 phase

 Title   : phase
 Usage   : $phase = $feat->phase()
           $feat->phase($phase)
 Function: get/set on start phase of predicted feature1
 Returns : [0,1,2]
 Args    : none if get, 0,1 or 2 if set. 

=cut

sub phase {
    my ($self) = @_;
    $self->throw("Method not impelemented here");
}

=head2 end_phase

 Title   : end_phase
 Usage   : $end_phase = $feat->end_phase()
           $feat->end_phase($end_phase)
 Function: get/set on end phase of predicted feature1
 Returns : [0,1,2]
 Args    : none if get, 0,1 or 2 if set. 

=cut

sub end_phase {
    my ($self) = @_;
    $self->throw("Method not impelemented here");
}

=head2 hphase

 Title   : hphase
 Usage   : $hphase = $fp->hphase()
           $fp->hphase($hphase)
 Function: get/set on start hphase of predicted feature2
 Returns : [0,1,2]
 Args    : none if get, 0,1 or 2 if set. 

=cut

sub hphase {
    my ($self) = @_;
    $self->throw("Method not impelemented here");
}

=head2 hend_phase

 Title   : hend_phase
 Usage   : $hend_phase = $feat->hend_phase()
           $feat->hend_phase($hend_phase)
 Function: get/set on end phase of predicted feature2
 Returns : [0,1,2]
 Args    : none if get, 0,1 or 2 if set. 

=cut

sub hend_phase {
    my ($self) = @_;
    $self->throw("Method not impelemented here");
}


#
# Bioperl 0.7 compliance
#

sub location {
    my ($self) = @_;
    return $self;
}

sub to_FTstring {
    my ($self) = @_;
    
    if( $self->strand == 1 ) {
	return $self->start."..".$self->end;
    } else {
	return "complement(".$self->start."..".$self->end.")";
    }
}

    
1;
