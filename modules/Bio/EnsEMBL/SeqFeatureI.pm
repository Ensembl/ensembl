
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

Bio::EnsEMBL::SeqFeatureI - Ensembl specific sequence feature interface.

=head1 SYNOPSIS

    my $feat = get_feature_from_somewhere;

    # Bio::SeqFeatureI methods can be used
    my $start = $feat->start;
    my $end   = $feat->end;

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


package Bio::EnsEMBL::SeqFeatureI;

use vars qw(@ISA);
use strict;
use Carp;

# Object preamble - inherits from Bio::Root::Object

use Bio::SeqFeatureI;

@ISA = qw(Bio::SeqFeatureI Exporter);


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

   $self->_abstractDeath;

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

   $self->_abstractDeath;
}


1;
