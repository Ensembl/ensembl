
#
# BioPerl module for ProteinFeature
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

ProteinFeature.pm - DESCRIPTION of Object

=head1 SYNOPSIS

my $feature = new Bio::EnsEMBL::ProteinFeature(-feature1 => $feat1,
					       -feature2 => $feat2,);

=head1 DESCRIPTION

This object inherits from Bio::EnsEMBL::FeaturePair. This extension has been implemented to work with the Protein object. Each Protein Feature should be stored in a Protein_FeaturePair object.

This object was formerly named Protein_FeaturePair.

=head1 CONTACT

mongin@ebi.ac.uk

=cut

package Bio::EnsEMBL::ProteinFeature;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::FeaturePair;

@ISA = qw(Bio::EnsEMBL::FeaturePair);

=head2 idesc

 Title   : idesc
 Usage   : $obj->idesc($newval)
 Function: 
 Returns : value of idesc
 Args    : newvalue (optional)


=cut

sub idesc{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'idesc'} = $value;
    }
    return $obj->{'idesc'};

}

=head2 interpro_ac

 Title   : interpro_ac
 Usage   : $obj->interpro_ac($newval)
 Function: 
 Returns : value of interpro_ac
 Args    : newvalue (optional)


=cut

sub interpro_ac{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'interpro_ac'} = $value;
    }
    return $obj->{'interpro_ac'};

}


1;
