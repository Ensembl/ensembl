
#
# Ensembl module for Bio::EnsEMBL::DB::WebExternalFeatureFactoryI
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DB::WebExternalFeatureFactoryI - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DB::WebExternalFeatureFactoryI;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::RootI

use Bio::EnsEMBL::DB::ExternalFeatureFactoryI;

          
@ISA = qw(Bio::EnsEMBL::DB::ExternalFeatureFactoryI);



=head2 get_Ensembl_SeqFeatures_clone_web

 Title   : get_Ensembl_SeqFeatures_clone_web
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Ensembl_SeqFeatures_clone_web{
   my ($self,@args) = @_;

   $self->warn("Object $self did not implement get_Ensembl_SeqFeatures_clone_web depiste having the correct interface!");
   return ();
}

1;
