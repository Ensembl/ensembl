
#
# Ensembl module for Bio::EnsEMBL::DBSQL::ExternalWrapper
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::ExternalWrapper - Makes a standard Ensembl database a ExternalFeatureFactoryI implementing object

=head1 SYNOPSIS

    # check out DB::ExternalFeatureFactoryI

=head1 DESCRIPTION

This class wraps a standard Ensembl database as if it is an
ExternalFeatureFactory interface, allowing it to serve up Genes and 
(assumming someone gets to write this as well) features.

The idea here is that this database will contain different data (eg, 
data on EMBL CDS from the original entry) which is updated at a different
cycle from the main stuff.

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::ExternalWrapper;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::RootI

use Bio::Root::RootI;
use Bio::EnsEMBL::DB::ExternalFeatureFactoryI;

@ISA = qw(Bio::EnsEMBL::DB::ExternalFeatureFactoryI Bio::Root::RootI);

sub new {
  my($class,$dbobj) = @_;
  

  my $self = {};
  bless $self,$class;


  if( !defined $dbobj || !$dbobj->isa('Bio::EnsEMBL::DBSQL::Obj') ) {
      $self->throw("No dbobj or not a dbobj [$dbobj]");
  }
  
  $self->dbobj($dbobj);

  return $self;
}


=head2 get_Ensembl_Genes_clone

 Title   : get_Ensembl_Genes_clone
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Ensembl_Genes_clone{
   my ($self,$cloneid) = @_;
   my $clone;
   eval {
       $clone = $self->dbobj->get_Clone($cloneid);
   };

   if( $@ ) {
       # return nothing
       return ();
   }

   return $clone->get_all_Genes();
}



=head2 dbobj

 Title   : dbobj
 Usage   : $obj->dbobj($newval)
 Function: 
 Returns : value of dbobj
 Args    : newvalue (optional)


=cut

sub dbobj{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'dbobj'} = $value;
    }
    return $obj->{'dbobj'};

}
