
#
# BioPerl module for Protein.pm
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Protein.pm - DESCRIPTION of Object

=head1 SYNOPSIS

my $protein = Bio::EnsEMBL::Protein->new ( -seq =>$sequence,
					      -accession_number  => $id,
					      -display_id => $id,
					      -primary_id => $id,
					      -id => $id,
					      -desc => $desc,
					      -moltype => $moltype,
					      );

=head1 DESCRIPTION

This object inherit from Bio::SeqI for the implementation and from PrimarySeq
A protein object hold the basic PrimarySeq requirements, hold annotation, DBlink, protein_features object.
It will also soon hold family and SNPs objects. The DB adaptor for this object is called Protein_Adaptor 

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Protein;

use vars qw(@ISA);
use strict;
use Bio::Root::Object;
use Bio::SeqI;
use Bio::Annotation::DBLink;
use Bio::EnsEMBL::Transcript;
use Bio::DBLinkContainerI;
use Bio::SeqIO;

# Object preamble - inheriets from Bio::Root::Object

@ISA = qw(Bio::PrimarySeq Bio::SeqI);




=head2 annotation

 Title   : annotation
 Usage   : $obj->annotation($newval)
 Function: 
 Returns : value of annotation
 Args    : newvalue (optional)


=cut

sub annotation{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'annotation'} = $value;
    }
    return $obj->{'annotation'};

}



=head2 top_SeqFeatures

 Title   : top_SeqFeatures
 Usage   : my @top_SeqFeatures = $protein->top_SeqFeatures
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub top_SeqFeatures{
   my ($self) = @_;
   
   my @seq_features = $self->each_Protein_feature;

   if (@seq_features) {
       return @seq_features;
   }
}

=head2 all_SeqFeature

 Title   : all_SeqFeature
 Usage   :
 Function:This method returns the same things than top_SeqFeature, sub SeqFeatures are not currently implemented in the protein object 
 Example :
 Returns : 
 Args    :


=cut

sub all_SeqFeature{
    
    my ($self) = @_;
   
    my @seq_feature = $self->each_Protein_feature;
    if (@seq_feature) {
	return @seq_feature;
    }

}


=head2 family

 Title   : family
 Usage   : $obj->family($newval)
 Function: This method contains the family object (Not used yet)
 Returns : value of family
 Args    : newvalue (optional)


=cut

sub family{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'family'} = $value;
    }
    return $obj->{'family'};

}


=head2 each_Protein_feature

 Title   : each_Protein_feature
 Usage   :my @features = $protein->each_Protein_feature
 Function:Retrieve an array of protein features (FeaturePair objects)
 Example :
 Returns : FeaturePair objects
 Args    :none


=cut

sub each_Protein_feature{
   my ($self,@args) = @_;

   if (defined ($self->{'_prot_feat'})) {
       return @{$self->{'_prot_feat'}};
   }

}


=head2 add_Protein_feature

 Title   : add_Protein_feature
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Protein_feature{
   my ($self,$value) = @_;

   if (!defined $value) {
     $self->throw("This [$value] is not a Protein feature");
   }

   push(@{$self->{'_prot_feat'}},$value);   

}


=head2 length

 Title   : length
 Usage   :my $lenght = $protein->length
 Function:
 Example :
 Returns :length of the aa sequence 
 Args    :none


=cut

sub length{
   my ($self) = @_;
   my $length = length($self->seq);


   return $length;

}

=head2 add_date

 Title   : add_date
 Usage   : $self->add_date($ref)
 Function: adds a date
 Example :
 Returns : 
 Args    :

=cut

sub add_date {
   my ($self) = shift;
   foreach my $dt ( @_ ) {
       push(@{$self->{'date'}},$dt);
   }
}

=head2 get_dates

 Title   : get_dates
 Usage   : foreach $dt ( $self->get_dates() )
 Function: gets an array of dates
 Example :
 Returns : 
 Args    :

=cut

sub get_dates {
   my ($self) = @_;
   return @{$self->{'date'}}; 
}

=head2 species

 Title   : species
 Usage   : $obj->species($newval)
 Function: 
 Returns : value of species
 Args    : newvalue (optional)


=cut

sub species{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'species'} = $value;
    }
    return $obj->{'species'};

}

=head2 primary_seq

 Title   : primary_seq
 Usage   : $obj->primary_seq($newval)
 Function: 
 Returns : value of primary_seq
 Args    : newvalue (optional)


=cut

sub primary_seq{
    my ($self) = @_;   
    return $self;
}




=head2 molecule

 Title   : molecule
 Usage   : $obj->molecule($newval)
 Function: 
 Returns : value of molecule
 Args    : newvalue (optional)


=cut

sub molecule{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'molecule'} = $value;
    }
    return $obj->{'molecule'};

}

=head2 moltype

 Title   : moltype
 Usage   : $obj->moltype($newval)
 Function: 
 Returns : value of moltype
 Args    : newvalue (optional)


=cut

sub moltype{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'moltype'} = $value;
    }
    return $obj->{'moltype'};

}


1;









