#
# BioPerl module for Bio::EnsEMBL::TimDB::Contig
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::TimDB::Contig - Perl wrapper over Tim's directories for Contigs

=head1 SYNOPSIS

    $contig = Bio::EnsEMBL::TimDB::Contig->new();
    
    # $sf is a Bio::SeqFeatureI type object. $seq is a Bio::Seq object

    $contig->add_SeqFeature($sf);
    $contig->seq($seq); 

=head1 DESCRIPTION

Tim's contigs

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...
package Bio::EnsEMBL::TimDB::Contig;
use vars qw($AUTOLOAD @ISA);
use strict;
use Bio::EnsEMBL::DB::ContigI;

# Object preamble - inheriets from Bio::Root::Object
use Bio::Root::Object;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::ContigI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  
  my $make = $self->SUPER::_initialize;

  $self->{'_sf_array'} = [];
  $self->{'_gene_array'} = [];
 
  # set stuff in self from @args
  return $make; # success - we hope!
}


=head2 get_all_SeqFeatures

 Title   : get_all_SeqFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SeqFeatures{
   my ($self) = @_;

   $self->throw("Tim has not reimplemented this function");

   return @{$self->{'_sf_array'}};
}

=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Genes{
    my ($self) = @_;
   $self->throw("Tim has not reimplemented this function");
    return @{$self->{'_gene_array'}};
}


=head2 add_SeqFeature

 Title   : add_SeqFeature
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub add_SeqFeature{
   my ($self,$sf) = @_;
   $self->throw("SeqFeatures cannot be added in TimDB");
}


=head2 add_Gene

 Title   : add_Gene
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub add_Gene{
   my ($self,$gene) = @_;
   $self->throw("Genes cannot be added to TimDB");
}


=head2 offset

 Title   : offset
 Usage   : $self->offset($newval)
 Function: 
 Returns : value of offset
 Args    : newvalue (optional)

=cut

sub offset{
    my $self = shift;

    $self->throw("Tim has not reimplemented this function");
    
    if( @_ ) {
	my $value = shift;
	$self->{'offset'} = $value;
    }
    return $self->{'offset'};
}


=head2 orientation

 Title   : orientation
 Usage   : $self->orientation($newval)
 Function: 
 Returns : value of orientation
 Args    : newvalue (optional)

=cut

sub orientation{
    my $self = shift;
    
    $self->throw("Tim has not reimplemented this function");
    
    if( @_ ) {
	my $value = shift;
       $self->{'orientation'} = $value;
    }
    return $self->{'orientation'};
}


=head2 seq

 Title   : seq
 Usage   : $self->seq($newval)
 Function: 
 Returns : value of seq
 Args    : newvalue (optional)

=cut

sub seq{
    my $self = shift;
    
    $self->throw("Tim has not reimplemented this function");
    
    if( @_ ) {
	my $value = shift;
	if(! $value->isa("Bio::Seq") ) {
	    $self->throw("$value is not a Bio::Seq!");
	}			
	
	$self->{'seq'} = $value;
    }
    return $self->{'seq'};
    
}


=head2 id

 Title   : id
 Usage   : $self->id($newval)
 Function: 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id{
   my $self = shift;
   if( @_ ) {
       my $value = shift;
       $self->{'id'} = $value;
   }
   return $self->{'id'};
}

1;
