
#
# BioPerl module for Bio::EnsEMBL::PerlDB::Clone
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::PerlDB::Clone - Pure Perl implementation of Clone

=head1 SYNOPSIS

    $clone = Bio::EnsEMBL::PerlDB::Clone->new();
 
    $clone->add_Contig($contig);
    

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::PerlDB::Clone;
use vars qw($AUTOLOAD @ISA);
use strict;
use Bio::EnsEMBL::DB::CloneI;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;


@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::CloneI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  $self->{'_contig_hash'} = {};
# set stuff in self from @args
 return $make; # success - we hope!
}

=head2 get_Contig

 Title   : get_Contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Contig{
   my ($self,$id) = @_;

   return $self->{'_contig_hash'}->{$id};
}


=head2 get_all_Contigs

 Title   : get_all_Contigs
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Contigs{
   my ($self) = @_;

   return values %{$self->{'_contig_hash'}};
}

=head2 add_Contig

 Title   : add_Contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Contig{
   my ($self,$contig) = @_;


   if( ! $contig->isa('Bio::EnsEMBL::DB::ContigI') ) {
       $self->warn("$contig is not a contigI object...");
   }

   $self->{'_contig_hash'}->{$contig->id()} = $contig;
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
   my %h;

   # read into a hash to make unique
   foreach my $contig ( $self->get_all_Contigs ) {
       foreach my $gene ( $contig->get_all_Genes ) {
	   $h{$gene->id()} = $gene;
       }

   }

   return values %h;

}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Example : 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id{
   my ($obj,$value) = @_;
   if( defined $value) {
       $obj->{'id'} = $value;
   }
   return $obj->{'id'};
   
}


=head2 embl_id

 Title   : embl_id
 Usage   : $obj->embl_id($newval)
 Function: 
 Returns : value of embl_id
 Args    : newvalue (optional)


=cut

sub embl_id{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'embl_id'} = $value;
    }
    return $obj->{'embl_id'};

}

=head2 embl_version

 Title   : embl_version
 Usage   : $obj->embl_version($newval)
 Function: 
 Returns : value of embl_version
 Args    : newvalue (optional)


=cut

sub embl_version{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'embl_version'} = $value;
    }
    return $obj->{'embl_version'};

}

=head2 seq_date

 Title   : seq_date
 Usage   : $obj->seq_date($newval)
 Function: 
 Returns : value of seq_date
 Args    : newvalue (optional)


=cut

sub seq_date{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'seq_date'} = $value;
    }
    return $obj->{'seq_date'};

}

=head2 version

 Title   : version
 Usage   : $obj->version($newval)
 Function: 
 Returns : value of version
 Args    : newvalue (optional)


=cut

sub version{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'version'} = $value;
    }
    return $obj->{'version'};

}
=head2 htg_phase

 Title   : htg_phase
 Usage   : $obj->htg_phase($newval)
 Function: 
 Returns : value of htg_phase
 Args    : newvalue (optional)


=cut

sub htg_phase{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'htg_phase'} = $value;
    }
    return $obj->{'htg_phase'};

}

1;









