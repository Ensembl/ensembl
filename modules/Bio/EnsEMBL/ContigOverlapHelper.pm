
#
# BioPerl module for ContigOverlapHelper
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::ContigOverlapHelper - Describes the overlap between two contigs

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

The contig overlap object is directional - you pick it up from a particular
contig and it tells you what position on that contig is equivalent to what
position on the other contig.

http://www.ncbi.nlm.nih.gov/genome/seq/Hs_Data/contig.xml

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the
object methods. Internal methods are usually
preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::ContigOverlapHelper;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::RootI;


@ISA = qw(Bio::Root::RootI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub new {
  my($pkg,@args) = @_;

  my $self = bless {}, $pkg;

  my ($sister,$sisterpos,$sisterpolarity,$selfposition,$distance,$source) = $self->_rearrange([qw( SISTER
											 SISTERPOSITION
											 SISTERPOLARITY
											 SELFPOSITION
											 DISTANCE
											 SOURCE
											 )], @args);
  $self->throw("No sister input")                 unless defined($sister);
  $self->throw("No sister position input")        unless defined($sisterpos);
  $self->throw("No sister polarity input")        unless defined($sisterpolarity);
  $self->throw("No source input")                 unless defined($source);
  $self->throw("No self position input")          unless defined($selfposition);

  $self->sister         ($sister);
  $self->sister_position($sisterpos);
  $self->sister_polarity($sisterpolarity);
  $self->self_position  ($selfposition);
  $self->distance       ($distance);
  $self->source         ($source);
  
  return $self;
}

=head2 sister

 Title   : sister
 Usage   : $obj->sister($newval)
 Function: Returns the RawContigI implementing object of the sister
           contig. ie the id of the sister object is
           $obj->sister->id();
 Returns : value of sister
 Args    : newvalue (optional)


=cut

sub sister{
   my $obj = shift;
   if( @_ ) {
       my $value = shift;
       if( !ref $value || ! $value->isa('Bio::EnsEMBL::DB::RawContigI') ) {
	   $obj->throw("Value [$value] is not a RawContigI. Problemo...");
       }
       $obj->{'sister_contig'} = $value;
   }
   return $obj->{'sister_contig'};
   
}

=head2 sister_position
    
 Title   : sister_position
 Usage   : $obj->sister_position($newval)
 Function: 
 Returns : value of sister_position
 Args    : newvalue (optional)


=cut

sub sister_position{
   my $obj = shift;
   if( @_ ) {
       my $value = shift;
       $obj->{'sister_position'} = $value;
   }
   return $obj->{'sister_position'};
   
}

=head2 sister_polarity

 Title   : sister_polarity
 Usage   : $obj->sister_polarity($newval)
 Function: 
 Returns : value of sister_polarity
 Args    : newvalue (optional)


=cut

sub sister_polarity {
   my $obj = shift;
   if( @_ ) {
       my $value = shift;
       $obj->{'sister_contig_polarity'} = $value;
   }
   return $obj->{'sister_contig_polarity'};
   
}

=head2 self_position

 Title   : self_position
 Usage   : $obj->self_position($newval)
 Function: 
 Returns : value of self_position
 Args    : newvalue (optional)


=cut

sub self_position{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'self_contig_position'} = $value;
  }
   return $obj->{'self_contig_position'};
   
}

=head2 distance

 Title   : distance
 Usage   : $obj->distance($newval)
 Function: 
 Example : 
 Returns : value of distance
 Args    : newvalue (optional)


=cut

sub distance{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'distance'} = $value;
    }
    return $obj->{'distance'};

}


=head2 source

 Title   : source
 Usage   : $obj->source
 Function: 
 Example : 
 Returns : source tag for oeverlap
 Args    : newvalue (optional)


=cut

sub source {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{source} = $value;
    }
    return $obj->{source};

}



1;
