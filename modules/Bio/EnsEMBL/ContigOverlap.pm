
#
# BioPerl module for ContigOverlap
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::ContigOverlap - Describes the overlap between two contigs

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

The contig overlap object is symmetrical.  It
describes the overlap between two  contigs.

http://www.ncbi.nlm.nih.gov/genome/seq/Hs_Data/contig.xml

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::ContigOverlap;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;


@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);

  my ($contiga,$contigb,$positiona,$positionb,$overlap_type,$source,$distance) 
      = $self->_rearrange([qw( CONTIGA
			       CONTIGB
			       POSITIONA
			       POSITIONB
			       OVERLAP_TYPE
			       SOURCE
			       DISTANCE
                               )], @args);
  unless ($contiga and $contigb
          and defined($positiona) and defined($positionb)
          and $overlap_type and $source) {
      $self->throw("You have to construct ContigOverlap objects with all six arguments (contiga,contigb,positiona,positionb,overlap_type,source) "
        ."Only got :". join(',', map "'$_'", ($contiga,$contigb,$positiona,$positionb,$overlap_type,$source)));
  }

  $self->contiga     ($contiga);
  $self->contigb     ($contigb);
  $self->positiona   ($positiona);
  $self->positionb   ($positionb);
  $self->overlap_type($overlap_type);
  $self->source      ($source);
  $self->distance    ($distance);

  # set stuff in self from @args
  return $make; # success - we hope!
}

=head2 contiga

 Title   : contiga
 Usage   : $obj->contiga($contig)
 Function: Get/Set for the 1st contig in the overlap
 Returns : Bio::EnsEMBL::DB::RawContigI
 Args    : newvalue (optional)


=cut

sub contiga {
   my $obj = shift;

   if( @_ ) {
       my $value = shift;
       if( !ref $value || ! $value->isa('Bio::EnsEMBL::DB::ContigI') ) {
	   $obj->throw("Value [$value] is not a RawContigI. Problemo...");
       }
       $obj->{'contiga'} = $value;
   }
   return $obj->{'contiga'};
   
}

=head2 contigb

 Title   : contigb
 Usage   : $obj->contigb($contig)
 Function: Get/Set for the 2nd contig in the overlap
 Returns : Bio::EnsEMBL::DB::RawContigI
 Args    : newvalue (optional)


=cut

sub contigb {
   my $obj = shift;

   if( @_ ) {
       my $value = shift;
       if( !ref $value || ! $value->isa('Bio::EnsEMBL::DB::ContigI') ) {
	   $obj->throw("Value [$value] is not a RawContigI. Problemo...");
       }
       $obj->{'contigb'} = $value;
   }
   return $obj->{'contigb'};
   
}

=head2 positiona
    
 Title   : positiona
 Usage   : $obj->positiona($newval)
 Function: 
 Returns : value of positiona
 Args    : newvalue (optional)


=cut

sub positiona {
   my $obj = shift;
   if( @_ ) {
       my $value = shift;
       $obj->{'positiona'} = $value;
   }
   return $obj->{'positiona'};
   
}

=head2 positionb
    
 Title   : positionb
 Usage   : $obj->positionb($newval)
 Function: 
 Returns : value of positionb
 Args    : newvalue (optional)


=cut

sub positionb {
   my $obj = shift;
   if( @_ ) {
       my $value = shift;
       $obj->{'positionb'} = $value;
   }
   return $obj->{'positionb'};
   
}

=head2 overlap_type

    $type = $self->overlap_type;
    $self->overlap_type($type);

Gets or sets the overlap_type, which is one of
the following four types:

=over 4

=item *

right2left

    5'--A------->3'
            5'--B------->3'

=item *

right2right

    5'--A------->3'
            3'<-------B--5'

=item *

left2left

    3'<-------A--5'
            5'--B------->3'

=item *

left2right

    3'<-------A--5'
            3'<-------B--5'

=back

=cut

{
    my %valid_type = map {$_, 1} qw(right2left right2right left2left left2right);
    
    sub overlap_type {
        my( $obj, $value ) = @_;
        if( $value ) {
            $obj->throw("invalid overlap type '$value'")
                unless $valid_type{$value};

            $obj->{'overlap_type'} = $value;
        }
        return $obj->{'overlap_type'};

    }
}

=head2 type

 Title   : source
 Usage   : $obj->source
 Function: String describing the source of the overlap
 Returns : value of source
 Args    : newvalue (optional)


=cut

sub source {
    my ($obj,$arg) = @_;

   if (defined($arg)) {
       $obj->{_source} = $arg;
   }
    return $obj->{_source};
    
}

=head2 distance

 Title   : distance
 Usage   : $obj->distance($dis)
 Function: 
 Returns : int
 Args    : newvalue (optional)


=cut

sub distance {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_distance} = $arg;
    }

    return $self->{_distance} || 1; 

}

=head2 invert

 Title   : invert
 Usage   : $obj->invert
 Function: 
 Returns : Reverses the sense of the overlap
 Args    : none


=cut

{
    my %type_complement = (
        'right2left'    => 'left2right',
        'left2right'    => 'right2left',
        'left2left'     => 'left2left',
        'right2right'   => 'right2right',
        );

    sub invert {
        my $self = shift;

        my $tmp = $self->contiga;
        $self->contiga($self->contigb);
        $self->contigb($tmp);

        $tmp = $self->positiona;
        $self->positiona($self->positionb);
        $self->positionb($tmp);

        my $oldtype = $self->overlap_type;
        my $newtype = $type_complement{$oldtype}
            or $self->throw("Invalid overlap type '$oldtype'");
        $self->overlap_type($newtype);
    }
}

1;
