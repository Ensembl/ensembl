
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
use Bio::EnsEMBL::ContigOverlapHelper;
use vars qw(@ISA);
use strict;

use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);

sub new {
  my($pkg,@args) = @_;

    my $self = bless {}, $pkg;

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

  return $self;
}

=head2 is_similar

 Title   : is_similar
 Usage   : $obj->is_similar($other_contig_overlap)
 Function: Finds out whether these two overlaps have the same 
           basic contigs connected at the same end
 Returns : 1 or 0
 Args    : a different contig_overlap object


=cut

sub is_similar {
    my ($self,$co ) = @_;

    if( !ref $co || !$co->isa('Bio::EnsEMBL::ContigOverlap') ) {
	$self->throw("Trying to test against a [$co] which is not good");
    }

    #warn $self->hash_string, "\n", $co->hash_string, "\n";

    if( $self->contiga->id eq $co->contiga->id ) {
	if( $self->contigb->id ne $co->contigb->id ) {
	    return 0;
	} 
	if( $self->overlap_type ne $co->overlap_type ) {
	    return 0;
	}
	return 1;
    } else {
	if( $self->contiga->id ne $co->contigb->id ) {
	    return 0;
	} 
	# yuk. inverted match
	if( $self->contigb->id ne $co->contiga->id ) {
	    return 0;
	} 
	my $tag = $self->_invert_overlap_type($co->overlap_type);
	if( $tag ne $self->overlap_type ) {
	    return 0;
	}
	return 1;
    }

    $self->throw("Should never reach here. Bad error!");
}

=head2 is_identical

 Title   : is_identical
 Usage   : $obj->is_identical($other_contig_overlap)
 Function: Finds out whether these two overlaps have the same 
           basic contigs connected at the same end with
           the same switch points
 Returns : 1 or 0
 Args    : a different contig_overlap object


=cut

sub is_identical {
    my ($self,$co ) = @_;

    if( !ref $co || !$co->isa('Bio::EnsEMBL::ContigOverlap') ) {
	$self->throw("Trying to test against a [$co] which is not good");
    }

    if( $self->contiga->id eq $co->contiga->id ) {
	if( $self->contigb->id ne $co->contigb->id ) {
	    return 0;
	} 
	if( $self->overlap_type ne $co->overlap_type ) {
	    return 0;
	}

	if( $self->positiona != $co->positiona || $self->positionb != $co->positionb || $self->distance != $co->distance ) {
	    return 0;
	}

	return 1;
    } else {
	if( $self->contiga->id ne $co->contigb->id ) {
	    return 0;
	} 
	# yuk. inverted match
	if( $self->contigb->id ne $co->contiga->id ) {
	    return 0;
	} 
	my $tag = $self->_invert_overlap_type($co->overlap_type);
	if( $tag ne $self->overlap_type ) {
	    return 0;
	}

	if( $self->positiona != $co->positionb || $self->positionb != $co->positiona || $self->distance != $co->distance ) {
	    return 0;
	}

	return 1;
    }

    $self->throw("Should never reach here. Bad error!");
}


sub _invert_overlap_type {
    my ($self,$type) = @_;

    if( $type eq 'right2left' ) { return 'left2right'; }
    if( $type eq 'left2right' ) { return 'right2left'; }

    return $type; # should be the two symetrical cases
}
	   
=head2 hash_string

 Title   : hash_string
 Usage   : $str = $co->hash_string
 Function: Returns contiga:contigb:positiona:positionb:overlap_type:distance
           Really a convience function for 
 Example : 
 Returns : 
 Args    :


=cut

sub hash_string{
   my ($self) = @_;

   return join(':',$self->contiga->id,$self->contigb->id,$self->positiona,$self->positionb,$self->overlap_type,$self->distance);

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

    return $self->{_distance} || 0; 

}

=head2 make_ContigOverlapHelper

    
    my($end, $helper) = $obj->make_ContigOverlapHelper('FOO');

For contig with id B<FOO> in the ContigOverlap
object, returns the end of contig B<FOO> which is
overlaped by the other contig (either "B<left>"
or "B<right>"), and a new
C<Bio::EnsEMBL::ContigOverlapHelper> for this
overlap.

=cut

{
    my %pol_a = (
        'right2left'   => ['right',  1],
        'right2right'  => ['right', -1],
        'left2right'   => ['left',   1],
        'left2left'    => ['left',  -1],
        );
    my %pol_b = (
        'right2left'   => ['left',   1],
        'right2right'  => ['right', -1],
        'left2right'   => ['right',  1],
        'left2left'    => ['left',  -1],
        );

    sub make_ContigOverlapHelper {
        my($self, $id) = @_;

        $self->throw("Can't return a ContigOverlapHelper without an id")
            unless $id;

        my $type = $self->overlap_type;

        my($end, $sis, $self_pos, $sister_pos, $sister_pol);
        if ($id eq $self->contiga->id) {
            $self_pos   = $self->positiona;
            $sis        = $self->contigb;
            $sister_pos = $self->positionb;
            ($end, $sister_pol) = @{$pol_a{$type}}
                or $self->throw("Illegal overlap type '$type'");
        }
        elsif ($id eq $self->contigb->id) {
            $self_pos   = $self->positionb;
            $sis        = $self->contiga;
            $sister_pos = $self->positiona;
            ($end, $sister_pol) = @{$pol_b{$type}}
                or $self->throw("Illegal overlap type '$type'");
        }
        else {
            $self->throw("ID $id not found in ContigOverlap object");
        }

        # Return new ContigOverlapHelper object
        my $helper = Bio::EnsEMBL::ContigOverlapHelper->new(
	    '-sister'           => $sis,
	    '-selfposition'     => $self_pos,
	    '-sisterposition'   => $sister_pos, 
	    '-sisterpolarity'   => $sister_pol,
	    '-distance'         => $self->distance,
	    '-source'           => $self->source,
	    );
        return($end, $helper);
    }
}

1;
