
#
# BioPerl module for Bio::EnsEMBL::Virtual::MapContig
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Virtual::MapContig - An internal object for VirtualContig

=head1 SYNOPSIS

    # get a MapContig from a Virtual::Contig
    my( $mapcontig,$position,$ori) = $vc->raw_contig_position(2000,1);

    # map contigs have-a raw contig
    print "Raw contig id is",$mapcontig->contig->id,"\n";

    # map contigs have a start/end in VC position
    print "At position",$mapcontig->start,":",$mapcontig->end,"\n";
    
    # map contigs have a start/end in RC position
    print "In RC ,"$mapcontig->rawcontig_start,":",$mapcontig->rawcontig_end,"\n";

    # they also have an orientation relative to VC
    print "Relative to VC this is in ",$mapcontig->orientation,"\n";

=head1 DESCRIPTION

This is an internal object to virtual contigs. A MapContig has-a
RawContig and contains the additional start/end points about what
region of the RawContig maps to the VC. It has start/end in both VC
and RC coordinates and the orientation of the RC relative to the
VC. As the length of the VC and RC coordinates have to the same, to
avoid messing this up, the end point of the RC coordinates is
calculated from start,end, and rawcontig_start.

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Virtual::MapContig;
use vars qw(@ISA);
use strict;
use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);

sub new {
    my ($class,@args) = @_;
    
    my $self = {};
    bless $self,$class;

    my ($contig,$start,$end,$rawcontig_start,$ori) = 
	$self->_rearrange([qw( 
			       RAWCONTIG
			       START 
			       END
			       RAWCONTIG_START 
			       ORIENTATION
			       )],@args);

    if( !defined $contig ||
	!defined $start  ||
	!defined $end  ||
	!defined $rawcontig_start  ||
	!defined $ori ) {
	$self->throw("Did not pass all arguments into MapContig");
    }

    $self->start($start);
    $self->end($end);
    $self->contig($contig);
    $self->orientation($ori);
    $self->rawcontig_start($rawcontig_start);

    if( $self->rawcontig_end > $self->contig->length ) {
	$self->throw("Attempting to build a mapcontig $start:$end starting at $rawcontig_start which is greater than the length of the contig ".$self->contig->length." on ".$self->contig->id);
    }

    return $self;
}

=head2 start

 Title   : start
 Usage   : $obj->start($newval)
 Function: 
 Example : 
 Returns : value of start
 Args    : newvalue (optional)


=cut

sub start{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'start'} = $value;
    }
    return $obj->{'start'};

}

=head2 end

 Title   : end
 Usage   : $obj->end($newval)
 Function: 
 Example : 
 Returns : value of end
 Args    : newvalue (optional)


=cut

sub end{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'end'} = $value;
    }
    return $obj->{'end'};

}

=head2 rawcontig_start

 Title   : rawcontig_start
 Usage   : $obj->rawcontig_start($newval)
 Function: 
 Example : 
 Returns : value of rawcontig_start
 Args    : newvalue (optional)


=cut

sub rawcontig_start{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'rawcontig_start'} = $value;
    }
    return $obj->{'rawcontig_start'};

}

=head2 rawcontig_end

 Title   : rawcontig_end
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub rawcontig_end{
   my ($self,@args) = @_;

   return $self->rawcontig_start + $self->_length - 1;
}


=head2 orientation

 Title   : orientation
 Usage   : $obj->orientation($newval)
 Function: 
 Example : 
 Returns : value of orientation
 Args    : newvalue (optional)


=cut

sub orientation{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'orientation'} = $value;
    }
    return $obj->{'orientation'};

}

=head2 contig

 Title   : contig
 Usage   : $obj->contig($newval)
 Function: 
 Example : 
 Returns : value of contig
 Args    : newvalue (optional)


=cut

sub contig{
   my ($obj,$value) = @_;
   if( defined $value) {
       if( !ref $value || !$value->isa('Bio::EnsEMBL::DB::RawContigI') ) {
	   $obj->throw("Value [$value] is not appropiate");
       }

      $obj->{'contig'} = $value;
    }
    return $obj->{'contig'};

}

=head2 start_in

 Title   : start_in
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub start_in{
   my ($self,@args) = @_;

   my ($p,$f,$l) = caller;
   $self->warn("$f:$l:start_in a deprecated method. use rawcontig_start instead");
   return $self->rawcontig_start();

}


=head2 end_in

 Title   : end_in
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub end_in{
   my ($self,@args) = @_;

   my ($p,$f,$l) = caller;
   $self->warn("$f:$l:end_in a deprecated method. use rawcontig_end instead");
   return $self->rawcontig_end();

}

=head2 _length

 Title   : _length
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _length{
   my ($self) = @_;

   return $self->end - $self->start +1;
}
    

1;

