# EnsEMBL sequence position 
#
# Copyright EMBL-EBI 2001
#
# Author: Arne Stabenau
# 
# Date : 26.02.2001
#


=head1 NAME

Bio::EnsEMBL::Range - Every object on the sequence either inherits from here or implements this functions.

=head1 SYNOPSIS


=head1 DESCRIPTION

A Range has start end strand as probably all sequence locations
have. But it can contain more than one of this. Additionally, for each
of this triples it contains "something" which represents the
context. Contexts have name and type and are supposed to be creatable
from this. (Unspecified how at the moment). So ranges can just have
name and type for context which makes them lightweigth (no additional
object created or stored).

Range is designed to be suitable for the RawContig VirtualContig
coordinate recalculation. Therefore they have the ability to consist
of more than one tuple. This will make it possible to represent
stretches on a Chromosome with a Range (in RawContig context).

Context might be required to return sequences (Implement the Bio::SeqI
interface) but for the moment we try first if the concept prooves
useful before we enforce more rules.

Main usefulness should come from ContextMapper objects which can transfer
Range objects in place to other contexts. 

As range objects can contain multiple simple ranges, it could be
feasable to represent partially successful mappings to other
contexts. Eg. a complete RawContig can map to three simple ranges (all
in one range object) in Chromosome context where the first and last
represent the beginning and end of the RawContig which did not map to
the golden path.

A more elaborate definition comes after prooving that the concept can
do all we want to do (not all that is out there) and that is easier
than the current concept. (They kill sticky Exons on business object
level). 




=head1 CONTACT

  Arne Stabenau: stabenau@ebi.ac.uk
  Ewan Birney  : birney@ebi.ac.uk

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Range;

use vars qw( @ISA );
use strict;


# not yet
# @ISA = qw( Bio::SplitLocationI );


# usually dont use the new here
sub new {
  my ( $class, @args ) = @_;
  my $self = bless {}, $class;

  return $self;
}


=head2 set_range

 Title   : set_range
 Usage   : $rangeObj->set_range( start, end, strand, context name, 
             context type )
 Function: Sets all positional things of this object. If you have more than one, this sets the first. add will set further ranges.
 Example :
 Returns :
 Args    :

=cut

sub set_range {
  my ( $self, @args ) = @_;
  $self->{_Range_simple} = \@args;
  $self->{_Range_simple_context} = undef;
  $self->{_Range_split} = undef;
}


=head2 set_range_context

 Title   : set_range_context
 Usage   : $rangeObj->set_range_context( start, end, strand, context object )
 Function: Sets all positional things of this object. If you have 
           more than one, this sets the first. add will set further 
           ranges. Instead of giving the name and type, give complete context
           object.
 Example :
 Returns :
 Args    :

=cut

sub set_range_context {
  my ( $self, @args ) = @_;
  $self->{_Range_simple} = \@args;
  $self->{_Range_simple_context} = $args[3];
  $self->{_Range_split} = undef;
}

=head2 add_range

 Title   : add_range
 Usage   : $rangeObj->add_range( start, end, strand, context name, 
             context type )
 Function: Adds another range to this range object. The continuous sequence
           represented by this object is now represented by more than one
           sequence range. 
 Example :
 Returns :
 Args    :

=cut

sub add_range {
  my ( $self, @args ) = @_;
  push( $self->{_Range_split}, \@args );
  push( $self->{_Range_split_context}, 0 ); 
}

=head2 add_range_context

 Title   : add_range_context
 Usage   : $rangeObj->add_range_context( start, end, strand, context object)

 Function: Adds another range to this range object. The continuous sequence
           represented by this object is now represented by more than one
           sequence range. 
 Example :
 Returns :
 Args    :

=cut


sub add_range_context {
  my ( $self, @args ) = @_;
  push( $self->{_Range_split}, \@args );
  push( $self->{_Range_split_context}, $args[3] );
}


=head2 get_context

 Title   : get_context
 Usage   : $rangeObj->get_context( )

 Function: Gives you the context object for the range. On split locations 
           its the context for the first one.
 Example :
 Returns :
 Args    :

=cut

sub get_context {
  my $self = shift;
  if( defined $self->{_Range_simple_context} ) {
    return $self->{_Range_simple_context};
  } else {
    # make a context from context type and context name
  }
}

=head2 get_start

 Title   : get_start
 Usage   : $rangeObj->get_start( )

 Function: Gives you the start basepair for the range. On split locations 
           its the start for the first one.
 Example :
 Returns :
 Args    :

=cut

sub get_start {
  my $self = shift;
  return $self->{_Range_simple}[0];
}

=head2 get_end

 Title   : get_end
 Usage   : $rangeObj->get_end( )

 Function: Gives you the end basepair for the range. On split locations 
           its the end for the first one.
 Example :
 Returns :
 Args    :

=cut

sub get_end {
  my $self = shift;
  return $self->{_Range_simple}[1];
}

=head2 get_strand

 Title   : get_start
 Usage   : $rangeObj->get_start( )

 Function: Gives you the strand for the range. On split locations 
           its the strand for the first one.
 Example :
 Returns : -1, 0, 1
 Args    :

=cut


sub get_strand {
  my $self = shift;
  return $self->{_Range_simple}[2];
}

=head2 get_context_name

 Title   : get_context_name
 Usage   : $rangeObj->get_context_name( )

 Function: Gives you the context name for the range. On split locations 
           its the context name for the first one.
 Example :
 Returns :
 Args    :

=cut


sub get_context_name {
  my $self = shift;
  if( defined $self->{_Range_simple_context} ) {
    return $self->{_Range_simple_context}->get_context_name();
  } else {
    return $self->{_Range_simple}[3];
  }
}

=head2 get_context_type

 Title   : get_context_type
 Usage   : $rangeObj->get_context_type( )

 Function: Gives you the context type for the range. On split locations 
           its the context type for the first one. With name and type
           Every range should be able to get its context object!
 Example :
 Returns :
 Args    :

=cut

sub get_context_type {
  my $self = shift;
  if( defined $self->{_Range_simple_context} ) {
    return $self->{_Range_simple_context}->get_context_type();
  } else {
    return $self->{_Range_simple}[4];
  }
}



=head2 get_context

 Title   : get_context
 Usage   : $rangeObj->get_context( )

 Function: Gives you the context for the range. On split locations 
           its the context for the first one. The context can deliver 
           sequence information. 

 Example :
 Returns :
 Args    :

=cut

sub get_context {
  my $self = shift;
  if( defined $self->{_Range_simple_context} ) {
    return $self->{_Range_simple_context};
  } else {
    # make context from [3] and [4]
  }
}

=head2 each_range

 Title   : each_range
 Usage   : $rangeObj->each_range( )
 Function: list of many [ start, end, strand, context name, context type ] 
 Example :
 Returns :
 Args    :

=cut

sub each_range {
  my $self = shift;
  my @res;
  my $context;

  if( defined $self->{_Range_simple_context} ) {
    push( @res, [ $self->{_Range_simple}[0..2], $self->{_Range_simple_context}->get_context_name(), $self->{_Range_simple_context}->get_context_type() ] );
  } else {
    push( @res, $self->{_Range_simple} );
  }
  
  for( my $i=0; $i<=$#self->{_Range_split}; $i++ ) {
    if( $self->{_Range_split_context}[$i] ) {
      $context = $self->{_Range_split_context}[$i];
      push( @res, $self->{_Range_split}[0..2], $context->get_context_name(), $context->get_context_type() );
    } else {
      push( @res, $self->{_Range_split}[$i] );
    }
  }
  return @res;
}

=head2 each_range_context

 Title   : each_range_context
 Usage   : $rangeObj->each_range_context( )

 Function: list of many [ start, end, strand, context ]
           This might be more expensive than just returning the
           context names and types. 
 Example :
 Returns :
 Args    :

=cut

sub each_range_context {
  my $self = shift;
  my @res;
  my $context;

  if( defined $self->{_Range_simple_context} ) {
    push( @res, [ $self->{_Range_simple}[0..2], $self->{_Range_simple_context}->get_context_name(), $self->{_Range_simple_context}->get_context_type() ] );
  } else {
    push( @res, $self->{_Range_simple} );
  }
  
  for( my $i=0; $i<=$#self->{_Range_split}; $i++ ) {
    if( $self->{_Range_split_context}[$i] ) {
      $context = $self->{_Range_split_context}[$i];
      push( @res, $self->{_Range_split}[0..2], $context->get_context_name(), $context->get_context_type() );
    } else {
      push( @res, $self->{_Range_split}[$i] );
    }
  }
  return @res;
}


