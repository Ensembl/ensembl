
#
# BioPerl module for Analysis::LogMessage
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Analysis::LogMessage - A single Log message for storing

=head1 SYNOPSIS

  # see Bio::EnsEMBL::Analysis::Log 

=head1 DESCRIPTION

Has actual fields for log messages etc.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Analysis::LogMessage;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;


@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my($author,$text) = 
      $self->_rearrange([qw(AUTHOR
			    TEXT
			    )],@args);

  
  my $make = $self->SUPER::_initialize;

  $author && $self->author($author);
  $text   && $self->text($text);

  # set caller stack and the date

  my $t = time();
  my $d = localtime($t);
  $self->date($t);
  $self->_set_stack();
  return $make; # success - we hope!
}

=head2 text

 Title   : text
 Usage   : $obj->text($newval)
 Function: Text is the actual textual message of this error
 Returns : value of text
 Args    : newvalue (optional)


=cut

sub text{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'text'} = $value;
    }
    return $obj->{'text'};

}

=head2 author

 Title   : author
 Usage   : $obj->author($newval)
 Function: This is the stored author for this analysis
 Returns : value of author
 Args    : newvalue (optional)


=cut

sub author{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'author'} = $value;
    }
    return $obj->{'author'};

}

=head2 date_string

 Title   : date_string
 Usage   : print "Message done on:", $mess->date_string, "\n"
 Function: Provides a string representation of the date
 Returns : scalar
 Args    : none


=cut

sub date_string{
   my ($self) = @_;
   my $d = localtime($self->date());
   return $d;
}

=head2 date

 Title   : date
 Usage   : $obj->date($newval)
 Function: Gives the unix time value for the date
 Returns : value of date
 Args    : newvalue (optional)


=cut

sub date{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'date'} = $value;
    }
    return $obj->{'date'};

}

=head2 stack

 Title   : stack
 Usage   : $obj->stack($newval)
 Function: Gives the stack as "file:line,file:line" type scalar.
           De-munge it with a double split if you want to, but probably
           only really want to give to humans anyway.
 Returns : value of stack
 Args    : newvalue (optional)

  Notes: Beware. This is different from ->stack_trace, which is
         inherited from Bio::Root::Object. 

=cut

sub stack{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'stack'} = $value;
    }
    return $obj->{'stack'};

}

=head2 _set_stack

 Title   : _set_stack
 Usage   : $self->_set_stack()
 Function: Sets the stack to the current, discarding the boring info
 Example : Should only be called from _initialize
 Returns : 
 Args    :


=cut

sub _set_stack{
   my ($self,@args) = @_;

   
   my $str;
   
   my @tmp = $self->stack_trace();
   # first 4 are boring, as they are always the same, being
   # this function, new, initialize and store...
   
   $str =shift @tmp;
   $str =shift @tmp;
   $str =shift @tmp;
   $str =shift @tmp;
   $str = "";

   foreach my $t ( @tmp ) {
       my @stack = @$t;
       
       if( $t ) {
	   
	   # Debugging
	   #  foreach my $tt ( @$t ) {
	   #      print STDERR "It is $tt\n";
	   #  }
	   
	   
	   # 1 is file, 2 is line number
	   
	   
	   $str .= join(':',$stack[1],$stack[2]);
	   $str .= ",";
       }
   }
   
   $self->stack($str);
   
}


1;






