
#
# BioPerl module for Analysis::Log
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Analysis::Log - Logs decisions made in an analysis run

=head1 SYNOPSIS

  my $obj = new Bio::EnsEMBL::Analysis::Log;

  $obj->store("This is a stored message");

  foreach my $mess ( $obj->each_LogMessage() ) {
     print "Message: ",$mess->text," by: ",$mess->author, " at:\n", $mess->stack, "\n\n"
  }

=head1 DESCRIPTION

Provides a Perl alternative to a log file for storing messages about what
decisions were made when. This gives us the ability to eventually store
them in a database or somewhere else and track them

For free you get the person (log-in) who made the message, the date, and the full
stack trace of this call.

NB - This is *not* the same as job tracking for Blast jobs. Nor should this be
used for exception handling, as that should be done by - guess what - the exception
handling stuff on Bio::Root::Object. 

This is for logging analysis decisions.

=head1 CONTACT

Ewan Birney <birney@sanger.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Analysis::Log;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::EnsEMBL::Analysis::LogMessage;

@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;
  my $a;
  my $make = $self->SUPER::_initialize;
  
  $a = `whoami`;
  chomp $a;

  $self->_author($a); # sets author as $a
  $self->{'_eal_mess'} = []; # an array of log messages
  
  $self->dump_to_error(0);

  return $make; # success - we hope!
}

=head2 store

 Title   : store
 Usage   : $log->store("A message to store")
 Function: Generates a LogMessage object, correctly
           initalised with stored author and routine etc.
 Returns : nothing
 Args    : scalar of message to store


=cut

sub store{
   my ($self,$text) = @_;
   my $mess;

   $mess = new Bio::EnsEMBL::Analysis::LogMessage;
   $mess->author($self->_author);
   $mess->text($text);
   $self->add_LogMessage($mess);

   if( $self->dump_to_error() ) {
       print STDERR "$text\n";
   }


}

=head2 add_LogMessage

 Title   : add_LogMessage
 Usage   : $log->add_LogMessage($mess)
 Function: Stores LogMessage object
 Returns : nothing
 Args    : A logmessage object


=cut

sub add_LogMessage{
   my ($self,$mess) = @_;

   # yup. lets get anal here

   if( ! $mess->isa("Bio::EnsEMBL::Analysis::LogMessage") ) {
       $self->throw("Will not add a non LogMessage object [$mess]!");
   }

   push(@{$self->{'_eal_mess'}},$mess);

}

=head2 each_LogMessage

 Title   : each_LogMessage
 Usage   : foreach $mess ( $log->each_LogMessage ) 
 Function: each of the message 
 Returns : An array of LogMessage objects
 Args    : none


=cut

sub each_LogMessage {
   my ($self) = @_;

   return @{$self->{'_eal_mess'}};
}

=head2 _author

 Title   : _author
 Usage   : $author = $log->_author,
           $log->_author('newauthor');
 Function: Gets the author. this is stored as otherwise
           we will spend too long in system calls for each
           message
 Returns : 
 Args    :


=cut

sub _author{
   my $self = shift;

   if( @_ ) {
       my $value = shift;
       $self->{'_eal_author'} = $value;
   }

   return $self->{'_eal_author'};

}

=head2 dump_to_error

 Title   : dump_to_error
 Usage   : $obj->dump_to_error($newval)
 Function: 
 Returns : value of dump_to_error
 Args    : newvalue (optional)


=cut

sub dump_to_error{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'dump_to_error'} = $value;
    }
    return $obj->{'dump_to_error'};

}


1;
