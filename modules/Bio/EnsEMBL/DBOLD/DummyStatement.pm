
#
# BioPerl module for Bio::EnsEMBL::DB::DummyStatement
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBOLD::DummyStatement - Dummy statement object for debugging

=head1 SYNOPSIS

 

=head1 DESCRIPTION

This object is made by the DB::Obj when the debugging level is > 10, to
provide a way of catching all the SQL correctly.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBOLD::DummyStatement;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;


@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

# set stuff in self from @args
  return $make; # success - we hope!
}

=head2 execute

 Title   : execute
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub execute {
   my ($self) = @_;
   my $fh = $self->_fileh();
   print $fh "SQL ", $self->_statement, "\n";

}


=head2 _statement

 Title   : _statement
 Usage   : $obj->_statement($newval)
 Function: 
 Returns : value of _statement
 Args    : newvalue (optional)


=cut

sub _statement{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_statement'} = $value;
    }
    return $obj->{'_statement'};

}

=head2 _fileh

 Title   : _fileh
 Usage   : $obj->_fileh($newval)
 Function: 
 Returns : value of _fileh
 Args    : newvalue (optional)


=cut

sub _fileh{
   my $obj = shift;
   if( @_ ) {
       my $value = shift;
       $obj->{'_fileh'} = $value;
   }
   return $obj->{'_fileh'};

}
