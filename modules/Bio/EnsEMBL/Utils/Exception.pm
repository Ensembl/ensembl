# EnsEMBL module for Bio::EnsEMBL::Utils::Exception
#
#

=head1 NAME

Bio::EnsEMBL::Utils::Exception - Utility functions for error handling

=head1 SYNOPSIS

    use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate verbose)

    eval {
      throw("this is an exception with a stack trace");
    );
    
    if($@) {
      print "Caught exception:\n$@";
    }

    #silence warnings
    verbose(-1);

    warning('this is a silent warning');

    verbose(0);
   
    warning('this is a warning');

    sub my_sub {
      deprecate('use other_sub() instead');
    }


=head1 DESCRIPTION

This is derived from the Bio::Root module in BioPerl.  Some formatting has
been changed and the deprecate function has been added.  Most notably the
object methods are now static class methods that can be called without 
inheriting from Bio::Root or Bio::EnsEMBL::Root.  This is especially useful
for throwing exceptions with stack traces outside of a blessed context.

The originaly implementations of these methods were by Steve Chervitz and
refactored by Ewan Birney.

It is recommended that these functions be used instead of inheriting 
unnecessarily from the Bio::EnsEMBL::Root or Bio::Root object.
The functions exported by this package provide a set of useful error handling
methods.

=head1 CONTACT

Post questions to the EnsEMBL development list: ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details exported static class methods. 

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Utils::Exception;

use Exporter;

use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(&throw &warning &stack_trace_dump 
                &stack_trace &verbose &deprecate);


=head2 throw

  Arg [1]    : string $msg
  Example    : use Bio::EnsEMBL::Utils::Exception qw(throw);
               throw('We have a problem');
  Description: Throws an exception which if not caught be an eval will
               provide a stack trace to STDERR and die.
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub throw {
  #for backwards compatibility with Bio::EnsEMBL::Root::throw
  #allow to be called as an object method as well as class method
  my $string = shift;
  $string = shift if(ref($string)); #skip object if one provided

  my $std = stack_trace_dump(3);

  my $out = "\n-------------------- EXCEPTION --------------------\n" .
            "MSG: $string\n" .
            "$std" .
            "---------------------------------------------------\n";
  die $out;
}



=head2 warning

  Arg [1]    : string warning(message);
  Example    : use Bio::EnsEMBL::Utils::Exception qw(warning)
               warning('This is a warning');
  Description: Places a warning. What happens now is down to the
               verbosity level which may be set through a call to verbose()
               verbosity 0 or not set => small warning
               verbosity -1 => no warning
               verbosity 1 => warning with stack trace
               verbosity 2 => converts warnings into throw
  Returntype : none
  Exceptions : warning every time
  Caller     : general

=cut

sub warning {
  my $string = shift;

  my $verbose = verbose() || 0;

  if ( $verbose == 2 ) {
    throw($string);
  } elsif ( $verbose == -1 ) {
    return;
  } elsif ( $verbose == 1 ) {
    my $out = "\n-------------------- WARNING ---------------------\n".
              "MSG: $string\n" .
              stack_trace_dump(3).
              "--------------------------------------------------\n";
	
    print STDERR $out;
    return;
  }

  my $out = "\n-------------------- WARNING ----------------------\n".
            "MSG: $string\n".
            "---------------------------------------------------\n";
  print STDERR $out;
}




=head2 verbose

  Arg [1]    : (optional) int 
  Example    : use Bio::EnsEMBL::Utils::Exception qw(verbose warning);
               verbose(-1);
               warning("No warning displayed");
               verbose(0);
               warning("Regular warning displayed");
  Description: Gets/Sets verbose level for how warning() behaves
               -1 = no warning
                0 = standard, small warning
                1 = warning with stack trace
                2 = warning becomes throw
  Returntype : int 
  Exceptions : none
  Caller     : general

=cut

my $VERBOSITY = 0;

sub verbose {
  $VERBOSITY = shift if(@_);

  return $VERBOSITY;
}



=head2 stack_trace_dump

  Arg [1]    : (optional) int $levels
               The number of levels to ignore from the top of the stack when
               creating the dump. This is useful when this is called internally
               from a warning or throw function when the immediate caller and 
               stack_trace_dump function calls are themselves uninteresting.
  Example    : use Bio::EnsEMBL::Utils::Exception qw(stack_trace_dump);
               print STDERR stack_trace_dump();
  Description: Returns a stack trace formatted as a string
  Returntype : string
  Exceptions : none
  Caller     : general, throw, warning

=cut

sub stack_trace_dump{
  my @stack = stack_trace();

  my $levels = 2; #default is 2 levels so stack_trace_dump call is not present
  $levels = shift if(@_);
  $levels = 1 if($levels < 1);
  
  while($levels) {
    $levels--;
    shift @stack;
  }

  my $out;
  my ($module,$function,$file,$position);


  foreach my $stack ( @stack) {
    ($module,$file,$position,$function) = @{$stack};
    $out .= "STACK $function $file:$position\n";
  }

  return $out;
}



=head2 stack_trace

  Arg [1]    : none
  Example    : use Bio::EnsEMBL::Utils::Exception qw(stack_trace)
  Description: Gives an array to a reference of arrays with stack trace info
               each coming from the caller(stack_number) call
  Returntype : array of listrefs of strings
  Exceptions : none
  Caller     : general, stack_trace_dump()

=cut

sub stack_trace {
  my $i = 0;
  my @out;
  my $prev;
  while ( my @call = caller($i++)) {

    # major annoyance that caller puts caller context as
    # function name. Hence some monkeying around...
    $prev->[3] = $call[3];
    push(@out,$prev);
    $prev = \@call;
  }
  $prev->[3] = 'toplevel';
  push(@out,$prev);
  return @out;
}


=head2 deprecate

  Arg [1]    : string $mesg
               A message describing why a method is deprecated
  Example    : use Bio::EnsEMBL::Utils::Exception qw(deprecate)
               sub old_sub {
                 deprecate('Please use new_sub() instead');
               }
  Description: Prints a warning to STDERR that the method which called 
               deprecate() is deprecated.  Also prints the line number and 
               file from which the deprecated method was called.  Deprecated
               warnings only appear once for each location the method was 
               called from
  Returntype : none
  Exceptions : warning every time
  Caller     : deprecated methods

=cut

my %DEPRECATED;

sub deprecate {
  my $mesg = shift;

  my @caller = caller(1);
  my $subname = $caller[3] ;
  my $line = $caller[2];

  #use only 2 subdirs for brevity when reporting the filename
  my $file;
  my @path = $caller[1];
  $file = pop(@path);
  my $i = 0;
  while(@path && $i < 2) {
    $i++;
    $file .= pop(@path);
  }

  #keep track of who called this method so that the warning is only displayed
  #once per deprecated call
  return if $DEPRECATED{"$line:$file"};

  if($VERBOSITY > -1) {
    print STDERR "\n------------------ DEPRECATED ---------------------\n" .
                 "Deprecated method call in file $file line $line.\n" .
                 "Method $subname is deprecated.\n" .
                 "$mesg\n" .
                 "---------------------------------------------------\n";
  }

  $DEPRECATED{"$line:$file"} = 1;
}

1;
