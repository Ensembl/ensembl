=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Utils::Exception - Utility functions for error handling

=head1 SYNOPSIS

  use Bio::EnsEMBL::Utils::Exception
    qw(throw warning deprecate verbose try catch);

  or to get all methods just

  use Bio::EnsEMBL::Utils::Exception;

  eval { throw("this is an exception with a stack trace") };
  if ($@) {
    print "Caught exception:\n$@";
  }

  # Or you can us the try/catch confortable syntax instead to deal with
  # throw or die.  Don't forget the ";" after the catch block.  With
  # this syntax, the original $@ is in $_ in the catch subroutine.

  try {
    throw("this is an exception with a stack trace");
  }
  catch { print "Caught exception:\n$_" };

  # silence warnings
  verbose('OFF');

  warning('this is a silent warning');

  #show deprecated and warning messages but not info
  verbose('DEPRECATE');

  warning('this is a warning');

  # show all messages
  verbose('ALL');

  info('this is an informational message');

  sub my_sub { deprecate('use other_sub() instead') }

  verbose('EXCEPTION');
  info( 'This is a high priority info message.', 1000 );

=head1 DESCRIPTION

This is derived from the Bio::Root module in BioPerl.  Some formatting
has been changed and the deprecate function has been added.  Most
notably the object methods are now static class methods that can be
called without inheriting from Bio::Root.  This is
especially useful for throwing exceptions with stack traces outside of a
blessed context.

The originaly implementations of these methods were by Steve Chervitz
and refactored by Ewan Birney.

It is recommended that these functions be used instead of inheriting
unnecessarily from the Bio::Root object.  The
functions exported by this package provide a set of useful error
handling methods.

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Exception;

use strict;
use warnings;

use Bio::EnsEMBL::ApiVersion;

use Exporter;

use vars qw(@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(throw warning stack_trace_dump 
             stack_trace verbose deprecate info try catch);

my $VERBOSITY         = 3000;
my $DEFAULT_INFO      = 4000;
my $DEFAULT_DEPRECATE = 3000;
my $DEFAULT_WARNING   = 2000;
my $DEFAULT_EXCEPTION = 1000;


=head2 throw

  Arg [1]    : string $msg
  Arg [2]    : (optional) int $level
               override the default level of exception throwing
  Example    : use Bio::EnsEMBL::Utils::Exception qw(throw);
               throw('We have a problem');
  Description: Throws an exception which if not caught by an eval will
               provide a stack trace to STDERR and die.  If the verbosity level
               is lower than the level of the throw, then no error message is
               displayed but the program will still die (unless the exception
               is caught).
  Returntype : none
  Exceptions : thrown every time
  Caller     : generally on error

=cut

sub throw {
  my $string = shift;

  # For backwards compatibility with Bio::EnsEMBL::Root::throw:  Allow
  # to be called as an object method as well as class method.  Root
  # function now deprecated so call will have the string instead.

  $string = shift if ( ref($string) );    # Skip object if one provided.
  $string = shift if ( $string eq "Bio::EnsEMBL::Utils::Exception" );

  my $level = shift;
  $level = $DEFAULT_EXCEPTION if ( !defined($level) );

  if ( $VERBOSITY < $level ) {
    die("\n");    # still die, but silently
  }

  my $std = stack_trace_dump(3);

  my $out = sprintf(
             "\n" .
             "-------------------- EXCEPTION --------------------\n" .
             "MSG: %s\n" .
             "%s" .
             "Date (localtime)    = %s\n" .
             "Ensembl API version = %s\n" .
             "---------------------------------------------------\n",
             $string, $std, scalar( localtime() ), software_version() );

  die($out);
} ## end sub throw



=head2 warning

  Arg [1]    : string warning(message);
  Arg [2]    : (optional) int level
               Override the default level of this warning changning the level
               of verbosity at which it is displayed.
  Example    : use Bio::EnsEMBL::Utils::Exception qw(warning)
               warning('This is a warning');
  Description: If the verbosity level is higher or equal to the level of this 
               warning then a warning message is printed to STDERR.  If the 
               verbosity lower then nothing is done.  Under the default
               levels of warning and verbosity warnings will be displayed.
  Returntype : none
  Exceptions : warning every time
  Caller     : general

=cut

sub warning {
  my $string = shift;

  # See throw() for this:
  $string = shift if ( ref($string) );    # Skip object if one provided.
  $string = shift if ( $string eq "Bio::EnsEMBL::Utils::Exception" );

  my $level = shift;
  $level = $DEFAULT_WARNING if ( !defined($level) );

  return if ( $VERBOSITY < $level );

  my @caller = caller;
  my $line = $caller[2] || '';

  # Use only two sub-dirs for brevity when reporting the file name.
  my $file;
  my @path = split( /\//, $caller[1] );
  $file = pop(@path);
  my $i = 0;
  while ( @path && $i < 2 ) {
    $i++;
    $file = pop(@path) . "/$file";
  }

  @caller = caller(1);
  my $caller_line;
  my $caller_file;
  $i = 0;
  if (@caller) {
    @path        = split( /\//, $caller[1] );
    $caller_line = $caller[2];
    $caller_file = pop(@path);
    while ( @path && $i < 2 ) {
      $i++;
      $caller_file = pop(@path) . "/$caller_file";
    }
  }

  my $out =
    sprintf( "\n" .
             "-------------------- WARNING ----------------------\n" .
             "MSG: %s\n" .
             "FILE: %s LINE: %d\n",
             $string, $file, $line );

  if ( defined($caller_file) ) {
    $out .= sprintf( "CALLED BY: %s  LINE: %d\n", $caller_file,
                     $caller_line );
  }
  $out .= sprintf(
           "Date (localtime)    = %s\n" .
           "Ensembl API version = %s\n" .
           "---------------------------------------------------\n",
           scalar( localtime() ), software_version() );

  warn($out);

} ## end sub warning



=head2 info

  Arg [1]    : string $string
               The message to be displayed
  Arg [2]    : (optional) int $level
               Override the default level of this message so it is displayed at
               a different level of verbosity than it normally would be.
  Example    : use Bio::EnsEMBL::Utils::Exception qw(verbose info)
  Description: This prints an info message to STDERR if verbosity is higher 
               than the level of the message.  By default info messages are not
               displayed.
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub info {
  my $string = shift;
  $string = shift if($string eq "Bio::EnsEMBL::Utils::Exception");
  my $level  = shift;

  $level = $DEFAULT_INFO if(!defined($level));

  return if($VERBOSITY < $level);

  print STDERR "INFO: $string\n";
}



=head2 verbose

  Arg [1]    : (optional) int 
  Example    : use Bio::EnsEMBL::Utils::Exception qw(verbose warning);
               #turn warnings and everything more important on (e.g. exception)
               verbose('WARNING'); 
               warning("Warning displayed");
               info("This won't be displayed");
               deprecate("This won't be diplayed"); 

               #turn exception messages on
               verbose('EXCEPTION'); 
               warning("This won't do anything");
               throw("Die with a message");

               #turn everying off
               verbose('OFF'); #same as verbose(0);               
               warning("This won't do anything");
               throw("Die silently without a message");

               #turn on all messages
               verbose('ALL');
               info("All messages are now displayed");

               if(verbose() > 3000) {
                 print "Verbosity is pretty high";
               }

  Description: Gets/Sets verbosity level which defines which messages are
               to be displayed.  An integer value may be passed or one of the
               following strings:
               'OFF'       (= 0)
               'EXCEPTION' (= 1000)
               'WARNING'   (= 2000)
               'DEPRECATE' (= 3000)
               'INFO'      (= 4000)
               'ALL'       (= 1000000)

  Returntype : int 
  Exceptions : none
  Caller     : general

=cut


sub verbose {
  if(@_) {
    my $verbosity = shift;
    $verbosity = shift if($verbosity eq "Bio::EnsEMBL::Utils::Exception");
    if($verbosity =~ /\d+/) { #check if verbosity is an integer
      $VERBOSITY = $verbosity;
    } else {
      $verbosity = uc($verbosity);
      if($verbosity eq 'OFF' || $verbosity eq 'NOTHING' || 
         $verbosity eq 'NONE') {
        $VERBOSITY = 0;
      } elsif($verbosity eq 'EXCEPTION' || $verbosity eq 'THROW') {
        $VERBOSITY = $DEFAULT_EXCEPTION;
      } elsif($verbosity eq 'WARNING' || $verbosity eq 'WARN') {
        $VERBOSITY = $DEFAULT_WARNING;
      } elsif($verbosity eq 'DEPRECATE' || $verbosity eq 'DEPRECATED') {
        $VERBOSITY = $DEFAULT_DEPRECATE;
      } elsif($verbosity eq 'INFO') {
        $VERBOSITY = $DEFAULT_INFO;
      } elsif($verbosity eq 'ON' || $verbosity eq 'ALL') {
        $VERBOSITY = 1e6;
      } else {
        $VERBOSITY = $DEFAULT_WARNING;
        warning("Unknown level of verbosity: $verbosity");
      }
    }
  }

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
  $levels = shift if($levels eq "Bio::EnsEMBL::Utils::Exception");
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
               called from.  No message is displayed if the level of verbosity
               is lower than the level of the warning.
  Returntype : none
  Exceptions : warning every time
  Caller     : deprecated methods

=cut

my %DEPRECATED;

sub deprecate {
  my $mesg = shift;
  $mesg = shift if($mesg eq "Bio::EnsEMBL::Utils::Exception"); #skip object if one provided

  my $level = shift;

  $level = $DEFAULT_DEPRECATE if(!defined($level));

  return if($VERBOSITY < $level);
                                 
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
  return if $DEPRECATED{"$line:$file:$subname"};

  if ( $VERBOSITY > -1 ) {
    print STDERR
      "\n------------------ DEPRECATED ---------------------\n"
      . "Deprecated method call in file $file line $line.\n"
      . "Method $subname is deprecated.\n"
      . "$mesg\n"
      . "Ensembl API version = "
      . software_version() . "\n"
      . "---------------------------------------------------\n";
  }

  $DEPRECATED{"$line:$file:$subname"} = 1;
}

=head2 try/catch

  Arg [1]    : anonymous subroutine
               the block to be tried
  Arg [2]    : return value of the catch function
  Example    : use Bio::EnsEMBL::Utils::Exception qw(throw try catch)
               The syntax is:
               try { block1 } catch { block2 };
               { block1 } is the 1st argument
               catch { block2 } is the 2nd argument
               e.g.
               try {
                 throw("this is an exception with a stack trace");
               } catch {
                 print "Caught exception:\n$_";
               };
               In block2, $_ is assigned the value of the first
               throw or die statement executed in block 1.

  Description: Replaces the classical syntax
               eval { block1 };
               if ($@) { block2 }
               by a more confortable one.
               In the try/catch syntax, the original $@ is in $_ in the catch subroutine.
               This try/catch implementation is a copy and paste from
               "Programming Perl" 3rd Edition, July 2000, by L.Wall, T. Christiansen
               & J. Orwant. p227, and is only possible because of subroutine prototypes.
  Returntype : depend on what is implemented the try or catch block
  Exceptions : none
  Caller     : general

=cut

## no critic
sub try (&$) {
  my ($try, $catch) = @_;
  eval { &$try };
  if ($@) {
    chop $@;
    local $_ = $@;
    &$catch;
  }
}

## no critic
sub catch (&) {
  shift;
}

1;
