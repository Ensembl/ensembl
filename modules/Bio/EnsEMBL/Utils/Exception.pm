package Bio::EnsEMBL::Utils::Exception;

use strict;
use vars qw($VERBOSITY);

$VERBOSITY = 0;

=head2 throw

  Title   : throw
  Usage   : $obj->throw("throwing exception message")
  Function: Throws an exception, which, if not caught with an eval brace
  will provide a nice stack trace to STDERR with the message
  Returns : nothing
  Args    : A string giving a descriptive error message


=cut

sub throw {

  my ($self,$string) = @_;

  my $std = $self->stack_trace_dump();

  my $out = "-------------------- EXCEPTION --------------------\n".
    "MSG: ".$string."\n".$std."-------------------------------------------\n";
  die $out;

}

=head2 warning

  Title   : warning
  Usage   : $object->warn("Warning message");
  Function: Places a warning. What happens now is down to the
            verbosity of the object  (value of $obj->verbose) 
            verbosity 0 or not set => small warning
            verbosity -1 => no warning
            verbosity 1 => warning with stack trace
            verbosity 2 => converts warnings into throw
  Example :
  Returns : 
  Args    :

=cut

sub warning {

  my ($self,$string) = @_;

  my $verbose = $self->verbose;
  $verbose = 0 unless defined $verbose;

  if ( $verbose == 2 ) {
    $self->throw($string);
  } elsif ( $verbose == -1 ) {
    return;
  } elsif ( $verbose == 1 ) {
    my $out = "-------------------- WARNING ---------------------\n".
      "MSG: ".$string."\n";
    $out .= $self->stack_trace_dump;
	
    print STDERR $out;
    return;
  }

  my $out = "-------------------- WARNING ---------------------\n".
    "MSG: ".$string."\n".
      "---------------------------------------------------\n";
  print STDERR $out;
}


=head2 throw

  Title   : throw
  Usage   : $obj->throw("throwing exception message")
  Function: Throws an exception, which, if not caught with an eval brace
            will provide a nice stack trace to STDERR with the message
  Returns : nothing
  Args    : A string giving a descriptive error message


=cut

sub throw {

  my ($self,$string) = @_;

  my $std = $self->stack_trace_dump();

  my $out = "-------------------- EXCEPTION --------------------\n".
  "MSG: ".$string."\n".$std."-------------------------------------------\n";

  die $out;

}

=head2 warning

  Title   : warning
  Usage   : $object->warning("Warning message");
  Function: Places a warning. What happens now is down to the
  verbosity of the object  (value of $obj->verbose) 
  verbosity 0 or not set => small warning
  verbosity -1 => no warning
  verbosity 1 => warning with stack trace
  verbosity 2 => converts warnings into throw
  Example :
  Returns : 
  Args    :

=cut

sub warning {

  my ($self,$string) = @_;

  my $verbose = $self->verbose;
  $verbose = 0 unless defined $verbose;


  if ( $verbose == 2 ) {
    $self->throw($string);
  } elsif ( $verbose == -1 ) {
    return;
  } elsif ( $verbose == 1 ) {
    my $out = "-------------------- WARNING ---------------------\n".
      "MSG: ".$string."\n";
    $out .= $self->stack_trace_dump;
	
    print STDERR $out;
    return;
  }    

  my $out = "-------------------- WARNING ---------------------\n".
    "MSG: ".$string."\n".
      "---------------------------------------------------\n";
  print STDERR $out;
}

=head2 verbose

  Title   : verbose
  Usage   : $self->verbose(1)
  Function: Sets verbose level for how ->warn behaves
            -1 = no warning
             0 = standard, small warning
             1 = warning with stack trace
             2 = warning becomes throw
  Returns : nothing
  Args    : -1,0,1 or 2


  =cut

sub verbose {

  my ($self,$value) = @_;

  if (ref($self) && (defined $value || ! defined $self->{'verbose'}) ) {
    $value = 0 unless defined $value;
    $self->{'verbose'} = $value;
  }
  return (ref($self) ? $self->{'_rootI_verbose'} : $VERBOSITY);
}

=head2 stack_trace_dump

 Title   : stack_trace_dump
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub stack_trace_dump{

  my ($self) = @_;

  my @stack = $self->stack_trace();

  shift @stack;
  shift @stack;
  shift @stack;

  my $out;
  my ($module,$function,$file,$position);


  foreach my $stack ( @stack) {
    ($module,$file,$position,$function) = @{$stack};
    $out .= "STACK $function $file:$position\n";
  }

  return $out;
}


=head2 stack_trace

 Title   : stack_trace
 Usage   : @stack_array_ref= $self->stack_trace
 Function: gives an array to a reference of arrays with stack trace info
           each coming from the caller(stack_number) call
 Returns : array containing a reference of arrays
 Args    : none


=cut

sub stack_trace {

  my ($self) = @_;

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

my %deprecate_cache;

sub deprecate {
  my $mesg = shift;

  my $subname = (caller( 1 ))[3] ;
  my $line_number = (caller( 1 ))[2];
  my ($file) = join( "/", (split(/\//, (caller( 1 ))[1]) )[-3,-2,-1] );
  
  return if $deprecate_cache{"$line_number:$file"};

  print STDERR "Deprecated method call in file $file line $line_number\n";
  print STDERR "Method $subname is deprecated.\n";
  print STDERR $mesg;

  $deprecate_cache{"$line_number:$file"} = 1;
}

1;
