=head1 NAME

Bio::EnsEMBL::Root

=head1 DESCRIPTION

Do not use Bio::EnsEMBL::Root anymore. It is included for backwards 
compatibility (every object in EnsEMBL used to inherit from this class) but 
will eventually be phased out. The replacement for the _rearrage method is the 
rearrange method which can be imported in the following way:

   use Bio::EnsEMBL::Utils::Argument qw(rearrange);

   #can now call rearrange as a class method (instead as object method)
   my ($start, $end) = rearrange(['START','END'], \@args);

If you want to use the throw or warn methods the replacement are the class
methods throw and warning from the Bio::EnsEMBL::Utils::Exception class:

   use Bio::EnsEMBL::Utils::Exception qw(throw warning);

   #can now call throw or warning even without bless reference
   warning('This is a warning');
   throw('This is an exception');


This module was stolen from BioPerl to avoid problems with moving to BioPerl 1
from 0.7

=head1 CONTACT

Functions originally from Steve Chervitz. Refactored by Ewan Birney.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Root;

use strict;
use vars qw($VERBOSITY);

$VERBOSITY = 0;

sub new{
  my($caller,@args) = @_;
  
  my $class = ref($caller) || $caller;
  return bless({}, $class);
}


=head2 throw

 Title   : throw
 Usage   : $obj->throw("throwing exception message")
 Function: Throws an exception, which, if not caught with an eval brace
           will provide a nice stack trace to STDERR with the message
 Returns : nothing
 Args    : A string giving a descriptive error message


=cut

sub throw{
   my ($self,$string) = @_;

   my $std = $self->stack_trace_dump();

   my $out = "-------------------- EXCEPTION --------------------\n".
     "MSG: ".$string."\n".$std."-------------------------------------------\n";
   die $out;

}

=head2 warn

 Title   : warn
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

sub warn{
    my ($self,$string) = @_;

    my $verbose = $self->verbose;
    $verbose = 0 unless defined $verbose;


    if( $verbose == 2 ) {
	$self->throw($string);
    } elsif( $verbose == -1 ) {
	return;
    } elsif( $verbose == 1 ) {
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

sub verbose{
   my ($self,$value) = @_;

   if(ref($self) && (defined $value || ! defined $self->{'verbose'}) ) {
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

sub stack_trace{
   my ($self) = @_;

   my $i = 0;
   my @out;
   my $prev;
   while( my @call = caller($i++)) {
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


=head2 _rearrange

 Usage     : $object->_rearrange( array_ref, list_of_arguments)
 Purpose   : Rearranges named parameters to requested order.
 Example   : $self->_rearrange([qw(SEQUENCE ID DESC)],@param);
           : Where @param = (-sequence => $s, 
	   :                 -id       => $i, 
	   :	             -desc     => $d);
 Returns   : @params - an array of parameters in the requested order.
           : The above example would return ($s, $i, $d)
 Argument  : $order : a reference to an array which describes the desired
           :          order of the named parameters.
           : @param : an array of parameters, either as a list (in
           :          which case the function simply returns the list),
           :          or as an associative array with hyphenated tags
           :          (in which case the function sorts the values 
           :          according to @{$order} and returns that new array.)
	   :	      The tags can be upper, lower, or mixed case
           :          but they must start with a hyphen (at least the
           :          first one should be hyphenated.)
 Source    : This function was taken from CGI.pm, written by Dr. Lincoln
           : Stein, and adapted for use in Bio::Seq by Richard Resnick and
           : then adapted for use in Bio::Root::Object.pm by Steve A. Chervitz.
 Comments  : (SAC)
           : This method may not be appropriate for method calls that are
           : within in an inner loop if efficiency is a concern.
           :
           : Parameters can be specified using any of these formats:
           :  @param = (-name=>'me', -color=>'blue');
           :  @param = (-NAME=>'me', -COLOR=>'blue');
           :  @param = (-Name=>'me', -Color=>'blue');
           :  @param = ('me', 'blue');  
           : A leading hyphenated argument is used by this function to 
           : indicate that named parameters are being used.
           : Therefore, the ('me', 'blue') list will be returned as-is.
           :
	   : Note that Perl will confuse unquoted, hyphenated tags as 
           : function calls if there is a function of the same name 
           : in the current namespace:
           :    -name => 'foo' is interpreted as -&name => 'foo'
	   :
           : For ultimate safety, put single quotes around the tag:
	   :    ('-name'=>'me', '-color' =>'blue');
           : This can be a bit cumbersome and I find not as readable
           : as using all uppercase, which is also fairly safe:
	   :    (-NAME=>'me', -COLOR =>'blue');
	   :
           : Personal note (SAC): I have found all uppercase tags to
           : be more managable: it involves less single-quoting,
           : the code is more readable, and there are no method naming conlicts.
           : Regardless of the style, it greatly helps to line
	   : the parameters up vertically for long/complex lists.

See Also   : L<_initialize>() 

=cut

#----------------'
sub _rearrange {
#----------------
    my($self,$order,@param) = @_;

    return unless @param;

    # If we've got parameters, we need to check to see whether
    # they are named or simply listed. If they are listed, we
    # can just return them. 

    return @param unless (defined($param[0]) && $param[0]=~/^-/); 

    # Now we've got to do some work on the named parameters.
    # The next few lines strip out the '-' characters which
    # preceed the keys, and capitalizes them.
    my $i;
    for ($i=0;$i<@param;$i+=2) {
			$param[$i]=~s/^\-//;
			$param[$i]=~tr/a-z/A-Z/;
    }

    # Now we'll convert the @params variable into an associative array.
    local($^W) = 0;  # prevent "odd number of elements" warning with -w.
    my(%param) = @param;

    my(@return_array);

    # What we intend to do is loop through the @{$order} variable,
    # and for each value, we use that as a key into our associative
    # array, pushing the value at that key onto our return array.
    my($key);

    foreach $key (@{$order}) {
			$key = uc($key);
			my($value) = $param{$key};
			delete $param{$key};
			push(@return_array,$value);
    }

    return (@return_array);
}


1;
