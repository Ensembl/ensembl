package Bio::EnsEMBL::Utils::Exception;

use vars qw(@ISA);
use strict;

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


sub rearrange {

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


