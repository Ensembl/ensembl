

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
