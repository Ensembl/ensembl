use strict;
use warnings;
 
use Getopt::Long;
use Cwd;
use XrefMapper::db;
 
use vars qw(@INC);
 
#effectively add this directory to the PERL5LIB automatically
my $dir = cwd() . '/' . __FILE__;
my @d = split(/\//, $dir);
pop(@d);
$dir = join('/', @d);
unshift @INC, $dir;
 
 
my $file;
my $verbose;
my $dumpcheck=undef;
my $maxdump=undef;
my $help;
 
GetOptions ('file=s'      => \$file,
            'verbose'     => \$verbose,
	    'dumpcheck'   => \$dumpcheck,
	    'maxdump=n'     => \$maxdump,
            'help'        => sub { &show_help(); exit 1;} );
 
usage("-file option is required")   if(!$file);
usage() if($help);
if(defined($dumpcheck) && defined($maxdump)){
  die "both dumpcheck and maxdump cannot both be specified\n";
}
 
open(FILE, $file) or die("Could not open input file '$file'");
 
my  @all_species;
my $xref;
#my $output=undef;
my $new=undef;
my $type;
while( my $line = <FILE> ) {
  chomp($line);
  next if $line =~ /^#/;
  next if !$line;

#  print $line."\n";
  my ($key, $value) = split("=",$line);

  if($key eq "species" || $key eq "xref"){
    if(defined($new)){ #save old one
      if($type eq "species"){
	push @all_species, $new;
      }
      else{
	$xref = $new;
      }
      $new = undef;
    }
    if($key eq "species"){
      $type = "species";
      eval "require XrefMapper::$value";
      my $module;
      if($@) {
	warn("Could not require mapper module XrefMapper::$value\n" .
	     "Using XrefMapper::BasicMapper instead:\n$@");
	require XrefMapper::BasicMapper;
	$module = "BasicMapper";
      }
      else{
	$module = $value;
      }
      {
	no strict 'refs';
	$new = "XrefMapper::$module"->new();
	$new->species($value);
      }
      if(defined($dumpcheck)){
	$new->dumpcheck("yes");
      }
      if(defined($maxdump)){
	$new->maxdump($maxdump);
      }
    }
    else{
      $type= "xref";
      $new = new XrefMapper::db();
    }
  }
  else{
    $new->$key($value);
  }
}

if(defined($new)){ #save last one
  if($type eq "species"){
    push @all_species, $new;
  }
  else{
    $xref = $new;
  }
  $new = undef;
}

for my $species ( @all_species ) {
  $species->xref($xref);
  $species->dump_seqs();
  $species->run_matching();
  $species->store();
}

 
print STDERR "*** All finished ***\n";
 
sub usage {
  my $msg = shift;
 
  print STDERR "$msg\n\n" if($msg);
  print STDERR <<EOF;
usage:   perl xref_mapper <options>
                                                                                    
options: -file <input_file>     input file with keyword pairs for  'species',
                                'host', 'port', 'dbname' ,'user', 'password' and 'directory'
                          
         -verbose               print out debug statements
 
         -maxdump <int>         dump out only int number of seqs.

         -dumpcheck             only dump if files do not exist.

         -help                  display this message
 
example: perl xref_mapper.pl -file mapper.input \\
              -verbose
 
EOF
 
  exit;
}
