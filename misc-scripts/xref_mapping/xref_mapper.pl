use strict;
use warnings;
 
use Getopt::Long;
use Cwd;
 
use vars qw(@INC);
 
#effectively add this directory to the PERL5LIB automatically
my $dir = cwd() . '/' . __FILE__;
my @d = split(/\//, $dir);
pop(@d);
$dir = join('/', @d);
unshift @INC, $dir;
 
 
my ($file, $verbose, $help);
 
GetOptions ('file=s'      => \$file,
            'verbose'     => \$verbose,
            'help'        => sub { &show_help(); exit 1;} );
 
usage("-file option is required")   if(!$file);
usage() if($help);
 
open(FILE, $file) or die("Could not open input file '$file'");
 
my  @all_species;

while( my $line = <FILE> ) {
  chomp($line);
  next if $line =~ /^#/;
  next if !$line;
 
  my ( $species, $host, $port, $dbname, $user, $password ,$dir) =
    split( "\t", $line );
 
  my $map;
  eval "require XrefMapper::$species";
  if($@) {
    warn("Could not require mapper module XrefMapper::$species\n" .
	 "Using XrefMapper::BasicMapper instead:\n$@");
    require XrefMapper::BasicMapper;
    $species = "BasicMapper";
  }
 
  {
    no strict 'refs';
    $map = "XrefMapper::$species"->new
      ( $species, $host, $port, $dbname, $user, $password ,$dir); 

  }
 
  push @all_species, $map;
 
}
 
for my $species ( @all_species ) {
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
                                                                                    
options: -file <input_file>     input file with tab delimited 'species',
                                'host', 'port', 'dbname' ,'user', 'password' and 'directory'
                                values on each line
                          
         -verbose               print out debug statements
 
         -help                  display this message
 
example: perl xref_mapper.pl -file mapper.input \\
              -verbose
 
EOF
 
  exit;
}
