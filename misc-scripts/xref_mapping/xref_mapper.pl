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
my $use_existing_mappings=undef;
my $maxdump=undef;
my $help;
my $upload = undef;
my $deleteexisting;;
my $location;

GetOptions ('file=s'              => \$file,
            'verbose'             => \$verbose,
	    'dumpcheck'           => \$dumpcheck,
	    'useexistingmappings' => \$use_existing_mappings,
	    'maxdump=n'           => \$maxdump,
	    'upload'              => \$upload,
	    'deleteexisting'      => \$deleteexisting,
	    'location=s'          => \$location,
            'help'                => sub { &show_help(); exit 1;} );
 
usage("-file option is required")   if(!$file);
usage() if($help);
if(defined($dumpcheck) && defined($maxdump)){
  die "both dumpcheck and maxdump cannot both be specified\n";
}
 
open(FILE, $file) or die("Could not open input file '$file'");
 
my  @all_species;
my $xref=undef;
my $species=undef;
my $type;

while( my $line = <FILE> ) {

  chomp($line);
  next if $line =~ /^#/;
  next if !$line;

  #  print $line."\n";
  my ($key, $value) = split("=",$line);
  if($key eq "species"){
    $type = "species";
    if(defined($species)){
      push @all_species, $species;
      $species = undef;
    }
    eval "require XrefMapper::$value";
    my $module;
    if($@) {
      warn("Did not find a specific mapping module XrefMapper::$value - using XrefMapper::BasicMapper instead\n");
      require XrefMapper::BasicMapper;
      $module = "BasicMapper";
    } else{
      $module = $value;
    }

    no strict 'refs';
    $species = "XrefMapper::$module"->new();
    $species->species($value);

    if(defined($dumpcheck)){
      $species->dumpcheck("yes");
    }
    if(defined($maxdump)){
      $species->maxdump($maxdump);
    }
    if(defined($use_existing_mappings)){
      $species->use_existing_mappings("yes");
    }

  }
  elsif($key eq "xref"){
    $type = "xref";
    $xref = new XrefMapper::db();
  }
  elsif($type eq "species"){ # processing species data
    $species->$key($value);
  }
  elsif($type eq "xref"){    # processing xref data
    $xref->$key($value);
  }
}
if(defined($species)){
  push @all_species, $species;
}

my $i = 1;

for my $species ( @all_species ) {

  $species->xref($xref); # attach xref object to species object

  print "\nDumping xref & Ensembl sequences" . info($i, @all_species) . "\n";
  $species->dump_seqs($location);

  print "\nRunning mapping" . info($i, @all_species) . "\n";
  $species->build_list_and_map();

  print "\nParsing mapping output" . info($i, @all_species) . "\n";
  $species->parse_mappings();

  print "\nUploading xrefs" . info($i, @all_species) . "\n" if ($upload);
  $species->do_upload($deleteexisting) if ($upload);

  $i++;

}

print STDERR "*** All finished ***\n";

sub info {

  my ($i, @all_species) = @_;

  return " for species $i of " . scalar(@all_species);

}


sub usage {
  my $msg = shift;

  print STDERR "$msg\n\n" if($msg);
  print STDERR <<EOF;
usage: perl xref_mapper <options>
                                                                                    
options:

 -file <input_file>     input file with keyword pairs for 'species','host',
                        'port', 'dbname' ,'user', 'password' and 'directory'

  -maxdump <int>        dump only <int> sequences.

  -dumpcheck            only dump if files do not exist.

  -location             only dump a subset of the genome. Format:
                          coord_system:version:name:start:end:strand
                          e.g.
                          chromosome:NCBI34:X:1000000:2000000:1
                        start, end, strand are optional
                        coord_system can also be 'seqlevel' or 'toplevel'

  -useexistingmapping   use existing *.map files

  -upload               upload xref, object_xref, identity_xref data, and set
                        display_xrefs for genes and transcripts. Data is written
                        to *.txt etc regardless of whether this option is used.

  -deleteexisting       delete existing data from xref, object_xref,
                        identity_xref and external synonym tables. Also set all
                        existing display_xref_id columns in gene and transcript
                        to null.

  -help                 display this message
 
example: perl xref_mapper.pl -file mapper.input 
EOF
 
  exit;
}
