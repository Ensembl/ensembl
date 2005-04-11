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
my $ensembl=undef;
my $mapper=undef;
my $type;

my %xref_hash=();
my %species_hash=();

while( my $line = <FILE> ) {

  chomp($line);
  next if $line =~ /^#/;
  next if !$line;

  #  print $line."\n";
  my ($key, $value) = split("=",$line);
  if($key eq "species"){
    $type = "species";
    $species_hash{'species'} = $value;
  }
  elsif($key eq "xref"){
    $type = "xref";
  }
  elsif($type eq "species"){ # processing species data
    $species_hash{lc($key)} = $value;
  }
  elsif($type eq "xref"){    # processing xref data
    $xref_hash{lc($key)} = $value;
  }
}


if(defined($xref_hash{host})){
  my ($host, $user, $dbname, $pass, $port);
  $host = $xref_hash{'host'};
  $user = $xref_hash{'user'};
  $dbname = $xref_hash{'dbname'};
  if(defined($xref_hash{'pass'})){
    $pass = $xref_hash{'pass'};
  }
  else{
    $pass = '';
  }
  if(defined($xref_hash{'port'})){
    $port = $xref_hash{'port'};
  }
  else{
    $port = 3306;
  }

  $xref = new XrefMapper::db(-host => $host,
			     -port => $port,
			     -user => $user, 
			     -pass => $pass,
			     -group   => 'core',
			     -dbname => $dbname);

  if(defined($xref_hash{'dir'})){
    $xref->dir($xref_hash{'dir'});
  }

}
else{
  die "No host name given for xref\n";
}

if(defined($species_hash{'species'})){
  my $value = $species_hash{'species'};
  if ($value !~ /_/) {
      print STDERR "\'$value\' is not a recognised species - please use full species name (e.g. homo_sapiens) in $file\n";
      exit(1);
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
  my ($host, $port, $user, $dbname, $pass);
  $host = $species_hash{'host'};
  $user = $species_hash{'user'};
  $dbname = $species_hash{'dbname'};
  if(defined($species_hash{'pass'})){
    $pass = $species_hash{'pass'};
  }
  else{
    $pass = '';
  }
  if(defined($species_hash{'port'})){
    $port = $species_hash{'port'};
  }
  else{
    $port = '';
  }
  
  
  $mapper = "XrefMapper::$module"->new();

  my $core = new XrefMapper::db(-host => $host,
			     -port => $port,
			     -user => $user, 
			     -pass => $pass,
			     -group   => 'core',
			     -dbname => $dbname);


  if(defined($species_hash{'dir'})){
    $core->dir($species_hash{'dir'});
  }

  $core->species($value);

  $mapper->core($core);
  
  if(defined($dumpcheck)){
    $mapper->dumpcheck("yes");
  }
  if(defined($maxdump)){
    $mapper->maxdump($maxdump);
  }
  if(defined($use_existing_mappings)){
    $mapper->use_existing_mappings("yes");
  }

  
}
else{
  die "No Species given\n";
}


$mapper->xref($xref); # attach xref object to mapper object

print "\nDumping xref & Ensembl sequences\n";
$mapper->dump_seqs($location);

print "\nChecking external_db table\n" if ($upload);
$mapper->upload_external_db() if ($upload);

print "\nRunning mapping\n";
$mapper->build_list_and_map();

print "\nParsing mapping output\n";
$mapper->parse_mappings();

print "\nUploading xrefs\n" if ($upload);
$mapper->do_upload($deleteexisting) if ($upload);


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
                        USE WITH CAUTION - MAY GIVE CONFLICTING RESULTS!

  -useexistingmapping   use existing *.map files

  -upload               upload xref, object_xref, identity_xref data, and set
                        display_xrefs for genes and transcripts. Data is written
                        to *.txt etc regardless of whether this option is used.
                        If external_db in core database is empty, it is populated
                        from ../external_db/external_dbs.txt

  -deleteexisting       delete existing data from xref, object_xref,
                        identity_xref and external synonym tables. Also set all
                        existing display_xref_id columns in gene and transcript
                        to null.

  -help                 display this message
 
example: perl xref_mapper.pl -file mapper.input 
EOF
 
  exit;
}
