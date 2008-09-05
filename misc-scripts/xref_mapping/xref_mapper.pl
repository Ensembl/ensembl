#!/usr/bin/perl

=head1 NAME

xref_mapper.pl - Create Ensembl gene model xref mappings

=head1 SYNOPSIS

=over 15

=item B<xref_mapper.pl>

[B<-help>]
[B<-man>]
[B<-file>S< >I<config_file>]
[B<-maxdump>S< >I<int>]
[B<-location>S< >I<string>]
[B<-logicname>S< >I<string>] 
[B<-useexistingmappings>]
[B<-upload>]
[B<-delete_external_db>]
[B<-notriage>]
[B<-recalc_display_xrefs>}
[B<-dumpcheck>]
[B<-external_db_file>]
[B<-nofarm>]

=back

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exit.

=item B<-man>

Print this command's manual page and exit.

=item B<-file> I<config_file>

Input file with keyword pairs for 'species','host', 'port', 'dbname',
'user', 'password' and 'directory'.

=item B<-maxdump> I<int>

dump only I<int> sequences.

=item B<-dumpcheck>

only dump if files do not exist.

=item B<-location> I<string>

only dump a subset of the genome. Format:

    coord_system:version:name:start:end:strand e.g.
    chromosome:NCBI34:X:1000000:2000000:1 

start, end, strand are optional. coord_system can also be 'seqlevel' or
'toplevel'

USE WITH CAUTION -MAY GIVE CONFLICTING RESULTS!

=item B<-logicname>

dump only the specified (analysis.logic_name) gene type from the core DB

=item B<-useexistingmappings>

use existing *.map files

=item B<-upload>

upload xref, object_xref, identity_xref data, and set display_xrefs
for genes and transcripts. Data is written to *.txt etc regardless of
whether this option is used. If external_db in core database is empty,
it is populated from ../external_db/external_dbs.txt

=item B<-delete_external_db>

deletes all entries of the external_db table if it contains any rows
and uploads new data into the table - you have to interactively
confirm the deletion before. Works only if option B<-upload> is used
as well.

=item B<-recalc_display_xrefs>

recalculate the display xrefs for all the genes and transcripts and
also recalculate the gene descriptions. This only reads the data
already stored in the core database and set these attributes based on
what is already stored. Useful if you have changed the prioritys for
the display xrefs etc or what sources are used in which order for gene
desriptions and merely want to recalc what should be
displayed. Generates .txt and .sql files into the species-specific
location specified by the --file. The --upload flag determines whether
these files are loaded into the target DB manually.

=item B<-notriage>

don't dump triage data

=item B<-nofarm>

run exonerate locally and do not use compute farm

=item B<EXAMPLE:>

perl xref_mapper.pl -file mapper.input

=back

=head1 DESCRIPTION

Creates xrefs to Ensembl gene models as follows;

1. Connects to the Ensembl core and xref databases using the settings 
   configured in B<-file>

2. Dumps cDNAs and peptides for the gene models in the core database.

3. Dumps the sequences for the xrefs in the xref database.

4. Uses the method specified in the XrefMapper::I<species>.pm file to
   associate xrefs with gene models.

=head1 SEE ALSO

xref_parser.pl

=head1 AUTHOR

Glenn and Ian

=cut


use strict;
use warnings;
 
use Getopt::Long;
use Pod::Usage;
use Cwd;
use XrefMapper::db;
 
use vars qw(@INC);
 
$| = 1;

my $file;
my $verbose;
my $dumpcheck=undef;
my $use_existing_mappings=undef;
my $maxdump=undef;
my $help;
my $man;
my $upload = undef;
my $delete_external_db ; 
my $location;
my $logic_name;
my $nofarm;
my $notriage=undef;
my $recalc_display_xrefs = undef;
my $external_db_file="../external_db/external_dbs.txt";


print "Options: ".join(" ",@ARGV)."\n";

GetOptions ('file=s'                    => \$file,
            'verbose'                   => \$verbose,
	    'dumpcheck'                 => \$dumpcheck,
	    'useexistingmappings'       => \$use_existing_mappings,
	    'maxdump=n'                 => \$maxdump,
	    'upload'                    => \$upload,
	    'delete_external_db'        => \$delete_external_db , 
	    'location=s'                => \$location,
            'logicname=s'               => \$logic_name,
	    'notriage'                  => \$notriage,
	    'recalc_display_xrefs_only' => \$recalc_display_xrefs,
	    'external_db_file=s'        => \$external_db_file,
            'help'                      => \$help,
	    'nofarm'                    => \$nofarm,
            'man'                       => \$man );
pod2usage(1)            if ($help);
pod2usage(VERBOSE => 2) if ($man);
#usage("-file option is required")   if(!$file);

pod2usage("\n[*DIE] -file option is required\n") if(!$file);

if(defined($dumpcheck) && defined($maxdump)){
  pod2usage( "\n[*DIE] Cannot specify both -dumpcheck and -maxdump\n" );
}

open(FILE, $file) or pod2usage("\nCannot open input file '$file':\n $!\n");
 
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
  if(defined($xref_hash{'password'})){
    $pass = $xref_hash{'password'};
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

  my $module;  
  my $class = "XrefMapper/$value.pm";
  eval {
    require $class;
  };
  if($@) {
    if ($@ =~ /Can\'t locate $class/) {
      warn("Did not find a specific mapping module XrefMapper::$value - using XrefMapper::BasicMapper instead\n");
      require XrefMapper::BasicMapper;
      $module = "BasicMapper";
    } else {
      die "$@";
    }

  } else{
    $module = $value;
  }

  no strict 'refs';
  my ($host, $port, $user, $dbname, $pass);
  $host = $species_hash{'host'};
  $user = $species_hash{'user'};
  $dbname = $species_hash{'dbname'};
  if(defined($species_hash{'password'})){
    $pass = $species_hash{'password'};
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

  if(defined($upload)){
    $mapper->upload("yes");
  }
  if(defined($external_db_file)){
    $mapper->external_db_file($external_db_file);
  }
  
  if(defined($dumpcheck)){
    $mapper->dumpcheck("yes");
  }
  if(defined($maxdump)){
    $mapper->maxdump($maxdump);
  }
  if(defined($nofarm)){
    $mapper->nofarm("yes");
  }
  if(defined($use_existing_mappings)){
    $mapper->use_existing_mappings("yes");
    $mapper->dumpcheck("yes");
  }
  if(defined($logic_name)){
    $mapper->logic_name($logic_name);
  }

  
}
else{
  die "No Species given\n";
}


$mapper->xref($xref); # attach xref object to mapper object



#$mapper->species_specific_cleanup();
#$mapper->species_specific_pre_attributes_set();
#exit;



if(defined($recalc_display_xrefs)){
  $mapper->genes_and_transcripts_attributes_set();
  print "Finished recalculating display xrefs and gene descriptions\n";
  exit;
}

#print "\nDumping xref & Ensembl sequences\n";
$mapper->dump_seqs($location);


print "\nChecking external_db table\n" if ($upload);
$mapper->upload_external_db($delete_external_db ) if $upload ;  

unless( $notriage ){
  unless( $mapper->count_unmapped_reasons ){
    die( "The unmapped_reason table in the target database is empty. ",
         "Either run this script using the -notriage flag, ",
         "or populate the table using e.g. ",
         "ensembl/misc-scripts/unmapped_reason/update_unmapped_reasons.pl\n" );
  }
}

$mapper->build_list_and_map();

$mapper->find_priority_sources();


print "\nParsing mapping output\n";
$mapper->parse_mappings($notriage);


if ($upload) {
  $mapper->do_upload();

  print "\nProcessing priority xrefs\n";
  $mapper->process_priority_xrefs();

  print "\nChecking pair data\n";
  $mapper->add_missing_pairs();

  if ($notriage) {
    $mapper->dump_xref_with_no_triage_data();
  } else {
    $mapper->dump_triage_data();
  }

  if ( !defined($notriage) ) {
    print "\nPriority unmapped xrefs sorting\n";
    $mapper->unmapped_data_for_prioritys();

    print "\nProcess DEPENDENT unmapped object data\n";
    $mapper->write_dependent_unmapped_objects();
  }

  print "\nChecking xrefs\n";
  $mapper->cleanup_database();

  # if special sources are set then make sure these ONLY have one
  # xref per transcript. This is based on % identity etc but if
  # two have the same and we cannot distinguish between them still keep
  # both.
  $mapper->check_special_sources();
} ## end if ($upload)


if ($upload) {
  $mapper->run_coordinatemapping($upload);
  $mapper->cleanup_database();
  $mapper->species_specific_pre_attributes_set();
  $mapper->genes_and_transcripts_attributes_set();
}
else {
  print "Gene descriptions, display_xrefs and status cannot be set until "
      . "xrefs have been uploaded;\n"
      . "  The -recalc_display_xrefs flag can be used to retrofit "
      . "these data after manual upload";
}

print "*** All finished ***\n";

sub info {
  my ( $i, @all_species ) = @_;
  return " for species $i of " . scalar(@all_species);
}
