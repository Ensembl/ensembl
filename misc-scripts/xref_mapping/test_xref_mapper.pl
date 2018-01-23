#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


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
[B<-useexistingmapping>]
[B<-upload>]
[B<-delete_external_db>]
[B<-notriage>]
[B<-no_recalc_displays>]
[B<-recalc_displays_only>]
[B<-delete_unmapped>]
[B<-dumpcheck>]

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

=item B<-no_recalc_displays>

Do NOT recalc the display xrefs

=item B<-location> I<string>

only dump a subset of the genome. Format:

    coord_system:version:name:start:end:strand e.g.
    chromosome:NCBI34:X:1000000:2000000:1 

start, end, strand are optional. coord_system can also be 'seqlevel' or
'toplevel'

USE WITH CAUTION -MAY GIVE CONFLICTING RESULTS!

=item B<-logicname>

dump only the specified (analysis.logic_name) gene type from the core DB

=item B<-useexistingmapping>

use existing *.map files

=item B<upload>

upload xref, object_xref, identity_xref data, and set display_xrefs
for genes and transcripts. Data is written to *.txt etc regardless of
whether this option is used. If external_db in core database is empty,
it is populated from ../external_db/external_dbs.txt

=item B<-delete_external_db>

deletes all entries of the external_db table if it contains any rows
and uploads new data into the table - you have to interactively
confirm the deletion before. Works only if option B<-upload> is used
as well.

=item B<-notriage>

don't dump triage data

=item B<-delete_unmapped>

deletes data from the unmapped_object table.

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
 
my $file;
my $verbose;
my $dumpcheck=undef;
my $no_recalc_displays=undef;
my $use_existing_mappings=undef;
my $maxdump=undef;
my $help;
my $man;
my $upload = undef;
my $delete_external_db ; 
my $location;
my $logic_name;
my $notriage=undef;
my $delete_unmapped = undef;

GetOptions ('file=s'                  => \$file,
            'verbose'                 => \$verbose,
	    'dumpcheck'               => \$dumpcheck,
	    'useexistingmappings'     => \$use_existing_mappings,
	    'maxdump=n'               => \$maxdump,
	    'upload'                  => \$upload,
	    'delete_external_db'      => \$delete_external_db , 
	    'location=s'              => \$location,
            'logicname=s'             => \$logic_name,
	    'notriage'                => \$notriage,
	    'no_recalc_displays'      => \$no_recalc_displays,
	    'delete_unmapped'         => \$delete_unmapped,
            'help'                    => \$help,
            'man'                     => \$man );
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
  
  if(defined($dumpcheck)){
    $mapper->dumpcheck("yes");
  }
  if(defined($no_recalc_displays)){
    $mapper->no_recalc_displays("yes");
  }
  if(defined($maxdump)){
    $mapper->maxdump($maxdump);
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


#### test bit


exit();

#### end test bit


print "\nDumping xref & Ensembl sequences\n";
$mapper->dump_seqs($location);


print "\nChecking external_db table\n" if ($upload);
$mapper->upload_external_db($delete_external_db ) if $upload ;  

$mapper->build_list_and_map();

print "\nParsing mapping output\n";
$mapper->parse_mappings($notriage);

$mapper->delete_unmapped() if ($delete_unmapped);

$mapper->do_upload() if ($upload);

print "\nChecking pair data\n" if($upload);
$mapper->add_missing_pairs() if($upload);

#$mapper->check_pairs() if(!$upload);

print "\nChecking xrefs\n" if ($upload);
$mapper->cleanup_database() if ($upload);

print  "*** All finished ***\n";

sub info {

  my ($i, @all_species) = @_;

  return " for species $i of " . scalar(@all_species);

}

