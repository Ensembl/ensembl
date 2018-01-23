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


my ($file, $user, $password, $verbose, $force, $help, $schema, $vega_schema, $limit);

GetOptions ('file=s'      => \$file,
            'schema=s'    => \$schema,
			'vega_schema=s' => \$vega_schema,
            'user=s'      => \$user,
            'password=s'  => \$password,
            'verbose'     => \$verbose,
            'force'       => \$force,
            'limit=s'     => \$limit,
            'help'        => sub { &show_help(); exit 1;} );

usage("-file option is required")   if(!$file);
usage("-schema option is required") if(!$schema);
usage() if($help);

open(FILE, $file) or die("Could not open input file '$file'");

my @all_species_converters;

while( my $line = <FILE> ) {
  chomp($line);
  next if $line =~ /^#/;
  next if !$line;

  my ( $species, $host, $source_db_name, $target_db_name ) = 
    split( "\t", $line );

  my $converter;
  if ($vega_schema) {
      $species = "vega::".$species;
      eval "require SeqStoreConverter::$species";
      if($@) {
	  warn("Could not require conversion module SeqStoreConverter::$species\n" .
	       "Using SeqStoreConverter::BasicConverter instead:\n$@");
	  require SeqStoreConverter::BasicConverter;
	  $species = "BasicConverter";
      }
      else {
	  warn "Using conversion module SeqStoreConverter::$species\n";
      }
  }
  
  else {
	eval "require SeqStoreConverter::$species";
	if($@) {
	  warn("Could not require conversion module SeqStoreConverter::$species\n" .
		   "Using SeqStoreConverter::BasicConverter instead:\n$@");
	  require SeqStoreConverter::BasicConverter;
	  $species = "BasicConverter";
	}
  }

  {
	  no strict 'refs';
	  $converter = "SeqStoreConverter::$species"->new
		  ( $user, $password, $host, $source_db_name, $target_db_name, 
			$schema, $vega_schema, $force, $verbose, $limit );
  }

  push @all_species_converters, $converter;

}

for my $converter ( @all_species_converters ) {
  $converter->debug( "\n\n*** converting " . $converter->source . " to " . 
                     $converter->target() . " ***");
  $converter->transfer_meta();
  $converter->create_coord_systems();
  $converter->create_seq_regions();
  $converter->create_assembly();
  $converter->create_attribs();
  $converter->set_top_level();

  $converter->transfer_dna();
  $converter->transfer_genes();
  $converter->transfer_prediction_transcripts();
  $converter->transfer_features();
  if ($vega_schema) {
      $converter->transfer_vega_stable_ids();
  } else {
      $converter->transfer_stable_ids();
  }
  $converter->copy_other_tables();
  $converter->copy_repeat_consensus();
  $converter->create_meta_coord();
  if ($vega_schema) {
      $converter->update_clone_info();
      $converter->remove_supercontigs();
      $converter->copy_internal_clone_names();
  }
}


print STDERR "*** All finished ***\n";

sub usage {
  my $msg = shift;

  print STDERR "$msg\n\n" if($msg);

  print STDERR <<EOF;
usage:   perl convert_seqstore <options>

options: -file <input_file>     input file with tab delimited 'species',
                                'host', 'source_db', 'target_db' values
                                on each line

         -schema <table_file>   file containing SQL schema definition

         -vega_schema <table_file> file containing SQL for additional Vega tables

         -user <user>           a mysql db user with read/write priveleges

         -password <password>   the mysql user's password

         -verbose               print out debug statements

         -force                 replace any target dbs that already exists

         -limit <num_rows>      limit the number of features transfered to
                                speed up testing

         -help                  display this message

example: perl convert_seqstore.pl -file converter.input \\
              -schema ../../sql/table.sql -user ensadmin -password secret \\
              -force -verbose

EOF
#'

  exit;
}
