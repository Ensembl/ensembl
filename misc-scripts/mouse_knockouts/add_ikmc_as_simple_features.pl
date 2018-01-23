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


#
# updates the external db tables on all of the core databases on a given host
#

use strict;
use warnings;


use Getopt::Long;
use DBI;
use IO::File;

my ( $host, $user, $pass,   $port, $dbname, $file);

GetOptions( "dbhost|host=s",        \$host,
            "dbuser|user=s",        \$user,
            "dbpass|pass=s",        \$pass,
            "dbport|port=i",        \$port,
            "file=s",        \$file,
            "dbname=s",    \$dbname);

$port ||= 3306;
$host ||= "ens-research";
$dbname ||= "ianl_mus_musculus_core_58_37k";

# connect to the database
my $dsn = "DBI:mysql:host=$host;port=$port;database=$dbname";

my $db = DBI->connect( $dsn, $user, $pass, {RaiseError => 1} );


my %analysis = ("No products available yet" => 0, 
                "ES cells available" => 0,
                "Mice available" => 0,
                "Vector available" => 0);

# if the 4 analysis do not exist create them
# if they do exist delete entrys in the simple feature table


my $find_sth = $db->prepare("SELECT  analysis_id FROM analysis where logic_name like ?")|| die "Could not prepare find_sth";

my $create_sth = $db->prepare("INSERT into analysis (logic_name) values (?)") || die "Could not prepare create_sth";

my $delete_sth = $db->prepare("DELETE from simple_feature where analysis_id = ?" )|| die "Could not prepare delete_sth";

foreach my $anal (keys %analysis){
  my $logic_name = "IKMC_".$anal;
  $logic_name =~ s/ /_/g;

  my $id =undef;
  $find_sth->execute($logic_name);
  $find_sth->bind_columns(\$id);
  $find_sth->fetch;
  if(defined($id)){ # delete existing simple features
    print STDERR "Analysis $logic_name already exists so clearing entrys for this ($id) in the simple_feature table\n";
    $delete_sth->execute($id) || die "Could not delete form simple table for id = $id";
  }
  else{             # create analysis
    $create_sth->execute($logic_name) || die "Could not create new analysis $logic_name";
    $id = $create_sth->{'mysql_insertid'};
    print STDERR "Creating new anlysis $logic_name ($id)\n";
  }
  $analysis{$anal} = $id;
}

$find_sth->finish();
$create_sth->finish();
$delete_sth->finish();


# make a hashes for gene 
# gene2seqregion{stable_id} = seq_region_id
# gene2start....
# gene2end  ....
my %gene2seqregion;
my %gene2start;
my %gene2end;
my %gene2strand;

my $gene_sth = $db->prepare("SELECT  stable_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand FROM gene") || die "Could not prepare gene_sth";

$gene_sth->execute();
my ($stable_id, $seq, $start, $end, $strand);
$gene_sth->bind_columns(\$stable_id, \$seq, \$start, \$end, \$strand);
while($gene_sth->fetch()){
  $gene2seqregion{$stable_id} = $seq;
  $gene2start{$stable_id}     = $start;
  $gene2end{$stable_id}       = $end;
  $gene2strand{$stable_id}    = $strand;
}
#process the file and add new simple features

my $insert_sth = $db->prepare("INSERT INTO simple_feature (seq_region_id, seq_region_start, seq_region_end, seq_region_strand, display_label, analysis_id ) VALUES(?, ?, ?, ?, ?, ?)") || die "Could not prepare insert_sth";

my $ikmc = get_filehandle($file);
my %count;
while ( $_ = $ikmc->getline() ) {

  chomp;
  my ($mgi, $label, $type, $stable_id) = split /\t/;

  if((defined($stable_id) and $stable_id) and defined($gene2seqregion{$stable_id})){
    if(!defined($analysis{$type})){
      print STDERR $_."\nUnknown type *$type*\n";
    }
    else{
      $insert_sth->execute($gene2seqregion{$stable_id}, 
                           $gene2start{$stable_id},
                           $gene2end{$stable_id},
                           $gene2strand{$stable_id}, 
                           $mgi,
                           $analysis{$type}) || die "Could not insert new values";
      $count{$type}++;
    }
  }
  else{
    if(!defined($stable_id) or !$stable_id){
      $count{"no stable_id"}++;
    }
    else{
      print STDERR "Could not find stable id $stable_id\n";
    }
  }
}

foreach my $key (keys %count){
  print $key."\t".$count{$key}."\n";
}


sub get_filehandle
{
    my ($file_name) = @_;

    my $io;

    my $alt_file_name = $file_name;
    $alt_file_name =~ s/\.(gz|Z)$//;

    if ( $alt_file_name eq $file_name ) {
        $alt_file_name .= '.gz';
    }

    if ( !-f $file_name ) {
        carp(   "File '$file_name' does not exist, "
              . "will try '$alt_file_name'" );
        $file_name = $alt_file_name;
    }

    if ( $file_name =~ /\.(gz|Z)$/ ) {
        # Read from zcat pipe
        $io = IO::File->new("zcat $file_name |")
          or carp("Can not open file '$file_name' with 'zcat'");
    } else {
        # Read file normally
        $io = IO::File->new($file_name)
          or carp("Can not open file '$file_name'");
    }

    if ( !defined $io ) { return undef }

    print "Reading from '$file_name'...\n";

    return $io;
}
