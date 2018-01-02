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

# $Id$

use strict;
use warnings;

# Finds all potential frameshifts (exons 1, 2, 4 or 5 bp apart)
# in a database and adds transcript attributes for them.
# Attribute value is intron number (first intron is 1, second 2 etc).

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Attribute;

use Getopt::Long;

my ($host, $port, $user, $pass, $dbpattern, $nostore, $nodelete, $print);

GetOptions('host|dbhost=s'      => \$host,
           'user|dbuser=s'      => \$user,
           'port|dbport=i'      => \$port,
           'pass|dbpass=s'      => \$pass,
           'dbpattern|dbname=s' => \$dbpattern,
       'nostore'     => \$nostore,
       'nodelete'    => \$nodelete,
       'print'       => \$print,
           'help'        => sub { usage(); exit(0); });

$port ||= 3306;

usage() if(!$user || !$dbpattern || !$host);

my $dsn = "DBI:mysql:host=$host";
$dsn .= ";port=$port" if ($port);

my $db = DBI->connect($dsn, $user, $pass);

my @dbnames = map {$_->[0] } @{ $db->selectall_arrayref( "show databases" ) };

for my $dbname ( @dbnames ) {

  next if ($dbname !~ /$dbpattern/);

  print $dbname . "\n";

  my $db_adaptor = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-host   => $host,
                               -user   => $user,
                               -pass   => $pass,
                               -dbname => $dbname,
                               -port   => $port);

  my $attribute_adaptor = $db_adaptor->get_AttributeAdaptor();
  my $transcript_adaptor = $db_adaptor->get_TranscriptAdaptor();
  my $gene_adaptor = $db_adaptor->get_GeneAdaptor();

  if (!$nodelete) {

    print STDERR "Deleting existing 'Frameshift' transcript attributes\n";
    my $dsth = $db_adaptor->dbc()->prepare("DELETE ta FROM transcript_attrib ta, attrib_type at WHERE at.attrib_type_id=ta.attrib_type_id AND at.code='Frameshift'");
    $dsth->execute();

  }

  print STDERR "Finding frameshifts in $dbname, creating transcript attributes ...\n";
  print STDERR "Attributes will not be stored in database\n" if ($nostore);

  my $count = 0;

  # get all transcripts then look at each of their introns in turn

  my @transcripts = @{$transcript_adaptor->fetch_all()};

  foreach my $transcript (@transcripts) {

    #print "Transcript " . $trans_no++ . " of " . scalar(@transcripts) . "\n";

    my $intron_number = 1;

    foreach my $intron (@{$transcript->get_all_Introns()}) {

    # only interested in the short ones
    if ($intron->length() < 6 && $intron->length() != 3) {

      print "Transcript " . $transcript->stable_id() . " intron $intron_number length " . $intron->length() . "\n"  if ($print);

      my $attribute = Bio::EnsEMBL::Attribute->new(-CODE => 'Frameshift',
                               -NAME => 'Frameshift',
                               -DESCRIPTION => 'Frameshift modelled as intron',
                               -VALUE => $intron_number);
    
      my @attribs = ($attribute);

      $attribute_adaptor->store_on_Transcript($transcript->dbID, \@attribs) if (!$nostore);

      $count++;

    }

    $intron_number++;

      } # foreach intron

    } # foreach transcript

  if ($count) {

    print "$count short intron attributes\n";
    print "Attributes not stored in database\n" if ($nostore);

  } else {

    print "No frameshift introns found!\n";

  }

}

# ----------------------------------------------------------------------

sub usage {

  print << "EOF";

  Finds all potential frameshifts (exons 1, 2 4 or 5 bp apart) in a database 
  and adds transcript  attributes for them. Attribute value is intron length.

  perl $0 {options}

 Options ([..] indicates optional):

   --host       The database server to connect to.

   [--port]     The port to use. Defaults to 3306.

   --user       Database username. Must allow writing.

   --pass       Password for user.

   --dbpattern  Regular expression to define which databases are affected.

  [--nostore]   Don't store the attributes, just print results.

  [--nodelete]  Don't delete any existing "Frameshift" attributes before creating new ones.

  [--print]     Print transcript stable ID, intron number and length.

  [--help]      This text.


EOF

  exit(0);

}

# ----------------------------------------------------------------------
