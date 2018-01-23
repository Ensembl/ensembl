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

# Figure out which seq_regions are non-redundant and set the appropriate
# attribute in seq_region_attrib

use strict;
use warnings;

use DBI;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;

my ($host, $port, $user, $password, $db, $verbose, $check);
$host = "127.0.0.1";
$port = 3306;
$password = "";
$user = "ensro";

GetOptions ('host=s'      => \$host,
            'user=s'      => \$user,
            'password=s'  => \$password,
            'port=s'      => \$port,
            'db=s'        => \$db,
            'verbose'     => \$verbose,
	    'check'       => \$check,
            'help'        => sub { &show_help(); exit 1;} );

die "Host must be specified"           unless $host;
die "Database must be specified"       unless $db;

my $dbi = DBI->connect("dbi:mysql:host=$host;port=$port;database=$db", "$user", "$password", {'RaiseError' => 1})  || die "Can't connect to target DB";

# ----------------------------------------

# check that there is an entry in the attrib_type table for nonredundant
# if there is, cache the ID for later; if not, make one
my $nr_attrib_type_id;
my $sth = $dbi->prepare("SELECT attrib_type_id FROM attrib_type WHERE code='nonredundant'");
$sth->execute();
while ( my @row= $sth->fetchrow_array()) {
  $nr_attrib_type_id = $row[0];
}
$sth->finish();

if ($nr_attrib_type_id) {

  debug("Attribute with name nonredundant already set in attrib_type with attrib_type_id " . $nr_attrib_type_id);

} else {

  debug("Attribute with name nonredundant not set in attrib_type; adding ...");

  $sth = $dbi->prepare("INSERT INTO attrib_type (code, name) VALUES ('nonredundant', 'Non-redundant sequence region')");
  $sth->execute();

  $nr_attrib_type_id = $sth->{mysql_insertid};

  debug("Added nonredundant attribute with ID " . $nr_attrib_type_id);

}

# ----------------------------------------

my $db_adaptor = new Bio::EnsEMBL::DBSQL::DBAdaptor(-user   => $user,
						    -dbname => $db,
						    -host   => $host,
						    -port   => $port,
						    -pass   => $password,
						    -driver => 'mysql' );

my $slice_adaptor = $db_adaptor->get_SliceAdaptor();

# Assume all entries in the top-level co-ordinate system are non-redundant
my @toplevel = @{$slice_adaptor->fetch_all('toplevel')};
debug("Got " . @toplevel . " sequence regions in the top-level co-ordinate system");
set_nr_attribute($slice_adaptor, $nr_attrib_type_id, $dbi, @toplevel);

# Rest of the co-ordinate systems, in "ascending" order
my @coord_systems = get_coord_systems_in_order($db_adaptor, $slice_adaptor);

debug("Starting pair-wise co-ordinate system comparison");

my @nr_slices; # will store non-redundant ones for later

for (my $lower_cs_idx = 0; $lower_cs_idx < @coord_systems; $lower_cs_idx++) {
  for (my $higher_cs_idx = $lower_cs_idx+1; $higher_cs_idx < @coord_systems; $higher_cs_idx++) {

    my $higher_cs = $coord_systems[$higher_cs_idx];
    my $lower_cs = $coord_systems[$lower_cs_idx];

    debug("$lower_cs:$higher_cs");

    # we are interested in the slices that do *not* project onto the "higher" coordinate system
    my @slices = @{$slice_adaptor->fetch_all($lower_cs)};

    my $projected_hit = 0;
    my $projected_miss = 0;

    foreach my $slice (@slices) {
      my @projected = @{$slice->project($higher_cs)};
      if (@projected > 0) {
        $projected_hit++;
      } else {
        $projected_miss++;
        push @nr_slices, $slice;
      }
      undef @projected;
    }
    debug ("Projecting " . $lower_cs . " onto " . $higher_cs . ": " .
	   $projected_hit . " hit, " . $projected_miss . " miss out of " . @slices . " total");

    undef @slices;

  }
}

set_nr_attribute($slice_adaptor, $nr_attrib_type_id, $dbi, @nr_slices);

#----------------------------------------

check_non_redundant($slice_adaptor) if $check;

$sth->finish();
$dbi->disconnect();

# ----------------------------------------
# Get all the co-ordinate systems in order
# Order is descending average length
# Return an array of co-ordinate system names

sub get_coord_systems_in_order() {

  my $db_adaptor = shift;
  my $slice_adaptor = shift;

  debug("Ordering co-ordinate systems by average length");

  my $cs_adaptor = $db_adaptor->get_CoordSystemAdaptor();
  my @coord_system_objs = @{$cs_adaptor->fetch_all()};

  # Calculate average lengths
  my %lengths;
  foreach my $cs (@coord_system_objs) {
    my @slices = @{$slice_adaptor->fetch_all($cs->name())};
    my $total_len = 0;
    foreach my $slice (@slices) {
      $total_len += $slice->length();
    }
    if ($total_len > 0) {
      $lengths{$cs->name()} = $total_len / scalar(@slices);
    } else {
      $lengths{$cs->name()} = 0;
      print "Warning - total length for " . $cs->name() . " is zero!\n";
    }
  }

  my @coord_systems = sort { $lengths{$a} <=> $lengths{$b} } keys %lengths;

  foreach my $cs_name (@coord_systems) {
    debug("Co-ord system: " . $cs_name . "  Average length: " . int $lengths{$cs_name});
  }

  debug("Got co-ordinate systems in order: " . join(', ', @coord_systems));

  return @coord_systems;

}

# ----------------------------------------------------------------------
# Misc / utility functions

sub show_help {

  print "Usage: perl set_non_redundant_attribs.pl {options}\n";
  print "Where options are:\n";
  print "  --host {hostname} The database host.\n";
  print "  --user {username} The database user. Must have write permissions\n";
  print "  --password {pass} The password for user, if required.\n";
  print "  --port {folder}   The database port to use.\n";
  print "  --db {schema}     The name of the database\n";
  print "  --check           Read back non-redundant slices from SliceAdaptor\n";
  print "  --verbose         Print extra output information\n";

}

# ----------------------------------------------------------------------

sub debug {

  my $str = shift;

  print $str . "\n" if $verbose;

}

# ----------------------------------------------------------------------

# Set the "nonredundant" attribute on a Slice or group of Slices
# arg 1: SliceAdaptor
# arg 2: internal ID of 'nonredundant' attrib_type
# arg 3: DB connection
# arg 3..n: Slices

sub set_nr_attribute {

  my ($slice_adaptor, $nr_attrib_type_id, $dbi, @targets) = @_;

  debug("Setting nonredundant attribute on " . @targets . " sequence regions");

  my $sth = $dbi->prepare("INSERT INTO seq_region_attrib (seq_region_id, attrib_type_id) VALUES (?,?)");

  foreach my $slice (@targets) {

    my $seq_region_id = $slice_adaptor->get_seq_region_id($slice);

    $sth->execute($seq_region_id, $nr_attrib_type_id);

  }

  $sth->finish();

}

# ----------------------------------------------------------------------

sub check_non_redundant {

  my $slice_adaptor = shift;

  my @all = @{$slice_adaptor->fetch_all_non_redundant()};

  print "Got " . @all . " non-redundant seq_regions from SliceAdaptor\n";

}

