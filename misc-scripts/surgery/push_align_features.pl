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

use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Getopt::Long;


#my $LIMIT = ' LIMIT 1000 ';
my $LIMIT = '';

my ($db, $dnadbname);

{
  my ($dbname, $host, $port, $user, $pass);

  GetOptions('dbname=s' => \$dbname,
             'dnadbname=s' => \$dnadbname,
             'user=s'   => \$user,
             'host=s'   => \$host,
             'port=i'   => \$port,
             'pass=s'   => \$pass);

  $port ||= 3306;

  usage() if(!$host || !$user || !$dbname);

  $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (-host => $host,
     -user => $user,
     -pass => $pass,
     -port => $port,
     -dbname => $dbname,
     -dnadbname => $dnadbname);

  if($dnadbname) {
    my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
      (-host => $host,
       -user => $user,
       -pass => $pass,
       -port => $port,
       -dbname => $dnadbname);
    $db->dnadb($dnadb);
  } else {
    $dnadbname = $dbname;
  }
}

my @daf_css = get_feat_coord_systems($db, 'dna_align_feature');
my @paf_css = get_feat_coord_systems($db, 'protein_align_feature');

### convert dna align features

my $sth = $db->dbc()->prepare("SHOW CREATE TABLE dna_align_feature");
$sth->execute();
my $create_table = $sth->fetchrow_arrayref()->[1];
$sth->finish();


$create_table =~ s/CREATE TABLE dna_align_feature/CREATE TABLE tmp_dna_align_feature/;
# difference b/w mysql 4 and mysql 3
$create_table =~ s/CREATE TABLE `dna_align_feature`/CREATE TABLE tmp_dna_align_feature/;

$sth = $db->dbc()->prepare($create_table);
$sth->execute();

$sth = $db->dbc()->prepare
    (qq{INSERT INTO tmp_dna_align_feature
        SELECT daf.dna_align_feature_id,
          a.asm_seq_region_id,
          if(a.ori = 1,
            (a.asm_start + daf.seq_region_start - a.cmp_start),
            (a.asm_start + a.cmp_end - daf.seq_region_end))
            as seq_region_start,
          if(a.ori = 1,
            (a.asm_start + daf.seq_region_end - a.cmp_start),
            (a.asm_end   + a.cmp_end - daf.seq_region_start))
            as seq_region_end,
          a.ori * daf.seq_region_strand as seq_region_strand,
          daf.hit_start, daf.hit_end, daf.hit_strand, daf.hit_name,
          daf.analysis_id, daf.score, daf.evalue, daf.perc_ident,
          daf.cigar_line
        FROM $dnadbname.assembly a,
             $dnadbname.seq_region asm_sr,
             $dnadbname.seq_region cmp_sr,
             dna_align_feature daf,
             $dnadbname.seq_region_attrib sra,
             $dnadbname.attrib_type at
        WHERE asm_sr.coord_system_id = ? AND cmp_sr.coord_system_id = ?
        AND   daf.seq_region_id = a.cmp_seq_region_id
        AND   a.asm_seq_region_id = asm_sr.seq_region_id
        AND   a.cmp_seq_region_id = cmp_sr.seq_region_id
        AND   daf.seq_region_start >= a.cmp_start
        AND   daf.seq_region_end <= a.cmp_end
        AND   asm_sr.seq_region_id = sra.seq_region_id
        AND   sra.attrib_type_id = at.attrib_type_id
        AND   at.code = 'toplevel'
        $LIMIT});

foreach my $cs_pair (@daf_css) {
  my ($top_cs, $feat_cs) = @$cs_pair;
  print STDERR "Converting dna align features between from coord system ",
    $feat_cs->name(), " to coord system ", $top_cs->name(), "\n";
  $sth->execute($top_cs->dbID(), $feat_cs->dbID());
}
$sth->finish();


### convert protein align features

$sth = $db->dbc()->prepare("SHOW CREATE TABLE protein_align_feature");
$sth->execute();
$create_table = $sth->fetchrow_arrayref()->[1];
$sth->finish();

$create_table =~ s/CREATE TABLE protein_align_feature/CREATE TABLE tmp_protein_align_feature/;
$create_table =~ s/CREATE TABLE `protein_align_feature`/CREATE TABLE tmp_protein_align_feature/;


$sth = $db->dbc()->prepare($create_table);
$sth->execute();

$sth = $db->dbc()->prepare
    (qq{INSERT INTO tmp_protein_align_feature
        SELECT paf.protein_align_feature_id,
          a.asm_seq_region_id,
          if(a.ori = 1,
            (a.asm_start + paf.seq_region_start - a.cmp_start),
            (a.asm_start + a.cmp_end - paf.seq_region_end))
            as seq_region_start,
          if(a.ori = 1,
            (a.asm_start + paf.seq_region_end - a.cmp_start),
            (a.asm_end   + a.cmp_end - paf.seq_region_start))
            as seq_region_end,
          a.ori * paf.seq_region_strand as seq_region_strand,
          paf.hit_start, paf.hit_end, paf.hit_name, paf.analysis_id,
          paf.score, paf.evalue, paf.perc_ident, paf.cigar_line
        FROM $dnadbname.assembly a,
             $dnadbname.seq_region asm_sr,
             $dnadbname.seq_region cmp_sr,
             protein_align_feature paf,
             $dnadbname.seq_region_attrib sra,
             $dnadbname.attrib_type at
        WHERE asm_sr.coord_system_id = ? AND cmp_sr.coord_system_id = ?
        AND   paf.seq_region_id = a.cmp_seq_region_id
        AND   a.asm_seq_region_id = asm_sr.seq_region_id
        AND   a.cmp_seq_region_id = cmp_sr.seq_region_id
        AND   paf.seq_region_start >= a.cmp_start
        AND   paf.seq_region_end <= a.cmp_end
        AND   asm_sr.seq_region_id = sra.seq_region_id
        AND   sra.attrib_type_id = at.attrib_type_id
        AND   at.code = 'toplevel'
        $LIMIT});

foreach my $cs_pair (@paf_css) {
  my ($top_cs, $feat_cs) = @$cs_pair;
  print STDERR "Converting protein align features between from coord system ",
    $feat_cs->name(), " to coord system ", $top_cs->name(), "\n";
  $sth->execute($top_cs->dbID(), $feat_cs->dbID());
}
$sth->finish();


print STDERR "Replacing existing align feature tables with new tables\n";

$db->dbc->do("drop table protein_align_feature");
$db->dbc->do("drop table dna_align_feature");
$db->dbc->do("alter table tmp_protein_align_feature rename protein_align_feature");
$db->dbc->do("alter table tmp_dna_align_feature rename dna_align_feature");


print STDERR "Updating meta_coord table\n";

$db->dbc->do("delete from meta_coord where table_name = 'dna_align_feature'");
$db->dbc->do("delete from meta_coord where table_name = 'protein_align_feature'");


$sth = $db->dbc->prepare("INSERT INTO meta_coord set table_name = ?, " .
                         "coord_system_id = ?");


my %seen = ();
foreach my $cs_pair (@daf_css) {
  my ($top_cs) = @$cs_pair;
  if(!$seen{$top_cs->dbID()}) {
    $sth->execute('dna_align_feature', $top_cs->dbID());
    $seen{$top_cs->dbID()} = 1;
  }
}

%seen = ();
foreach my $cs_pair (@paf_css) {
  my ($top_cs) = @$cs_pair;
  if(!$seen{$top_cs->dbID()}) {
    $sth->execute('protein_align_feature', $top_cs->dbID());
    $seen{$top_cs->dbID()} = 1;
  }
}

$sth->finish();



sub get_feat_coord_systems {
  my $db = shift;
  my $type = shift;

  my $csa = $db->get_CoordSystemAdaptor();
  my $mcc = $db->get_MetaCoordContainer();

  # determine what coordinate systems have top-level seq_regions

  my $sth = $db->dbc->prepare
    (qq{SELECT distinct(sr.coord_system_id)
        FROM   seq_region sr, seq_region_attrib sra, attrib_type at
        WHERE  sr.seq_region_id = sra.seq_region_id
        AND    sra.attrib_type_id = at.attrib_type_id
        AND    at.code = 'toplevel'});

  $sth->execute();

  my @top_css = map {$csa->fetch_by_dbID($_->[0])}
                @{$sth->fetchall_arrayref()};


  $sth->finish();

  # determine what coord systems features are stored in

  my @feat_css = @{$mcc->fetch_all_CoordSystems_by_feature_type($type)};

  my @css;

  foreach my $feat_cs (@feat_css) {
    foreach my $top_cs (@top_css) {
      if(!$feat_cs->equals($top_cs)) {
        if($feat_cs->rank() > $top_cs->rank()) {
          my @mp = @{$csa->get_mapping_path($top_cs, $feat_cs)};
          if(@mp == 2) {
            if($mp[0]->equals($top_cs)) {
              push @css, \@mp;
            } else {
              die("Unexpected: ". $top_cs->name() .
                  " coord system should not be component of " .
                  $feat_cs->name(), " coord system.");
            }
          }
        } else {
          die("There is no 1 step mapping path between coord systems " .
              $feat_cs->name()." and ".$top_cs->name().". This script " .
              "requires one step mapping paths.");
        }
      }
    }
  }


  if(!@css) {
    die("Nothing to do for ${type}s (or missing meta/assembly info)\n");
  }

  return @css;
}



sub usage {
  print STDERR
qq{
This program will convert the coordinates of alignment features from
their current cooridinate systems to coordinates on toplevel sequence
regions.  This can speed up the retrieval of these features dramatically.

Features which are partially or entirely non-golden portions of
sequence regions will be discarded.

This program is only capable of mapping if there is a direct relationship
between the coordinate systems defined in the assembly and meta tables.

This program can use a different database to obtain the assembly information
but the database must be on the same MySQL instance.  This is useful for
satellite databases such as the est database.

usage:
  perl push_align_features -dbname <dbname> -user <user> -host <host> \
                           [-port <port>] [-pass <pass>] [-dnadbname <dbname>]
};

  exit;
}
