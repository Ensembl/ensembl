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

##!/usr/local/ensembl/bin/perl

=head1 NAME

load_patch_fix_to_ref.pl - load patch fixes in the old
database (database with alternative assembly), to the new assembly
The patch fixes should have been integrated into the reference assembly
If the patch shares contigs with the new ref, we have a direct mapping


=head1 SYNOPSIS

load_patch_fix_to_ref.pl [arguments]

Required arguments:

  --host=HOST                 new core db host HOST
  --port=PORT                 new core db port PORT
  --user=USER                 new core db username USER
  --pass=PASS                 new core db passwort PASS
  --dbname=NAME               new core db name NAME
  --altdbname=NAME            old core db name NAME  
  --assembly=NAME             new assembly NAME
  --altassembly=NAME          old assembly NAME
  
Optional arguments:

  --conffile=filename     read parameters from FILE
                                        (default: conf/Conversion.ini)

                          (if different from --host, --port, --user, --pass):
  --althost=hOST          old core db host HOST
  --altport=PORT          old core db port PORT
  --altuser=USER          old core db username USER
  --altpass=PASS          old core db passwort PASS

  --logfile, --log=FILE               log to FILE (default: *STDOUT)
  --logpatch=PATH                      write logfile to PATH (default: .)
  --logappend, --log_append           append to logfile (default: truncate)

  -v, --verbose=0|1                   verbose logging (default: false)
  -i, --interactive=0|1               run script interactively (default: true)
  -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

This script will map patch fixes from the old assembly to the reference in the new assembly
It will load the mappings as new entries in the assembly table

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin);
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use vars qw($SERVERROOT);

BEGIN {
    $SERVERROOT = "$Bin/../";
}


$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'althost=s',
    'altport=n',
    'altuser=s',
    'altpass=s',
    'altdbname=s',
    'assembly=s',
    'altassembly=s',
);
$support->allowed_params(
    $support->get_common_params,
    'althost',
    'altport',
    'altuser',
    'altdbname',
    'assembly',
    'altassembly',
);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params(
  'altdbname',
  'assembly',
  'altassembly',
  'dbname',
  'host',
  'user'
);


#####
# connect to database and get adaptors
#
my ($dba, $dbh, $sth);

if ( !defined($support->param('pass')) ) {
    $support->param('pass', '');
}


# first set connection parameters for alternative db and test db
if ( !defined($support->param('althost')) ) { $support->param('althost',$support->param('host')); }
if ( !defined($support->param('altport')) ) { $support->param('altport',$support->param('port')); }
if ( !defined($support->param('altuser')) ) { $support->param('altuser',$support->param('user')); }
if ( !defined($support->param('altpass')) ) { $support->param('altpass',$support->param('pass')); }

# reference database
$dba->{'ref'} = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host => $support->param('host'),
                                            -user   => $support->param('user'),
                                            -pass   => $support->param('pass'),
                                            -port   => $support->param('port'),
                                            -dbname => $support->param('dbname')  );

$dbh->{'ref'} = $dba->{'ref'}->dbc->db_handle;
my $ref_helper = $dba->{'ref'}->dbc->sql_helper();

# database containing the alternative assembly
$dba->{'alt'} = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host => $support->param('althost'),
                                            -user   => $support->param('altuser'),
                                            -pass   => $support->param('altpass'),
                                            -port   => $support->param('altport'),
                                            -dbname => $support->param('altdbname')  );

$dbh->{'alt'} = $dba->{'alt'}->dbc->db_handle;
my $alt_helper = $dba->{'alt'}->dbc->sql_helper();


my (%alt_names_hash, %ref_names_hash);

my $region_sql = qq(
       SELECT s.name, seq_region_id
         FROM seq_region s, coord_system cs
        WHERE s.coord_system_id = cs.coord_system_id
          AND cs.name = 'chromosome'
          AND cs.version = ?
);

my $alt_regions = $ref_helper->execute(-SQL => $region_sql, -PARAMS => [$support->param('altassembly')]);
foreach my $alt_region (@$alt_regions) {
  $alt_names_hash{$alt_region->[0]} = $alt_region->[1];
}

my $ref_regions = $ref_helper->execute(-SQL => $region_sql, -PARAMS => [$support->param('assembly')]);
foreach my $ref_region (@$ref_regions) {
  $ref_names_hash{$ref_region->[0]} = $ref_region->[1];
}

my ($c, $patch_name, $contig_name, $chr_name, $asm_start, $asm_end, $ori);
my ($alt_asm_start, $alt_asm_end, $alt_cmp_start, $alt_cmp_end, $alt_ori);
my ($ref_asm_start, $ref_asm_end, $ref_cmp_start, $ref_cmp_end, $ref_ori);

## Get patch to contig mapping from old assembly

my $alt_sql = qq(
       SELECT s1.name, s2.name, s3.name, asm_start, asm_end, cmp_start, cmp_end, asm.ori
       FROM assembly asm, seq_region s1, coord_system cs1, seq_region s3,
            seq_region s2, coord_system cs2, assembly_exception ax
       WHERE ax.exc_type = "PATCH_FIX" AND ax.seq_region_id = s1.seq_region_id
         AND s1.seq_region_id = asm_seq_region_id
         AND s2.seq_region_id = cmp_seq_region_id
         AND s1.coord_system_id = cs1.coord_system_id AND cs1.name = "chromosome"
         AND s2.coord_system_id = cs2.coord_system_id AND cs2.name = "contig"
         AND s3.seq_region_id = ax.exc_seq_region_id;
);

## Get chromosome to contig mapping for the equivalent contig on the new assembly

my $ref_sql = qq(
        SELECT asm_start, asm_end, cmp_start, cmp_end, ori
        FROM assembly asm, seq_region s1, coord_system cs1,
            seq_region s2, coord_system cs2
        WHERE s1.seq_region_id = asm_seq_region_id
         AND s2.seq_region_id = cmp_seq_region_id
         AND s1.coord_system_id = cs1.coord_system_id AND cs1.name = "chromosome"
         AND s2.coord_system_id = cs2.coord_system_id AND cs2.name = "contig"
         AND s1.name = ? and s2.name = ?
         AND (cmp_start = ? OR cmp_end = ?);
);

my $update_sql = qq(
         INSERT IGNORE INTO assembly VALUES(?,?,?,?,?,?,?);
);


my $patch_mappings = $alt_helper->execute(-SQL => $alt_sql);

foreach my $mapping (@$patch_mappings) {

  my @mapping = @$mapping;
  $patch_name = $mapping->[0];
  $contig_name = $mapping->[1];
  $chr_name = $mapping->[2];
  $alt_asm_start = $mapping->[3];
  $alt_asm_end = $mapping->[4];
  $alt_cmp_start = $mapping->[5];
  $alt_cmp_end = $mapping->[6];
  $alt_ori = $mapping->[7];

  my $ref_mappings = $ref_helper->execute(-SQL => $ref_sql, -PARAMS => [$chr_name, $contig_name, $alt_cmp_start, $alt_cmp_end]);

  if (scalar(@$ref_mappings) == 0) {
    $support->log_stamped("Could not map $patch_name. Found no mapping in reference for $chr_name and $contig_name\n", 1);
    next;
  } else {
    foreach my $ref_mapping(@$ref_mappings) {
      $ref_asm_start = $ref_mapping->[0];
      $ref_asm_end = $ref_mapping->[1];
      $ref_cmp_start = $ref_mapping->[2];
      $ref_cmp_end = $ref_mapping->[3];
      $ref_ori = $ref_mapping->[4];

      if ($ref_ori == $alt_ori) {
        $ori = 1;
      } else {
        $ori = -1;
      }

      if ($ref_cmp_start == $alt_cmp_start && $ref_cmp_end == $alt_cmp_end) {

        $support->log_stamped("About to add direct assembly entry for $chr_name: " . $ref_names_hash{$chr_name} . ", $patch_name: " . $alt_names_hash{$patch_name} . " on coordinates $ref_asm_start, $ref_asm_end, $alt_asm_start-$alt_asm_end\n", 1);
        $c += $ref_helper->execute_update(-SQL => $update_sql, -PARAMS => [$ref_names_hash{$chr_name}, $alt_names_hash{$patch_name}, $ref_asm_start, $ref_asm_end, $alt_asm_start, $alt_asm_end, $ori]);

      } elsif ($ref_cmp_start == $alt_cmp_start) {

        if ($alt_cmp_end > $ref_cmp_end) {
          $asm_end = $alt_asm_end - ($alt_cmp_end - $ref_cmp_end);
          if (($asm_end - $alt_asm_start) != ($ref_asm_end - $ref_asm_start)) {
            $support->log_stamped("$patch_name does not have a correct mapping to $chr_name with $contig_name. $alt_asm_start-$asm_end has different length from $ref_asm_start-$ref_asm_end\n", 1);
            last;
          }
          $support->log_stamped("About to add adjusted alt end assembly entry for $chr_name: " . $ref_names_hash{$chr_name} . ", $patch_name: " . $alt_names_hash{$patch_name} . " on coordinates $ref_asm_start, $ref_asm_end, $alt_asm_start-$asm_end\n", 1);
          $c += $ref_helper->execute_update(-SQL => $update_sql, -PARAMS => [$ref_names_hash{$chr_name}, $alt_names_hash{$patch_name}, $ref_asm_start, $ref_asm_end, $alt_asm_start, $asm_end, $ori]);
        } else {
          $asm_end = $ref_asm_end - ($ref_cmp_end - $alt_cmp_end);
          if (($asm_end - $ref_asm_start) != ($alt_asm_end - $alt_asm_start)) {
            $support->log_stamped("$patch_name does not have a correct mapping to $chr_name with $contig_name. $alt_asm_start-$alt_asm_end has different length from $ref_asm_start-$asm_end\n", 1);
            last;
          }
          $support->log_stamped("About to add adjusted ref end assembly entry for $chr_name: " . $ref_names_hash{$chr_name} . ", $patch_name: " . $alt_names_hash{$patch_name} . " on coordinates $ref_asm_start, $asm_end, $alt_asm_start-$alt_asm_end\n", 1);
          $c += $ref_helper->execute_update(-SQL => $update_sql, -PARAMS => [$ref_names_hash{$chr_name}, $alt_names_hash{$patch_name}, $ref_asm_start, $asm_end, $alt_asm_start, $alt_asm_end, $ori]);
        }

      } elsif ($ref_cmp_end == $alt_cmp_end) {

        if ($alt_cmp_start < $ref_cmp_start) {
          $asm_start = $alt_asm_start - $alt_cmp_start + $ref_cmp_start;
          if (($alt_asm_end - $asm_start) != ($ref_asm_end - $ref_asm_start)) {
            $support->log_stamped("$patch_name does not have a correct mapping to $chr_name with $contig_name. $asm_start-$alt_asm_end has different length from $ref_asm_start-$ref_asm_end\n", 1);
            last;
          }
          $support->log_stamped("About to add adjusted alt start assembly entry for $chr_name: " . $ref_names_hash{$chr_name} . ", $patch_name: " . $alt_names_hash{$patch_name} . " on coordinates $ref_asm_start, $ref_asm_end, $asm_start-$alt_asm_end\n", 1);
          $c += $ref_helper->execute_update(-SQL => $update_sql, -PARAMS => [$ref_names_hash{$chr_name}, $alt_names_hash{$patch_name}, $ref_asm_start, $ref_asm_end, $asm_start, $alt_asm_end, $ori]);
        } else {
          $asm_start = $ref_asm_start - $alt_cmp_start + $ref_cmp_start;
          if (($ref_asm_end - $asm_start) != ($alt_asm_end - $alt_asm_start)) {
            $support->log_stamped("$patch_name does not have a correct mapping to $chr_name with $contig_name. $alt_asm_start-$alt_asm_end has different length from $asm_start-$ref_asm_end\n", 1);
            last;
          }
          $support->log_stamped("About to add adjusted ref start assembly entry for $chr_name: " . $ref_names_hash{$chr_name} . ", $patch_name: " . $alt_names_hash{$patch_name} . " on coordinates $asm_start, $ref_asm_end, $alt_asm_start-$alt_asm_end\n", 1);
          $c += $ref_helper->execute_update(-SQL => $update_sql, -PARAMS => [$ref_names_hash{$chr_name}, $alt_names_hash{$patch_name}, $asm_start, $ref_asm_end, $alt_asm_start, $alt_asm_end, $ori]);
        }

      } else {
        $support->log_stamped("Could not map $patch_name. Found no mapping in reference for $chr_name and $contig_name for $alt_asm_start-$alt_asm_end and $ref_asm_start-$ref_asm_end\n", 1);
      }
    }
  }
}

# finish logfile
$support->finish_log;

