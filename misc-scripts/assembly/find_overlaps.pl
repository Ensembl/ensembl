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

find_overlaps.pl - find assembly entries with overlapping regions

=head1 SYNOPSIS

find_overlaps.pl [arguments]

Required arguments:

  --dbname, db_name=NAME              database name NAME
  --host, --dbhost, --db_host=HOST    database host HOST
  --port, --dbport, --db_port=PORT    database port PORT
  --user, --dbuser, --db_user=USER    database username USER
  --pass, --dbpass, --db_pass=PASS    database passwort PASS

Optional arguments:

  --conffile, --conf=FILE             read parameters from FILE
                                      (default: conf/Conversion.ini)

  --logfile, --log=FILE               log to FILE (default: *STDOUT)
  --logpath=PATH                      write logfile to PATH (default: .)
  --logappend, --log_append           append to logfile (default: truncate)

  -v, --verbose=0|1                   verbose logging (default: false)
  -i, --interactive=0|1               run script interactively (default: true)
  -n, --dry_run, --dry=0|1            don't write results to database
  -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

This script looks for overlapping assembly and component seq_regions in the assembly table of
an Ensembl core-style database. It prints the number of overlapping entries per
coordinate system and optionally (when run with the --verbose option) a list of
the overlapping assembly entries.

The script is intended to detect and debug problems in the assembly. It
complements the AssemblyMultipleOverlap healthcheck.


=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<http://lists.ensembl.org/mailman/listinfo/dev>

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin);
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport("$Bin/../../..");

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
);
$support->allowed_params(
    $support->get_common_params,
);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

# database connection
my $dba = $support->get_database('ensembl');
my $dbh = $dba->dbc->db_handle;

my $sql = qq(
  SELECT cs1.name AS asm_name, cs1.version AS asm_version,
         cs2.name AS cmp_name, cs2.version AS cmp_version
  FROM assembly a, seq_region sr1, seq_region sr2,
       coord_system cs1, coord_system cs2
  WHERE a.asm_seq_region_id = sr1.seq_region_id
  AND sr1.coord_system_id = cs1.coord_system_id
  AND a.cmp_seq_region_id = sr2.seq_region_id
  AND sr2.coord_system_id = cs2.coord_system_id
  GROUP BY asm_name, asm_version, cmp_name, cmp_version
);


my $fmt1 = "%6s %6s %10s %10s %10s %10s %3s %3s %3s\n";


foreach  my $type ( qw(asm cmp)){
my $sth = $dbh->prepare($sql);
$sth->execute;
while (my ($asm_name, $asm_version, $cmp_name, $cmp_version) =
  $sth->fetchrow_array) {

  my $sql1 = qq(
    SELECT a.*, sr1.name as asm_sr_name, sr2.name as cmp_sr_name
    FROM assembly a, seq_region sr1, seq_region sr2,
         coord_system cs1, coord_system cs2
    WHERE a.asm_seq_region_id = sr1.seq_region_id
    AND a.cmp_seq_region_id = sr2.seq_region_id
    AND sr1.coord_system_id = cs1.coord_system_id
    AND sr2.coord_system_id = cs2.coord_system_id
    AND cs1.name = '$asm_name'
    AND cs2.name = '$cmp_name'
  );

  if ($asm_version) {
    $sql1 .= " AND cs1.version = '$asm_version'";
  } else {
    $sql1 .= " AND cs1.version IS NULL";
    $asm_version = 'NULL';
  }
  
  if ($cmp_version) {
    $sql1 .= " AND cs2.version = '$cmp_version'";
  } else {
    $sql1 .= " AND cs2.version IS NULL";
    $cmp_version = 'NULL';
  }
  

  $sql1 .= " ORDER BY a.".$type."_seq_region_id, a.".$type."_start ASC, a.".$type."_end";

  $support->log_stamped("$asm_name.$asm_version $cmp_name.$cmp_version\n");

  my $sth1 = $dbh->prepare($sql1);
  $sth1->execute;

  # do an initial fetch
  my $last = $sth1->fetchrow_hashref;
#  foreach my $key (qw(asm_seq_region_id asm_sr_name cmp_seq_region_id cmp_sr_name asm_start asm_end cmp_start cmp_end ori)){
#    print $key." = ".$last->{$key}."\t";
#  }
#  print "\n";  
  my $i = 0;
  
  $support->log_verbose($type . ' overlaps:\n');
  while ($last and (my $r = $sth1->fetchrow_hashref)) {

#    foreach my $key (qw(asm_seq_region_id asm_sr_name cmp_seq_region_id cmp_sr_name asm_start asm_end cmp_start cmp_end ori)){
#      print $key." = ".$r->{$key}."\t";
#    }
#    print "\n";
        
    # look for overlaps with last segment
    if ($last->{$type.'_seq_region_id'} == $r->{$type.'_seq_region_id'} and
        $last->{$type.'_end'} >= $r->{$type.'_start'}) {

      $i++;
      
      # debug warnings
      $support->log_verbose('last:'.sprintf($fmt1, map { $last->{$_} }
        qw(asm_seq_region_id asm_sr_name cmp_seq_region_id cmp_sr_name
        asm_start asm_end cmp_start cmp_end ori)), 1);
      $support->log_verbose('this:'.sprintf($fmt1, map { $r->{$_} }
        qw(asm_seq_region_id asm_sr_name cmp_seq_region_id cmp_sr_name
        asm_start asm_end cmp_start cmp_end ori))."\n", 1);

    }

    $last = $r;
  }
    
  $sth1->finish;

  $support->log("\nFound $i overlaps.\n\n", 1);
}
$sth->finish;
}

# finish logfile
$support->finish_log;

