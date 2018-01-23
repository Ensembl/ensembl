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

merge_cmp_gaps.pl - merge small gaps between mappings and bridge mappings without gaps

=head1 SYNOPSIS

fix_overlaps.pl [arguments]

Required arguments:

  --dbname, db_name=NAME              database name NAME
  --host, --dbhost, --db_host=HOST    database host HOST
  --port, --dbport, --db_port=PORT    database port PORT
  --user, --dbuser, --db_user=USER    database username USER
  --pass, --dbpass, --db_pass=PASS    database passwort PASS
  --assembly=ASSEMBLY                 assembly version ASSEMBLY

  --altassembly=ASSEMBLY              alternative assembly version ASSEMBLY
  --chromosomes, --chr=LIST           only process LIST seq_regions

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
    
The script merges adjacent assembly segments which can result from alternating
alignments from clone identity and blastz alignment.


=head1 AUTHOR

Monika Komorowska <monika@ebi.ac.uk>, Ensembl core API team

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
    'assembly=s',
    'altassembly=s',
    'chromosomes|chr=s@',
);
$support->allowed_params(
    $support->get_common_params,
    'assembly',
    'altassembly',
    'chromosomes',
);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

$support->comma_to_list('chromosomes');

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params(
    'assembly',
    'altassembly',
    'chromosomes'
);

# database connection
my $dba = $support->get_database('ensembl');
my $dbh = $dba->dbc->db_handle;

my $assembly = $support->param('assembly');
my $altassembly = $support->param('altassembly');

my $sql = qq(
  SELECT a.*
  FROM assembly a, seq_region sr1, seq_region sr2,
       coord_system cs1, coord_system cs2
  WHERE a.asm_seq_region_id = sr1.seq_region_id
  AND a.cmp_seq_region_id = sr2.seq_region_id
  AND sr1.coord_system_id = cs1.coord_system_id
  AND sr2.coord_system_id = cs2.coord_system_id
  AND cs1.version = '$assembly'
  AND cs2.version = '$altassembly'
  AND sr2.name = ?
  ORDER BY a.ori, a.cmp_start
);

my $sth = $dbh->prepare($sql);

my $fmt1 = "%10s %10s %10s %10s %3s\n";


foreach my $chr ($support->param('chromosomes')) {
  

  $support->log_stamped("\nToplevel seq_region $chr...\n");

  $sth->execute($chr);

  my @rows = ();

  # do an initial fetch
  my $last = $sth->fetchrow_hashref;

  # skip seq_regions for which we don't have data
  unless ($last) {
    $support->log("No mappings found. Skipping.\n", 1);
    next;
  }
  
  push @rows, $last;
  
  my $i = 0;
  my $j = 0;
  my $k = 0;
  
  while ($last and (my $r = $sth->fetchrow_hashref)) {
    # merge adjacent assembly segments (these can result from alternating
    # alignments from clone identity and blastz alignment)
    if ($last->{'asm_end'} == ($r->{'asm_start'} - 1) and
        $last->{'cmp_end'} == ($r->{'cmp_start'} - 1)) {

      $j++;

      # debug warnings
      $support->log_verbose('merging - last: '.sprintf($fmt1,
        map { $last->{$_} } qw(asm_start asm_end cmp_start cmp_end ori)), 1);
      $support->log_verbose('this:           '.sprintf($fmt1, map { $r->{$_} }
        qw(asm_start asm_end cmp_start cmp_end ori)), 1);

      # remove last row
      pop(@rows);

      # merge segments and add new row
      $last->{'asm_end'} = $r->{'asm_end'};
      $last->{'cmp_end'} = $r->{'cmp_end'};
      push @rows, $last;

      next;
    }
    
    # bridge small gaps (again, these can result from alternating alignments
    # from clone identity and blastz alignment). A maximum gap size of 10bp is
    # allowed
    my $asm_gap = $r->{'asm_start'} - $last->{'asm_end'} - 1;
    my $cmp_gap = $r->{'cmp_start'} - $last->{'cmp_end'} - 1;

    if ($asm_gap == $cmp_gap and $asm_gap <= 10 and $asm_gap > 0) {

      $k++;

      # debug warnings
      $support->log_verbose('bridging - last: '.sprintf($fmt1,
        map { $last->{$_} } qw(asm_start asm_end cmp_start cmp_end ori)), 1);
      $support->log_verbose('this:            '.sprintf($fmt1, map { $r->{$_} }
        qw(asm_start asm_end cmp_start cmp_end ori)), 1);

      # remove last row
      pop(@rows);

      # merge segments and add new row
      $last->{'asm_end'} = $r->{'asm_end'};
      $last->{'cmp_end'} = $r->{'cmp_end'};
      push @rows, $last;

      next;
    }

    push @rows, $r;
    $last = $r;
  }

  $support->log("Merged $j mappings.\n", 1);
  $support->log("Bridged $k gaps.\n", 1);

  if ((!$support->param('dry_run')) && ($j + $k > 0) ) {

      # delete all current mappings from the db and insert the corrected ones
      my $c = $dbh->do(qq(
    DELETE a
    FROM assembly a, seq_region sr1, seq_region sr2,
         coord_system cs1, coord_system cs2
    WHERE a.asm_seq_region_id = sr1.seq_region_id
    AND a.cmp_seq_region_id = sr2.seq_region_id
    AND sr1.coord_system_id = cs1.coord_system_id
    AND sr2.coord_system_id = cs2.coord_system_id
    AND cs1.version = '$assembly'
    AND cs2.version = '$altassembly'
    AND sr2.name = '$chr'
  ));

      $support->log("\nDeleted $c entries from the assembly table.\n");

      # now insert the fixed entries
      $sql = qq(INSERT IGNORE INTO assembly VALUES (?, ?, ?, ?, ?, ?, ?));
      my $sth1 = $dbh->prepare($sql);
  
      foreach my $r (@rows) {
	  $sth1->execute(map { $r->{$_} } qw(asm_seq_region_id cmp_seq_region_id asm_start asm_end cmp_start cmp_end ori));
      }

      $support->log("Added ".scalar(@rows)." fixed entries to the assembly table.\n");

  }


}

$sth->finish;


# finish logfile
$support->finish_log;

