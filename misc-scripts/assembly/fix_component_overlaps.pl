#!/usr/local/ensembl/bin/perl

=head1 NAME

fix_overlaps.pl - remove overlapping mappings between two closely related
assemblies

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

Optional arguments:

  --chromosomes, --chr=LIST           only process LIST toplevel seq_regions

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

This script removes overlapping component mappings from the assembly table.

=head1 RELATED FILES

The whole process of creating a whole genome alignment between two assemblies
is done by a series of scripts. Please see

  ensembl/misc-scripts/assembly/README

for a high-level description of this process, and POD in the individual scripts
for the details.


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
    'altassembly'
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
  ORDER BY a.cmp_start
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
    
    # look for overlaps with last segment
    if ($last->{'cmp_end'} >= $r->{'cmp_start'}) {

      $i++;
      
      # debug warnings
      $support->log_verbose('last:   '.sprintf($fmt1, map { $last->{$_} }
        qw(asm_start asm_end cmp_start cmp_end ori)), 1);
      $support->log_verbose('before: '.sprintf($fmt1, map { $r->{$_} }
        qw(asm_start asm_end cmp_start cmp_end ori)), 1);

      # skip if this segment ends before the last one
      if ($r->{'cmp_end'} <= $last->{'cmp_end'}) {
        $support->log_verbose("skipped\n\n", 1);
        next;
      }
    
      my $overlap = $last->{'cmp_end'} - $r->{'cmp_start'} + 1;

      $r->{'cmp_start'} += $overlap;

      if ($r->{'ori'} == -1) {
        $r->{'asm_end'} -= $overlap;
      } else {
        $r->{'asm_start'} += $overlap;
      }
      
      $support->log_verbose('after:  '.sprintf($fmt1, map { $r->{$_} }
        qw(asm_start asm_end cmp_start cmp_end ori))."\n", 1);
    }

    push @rows, $r;
    $last = $r;
  }

  $support->log("Fixed $i mappings.\n", 1);


  if (!$support->param('dry_run') && $i > 0  ) {

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
      $sql = qq(INSERT INTO assembly VALUES (?, ?, ?, ?, ?, ?, ?));
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

