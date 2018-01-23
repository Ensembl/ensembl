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

align_by_clone_identity.pl - create a whole genome alignment between two closely
related assemblies, step 1

=head1 SYNOPSIS

align_by_clone_identity.pl [arguments]

Required arguments:

    --dbname, db_name=NAME              database name NAME
    --host, --dbhost, --db_host=HOST    database host HOST
    --port, --dbport, --db_port=PORT    database port PORT
    --user, --dbuser, --db_user=USER    database username USER
    --pass, --dbpass, --db_pass=PASS    database passwort PASS
    --assembly=ASSEMBLY                 assembly version ASSEMBLY

    --altdbname=NAME                    alternative database NAME
    --altassembly=ASSEMBLY              alternative assembly version ASSEMBLY

Optional arguments:

    --chromosomes, --chr=LIST           only process LIST toplevel seq_regions
    --skipclones=FILE                   read list of clones to skip from FILE

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

This script is part of a series of scripts to create a mapping between two
assemblies. It assembles the toplevel coordinate systems of two different
assemblies of a genome by creating a whole genome alignment between the two.

The process assumes that the two assemblies are reasonably similar, i.e. there
are no major rearrangements or clones moved from one toplevel seq_region to
another.

See "Related scripts" below for an overview of the whole process.

This particular script creates a whole genome alignment between two closely
related assemblies. You will need a database containing the reference assembly
and the alternative toplevel seq_regions which can be created using
load_alternative_assembly.pl.

The alignment is created in two steps:

    1. Match clones with same name and version directly and create alignment
       blocks for these regions. Clones can be tagged manually to be excluded
       from these direct matches by listing them in a file of clones to skip
       (--skipclones argument). This can be useful to get better results in
       regions with major assembly differences.
       
       The result is stored in the assembly table as an assembly between the
       toplevel seq_regions of both genome assemblies.

    2. Store non-aligned blocks in a temporary table (tmp_align). They can
       later be aligned using blastz by align_nonident_regions.pl.

=head1 RELATED FILES

The whole process of creating a whole genome alignment between two assemblies
is done by a series of scripts. Please see

  ensembl/misc-scripts/assembly/README

for a high-level description of this process, and POD in the individual scripts
for the details.


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
use vars qw($SERVERROOT);

BEGIN {
    $SERVERROOT = "$Bin/../../..";
    unshift(@INC, "$SERVERROOT/ensembl/modules");
    unshift(@INC, "$SERVERROOT/bioperl-live");
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::Attribute;

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'assembly=s',
    'althost=s',
    'altport=i',
    'altuser=s',
    'altpass=s',
    'altdbname=s',
    'altassembly=s',
    'chromosomes|chr=s@',
    'skipclones|skip_clones=s',
);
$support->allowed_params(
    $support->get_common_params,
    'assembly',
    'althost',
    'altport',
    'altuser',
    'altpass',
    'altdbname',
    'altassembly',
    'chromosomes',
    'skipclones',
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
    'altdbname',
    'altassembly'
);

# suggest to run non-verbosely
my $txt = qq(Running this script with the --verbose option will create a lot of output.
    It is recommended to do this only for debug purposes.
    Shall I switch to non-verbose logging for you?);

if ($support->param('verbose') and $support->user_proceed($txt)) {
  $support->param('verbose', 0);
}

#####
# connect to database and get adaptors
#
my ($dba, $dbh, $sql, $sth);

# first set connection parameters for alternative db if not different from
# reference db
map { $support->param("alt$_", $support->param($_)) unless ($support->param("alt$_")) } qw(host port user);

# reference database
my $R_dba = $support->get_database('ensembl');
my $R_dbh = $R_dba->dbc->db_handle;
my $R_sa = $R_dba->get_SliceAdaptor;

# database containing the alternative assembly
my $A_dba = $support->get_database('core', 'alt');
my $A_sa = $A_dba->get_SliceAdaptor;

#####
# create temporary table for storing non-aligned blocks
#
unless ($support->param('dry_run')) {
  $R_dbh->do(qq(
        CREATE TABLE IF NOT EXISTS tmp_align (
          tmp_align_id int(10) unsigned NOT NULL auto_increment,
          alt_seq_region_name varchar(20) NOT NULL,
          alt_start int(10) UNSIGNED NOT NULL,
          alt_end int(10) UNSIGNED NOT NULL,
          ref_seq_region_name varchar(20) NOT NULL,
          ref_start int(10) UNSIGNED NOT NULL,
          ref_end int(10) UNSIGNED NOT NULL,

          PRIMARY KEY (tmp_align_id)
          )
  ));

  # clear tmp_align table of entries from previous runs
  $R_dbh->do(qq(DELETE FROM tmp_align));
}

#####
# get reference and alternative toplevel seq_regions
#
my $R_chrlength = $support->get_chrlength($R_dba, $support->param('assembly'), 'toplevel');
my $A_chrlength = $support->get_chrlength($R_dba, $support->param('altassembly'), 'toplevel');

#####
# read list of clones to skip from file
#
$support->log("Reading list of clones to skip from file...\n");
my %skip = ();
if ($support->param('skipclones')) {
  my $infh = $support->filehandle('<', $support->param('skipclones'));
  while (<$infh>) {
    chomp;
    $skip{$_} = 1;
  }
}
$support->log("Done.\n\n");

#####
# loop over toplevel seq_regions
#
$support->log_stamped("Looping over toplevel seq_regions...\n");

my $match = {};
my $nomatch = {};
my %stats_total;
my @block_length;

my $fmt1 = "%-40s%10.0f\n";
my $fmt2 = "%-40s%9.2f%%\n";
my $fmt3 = "%-12s%-12s%-12s%-12s%-12s%-9s\n";
my $fmt4 = "%10.0f  %10.0f    %7.0f   %10.0f  %10.0f  %7.0f\n";
my $fmt5 = "%-40s%10s\n";
my $fmt6 = "%-10s%-12s%-10s%-12s\n";

my $sth1 = $R_dbh->prepare(qq(
    INSERT IGNORE INTO assembly (asm_seq_region_id, cmp_seq_region_id,
      asm_start, asm_end, cmp_start, cmp_end, ori)
    VALUES (?, ?, ?, ?, ?, ?, 1)
));
my $sth2 = $R_dbh->prepare(qq(
    INSERT INTO tmp_align values(NULL, ?, ?, ?, ?, ?, ?
)));

foreach my $R_chr ($support->sort_chromosomes($R_chrlength)) {

  $support->log_stamped("Toplevel seq_region $R_chr...\n", 1);

  my $A_chr = $R_chr;

  # fetch toplevel seq_region slices
  my $R_slice = $R_sa->fetch_by_region('toplevel', $R_chr, undef, undef, undef, $support->param('assembly'));
  my $A_slice = $A_sa->fetch_by_region('toplevel', $A_chr, undef, undef, undef, $support->param('altassembly'));

  unless ($R_slice and $A_slice) {
    $support->log("Seq_region not found in either ref or alt assembly. Skipping.\n", 2);
    next;
  }

  # we need to fetch the alternative slice from the reference db explicitely by
  # coord_system, since toplevel attribute is not set there
  my $cs_name = $A_slice->coord_system_name;
  my $A_slice_ref = $R_sa->fetch_by_region($cs_name, $A_chr, undef, undef, undef, $support->param('altassembly'));

  # project to contigs (not to clones because for WGS assembly regions there
  # are no clones)
  my @R_contigs = @{ $R_slice->project('contig') };
  my @A_contigs = @{ $A_slice->project('contig') };

  # loop over alternative clontigs
  my $last = 0;
  my $j = 0;
  my $match_flag = 0;
  my $last_A_seg;
  my %stats_chr;

  foreach my $A_seg (@A_contigs) {

    my $A_contig = $A_seg->to_Slice;

    # project contig to clone for clone name comparison
    my @A_clones = @{ $A_contig->project('clone') };

    my $A_clone;
    $A_clone = $A_clones[0]->to_Slice if (@A_clones);

    $support->log_verbose("Alternative contig ($j) ".$A_contig->seq_region_name.":".$A_contig->start."-".$A_contig->end.":".$A_contig->strand." $A_chr:".$A_seg->from_start."-".$A_seg->from_end."\n", 2);

    # walk reference contigs
    REFCLONES:
    for (my $i = $last; $i < scalar(@R_contigs); $i++) {

      my $R_contig = $R_contigs[$i]->to_Slice;

      # project contig to clone for clone name comparison
      my @R_clones = @{ $R_contig->project('clone') };

      my $R_clone;
      $R_clone = $R_clones[0]->to_Slice if (@R_clones);

      # same clone name.version and contig and clone strand found
      if ($A_clone and $R_clone and
          $A_clone->seq_region_name eq $R_clone->seq_region_name and
          $A_clone->strand == $R_clone->strand and
          $A_contig->strand == $R_contig->strand) {

        # same clone start/end -> identical assembly
        if ($A_clone->start == $R_clone->start and $A_clone->end == $R_clone->end) {
          # check if clone is tagged to be skipped
          # this can be used to resolve some odd assembly differences
          if ($skip{$A_clone->seq_region_name}) {

            $support->log_verbose("Skipping matching reference clone ($i)".
                $R_clone->seq_region_name.":".$R_clone->start."-".
                $R_clone->end.":".$R_clone->strand."$R_chr:".
                $R_contigs[$i]->from_start."-".$R_contigs[$i]->from_end.
                "\n", 2);

            &found_nomatch($R_chr, $A_chr, $match, $nomatch, $A_seg,
                $last_A_seg, $R_contigs[$i], $R_contigs[$i-1],
                $match_flag, $i, $j
            );

            $stats_chr{'skipped'}++;
            $match_flag = 0;

          } else {

            $support->log_verbose("Found matching reference clone ($i)".
                $R_clone->seq_region_name.":".$R_clone->start."-".
                $R_clone->end.":".$R_clone->strand."$R_chr:".
                $R_contigs[$i]->from_start."-".$R_contigs[$i]->from_end.
                "\n", 2);

            &found_match($R_chr, $A_chr, $match, $nomatch, $A_seg,
                $last_A_seg, $R_contigs[$i], $R_contigs[$i-1],
                $match_flag, $i, $j
            );

            $stats_chr{'identical'}++;
            $match_flag = 1;
          }

        # start/end mismatch
        } else {

          $support->log_verbose("Start/end mismatch for clone ($i) ".$R_contig->seq_region_name.":".$R_contig->start."-".$R_contig->end.":".$R_contig->strand." $R_chr:".$R_contigs[$i]->from_start."-".$R_contigs[$i]->from_end."\n", 2);

          &found_nomatch(
              $R_chr, $A_chr, $match, $nomatch, $A_seg, $last_A_seg,
              $R_contigs[$i], $R_contigs[$i-1], $match_flag, $i, $j
          );

          $stats_chr{'mismatch'}++;
          $match_flag = 0;
        }
        $i++;
        $last = $i;
        last REFCLONES;

      # different clones or no clones found
      } else {

        $support->log_verbose("Skipping clone ($i)".$R_contig->seq_region_name.":".$R_contig->start."-".$R_contig->end.":".$R_contig->strand." $R_chr:".$R_contigs[$i]->from_start."-".$R_contigs[$i]->from_end."\n", 2);

        &found_nomatch($R_chr, $A_chr, $match, $nomatch, $A_seg, $last_A_seg, $R_contigs[$i], $R_contigs[$i-1], $match_flag, $i, $j);

        $match_flag = 0;

      }
    }

    $last_A_seg = $A_seg;
    $j++;
  }

  # adjust the final clone count
  if ($match_flag) {

    # last clone was a match, adjust matching clone count
    if ($match->{$R_chr}) {

      my $c = scalar(@{ $match->{$R_chr} }) - 1;
      $match->{$R_chr}->[$c]->[2] = scalar(@A_contigs) - $match->{$R_chr}->[$c]->[2];
      $match->{$R_chr}->[$c]->[5] = scalar(@R_contigs) - $match->{$R_chr}->[$c]->[5];

    }

  } else {

    # last clone was a non-match, adjust non-matching clone count
    if ($nomatch->{$R_chr}) {

      my $c = scalar(@{ $nomatch->{$R_chr} }) - 1;
      $nomatch->{$R_chr}->[$c]->[2] = scalar(@A_contigs) - $nomatch->{$R_chr}->[$c]->[2];
      $nomatch->{$R_chr}->[$c]->[5] = scalar(@R_contigs) - $nomatch->{$R_chr}->[$c]->[5];

    }

  }

  # filter single assembly inserts from non-aligned blocks (i.e. cases where 
  # a block has clones only in one assembly, not in the other) - there is
  # nothing to be done with them
  @{ $nomatch->{$R_chr} } = grep { $_->[2] > 0 and $_->[5] > 0 }
    @{ $nomatch->{$R_chr} } if ($nomatch->{$R_chr});

  # store directly aligned blocks in assembly table
  unless ($support->param('dry_run')) {

    $support->log("Adding assembly entries for directly aligned blocks...\n", 1);
    my $c;

    for ($c = 0; $c < scalar(@{ $match->{$R_chr} || [] }); $c++) {
      $sth1->execute(
          $R_sa->get_seq_region_id($R_slice),
          $R_sa->get_seq_region_id($A_slice_ref),
          $match->{$R_chr}->[$c]->[3],
          $match->{$R_chr}->[$c]->[4],
          $match->{$R_chr}->[$c]->[0],
          $match->{$R_chr}->[$c]->[1]
      );
    }

    $support->log("Done inserting $c entries.\n", 1);
  }

  # store non-aligned blocks in tmp_align table
  unless ($support->param('dry_run')) {
    
    if ($nomatch->{$R_chr}) {
      
      $support->log("Storing non-aligned blocks in tmp_align table...\n", 1);
      my $c;
      
      for ($c = 0; $c < scalar(@{ $nomatch->{$R_chr} }); $c++) {
        $sth2->execute(
            $nomatch->{$R_chr}->[$c]->[6],
            $nomatch->{$R_chr}->[$c]->[0],
            $nomatch->{$R_chr}->[$c]->[1],
            $R_chr,
            $nomatch->{$R_chr}->[$c]->[3],
            $nomatch->{$R_chr}->[$c]->[4],
        );
      }

      $support->log("Done inserting $c entries.\n", 1);
    }
  }

  # stats for this toplevel seq_region
  $stats_chr{'A_only'} = scalar(@A_contigs) - $stats_chr{'identical'} - $stats_chr{'mismatch'};
  $stats_chr{'R_only'} = scalar(@R_contigs) - $stats_chr{'identical'} - $stats_chr{'mismatch'};
  for (my $c = 0; $c < scalar(@{ $match->{$R_chr} || [] }); $c++) {
    $stats_chr{'A_matchlength'} += $match->{$R_chr}->[$c]->[1] - $match->{$R_chr}->[$c]->[0];
    $stats_chr{'R_matchlength'} += $match->{$R_chr}->[$c]->[4] - $match->{$R_chr}->[$c]->[3];
  }
  $stats_chr{'A_coverage'} = 100 * $stats_chr{'A_matchlength'} /  $A_slice->length;
  $stats_chr{'R_coverage'} = 100 * $stats_chr{'R_matchlength'} /  $R_slice->length;
  map { $stats_total{$_} += $stats_chr{$_} } keys %stats_chr;

  $support->log("\nStats for toplevel seq_region $R_chr:\n\n", 1);
  $support->log(sprintf($fmt5, "Alternative toplevel seq_region name:", $A_chr), 2);
  $support->log(sprintf($fmt1, "Length (alternative):", $A_slice->length), 2);
  $support->log(sprintf($fmt1, "Length (reference):", $R_slice->length), 2);
  $support->log(sprintf($fmt1, "Identical clones:", $stats_chr{'identical'}), 2);
  $support->log(sprintf($fmt1, "Identical clones that were skipped:", $stats_chr{'skipped'}), 2);
  $support->log(sprintf($fmt1, "Clones with start/end mismatch:", $stats_chr{'mismatch'}), 2);
  $support->log(sprintf($fmt1, "Clones only in alternative assembly:", $stats_chr{'A_only'}), 2);
  $support->log(sprintf($fmt1, "Clones only in refernce assembly:", $stats_chr{'R_only'}), 2);
  $support->log(sprintf($fmt2, "Direct match coverage (alternative):", $stats_chr{'A_coverage'}), 2);
  $support->log(sprintf($fmt2, "Direct match coverage (reference):", $stats_chr{'R_coverage'}), 2);

  # Aligned blocks
  if ($match->{$R_chr}) {
    
    $support->log("\nDirectly aligned blocks:\n\n", 1);
    $support->log(sprintf($fmt3, qw(ALT_START ALT_END ALT_CLONES REF_START REF_END REF_CLONES)), 2);
    $support->log(('-'x71)."\n", 2);
    
    for (my $c = 0; $c < scalar(@{ $match->{$R_chr} }); $c++) {
      
      $support->log(sprintf($fmt4, @{ $match->{$R_chr}->[$c] }), 2);
      
      # sanity check: aligned region pairs must have same length
      my $e_len = $match->{$R_chr}->[$c]->[1] - $match->{$R_chr}->[$c]->[0] + 1;
      my $v_len = $match->{$R_chr}->[$c]->[4] - $match->{$R_chr}->[$c]->[3] + 1;
      
      $support->log_warning("Length mismatch: $e_len <> $v_len\n", 2) unless ($e_len == $v_len);

    }
  }

  # Non-aligned blocks
  if ($nomatch->{$R_chr}) {
    
    $support->log("\nNon-aligned blocks:\n\n", 1);
    $support->log(sprintf($fmt3, qw(ALT_START ALT_END ALT_CLONES REF_START REF_END REF_CLONES)), 2);
    $support->log(('-'x71)."\n", 2);
    
    for (my $c = 0; $c < scalar(@{ $nomatch->{$R_chr} }); $c++) {
      
      $support->log(sprintf($fmt4, @{ $nomatch->{$R_chr}->[$c] }), 2);

      # find longest non-aligned block
      my $A_length = $nomatch->{$R_chr}->[$c]->[1] - $nomatch->{$R_chr}->[$c]->[0] + 1;
      my $R_length = $nomatch->{$R_chr}->[$c]->[4] - $nomatch->{$R_chr}->[$c]->[3] + 1;
      push @block_length, [$A_chr, $A_length, $R_chr, $R_length];
    
    }
  }

  $support->log_stamped("\nDone with toplevel seq_region $R_chr.\n", 1);
}

# overall stats
$support->log("\nOverall stats:\n");
$support->log(sprintf($fmt1, "Identical clones:", $stats_total{'identical'}), 1);
$support->log(sprintf($fmt1, "Identical clones that were skipped:", $stats_total{'skipped'}), 1);
$support->log(sprintf($fmt1, "Clones with start/end mismatch:", $stats_total{'mismatch'}), 1);
$support->log(sprintf($fmt1, "Clones only in alternative assembly:", $stats_total{'A_only'}), 1);
$support->log(sprintf($fmt1, "Clones only in reference assembly:", $stats_total{'R_only'}), 1);

$support->log("\nNon-match block lengths:\n");
$support->log(sprintf($fmt6, qw(ALT_CHR ALT_LENGTH REF_CHR REF_LENGTH)), 1);
$support->log(('-'x42)."\n", 1);
foreach my $block (sort { $a->[1] <=> $b->[1] } @block_length) {
  $support->log(sprintf("%-10s%10.0f  %-10s%10.0f\n", @{ $block }), 1);
}

$support->log_stamped("\nDone.\n");

# finish logfile
$support->finish_log;


### end main


=head2 found_match

  Arg[1]      : String $R_chr - reference toplevel seq_region name 
  Arg[2]      : String $A_chr - alternative toplevel seq_region name 
  Arg[3]      : Hashref $match - datastructure to store aligned blocks
  Arg[4]      : Hashref $nomatch - datastructure to store non-aligned blocks
  Arg[5]      : Bio::EnsEMBL::ProjectionSegment $A_seg - current alternative
                segment
  Arg[6]      : Bio::EnsEMBL::ProjectionSegment $last_A_seg - last alternative
                segment
  Arg[7]      : Bio::EnsEMBL::ProjectionSegment $R_seg - current reference
                segment
  Arg[8]      : Bio::EnsEMBL::ProjectionSegment $last_R_seg - last reference
                segment
  Arg[9]      : Boolean $match_flag - flag indicating if last clone was a match
  Arg[10]     : Int $i - reference clone count
  Arg[11]     : Int $j - alternative clone count
  Description : This function is called when two clones match (i.e. have the
                same name.version in both assemblies). Depending on the state
                of the last clone (match or nomatch), it extends aligned blocks
                or finishes the non-aligned block and creates a new aligned
                block.
  Return type : none
  Exceptions  : none
  Caller      : internal

=cut

sub found_match {
  my ($R_chr, $A_chr, $match, $nomatch, $A_seg, $last_A_seg, $R_seg, $last_R_seg, $match_flag, $i, $j) = @_;

  # last clone was a match
  if ($match_flag) {

    # adjust align block end
    if ($match->{$R_chr}) {

      my $c = scalar(@{ $match->{$R_chr} }) - 1;

      # if the gaps between this clone and the last are different, start
      # a new block
      if (($A_seg->from_start - $match->{$R_chr}->[$c]->[1]) !=
          ($R_seg->from_start - $match->{$R_chr}->[$c]->[4])) {

        $support->log("Gap size mismatch at A:$A_chr:".$match->{$R_chr}->[$c]->[1].'-'.$A_seg->from_start.", R:$R_chr:".$match->{$R_chr}->[$c]->[4].'-'.$R_seg->from_start."\n", 2);

        # finish the last align block
        $match->{$R_chr}->[$c]->[1] = $last_A_seg->from_end;
        $match->{$R_chr}->[$c]->[2] = $j - $match->{$R_chr}->[$c]->[2];
        $match->{$R_chr}->[$c]->[4] = $last_R_seg->from_end;
        $match->{$R_chr}->[$c]->[5] = $i - $match->{$R_chr}->[$c]->[5];

        # start a new align block
        push @{ $match->{$R_chr} }, [
          $A_seg->from_start,
          $A_seg->from_end,
          $j,
          $R_seg->from_start,
          $R_seg->from_end,
          $i,
          $A_chr,
        ];

      # adjust align block end
      } else {
        $match->{$R_chr}->[$c]->[1] = $A_seg->from_end;
        $match->{$R_chr}->[$c]->[4] = $R_seg->from_end;
      }
    }

  # last clone was a non-match
  } else {
    
    # start a new align block
    push @{ $match->{$R_chr} }, [
      $A_seg->from_start,
      $A_seg->from_end,
      $j,
      $R_seg->from_start,
      $R_seg->from_end,
      $i,
      $A_chr,
    ];

    # finish the last non-align block
    if ($nomatch->{$R_chr} and $last_A_seg) {
      my $c = scalar(@{ $nomatch->{$R_chr} }) - 1;
      $nomatch->{$R_chr}->[$c]->[1] = $last_A_seg->from_end;
      $nomatch->{$R_chr}->[$c]->[2] = $j - $nomatch->{$R_chr}->[$c]->[2];
      $nomatch->{$R_chr}->[$c]->[4] = $last_R_seg->from_end;
      $nomatch->{$R_chr}->[$c]->[5] = $i - $nomatch->{$R_chr}->[$c]->[5];
    }
  }
}

=head2 found_nomatch

  Arg[1]      : String $R_chr - reference toplevel seq_region name 
  Arg[2]      : String $A_chr - alternative toplevel seq_region name 
  Arg[3]      : Hashref $match - datastructure to store aligned blocks
  Arg[4]      : Hashref $nomatch - datastructure to store non-aligned blocks
  Arg[5]      : Bio::EnsEMBL::ProjectionSegment $A_seg - current alternative
                segment
  Arg[6]      : Bio::EnsEMBL::ProjectionSegment $last_A_seg - last alternative
                segment
  Arg[7]      : Bio::EnsEMBL::ProjectionSegment $R_seg - current reference 
                segment
  Arg[8]      : Bio::EnsEMBL::ProjectionSegment $last_R_seg - last reference
                segment
  Arg[9]      : Boolean $match_flag - flag indicating if last clone was a match
  Arg[10]     : Int $i - reference clone count
  Arg[11]     : Int $j - alternative clone count
  Description : This function is called when two clones don't match (either
                different name.version or length mismatch in the two
                assemblies). Depending on the state of the last clone (nomatch
                or match), it extends non-aligned blocks or finishes the
                aligned block and creates a new non-aligned block.
  Return type : none
  Exceptions  : none
  Caller      : internal

=cut

sub found_nomatch {
  my ($R_chr, $A_chr, $match, $nomatch, $A_seg, $last_A_seg, $R_seg, $last_R_seg, $match_flag, $i, $j) = @_;

  # last clone was a match
  if ($match_flag) {
    
    # start a new non-align block
    push @{ $nomatch->{$R_chr} }, [
      $A_seg->from_start,
      $A_seg->from_end,
      $j,
      $R_seg->from_start,
      $R_seg->from_end,
      $i,
      $A_chr,
    ];

    # finish the last align block
    if ($match->{$R_chr}) {
      my $c = scalar(@{ $match->{$R_chr} }) - 1;
      $match->{$R_chr}->[$c]->[1] = $last_A_seg->from_end;
      $match->{$R_chr}->[$c]->[2] = $j - $match->{$R_chr}->[$c]->[2];
      $match->{$R_chr}->[$c]->[4] = $last_R_seg->from_end;
      $match->{$R_chr}->[$c]->[5] = $i - $match->{$R_chr}->[$c]->[5];
    }

  # last clone was a non-match
  } else {
    
    # adjust non-align block end
    if ($nomatch->{$R_chr}) {
      my $c = scalar(@{ $nomatch->{$R_chr} }) - 1;
      $nomatch->{$R_chr}->[$c]->[1] = $A_seg->from_end;
      $nomatch->{$R_chr}->[$c]->[4] = $R_seg->from_end;

    # we're at the beginning of a seq_region, so start a new non-align block
    } else {
      push @{ $nomatch->{$R_chr} }, [
        $A_seg->from_start,
        $A_seg->from_end,
        $j,
        $R_seg->from_start,
        $R_seg->from_end,
        $i,
        $A_chr,
      ];
    }
  }
}

