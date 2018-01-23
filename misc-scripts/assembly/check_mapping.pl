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

check_mapping.pl - script to check whole genome alignment between two
assemblies.

=head1 SYNOPSIS

check_mapping.pl [arguments]

Required arguments:

  --dbname, db_name=NAME              database name NAME
  --host, --dbhost, --db_host=HOST    database host HOST
  --port, --dbport, --db_port=PORT    database port PORT
  --user, --dbuser, --db_user=USER    database username USER
  --pass, --dbpass, --db_pass=PASS    database passwort PASS
  --assembly=ASSEMBLY                 assembly version ASSEMBLY

  --altdbname=NAME                    alternative database NAME
  --altassembly=ASSEMBLY              alternative assembly version ASSEMBLY

  --nolog                             outside of a webserver environment, 
                                      logging needs to be handled manually.
                                      Use in conjunction with shell redirects.

Optional arguments:

  --althost=HOST                      alternative databases host HOST
  --altport=PORT                      alternative database port PORT
  --altuser=USER                      alternative database username USER
  --altpass=PASS                      alternative database password PASS

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

This script checks if the whole genome alignment between two assemblies is
correct. It does so by comparing the sequence in the reference database with
the sequence of the projected fragments in the alternative database.

=head1 RELATED FILES

The whole process of creating a whole genome alignment between two assemblies
is done by a series of scripts. Please see

  ensembl/misc-scripts/assembly/README

for a high-level description of this process, and POD in the individual scripts
for the details.


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

use Algorithm::Diff qw(diff);

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'assembly=s',
    'altdbname=s',
    'altassembly=s',
    'althost=s',
    'altport=n',
    'altpass=s',
    'chromosomes|chr=s@',
);
$support->allowed_params(
    $support->get_common_params,
    'assembly',
    'altdbname',
    'altassembly',
    'althost',
    'altport',
    'altpass',
    'chromosomes',
);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

$support->comma_to_list('chromosomes');

# ask user to confirm parameters to proceed
#$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params(
    'assembly',
    'altdbname',
    'altassembly'
);

# first set connection parameters for alternative db if not different from
# reference db
map { $support->param("alt$_", $support->param($_)) unless ($support->param("alt$_")) } qw(host port user);

# reference database
my $R_dba = $support->get_database('ensembl');
my $R_sa = $R_dba->get_SliceAdaptor;

# database containing the alternative assembly
my $A_dba = $support->get_database('core', 'alt');
my $A_sa = $A_dba->get_SliceAdaptor;

$support->log("Looping over toplevel seq_regions...\n\n");

my @global_diff_bins;

foreach my $chr ($support->sort_chromosomes) {
  $support->log_stamped("Toplevel seq_region $chr...\n", 1);

  my $R_slice = $R_sa->fetch_by_region('toplevel', $chr);
  my $A_slice = $A_sa->fetch_by_region('toplevel', $chr);

  unless ($A_slice) {
    $support->log("Not found in alternative db. Skipping.\n", 2);
    next;
  }
  
  my $cs_name = $A_slice->coord_system_name;

  # compare reference and alternative sequence
  my @segments = @{ $R_slice->project($cs_name, $support->param('altassembly')) };
  
  my $i;
  my $k;

  foreach my $seg (@segments) {
    # reference sequence
    my $R_sub_slice = $R_slice->sub_Slice($seg->from_start, $seg->from_end);
    my $R_seq = $R_sub_slice->seq;
    
    # alternative sequence
    my $A_proj_slice = $seg->to_Slice;
    
    # ignore PAR region (i.e. we project to the symlinked seq_region)
    next if ($A_proj_slice->seq_region_name ne $chr);
    
    my $A_sub_slice = $A_slice->sub_Slice($A_proj_slice->start, $A_proj_slice->end, $A_proj_slice->strand);
    my $A_seq = $A_sub_slice->seq;

    # compare
    if ($R_seq eq $A_seq) {
      # sequences are identical -> ok
      $support->log_verbose("Sequence match at ".$R_sub_slice->name."\n", 2);

    } else {
      # not identical -> something is wrong
      $support->log("Sequence mismatch at ".$R_sub_slice->name."\n", 2);

      my $R_sub_seq;
      my $A_sub_seq;
      
      if ($R_sub_slice->length > 20) {
        $R_sub_seq = substr($R_seq, 0, 10)."...".substr($R_seq, -10, 10);
        $A_sub_seq = substr($A_seq, 0, 10)."...".substr($A_seq, -10, 10);
      } else {
        $R_sub_seq = substr($R_seq, 0, 20);
        $A_sub_seq = substr($A_seq, 0, 20);
      }
      
      $support->log("Ref: $R_sub_seq\n", 3);
      $support->log("Alt: $A_sub_seq\n\n", 3);

      $i++;
      
      
      if (length($R_seq) == length($A_seq)){
          #my $diffs = ($R_seq ^ $A_seq) =~ tr/\0//c;  # A concatenation of differences
          # this approach is x10 faster than relying on Algorithm::Diff, as long as there
          # are no InDels, and the lengths are comparable.
          my $mask = ($R_seq ^ $A_seq);
          my @diffs = split (//,$mask);
          my ($in_change,$change_start,$change_end);
          for (my $x=0; $x<scalar(@diffs); $x++) {
              if ($in_change) {
                  if ($diffs[$x] eq "\0") {
                      $in_change = 0;
                      $change_end = $x;
                      
                      my $length = $change_end - $change_start + 1;
                      $global_diff_bins[$length]++;
                  }
                  else {
                      next;
                  }                  
              } elsif ($diffs[$x] ne "\0") {
                  $in_change = 1;
                  $change_start = $x;
              }
              
          }
      } else {
          my @Ref = split(//,$R_seq);
          my @Alt = split(//,$A_seq);
          my @diffs = diff( \@Ref, \@Alt );
          foreach (@diffs) {
              my $length = 0;
              foreach my $desc (@{$_}) {
                  if ($desc->[0] eq '+') {$length++;}
                  if ($desc->[0] eq '-') {$length--;}
              };
              $global_diff_bins[$length]++;
          }
      }
      
    }

    $k++;
  }

  if ($i) {
    $support->log("Total: $i (of $k) alignments contain sequence mismatches.\n", 2);
  } else {
    $support->log("All $k alignments ok.\n", 2);
  }

  $support->log_stamped("Done.\n\n", 1);
}

$support->log("Summary of changes across all chromosomes\n\n",1);
$support->log("|Bin  |Frequency\n",2);
for (my $bin = 0; $bin < scalar(@global_diff_bins); $bin++) {
    if (defined ($global_diff_bins[$bin])) {
        $support->log("|$bin  |$global_diff_bins[$bin]\n");
    }
}

# finish logfile
$support->finish_log;

