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

compare_assemblies.pl - compare two assemblies

=head1 SYNOPSIS

compare_assemblies.pl [arguments]

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

This script compares two assemblies. At the moment all it does is to list
clones that are on different toplevel seq_regions in the two assemblies.

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
    'chromosomes|chr=s@',
);
$support->allowed_params(
    $support->get_common_params,
    'assembly',
    'altdbname',
    'altassembly',
    'althost',
    'altport',
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
    'altdbname',
    'altassembly'
);

# first set connection parameters for alternative db if not different from
# reference db
map { $support->param("alt$_", $support->param($_)) unless ($support->param("alt$_")) } qw(host port user pass);

# reference database
my $R_dba = $support->get_database('ensembl');
my $R_sa = $R_dba->get_SliceAdaptor;

# database containing the alternative assembly
my $A_dba = $support->get_database('core', 'alt');
my $A_sa = $A_dba->get_SliceAdaptor;

my $fmt1 = "%-20s%-12s%-12s\n";
my $fmt2 = "%-35s%10.0f\n";
my @diff_total;
my %stats_total;

$support->log_stamped("Looping over toplevel seq_regions...\n");

foreach my $R_chr ($support->sort_chromosomes($support->get_chrlength)) {
    $support->log_stamped("\nChromosome $R_chr...\n", 1);

    # fetch toplevel seq_region slice and project to clones
    my $R_slice = $R_sa->fetch_by_region('toplevel', $R_chr);
    my @R_clones = @{ $R_slice->project('clone') };

    # loop over reference clones
    my @diff;
    my %stats;
    foreach my $R_seg (@R_clones) {
        $stats{'num_clones'}++;
        my $R_clone = $R_seg->to_Slice;
        my $R_clone_name = $R_clone->seq_region_name;

        # fetch clone from alternative db and project to toplevel seq_region
        my $A_clone = $A_sa->fetch_by_region('clone', $R_clone_name);
        if ($A_clone) {
            my ($A_seg) = @{ $A_clone->project('toplevel') };
            if ($A_seg) {
                my $A_slice = $A_seg->to_Slice;

                # compare toplevel seq_region names
                my $A_chr = $A_slice->seq_region_name;
                unless ($R_chr eq $A_chr) {
                    push @diff, [$R_clone_name, $R_chr, $A_chr];
                }
            } else {
                $stats{'does_not_project'}++;
                $support->log_verbose("Clone $R_clone_name doesn't project to toplevel seq_region.\n", 2);
            }
        } else {
            $stats{'not_in_alt'}++;
            $support->log_verbose("Clone $R_clone_name not in alternative db.\n", 2);
        }
    }
    push @diff_total, @diff;
    map { $stats_total{$_} += $stats{$_} }
        qw(num_clones does_not_project not_in_alt);

    # print results for this toplevel seq_region
    $support->log("\nStats for toplevel seq_region $R_chr:\n\n", 2);
    $support->log(sprintf($fmt2, "Total number of clones:", $stats{'num_clones'}), 3);
    $support->log(sprintf($fmt2, "Clones not in alternative db:", $stats{'not_in_alt'}), 3);
    $support->log(sprintf($fmt2, "Clones not in alternative assembly:", $stats{'does_not_project'}), 3);
    if (@diff) {
        $support->log("\nClones on different toplevel seq_regions in the 2 assemblies:\n\n", 3);
        $support->log(sprintf($fmt1, qw(CLONE R_CHR A_CHR)), 4);
        $support->log(('-'x44)."\n", 4);
        foreach my $d (@diff) {
            $support->log(sprintf($fmt1, @{ $d }), 4);
        }
    } else {
        $support->log("\nAll matching clones on same toplevel seq_region in the 2 assemblies.\n", 3);
    }
}

$support->log_stamped("Done.\n\n");

# print overall results
$support->log("\nOverall stats:\n");
$support->log(sprintf($fmt2, "Total number of clones:", $stats_total{'num_clones'}), 1);
$support->log(sprintf($fmt2, "Clones not in alternative db:", $stats_total{'not_in_alt'}), 1);
$support->log(sprintf($fmt2, "Clones not in alternative assembly:", $stats_total{'does_not_project'}), 1);
if (@diff_total) {
    $support->log("\nClones on different toplevel seq_regions in the 2 assemblies:\n\n", 1);
    $support->log(sprintf($fmt1, qw(CLONE R_CHR A_CHR)), 2);
    $support->log(('-'x44)."\n", 2);
    foreach my $d (@diff_total) {
        $support->log(sprintf($fmt1, @{ $d }), 2);
    }
} else {
    $support->log("\nAll clones on same toplevel seq_region in the 2 assemblies.\n", 2);
}
$support->log("\n");

# finish logfile
$support->finish_log;

