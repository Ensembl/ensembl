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

align_nonident_regions_wrapper.pl - lsf wrapper for align_nonident_regions.pl

=head1 SYNOPSIS

align_nonident_regions_wrapper.pl [arguments]

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
    --bindir=DIR                        look for program binaries in DIR
    --tmpfir=DIR                        use DIR for temporary files (useful for
                                        re-runs after failure)

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

This script is a wrapper around align_nonident_regions.pl to run one toplevel seq_region
at a time via lsf. See documentation in align_nonident_regions.pl for details.

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
    unshift(@INC, ".");
    unshift(@INC, "$SERVERROOT/ensembl/modules");
    unshift(@INC, "$SERVERROOT/bioperl-live");
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use AssemblyMapper::BlastzAligner;

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
    'bindir=s',
    'tmpdir=s',
    'chromosomes|chr=s@',
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
    'bindir',
    'tmpdir',
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

#####
# connect to database and get adaptors
#
my $R_dba = $support->get_database('ensembl');

# loop over toplevel seq_regions
$support->log_stamped("Looping over toplevel seq_regions...\n");

foreach my $chr ($support->sort_chromosomes($support->get_chrlength)) {
    
    $support->log_stamped("Toplevel seq_region $chr...\n", 1);

    # run align_nonident_regions.pl via lsf
    my $options = $support->create_commandline_options({
        'allowed_params' => 1,
        'replace' => {
            logfile     => "align_nonident_regions.chr.$chr.log",
            chromosomes => $chr,
            interactive => 0,
        },
    });
    $support->log("Running align_nonident_regions.pl via lsf...\n", 2);

    my $logpath = $support->param('logpath').'/lsf';
    unless (-d $logpath) {
        system("mkdir $logpath") == 0 or
            $support->log_error("Can't create lsf log dir $logpath: $!\n");
    }

    my $cmd = qq(bsub -o $logpath/align_nonident_regions.$chr.%J.out -e $logpath/align_nonident_regions.$chr.%J.err ./align_nonident_regions.pl $options
    );

    system("$cmd") == 0
        or $support->log_error("Error running align_nonident_regions.pl: $!i\n");

    $support->log_stamped("Done.\n", 2);

    $support->log_stamped("Done with toplevel seq_region $chr.\n", 1);
}

$support->log_stamped("Done.\n");

# finish logfile
$support->finish_log;


