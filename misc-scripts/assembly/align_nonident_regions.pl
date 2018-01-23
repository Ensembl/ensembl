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

align_nonident_regions.pl - create whole genome alignment between two closely
related assemblies for non-identical regions

=head1 SYNOPSIS

align_nonident_regions.pl [arguments]

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

    --chromosomes, --chr=LIST           only process LIST chromosomes
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

This script is part of a series of scripts to create a mapping between two
assemblies. It assembles the chromosome coordinate systems of two different
assemblies of a genome by creating a whole genome alignment between the two.

The process assumes that the two assemblies are reasonably similar, i.e. there
are no major rearrangements or clones moved from one chromosome to another.

See "Related files" below for an overview of the whole process.

This particular script creates a whole genome alignment between two closely
related assemblies for non-identical regions. These regions are identified by
another script (align_by_clone_identity.pl) and stored in a temporary database
table (tmp_align).

Alignments are calculated by this algorithm:

    1. fetch region from tmp_align
    2. write soft-masked sequences to temporary files
    3. align using blastz
    4. filter best hits (for query sequences, i.e. alternative regions) using
       axtBest
    5. parse blastz output to create blocks of exact matches only
    7. write alignments to assembly table

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

# create BlastzAligner object
my $aligner = AssemblyMapper::BlastzAligner->new(-SUPPORT => $support);

# create tmpdir to store input and output
$aligner->create_tempdir($support->param('tmpdir'));

# loop over non-aligned regions in tmp_align table
$support->log_stamped("Looping over non-aligned blocks...\n");

$sql = qq(SELECT * FROM tmp_align);
if ($support->param('chromosomes')) {
  my $chr_string = join("', '", $support->param('chromosomes'));
  $sql .= " WHERE ref_seq_region_name IN ('$chr_string')";
}
$sth = $R_dbh->prepare($sql);
$sth->execute;

while (my $row = $sth->fetchrow_hashref) {

  my $id = $row->{'tmp_align_id'};
  $aligner->id($id);
  $aligner->seq_region_name($row->{'ref_seq_region_name'});

  $support->log_stamped("Block with tmp_align_id = $id\n", 1);
  
  my $A_slice = $A_sa->fetch_by_region(
      'toplevel',
      $row->{'alt_seq_region_name'},
      $row->{'alt_start'},
      $row->{'alt_end'},
      1,
      $support->param('altassembly'),
  );
  
  my $R_slice = $R_sa->fetch_by_region(
      'toplevel',
      $row->{'ref_seq_region_name'},
      $row->{'ref_start'},
      $row->{'ref_end'},
      1,
      $support->param('assembly'),
  );

  # write sequences to file, and convert sequence files from fasta to nib
  # format (needed for lavToAxt)
  my $A_basename = "alt_seq.$id";
  my $R_basename = "ref_seq.$id";
  
  $support->log("Writing sequences to fasta and nib files...\n", 2);
  
  $aligner->write_sequence(
      $A_slice,
      $support->param('altassembly'),
      $A_basename
  );
  
  $aligner->write_sequence(
      $R_slice,
      $support->param('assembly'),
      $R_basename
  );
  
  $support->log("Done.\n", 2);

  # align using blastz
  $support->log("Running blastz...\n", 2);
  $aligner->run_blastz($A_basename, $R_basename);
  $support->log("Done.\n", 2);

  # convert blastz output from lav to axt format
  $support->log("Converting blastz output from lav to axt format...\n", 2);
  $aligner->lav_to_axt;
  $support->log("Done.\n", 2);

  # find best alignment with axtBest
  $support->log("Finding best alignment with axtBest...\n", 2);
  $aligner->find_best_alignment;
  $support->log("Done.\n", 2);

  # parse blastz output, and convert relative alignment coordinates to
  # chromosomal coords
  $support->log("Parsing blastz output...\n", 2);
  
  $aligner->parse_blastz_output;
  
  $aligner->adjust_coords(
      $row->{'alt_start'},
      $row->{'alt_end'},
      { $id => [ $row->{'ref_start'}, $row->{'ref_end'} ] }
  );
  
  $support->log("Done.\n", 2);

  # cleanup temp files
  $support->log("Cleaning up temp files...\n", 2);
  $aligner->cleanup_tmpfiles(
    "$A_basename.fa",
    "$A_basename.nib",
    "$R_basename.fa",
    "$R_basename.nib",
  );
  $support->log("Done.\n", 2);

  # log alignment stats
  $aligner->log_block_stats(2);

  $support->log_stamped("Done with block $id.\n", 1);

}

$support->log_stamped("Done.\n");

# write alignments to assembly table
unless ($support->param('dry_run')) {
  $aligner->write_assembly($R_dba);
}

# cleanup
$support->log_stamped("\nRemoving tmpdir...\n");
$aligner->remove_tempdir;
$support->log_stamped("Done.\n\n");

# overall stats
$aligner->log_overall_stats;

# remind to drop tmp_align
$support->log("\nDon't forget to drop the tmp_align table when all is done!\n\n");

# finish logfile
$support->finish_log;


