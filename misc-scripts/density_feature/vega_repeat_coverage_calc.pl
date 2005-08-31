#!/usr/local/bin/perl

=head1 NAME

vega_repeat_coverage_calc.pl - calculate the repeat coverage

=head1 SYNOPSIS

vega_repeat_coverage_calc.pl [options]

General options:
    --conffile, --conf=FILE             read parameters from FILE
                                        (default: conf/Conversion.ini)

    --dbname, db_name=NAME              use database NAME
    --host, --dbhost, --db_host=HOST    use database host HOST
    --port, --dbport, --db_port=PORT    use database port PORT
    --user, --dbuser, --db_user=USER    use database username USER
    --pass, --dbpass, --db_pass=PASS    use database passwort PASS
    --logfile, --log=FILE               log to FILE (default: *STDOUT)
    --logpath=PATH                      write logfile to PATH (default: .)
    --logappend, --log_append           append to logfile (default: truncate)
    -v, --verbose                       verbose logging (default: false)
    -i, --interactive=0|1               run script interactively (default: true)
    -n, --dry_run, --dry=0|1            don't write results to database
    -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

This script calculates the repeat coverage for given database.

The block size is determined so that you have 150 bins for the smallest
chromosome over 5 Mb in length. For chromosomes smaller than 5 Mb, an
additional smaller block size is used to yield 150 bins for the overall
smallest chromosome. This will result in reasonable resolution for small
chromosomes and high performance for big ones.

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Patrick Meidl <pm2@sanger.ac.uk>

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin);
use vars qw($SERVERROOT);

BEGIN {
    $SERVERROOT = "$Bin/../../..";
    unshift(@INC, "$SERVERROOT/ensembl-otter/modules");
    unshift(@INC, "$SERVERROOT/ensembl/modules");
    unshift(@INC, "$SERVERROOT/bioperl-live");
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;
use Bio::EnsEMBL::Mapper::RangeRegistry;
use POSIX;

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->allowed_params($support->get_common_params);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

# connect to database and get adaptors
my $dba = $support->get_database('ensembl');
my $dfa = $dba->get_DensityFeatureAdaptor;
my $dta = $dba->get_DensityTypeAdaptor;
my $aa = $dba->get_AnalysisAdaptor;

# Create Analysis object
my $analysis = new Bio::EnsEMBL::Analysis (
        -program     => "repeat_coverage_calc.pl",
        -database    => "ensembl",
        -gff_source  => "repeat_coverage_calc.pl",
        -gff_feature => "density",
        -logic_name  => "PercentageRepeat",
);
$aa->store($analysis) unless ($support->param('dry_run'));

# split chromosomes by size and determine block size
my $chr_slices = $support->split_chromosomes_by_size(5000000);

# loop over block sizes
foreach my $block_size (keys %{ $chr_slices }) {
    $support->log("Available chromosomes using block size of $block_size:\n    ");
    $support->log(join("\n    ", map { $_->seq_region_name } @{ $chr_slices->{$block_size} })."\n");

    # create DensityType objects
    my $density_type = Bio::EnsEMBL::DensityType->new(
         -analysis   => $analysis,
         -block_size => $block_size,
         -value_type => 'ratio',
    );
    $dta->store($density_type) unless ($support->param('dry_run'));

    # loop over chromosomes
    $support->log_stamped("Looping over chromosomes...\n");
    my ($current_start, $current_end);
    foreach my $slice (@{ $chr_slices->{$block_size} }) {
        $current_start = 1;
        my $chr = $slice->seq_region_name;
        my $i;
        my $bins = POSIX::ceil($slice->end/$block_size);

        $support->log_stamped("Chromosome $chr with block size $block_size...\n", 1);

        # loop over blocks
        while($current_start <= $slice->end) {
            $i++;
            $current_end = $current_start + $block_size - 1;
            if ($current_end > $slice->end) {
                $current_end = $slice->end;
            }
            my $this_block_size = $current_end - $current_start + 1;
            my $sub_slice = $slice->sub_Slice($current_start, $current_end);
            my $rr = Bio::EnsEMBL::Mapper::RangeRegistry->new or die;
            foreach my $repeat (@{ $sub_slice->get_all_RepeatFeatures }) {
                $rr->check_and_register("1", $repeat->start, $repeat->end);
            }
            my $count = 0;
            my $non_repeats = $rr->check_and_register("1", 1, $this_block_size);
            if (defined $non_repeats) {
                foreach my $non_repeat (@{ $non_repeats }) {
                    $count += ($non_repeat->[1] - $non_repeat->[0]) + 1;
                }
            }
            my $percentage_repeat = (($this_block_size-$count)/$this_block_size)*100;
            $support->log_verbose("Chr: $chr | Bin: $i/$bins | ", 2);
            $support->log_verbose("\%repeat ".sprintf("%.2f", $percentage_repeat)."\n");
            my $df = Bio::EnsEMBL::DensityFeature->new(
                -seq_region    => $slice,
                -start         => $current_start,
                -end           => $current_end,
                -density_type  => $density_type,
                -density_value => $percentage_repeat,
            );
            $dfa->store($df) unless ($support->param('dry_run'));
            $current_start = $current_end + 1;
        }
        $support->log_stamped("Done.\n", 1);
    }
    $support->log_stamped("Done.\n");
}

# finish logfile
$support->finish_log;

