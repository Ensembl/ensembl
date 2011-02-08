#!/usr/local/ensembl/bin/perl

=head1 NAME

glovar_snp_density.pl - Script to calculate glovar SNP density and stats (and
optioanlly prepare AV index dumps) for Vega.

=head1 SYNOPSIS

glovar_snp_density.pl [options]

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

Specific options:

    --chromosomes, --chr=LIST           only process LIST chromosomes
    --avdump=0|1                        create AV dump
    --glovardbname=NAME                 use Glovar database NAME
    --glovarhost=HOST                   use Glovar database host HOST
    --glovarport=PORT                   use Glovar database port PORT
    --glovaruser=USER                   use Glovar database username USER
    --glovarpass=PASS                   use Glovar database passwort PASS
    --oracle_home=PATH                  set $ORACLE_HOME env variable to PATH
    --ld_library_path=PATH              set $LD_LIBRARY_PATH env variable to
                                        PATH
    --glovar_snp_consequence_exp=NUM    use NUM glovar SNP consequence
                                        experiment

=head1 DESCRIPTION

This script calculates Glovar SNP densities and total numbers per chromosome
for use in mapview. Can be run for individual chromosomes if desired (default:
all chromosomes). Since it uses a lot of memory, there is a wrapper script
which runs this script one chromosome at a time (glovar_snp_wrapper.pl). It
optionally also dumps SNPs into a file for generating the AV search index.

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

Post questions to the EnsEMBL development list dev@ensembl.org

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin);
use vars qw($SERVERROOT);

BEGIN {
    $SERVERROOT = "$Bin/../../..";
    unshift(@INC, "$SERVERROOT/ensembl/modules");
    unshift(@INC, "$SERVERROOT/ensembl-external/modules");
    unshift(@INC, "$SERVERROOT/ensembl-variation/modules");
    unshift(@INC, "$SERVERROOT/bioperl-live");
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;
use POSIX;

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'chromosomes|chr=s@',
    'avdump=s',
    'glovarhost=s',
    'glovarport=s',
    'glovaruser=s',
    'glovarpass=s',
    'glovardbname=s',
    'oracle_home=s',
    'ld_library_path=s',
    'glovar_snp_consequence_exp=n',
);
$support->allowed_params($support->get_common_params,
    'chromosomes',
    'avdump',
    'glovarhost',
    'glovarport',
    'glovaruser',
    'glovarpass',
    'glovardbname',
    'oracle_home',
    'ld_library_path',
    'glovar_snp_consequence_exp',
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

# connect to database and get adaptors
my $dba = $support->get_database('ensembl');
my $dba_glovar = $support->get_glovar_database;
my $dfa = $dba->get_DensityFeatureAdaptor;
my $dta = $dba->get_DensityTypeAdaptor;
my $aa  = $dba->get_AnalysisAdaptor;
my $attrib_adaptor = $dba->get_AttributeAdaptor;

# split chromosomes by size and determine block size
my $chr_slices = $support->split_chromosomes_by_size(5000000);

# create Analysis object
my $analysis = Bio::EnsEMBL::Analysis->new(
        -program     => "glovar_snp_density.pl",
        -database    => "vega",
        -gff_source  => "glovar_snp_density.pl",
        -gff_feature => "density",
        -logic_name  => "snpdensity",
);
$aa->store( $analysis ) unless ($support->param('dry_run'));

# settings for AV index dump
use constant SNP_LINE => join("\t",
    'Glovar SNP', '%s', '/%s/snpview?snp=%s&source=glovar', '%s',
    "Single nucleotide polymorphism (SNP) %s [Alleles: %s]. Alternative IDs: %s.\n"
);
my $species = $support->species;
my $fh;
if ($support->param('avdump')) {
    $support->log("Preparing directories for AV index dump...\n");
    my $dumpdir = $support->serverroot."/sanger-plugins/vega/utils/indexing/input";
    unless (-e $dumpdir) {
        mkdir $dumpdir, 0777 or
            $support->log_error("Could not creat directory $dumpdir: $!\n");
    }
    unless (-e "$dumpdir/$species") {
        mkdir "$dumpdir/$species", 0777 or
            $support->log_error("Could not creat directory $dumpdir/$species: $!\n");
    }
    $fh = $support->filehandle('>>', "$dumpdir/$species/SNP.txt");
    $support->log("Done.\n");
}

# loop over block sizes
my %av_done;
foreach my $block_size (keys %{ $chr_slices }) {
    $support->log("Available chromosomes using block size of $block_size:\n    ");
    $support->log(join("\n    ", map { $_->seq_region_name } @{ $chr_slices->{$block_size} })."\n");

    # create DensityType objects
    my $dt = Bio::EnsEMBL::DensityType->new(
            -analysis   => $analysis,
            -block_size => $block_size,
            -value_type => 'sum',
    );
    $dta->store($dt) unless ($support->param('dry_run'));

    # looping over chromosomes
    $support->log_stamped("Looping over chromosomes...\n");
    my ($current_start, $current_end);
    foreach my $slice (@{ $chr_slices->{$block_size} }) {
        $current_start = 1;
        my $chr = $slice->seq_region_name;
        my ($total, $i);
        my $bins = POSIX::ceil($slice->end/$block_size);
        
        $support->log_stamped("Chromosome $chr with block size $block_size...\n", 1);

        # loop over blocks
        while ($current_start <= $slice->end) {
            $i++;
            $current_end = $current_start + $block_size - 1;
            if ($current_end > $slice->end) {
                $current_end = $slice->end;
            }
            my $sub_slice = $slice->sub_Slice($current_start, $current_end);
            my $count = 0;

            my $varfeats;
            eval { $varfeats = $sub_slice->get_all_ExternalFeatures('GlovarSNP'); };
            if ($@) {
                $support->log_warning($@);
                $current_start = $current_end + 1;
                next;
            }

            # only count varfeats that don't overlap slice start
            # also, avoid duplicate counting
            my %varfeats = map { $_->start > 0 ? ($_->variation_name => 1) : () } @{$varfeats};
            $count = scalar(keys %varfeats);

            # AV index dump
            if ($support->param('avdump') && (! $av_done{$chr})) {
                foreach my $varfeat (@{$varfeats}) {
                    next if ($varfeat->start < 1);
                    my $snpid = $varfeat->variation_name;

                    # dblinks
                    my @sources = @{ $varfeat->variation->get_all_synonym_sources };
                    my (@IDs, @desc);
                    foreach my $source (@sources) {
                        my @extIDs = @{ $varfeat->variation->get_all_synonyms($source) };
                        push @IDs, @extIDs;
                        push @desc, "$source: @extIDs";
                    }
                    
                    print $fh sprintf SNP_LINE,
                        $snpid,
                        $species,
                        $snpid,
                        join(" ", @IDs),
                        $snpid,
                        $varfeat->allele_string,
                        join(", ", @desc)
                    ; 
                }
            }

            # density
            my $df = Bio::EnsEMBL::DensityFeature->new(
                -seq_region    => $slice,
                -start         => $current_start,
                -end           => $current_end,
                -density_type  => $dt,
                -density_value => $count,
            );
            $current_start = $current_end + 1;
            $dfa->store($df) unless ($support->param('dry_run'));
            $total += $count;

            # logging
            $support->log_verbose("Chr: $chr | Bin: $i/$bins | Count: $count | ".$support->date_and_mem."\n", 2);
        }

        # set flag to do AV dump only once for each chromosome
        $av_done{$chr} = 1;
                
        # stats
        $support->log("Total for chr $chr: $total\n", 1);
        my $stat = Bio::EnsEMBL::Attribute->new(
            -NAME => 'SNPs',
            -CODE => 'SNPCount',
            -VALUE => $total,
            -DESCRIPTION => 'Total Number of SNPs',
        );
        $attrib_adaptor->store_on_Slice($slice, [$stat]) unless ($support->param('dry_run'));

        $support->log_stamped("Done.\n", 1);
    }
    $support->log_stamped("Done.\n");
}

# finish logfile
$support->finish_log;

