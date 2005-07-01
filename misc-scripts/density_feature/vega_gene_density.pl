#!/usr/local/bin/perl

=head1 NAME

vega_gene_density.pl - script to calculate gene densities and stats in Vega

=head1 SYNOPSIS

vega_gene_density.pl [options]

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

=head1 DESCRIPTION

This script calculates Vega gene densities and total numbers per chromosome
for use in mapview.

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

based on code by
    Graham McVicer <mcvicker@ebi.ac.uk>
    Steve Trevanion <st3@sanger.ac.uk>

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
use POSIX;

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options('chromosomes|chr=s@');
$support->allowed_params($support->get_common_params, 'chromosomes');

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
my $slice_adaptor = $dba->get_SliceAdaptor;
my $dfa = $dba->get_DensityFeatureAdaptor;
my $dta = $dba->get_DensityTypeAdaptor;
my $aa = $dba->get_AnalysisAdaptor;
my $attrib_adaptor = $dba->get_AttributeAdaptor;

# split chromosomes by size and determine block size
my $chr_slices = $support->split_chromosomes_by_size(5000000);

# create Analysis and DensityType objects
my %gene_types = (
    "Known" => "knownGeneDensity",
    "Novel_CDS" => "novelCDSDensity",
    "Novel_Transcript" => "novelTransDensity",
    "Putative" => "putativeTransDensity",
    "Predicted_Gene" => "predictedTransDensity",
    "Ig_Pseudogene_Segment" => "IgPseudoSegDensity",
    "Ig_Segment" => "IgSegDensity",
    "Processed_pseudogene" => "pseudoGeneDensity",
    "Pseudogene" => "pseudoGeneDensity",
    "Unprocessed_pseudogene" => "pseudoGeneDensity",
    "Known_in_progress" => "knownGeneDensity",
    "Novel_CDS_in_progress" => "novelCDSDensity", 
);
$support->log("Available gene types:\n    ");
$support->log(join("\n    ", sort keys %gene_types)."\n");

my (%density_types, $dtcache);
foreach my $type (keys %gene_types) {
    $density_types{$gene_types{$type}} = 1;
    my $analysis = new Bio::EnsEMBL::Analysis (
        -program     => "vega_gene_density.pl",
        -database    => "ensembl",
        -gff_source  => "vega_gene_density.pl",
        -gff_feature => "density",
        -logic_name  => $gene_types{$type}
    );
    $aa->store($analysis) unless ($support->param('dry_run'));
    foreach my $block_size (keys %{ $chr_slices }) {
        my $dt = Bio::EnsEMBL::DensityType->new(
            -analysis   => $analysis,
            -block_size => $block_size,
            -value_type => 'sum'
        );
        $dta->store($dt) unless ($support->param('dry_run'));
        $dtcache->{$block_size}->{$gene_types{$type}} = $dt;
    }
}

# loop over block sizes
foreach my $block_size (keys %{ $chr_slices }) {
    $support->log("Available chromosomes using block size of $block_size:\n    ");
    $support->log(join("\n    ", map { $_->seq_region_name } @{ $chr_slices->{$block_size} })."\n");

    # looping over chromosomes
    $support->log_stamped("Looping over chromosomes...\n");
    my ($current_start, $current_end);
    foreach my $slice (@{ $chr_slices->{$block_size} }) {
        $current_start = 1;
        my $chr = $slice->seq_region_name;
        my (%total, $i, %gene_names);
        my $bins = POSIX::ceil($slice->end/$block_size);
        
        $support->log_stamped("Chromosome $chr with block size $block_size...\n", 1);
        
        # loop over blocks
        my @density_features;
        while($current_start <= $slice->end) {
            $i++;
            $current_end = $current_start + $block_size - 1;
            if ($current_end > $slice->end) {
                $current_end = $slice->end;
            }
            my $sub_slice = $slice->sub_Slice($current_start, $current_end);
            my %num = ();
            
            # count genes by type
            my $genes;
            eval { $genes = $sub_slice->get_all_Genes; };
            if ($@) {
                $support->log_warning("$@");
                $current_start = $current_end + 1;
                next;
            }
            foreach my $gene (@{$genes}) {
                # only count genes that don't overlap the subslice start
                # (since these were already counted in the last bin)
                if ($gene->start >= 1) {
                    $total{$gene->type}++;
                }
                $num{$gene_types{$gene->type}}++;
            }
            
            # create DensityFeature objects for each type
            foreach my $type (keys %density_types) {
                push @density_features, Bio::EnsEMBL::DensityFeature->new(
                    -seq_region    => $slice,
                    -start         => $current_start,
                    -end           => $current_end,
                    -density_type  => $dtcache->{$block_size}->{$type},
                    -density_value => $num{$type} || 0
                );
            }
            $current_start = $current_end + 1;
            
            # logging
            $support->log_verbose("Chr: $chr | Bin: $i/$bins | Counts: ".
                join(",", map { $num{$gene_types{$_}} || 0 }
                    sort keys %gene_types)."\n", 2);
        }
        
        # store DensityFeatures for the chromosome
        $dfa->store(@density_features) unless ($support->param('dry_run'));
        
        # stats
        my @attribs;
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => '12:Known genes',
            -CODE => 'KnownGeneCount',
            -VALUE => $total{'Known'} || 0,
            -DESCRIPTION => 'Total Number of Known genes'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => '14:Novel CDS',
            -CODE => 'NovelCDSCount',
            -VALUE => $total{'Novel_CDS'} || 0,
            -DESCRIPTION => 'Total Number of Novel CDSs'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => '16:Novel transcripts',
            -CODE => 'NovelTransCount',
            -VALUE => $total{'Novel_Transcript'} || 0,
            -DESCRIPTION => 'Total Number of Novel transcripts'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => '24:Putative',
            -CODE => 'PutTransCount',
            -VALUE => $total{'Putative'} || 0,
            -DESCRIPTION => 'Total Number of Putative transcripts'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => '30:Predicted transcripts',
            -CODE => 'PredTransCount',
            -VALUE => $total{'Predicted_Gene'} || 0,
            -DESCRIPTION => 'Total Number of Predicted transcripts'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => '26:Ig segments',
            -CODE => 'IgSegCount',
            -VALUE => $total{'Ig_Segment'} || 0,
            -DESCRIPTION => 'Total Number of Ig Segments'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => '28:Ig pseudogene segments',
            -CODE => 'IgPsSegCount',
            -VALUE => $total{'Ig_Pseudogene_Segment'} || 0,
            -DESCRIPTION => 'Total Number of Ig Pseudogene Segments'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => '18:Total pseudogenes',
            -CODE => 'TotPsCount',
            -VALUE => ($total{'Pseudogene'}
                        + $total{'Processed_pseudogene'}
                        + $total{'Unprocessed_pseudogene'}) || 0,
            -DESCRIPTION => 'Total Number of Pseudogenes'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => '20:Processed pseudogenes',
            -CODE => 'ProcPsCount',
            -VALUE => $total{'Processed_pseudogene'} || 0,
            -DESCRIPTION => 'Number of Processed pseudogenes'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => '22:Unprocessed pseudogenes',
            -CODE => 'UnprocPsCount',
            -VALUE => $total{'Unprocessed_pseudogene'} || 0,
            -DESCRIPTION => 'Number of Unprocessed pseudogenes'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => '13:Known genes (in progress)',
            -CODE => 'KnwnprogCount',
            -VALUE => $total{'Known_in_progress'} || 0,
            -DESCRIPTION => 'Number of Known Genes in progress'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => '15:Novel CDS (in progress)',
            -CODE => 'NovCDSprogCount',
            -VALUE => $total{'Novel_CDS_in_progress'} || 0,
            -DESCRIPTION => 'Number of novel CDS in progress'
        );
        
        # only store unclassified pseudogenes if there are no processed and
        # unprocessed pseudos, ie if total pseudos eq pseudos
        unless ($total{'Unprocessed_pseudogene'} == 0 && $total{'Processed_pseudogene'} == 0) {  
            push @attribs, Bio::EnsEMBL::Attribute->new
                (-NAME => '23:Unclassified pseudogenes',
                 -CODE => 'UnclassPsCount',
                 -VALUE => $total{'Pseudogene'} || 0,
                 -DESCRIPTION => 'Number of Unclassified pseudogenes');
        }
        
        $attrib_adaptor->store_on_Slice($slice, \@attribs) unless ($support->param('dry_run'));
        
        # log stats
        $support->log("\n");
        $support->log("Totals for chr $chr:\n", 1);
        my $fmt = "%-30s%-40s\n";
        $support->log(sprintf($fmt, qw(TYPE COUNT)), 2);
        $support->log(('-'x70)."\n", 2);
        map { $support->log(sprintf($fmt, $_, $total{$_}), 2) } sort keys %total;
        $support->log_stamped("Done.\n", 1);
        $support->log("\n");
    }
    $support->log_stamped("Done.\n");
}

# finish logfile
$support->finish_log;

