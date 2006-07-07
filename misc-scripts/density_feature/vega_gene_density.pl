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
for use in mapview. It also checks for new biotype/status pairs and warns to
adapt the appropriate modules to deal with them.

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
my $dbh = $dba->dbc->db_handle;

# split chromosomes by size and determine block size
my $chr_slices = $support->split_chromosomes_by_size(5000000);

# known biotype/status pairs and associated density type
my %gene_types = (
    "protein_coding_KNOWN"             => "knownPCodDensity",
    "protein_coding_in_progress_KNOWN" => "knownPCodDensity",
	"processed_transcript_KNOWN"       => "knownPTransDensity",
    "protein_coding_NOVEL"             => "novelPCodDensity",
    "protein_coding_in_progress_NOVEL" => "novelPCodDensity",
    "processed_transcript_NOVEL"       => "novelPTransDensity",
    "processed_transcript_PUTATIVE"    => "putativePTransDensity",
    "protein_coding_PREDICTED"         => "predictedPCodDensity",
    "Ig_pseudogene_segment_UNKNOWN"    => "IgPseudoSegDensity",
    "Ig_segment_NOVEL"                 => "IgSegDensity",
	"Ig_segment_KNOWN"                 => "IgSegDensity",
    "pseudogene_UNKNOWN"               => "pseudoGeneDensity",
    "processed_pseudogene_UNKNOWN"     => "pseudoGeneDensity",
    "unprocessed_pseudogene_UNKNOWN"   => "pseudoGeneDensity",
);

# check for new biotype/status pairs
my $sql = qq(
    SELECT biotype, status
    FROM gene
    GROUP by biotype, status
);
my $sth = $dbh->prepare($sql);
$sth->execute;
my (%type_status, $new);
while (my ($biotype, $status) = $sth->fetchrow_array) {
	my $type = $biotype.'_'.$status;
    if ($gene_types{$type}) {
        $type_status{$type} = 'no';
    } else {
        $type_status{$type} = 'YES';
        $new = 1;
    }
}
my $FMT = "%-50s%-20s\n";
$support->log("Checking for new biotype/status pairs...\n\n");
$support->log(sprintf($FMT, qw(BIOTYPE/STATUS NEW)), 1);
$support->log(('-'x70)."\n", 1);
map { $support->log(sprintf($FMT, $_, $type_status{$_}), 1) }
    sort keys %type_status;
$support->log("\n");
if ($new) {
    $support->log_warning("There are new biotype/status pairs! You might need to adapt Bio::EnsEMBL::ColourMap, EnsEMBL::Sanger_vega::Component::Chromosome and configure mapview to show them.\n\n");
}


# create Analysis and DensityType objects
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
    $support->log_stamped("\nLooping over chromosomes...\n");
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
                my $gene_type = $gene->biotype . '_' . $gene->status;
                if ($gene->start >= 1) {
                    $total{$gene_type}++;
                }
                $num{$gene_types{$gene_type}}++;
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
            -NAME => 'protein_coding_KNOWN',
            -CODE => 'KnownPCCount',
            -VALUE => $total{'protein_coding_KNOWN'} || 0,
            -DESCRIPTION => 'Number of Known Protein Coding',
        );

        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => 'processed_transcript_KNOWN',
            -CODE => 'KnownPTCount',
            -VALUE => $total{'processed_transcript_KNOWN'} || 0,
            -DESCRIPTION => 'Number of Known Processed Transcripts',
        );

        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => 'protein_coding_NOVEL',
            -CODE => 'NovelPCCount',
            -VALUE => $total{'protein_coding_NOVEL'} || 0,
            -DESCRIPTION => 'Number of Novel Protein Coding'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => 'processed_transcript_NOVEL',
            -CODE => 'NovelPTCount',
            -VALUE => $total{'processed_transcript_NOVEL'} || 0,
            -DESCRIPTION => 'Number of Novel transcripts'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => 'processed_transcript_PUTATIVE',
            -CODE => 'PutPTCount',
            -VALUE => $total{'processed_transcript_PUTATIVE'} || 0,
            -DESCRIPTION => 'Number of Putative processed transcripts'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => 'protein_coding_PREDICTED',
            -CODE => 'PredPCCount',
            -VALUE => $total{'protein_coding_PREDICTED'} || 0,
            -DESCRIPTION => 'Number of Predicted transcripts'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => 'total_Ig_segment_',
            -CODE => 'IgSegCount',
            -VALUE => $total{'Ig_segment_NOVEL'}
						+ $total{'Ig_segment_KNOWN'} || 0,
            -DESCRIPTION => 'Number of Ig Segments'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => 'Ig_pseudogene_segment_',
            -CODE => 'IgPsSegCount',
            -VALUE => $total{'Ig_pseudogene_segment_'} || 0,
            -DESCRIPTION => 'Number of Ig Pseudogene Segments'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => 'total_pseudogene_',
            -CODE => 'TotPsCount',
            -VALUE => ($total{'pseudogene_'}
                        + $total{'processed_pseudogene_'}
                        + $total{'unprocessed_pseudogene_'}) || 0,
            -DESCRIPTION => 'Total Number of Pseudogenes'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => 'processed_pseudogene_',
            -CODE => 'ProcPsCount',
            -VALUE => $total{'processed_pseudogene_'} || 0,
            -DESCRIPTION => 'Number of Processed pseudogenes'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => 'unprocessed_pseudogene_',
            -CODE => 'UnprocPsCount',
            -VALUE => $total{'unprocessed_pseudogene_'} || 0,
            -DESCRIPTION => 'Number of Unprocessed pseudogenes'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => 'protein_coding_in_progress_KNOWN',
            -CODE => 'KnwnPCProgCount',
            -VALUE => $total{'protein_coding_in_progress_KNOWN'} || 0,
            -DESCRIPTION => 'Number of Known Protein coding genes in progress'
        );
        
        push @attribs, Bio::EnsEMBL::Attribute->new(
            -NAME => 'protein_coding_in_progress_NOVEL',
            -CODE => 'NovPCProgCount',
            -VALUE => $total{'protein_coding_in_progress_NOVEL'} || 0,
            -DESCRIPTION => 'Number of Novel Protein coding genes in progress'
        );

        $attrib_adaptor->store_on_Slice($slice, \@attribs) unless ($support->param('dry_run'));
        
        # log stats
        $support->log("\n");
        $support->log("Totals for chr $chr:\n", 1);
        $support->log(sprintf($FMT, qw(TYPE COUNT)), 2);
        $support->log(('-'x70)."\n", 2);
        map { $support->log(sprintf($FMT, $_, $total{$_}), 2) } sort keys %total;
        $support->log_stamped("\nDone.\n\n", 1);
    }
    $support->log_stamped("Done.\n");
}

# finish logfile
$support->finish_log;

