#!/usr/local/bin/perl

=head1 NAME

vega_seq_region_stats.pl -
script to calculate gene and snp numbers for mapview stats

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Patrick Meidl <pm2@sanger.ac.uk>

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=cut

use strict;
BEGIN {
    $ENV{'ENSEMBL_SERVERROOT'} = "../../..";
    unshift(@INC,"$ENV{'ENSEMBL_SERVERROOT'}/conf");
    unshift(@INC,"$ENV{'ENSEMBL_SERVERROOT'}/ensembl-compara/modules");
    unshift(@INC,"$ENV{'ENSEMBL_SERVERROOT'}/ensembl-draw/modules");
    unshift(@INC,"$ENV{'ENSEMBL_SERVERROOT'}/ensembl-external/modules");
    unshift(@INC,"$ENV{'ENSEMBL_SERVERROOT'}/ensembl-otter/modules");
    unshift(@INC,"$ENV{'ENSEMBL_SERVERROOT'}/modules");
    unshift(@INC,"$ENV{'ENSEMBL_SERVERROOT'}/ensembl/modules");
    unshift(@INC,"$ENV{'ENSEMBL_SERVERROOT'}/bioperl-live");
}

use SiteDefs;
use EnsWeb;
use EnsEMBL::DB::Core;
use Getopt::Long;

my ($species, $run, $help);
&GetOptions(
    "species=s" => \$species,
    "run=s"	=> \$run,
    "help"      => \$help,
    "h"         => \$help,
);

if($help || !($species && $run)){
    print qq(Usage:
    ./vega_seq_region_stats.pl
        --species=Homo_sapiens
        --run=gene,snp
        [--help]\n\n);
    exit;
}

$ENV{'ENSEMBL_SPECIES'} = $species;

## set db user/pass to allow write access
my $db_ref = $EnsWeb::species_defs->databases;
$db_ref->{'ENSEMBL_DB'}{'USER'} = $EnsWeb::species_defs->ENSEMBL_WRITE_USER;
$db_ref->{'ENSEMBL_DB'}{'PASS'} = $EnsWeb::species_defs->ENSEMBL_WRITE_PASS;

## connect to databases
my %run;
@run{split(/,/, $run)} = 1;
my @dbs = qw(core);
push @dbs, 'glovar' if ($run{'snp'});
my $databases = &EnsEMBL::DB::Core::get_databases(@dbs);

die "Problem connecting to databases: $databases->{'error'}\n"
    if  $databases->{'error'} ;
warn "Database error: $databases->{'non_fatal_error'}\n"
    if $databases->{'non_fatal_error'};

## get slice and attribute adaptors, loop over all toplevel slices
my $slice_adaptor = $databases->{'core'}->get_SliceAdaptor();
my $attrib_adaptor = $databases->{'core'}->get_AttributeAdaptor();
my $top_slices = $slice_adaptor->fetch_all("toplevel");

foreach my $slice (@$top_slices) {
    my %num;
    my @attribs;
    
    print STDERR "Processing seq_region ", $slice->seq_region_name(), ":\n";

    ## genes
    if ($run{'gene'}) {
        print STDERR "\tGenes...\n";
        my @genes = @{$slice->get_all_Genes()};
        foreach my $gene (@genes) {
            my $type = $gene->type;
            $num{$type}++;
            if ($type =~ /seudogene$/) {
                $num{'Total_Pseudogenes'}++;
            }
        }

        print STDERR "Slice", $slice->seq_region_name(),
            " has the following features:\n\n";
        foreach my $type (keys %num) {
            print STDERR "$type = $num{$type}\n";
        }

        push @attribs, Bio::EnsEMBL::Attribute->new
        (-NAME => 'Known genes',
         -CODE => 'KnownGeneCount',
         -VALUE => $num{'Known'},
         -DESCRIPTION => 'Total Number of Known genes');

        push @attribs, Bio::EnsEMBL::Attribute->new
        (-NAME => 'Novel CDS',
         -CODE => 'NovelCDSCount',
         -VALUE => $num{'Novel_CDS'},
         -DESCRIPTION => 'Total Number of Novel CDSs');

        push @attribs, Bio::EnsEMBL::Attribute->new
        (-NAME => 'Novel transcripts',
         -CODE => 'NovelTransCount',
         -VALUE => $num{'Novel_Transcript'},
         -DESCRIPTION => 'Total Number of Novel transcripts');

        push @attribs, Bio::EnsEMBL::Attribute->new
        (-NAME => 'Putative transcripts',
         -CODE => 'PutTransCount',
         -VALUE => $num{'Putative'},
         -DESCRIPTION => 'Total Number of Putative transcripts');

        push @attribs, Bio::EnsEMBL::Attribute->new
        (-NAME => 'Predicted transcripts',
         -CODE => 'PredTransCount',
         -VALUE => $num{'Predicted_Gene'},
         -DESCRIPTION => 'Total Number of Predicted transcripts');

        push @attribs, Bio::EnsEMBL::Attribute->new
        (-NAME => 'Total Pseudogenes',
         -CODE => 'TotPsCount',
         -VALUE => $num{'Total_Pseudogenes'},
         -DESCRIPTION => 'Total Number of Pseudogenes');

        push @attribs, Bio::EnsEMBL::Attribute->new
        (-NAME => 'Unclassified pseudogenes',
         -CODE => 'UnclassPsCount',
         -VALUE => $num{'Pseudogene'},
         -DESCRIPTION => 'Number of Unclassified pseudogenes');

        push @attribs, Bio::EnsEMBL::Attribute->new
        (-NAME => 'Processed pseudogenes',
         -CODE => 'ProcPsCount',
         -VALUE => $num{'Processed_pseudogene'},
         -DESCRIPTION => 'Number of Processed pseudogenes');

        push @attribs, Bio::EnsEMBL::Attribute->new
        (-NAME => 'Unprocessed pseudogenes',
         -CODE => 'UnprocPsCount',
         -VALUE => $num{'Unprocessed_pseudogene'},
         -DESCRIPTION => 'Number of Unprocessed pseudogenes');

        push @attribs, Bio::EnsEMBL::Attribute->new
        (-NAME => 'Ig Segments',
         -CODE => 'IgSegCount',
         -VALUE => $num{'Ig_Segment'},
         -DESCRIPTION => 'Total Number of Ig Segments');

        push @attribs, Bio::EnsEMBL::Attribute->new
        (-NAME => 'Ig Pseudogene Segments',
         -CODE => 'IgPsSegCount',
         -VALUE => $num{'Ig_Pseudogene_Segment'},
         -DESCRIPTION => 'Total Number of Ig Pseudogene Segments');
    }

    ## snps
    if ($run{'snp'}) {
        print STDERR "\tSNPs...\n";
	my $snps = $slice->get_all_ExternalFeatures('GlovarSNP');
	push @attribs, Bio::EnsEMBL::Attribute->new
	(-NAME => 'SNPs',
	 -CODE => 'SNPCount',
	 -VALUE => scalar( @$snps ),
	 -DESCRIPTION => 'Total Number of SNPs');
    }
    
    $attrib_adaptor->store_on_Slice($slice, \@attribs);
    #  print_chromo_stats([$slice]);
}



sub print_chromo_stats {
    my $chromosomes = shift;
    
    foreach my $chr (@$chromosomes) {
	print "\nchromosome: ",$chr->seq_region_name(),"\n";
	foreach my $attrib (@{$chr->get_all_Attributes()}) {
	    print "  ", $attrib->name(), ": ", $attrib->value(), "\n";
	}
    }
}

1;


