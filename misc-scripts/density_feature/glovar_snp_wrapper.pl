#!/usr/local/bin/perl

=head1 NAME

glovar_snp_wrapper.pl -
Wrapper for glovar_snp_density.pl

=head1 SYNOPSIS

    ./glovar_snp_density.pl
        --species=Homo_sapiens
        [--dry-run|-n]
        [--avdump|-a]

=head1 DESCRIPTION

Wrapper for glovar_snp_density.pl to run it chromosome by chromosome. This is
an attempt to avoid high memory footprints caused by a memory leak somerwhere
in the API. See glovar_snp_density.pl for more detailed documentation.

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
use Getopt::Long;

my ($species, $dry, $avdump);
&GetOptions(
    "species=s" => \$species,
    "dry-run"   => \$dry,
    "n"         => \$dry,
    "avdump"    => \$avdump,
    "a"         => \$avdump,
);

unless ($species) {
    print qq(Usage:
    ./glovar_snp_density.pl
        --species=Homo_sapiens
        [--avdump|-a]
        [--dry-run|-n]\n\n);
    exit;
}

$ENV{'ENSEMBL_SPECIES'} = $species;

## run glovar_snp_density.pl for each chromsome in this species
my $command = "./glovar_snp_density.pl --species=$species";
$command .= " -n" if ($dry);
$command .= " -a" if ($avdump);
foreach my $chr (@{$EnsWeb::species_defs->ENSEMBL_CHROMOSOMES}) {
    warn "$command --chr=$chr\n";
    system("$command --chr=$chr") == 0 or die "$command --chr=$chr failed: $!\n";
}

