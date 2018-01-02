=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package XrefParser::TAIRIDParser;

=pod

=head1 NAME

XrefParser::TAIRIDParser

=head1 DESCRIPTION

Takes the TAIR cDNA FASTA file and uses the description lines to
create direct xrefs.

Sample data:

>AT5G55930.1 | Symbols: ATOPT1, OPT1 | oligopeptide transporter 1 | chr5:22652715-22656106 FORWARD LENGTH=2820
CAAAATTCATGTGGTGTAAATTGTCTAAAGTCTGTATTTTTTTTTATTGACATCCATTTTTTTTGTGTCGAAAGTCTAT
...

=head1 AUTHOR

Ken Youens-Clark E<lt>kclark@cshl.eduE<gt>.

=cut

use strict;

use base qw( XrefParser::BaseParser );

use Data::Dumper;
use Readonly;

Readonly my $EMPTY_STR        => q{};
Readonly my $DIRECT           => 'DIRECT';
Readonly my $GENE             => 'Gene';
Readonly my $NASC_GENE_ID     => 'NASC_GENE_ID';
Readonly my $TAIR_LOCUS_MODEL => 'TAIR_LOCUS_MODEL';
Readonly my $TAIR_LOCUS       => 'TAIR_LOCUS';
Readonly my $TAIR_SYMBOL      => 'TAIR_SYMBOL';
Readonly my $TAIR_TRANSLATION => 'TAIR_TRANSLATION';
Readonly my $TRANSCRIPT       => 'Transcript';
Readonly my $TRANSLATION      => 'Translation';


# -----------------------------------------------------------------
sub run {
    my ( $self, $args ) = @_;

    my $notify = sub { print @_, "\n" if $args->{'verbose'} };

    my $files = $args->{'files'};
    my $file  = ref $files eq 'ARRAY' ? shift @$files : $EMPTY_STR;

    if ( $file ) {
        $notify->(sprintf "%s Processing file '%s'", __PACKAGE__, $file); 
    }
    else {
        printf STDERR "%s called without a 'files' argument\n%s", 
            __PACKAGE__, Dumper($args);
        return 1; # error
    }

    my $tair_io = $self->get_filehandle($file);

    if ( !defined $tair_io ) {
        print STDERR "ERROR: Could not open $file\n";
        return 1; # 1 is an error
    }

    my $source_id  = $args->{'source_id'}  ||
                     XrefParser::BaseParser->get_source_id_for_filename($file);
    my $species_id = $args->{'species_id'} ||
                     XrefParser::BaseParser->get_species_id_for_filename($file);

    my $tairg_source_id = $self->get_source_id_for_source_name($TAIR_LOCUS);
    my $tairs_source_id = $self->get_source_id_for_source_name($TAIR_SYMBOL);
    my $tairl_source_id = 
        $self->get_source_id_for_source_name($TAIR_LOCUS_MODEL);
    my $tairt_source_id =
        $self->get_source_id_for_source_name($TAIR_TRANSLATION);
    my $nascg_source_id = $self->get_source_id_for_source_name($NASC_GENE_ID);

    my $line_num = 0;
    my %xrefs_added;
    while ( my $line = $tair_io->getline() ) {
        # Only process FASTA header lines
        next unless $line =~ s/^>//; 

        chomp $line;

        my ( $gene_stable_id, $symbol_str, $desc ) = split /\s*\|\s*/, $line;

        next unless $gene_stable_id;

        $desc ||= $EMPTY_STR;

        if ( $args->{'verbose'} ) {
            printf "%-70s\r",
                sprintf( '%9d: Processing %s', ++$line_num, $gene_stable_id )
            ;
        }

        #
        # Transcript, e.g., "AT5G55930.1"
        #
        if ( $gene_stable_id =~ /^([A-Z0-9]+) \. (\d+)$/xms ) {
            my $transcript_id = $gene_stable_id;
            $gene_stable_id   = $1;

            my $tairl_xref_id = $self->add_xref({ 
                source_id     => $tairl_source_id, 
                species_id    => $species_id, 
                info_type     => $DIRECT,
                acc           => $transcript_id, 
                label         => $transcript_id, 
                desc          => $desc,
            });

            $self->add_direct_xref( 
                $tairl_xref_id, $transcript_id, $TRANSCRIPT, $DIRECT 
            );

            $xrefs_added{ $TAIR_LOCUS_MODEL }++;

            my $tairt_xref_id =  $self->add_xref({
                source_id     => $tairt_source_id,
                species_id    => $species_id,
                info_type     => $DIRECT,
                acc           => $gene_stable_id,
                label         => $gene_stable_id,
                desc          => $EMPTY_STR,
            });

            $self->add_direct_xref( 
                $tairt_xref_id, $transcript_id, $TRANSLATION, $DIRECT 
            );

            $xrefs_added{ $TAIR_TRANSLATION }++;
        }

        #
        # Gene IDs for TAIR and NASC
        #
        my $tairg_xref_id =  $self->add_xref({
            source_id     => $tairg_source_id,
            species_id    => $species_id,
            info_type     => $DIRECT,
            acc           => $gene_stable_id,
            label         => $gene_stable_id,
            desc          => $desc,
        });

        $self->add_direct_xref( 
            $tairg_xref_id, $gene_stable_id, $GENE, $DIRECT 
        );

        $xrefs_added{ $TAIR_LOCUS }++;

        my $nascg_xref_id =  $self->add_xref({
            source_id     => $nascg_source_id,
            species_id    => $species_id,
            info_type     => $DIRECT,
            acc           => $gene_stable_id,
            label         => $gene_stable_id . '-TAIR-G',
            desc          => $desc,
        });

        $self->add_direct_xref( 
            $nascg_xref_id, $gene_stable_id, $GENE, $DIRECT 
        );

        $xrefs_added{ $NASC_GENE_ID }++;

        #
        # Symbols, e.g., "ATOPT1, OPT1"
        #
        if ( $symbol_str ) {
            $symbol_str =~ s/^\s*Symbols:\s*//;
            if ( my @symbols = map { $_ || () } split /\s*,\s*/, $symbol_str ) {
                if ( my $main_sym = shift @symbols ) {
                    my $sym_xref_id = $self->add_xref({
                        source_id   => $tairs_source_id,
                        species_id  => $species_id, 
                        info_type   => $DIRECT,
                        acc         => $main_sym,
                        label       => $main_sym,
                        desc        => $EMPTY_STR, 
                    });

                    #
                    # Add only first symbol to the gene
                    #
                    $self->add_direct_xref( 
                        $sym_xref_id, $gene_stable_id, $GENE, $DIRECT 
                    );

                    $xrefs_added{ $TAIR_SYMBOL }++;

                    #
                    # Add the remainder as "external_synonym"
                    #
                    for my $symbol ( @symbols ) {
                        $self->add_to_syn(
                            $main_sym, $tairs_source_id, 
                            $symbol, $species_id
                        );

                        $xrefs_added{'SYNONYMS'}++;
                    }
                }
            }
        }
    }

    $tair_io->close();

    $notify->(
        join("\n",
            $EMPTY_STR,
            map { 
                sprintf " - Added %9d %s xrefs", $xrefs_added{ $_ }, $_ 
            }
            sort keys %xrefs_added
        )
    );

    return 0; # successful
}

1;
