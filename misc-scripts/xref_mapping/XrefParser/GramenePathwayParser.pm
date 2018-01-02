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

package XrefParser::GramenePathwayParser;

=pod

=head1 NAME

XrefParser::GramenePathwayParser

=head1 DESCRIPTION

Parse pathway dumps from Gramene.  File format (and example data):

 gene_name                  AT1G66030
 enzyme_name                fatty acid (omega-1)-hydroxylase
 reaction_id                RXN-7796
 reaction_name              
 ec                         2.7.7.-
 pathway_id                 PWY-5129
 pathway_name               sphingolipid biosynthesis (plants)

=head1 AUTHOR

Ken Youens-Clark E<lt>kclark@cshl.eduE<gt>.

=cut

use strict;
use Text::RecordParser::Tab;
use base 'XrefParser::BaseParser';

sub run {
    my ($self, $args) = @_;
    my $source_id     = $args->{'source_id'};
    my $species_id    = $args->{'species_id'};
    my $files         = $args->{'files'};
    my $release_file  = $args->{'rel_file'};
    my $verbose       = $args->{'verbose'};
    my $file          = ref $files eq 'ARRAY' ? shift @$files : '';

    if ( !$file ) {
        printf STDERR "%s called without a 'files' argument\n%s", 
            __PACKAGE__, Dumper($args);
        return 1; # error
    }

    my $p = Text::RecordParser::Tab->new( $file );

    my $direct_xref_count = 0;
    while ( my $rec = $p->fetchrow_hashref ) {
        my $gene = $rec->{'gene_name'} or next;

        if ( my $ec = $rec->{'ec'} ) {
            my $ec_xref_id = $self->add_xref({
                source_id  => $source_id,
                species_id => $species_id,
                acc        => $ec,
                label      => '',
                desc       => '',
                info_type  => 'DIRECT',
            });

            $self->add_direct_xref( $ec_xref_id, $gene, 'Gene', 'DIRECT' );
            $direct_xref_count++;
        }

        if ( my $pathway_id = $rec->{'pathway_id'} ) {
            my $pathway_xref_id = $self->add_xref({
                source_id  => $source_id,
                species_id => $species_id,
                acc        => $pathway_id,
                label      => $rec->{'pathway_name'},
                desc       => '',
                info_type  => 'DIRECT'
            });

            $self->add_direct_xref( $pathway_xref_id, $gene, 'Gene', 'DIRECT' );
            $direct_xref_count++;
        }
    }

    printf "Parsed pathway Ids from file '%s,' added %s direct_xrefs\n",
        $file, $direct_xref_count;

    return 0; # success
}

1;
