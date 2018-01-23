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

package XrefParser::ECParser;

=pod 

=head1 NAME

XrefParser::ECParser

=head1 DESCRIPTION

Parse the gene to EC text file
The format is;
<gene_loc>\t<EC>
multiple ECs maybe separated by ', '
the file got from parsing pathway pathologic files. 
the annotation got from GO annotation and ec2go and protein EC.

NB: This came from the Gramene project.  Updated by Ken to work 
at EBI, 2/2012.

=head1 AUTHOR

Ken Youens-Clark E<lt>kclark@cshl.eduE<gt>.

=cut

# ------------------------------------------------------------------

use strict;
use base qw( XrefParser::BaseParser );

sub run {
    my ($self, $args) = @_;
    my $source_id     = $args->{'source_id'};
    my $species_id    = $args->{'species_id'};
    my $files         = $args->{'files'};
    my $release_file  = $args->{'rel_file'};
    my $verbose       = $args->{'verbose'};

    my $file = ref $files eq 'ARRAY' ? shift @$files : '';

    return unless $file;

    my $file_io = $self->get_filehandle($file);

    my @xrefs;
    my ( $xref_count, $direct_xref_count ) = ( 0, 0 );
    while ( my $line = $file_io->getline() ) {
        chomp $line;

        my ( $stable_id, $ecs ) = split( /\t/, $line );

        for my $acc (  split /\s*,\s*/, $ecs ) {
            my $xref_id = $self->add_xref({
                acc        => $acc,
                label      => $acc,
                desc       => '',
                source_id  => $source_id,
                species_id => $species_id,
                info_type  => 'DIRECT'
            });

            $xref_count++;

            $self->add_direct_xref( $xref_id, $stable_id, 'Gene', '' );

            $direct_xref_count++;
        }
    }

    printf "Parsed EC IDs from %s\nadded %sxrefs and %s direct_xrefs\n",
        $file, $xref_count, $direct_xref_count,
    ;

    return 0;
}

1;
