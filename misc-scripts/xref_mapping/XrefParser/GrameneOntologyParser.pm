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

package XrefParser::GrameneOntologyParser;

=pod

=head1 NAME

XrefParser::GrameneOntologyParser

=head1 DESCRIPTION

Parse ontology files.  File format (and example data):

 1: DB, database contributing the file 
 2: DB_Object_ID (gene stable id)
 3: DB_Object_Symbol, see below
 4: Qualifier (optional), one or more of 'NOT', 'contributes_to',
    'colocalizes_with' as qualifier(s) for a GO annotation, when needed,
    multiples separated by pipe (|)
 5: Ontology accession, e.g., "PO:0007633"
 6: DB:Reference(|DB:Reference), the reference associated with the GO
    annotation
 7: Evidence, the evidence code for the GO annotation
 8: With (or) From (optional), any With or From qualifier for the GO
    annotation
 9: Aspect, which ontology the GO term belongs (Function, Process or
    Component)
10: DB_Object_Name(|Name) (optional), a name for the gene product in
    words, e.g. 'acid phosphatase'
11: DB_Object_Synonym(|Synonym) (optional), see below
12: DB_Object_Type, type of object annotated, e.g. gene, protein, etc.
13: taxon(|taxon), taxonomic identifier of species encoding gene
    product
14: Date, date GO annotation was made in the format
15: Assigned_by, source of the annotation (either "TAIR" or "TIGR")

TAIR
locus:2007106
EXPA1

PO:0000293
TAIR:Publication:501680952
NAS

S
AT1G69530
AT1G69530|ATEXPA1|EXP1|AT-EXP1|ATEXP1|ATHEXP ALPHA 1.2|EXPA1|expansin A1|EXPANSIN 1|F10D13.18|F10D13_18
protein
taxon:3702
20040428
TAIR

TAIR:locus:2007106

=head1 AUTHOR

Ken Youens-Clark E<lt>kclark@cshl.eduE<gt>.

=cut

use strict;
use POSIX qw(strftime);
use File::Basename;

use base qw( XrefParser::BaseParser );

# Creates ontology xrefs from a GAF annotation file.
# Also requires an OBO file mapping IDs to terms.
# The ontology xrefs are dependent on xrefs 

#----------------------------------------------------------------------
sub run {
    my ($self, $args) = @_;
    my $source_id     = $args->{'source_id'};
    my $species_id    = $args->{'species_id'};
    my $files         = $args->{'files'};
    my $release_file  = $args->{'rel_file'};
    my $notify        = sub { print @_ if $args->{'verbose'} };
    my $file          = ref $files eq 'ARRAY' ? shift @$files : '';

    if ( $file ) {
        $notify->(sprintf "%s parsing file '$file'\n", __PACKAGE__);
    }
    else {
        printf STDERR "%s called without a 'files' argument\n%s", 
            __PACKAGE__, Dumper($args);
        return 1; # error
    }


    # get the 'PO' or 'GO' source_name
    my $source_name = $self->get_source_name_for_source_id($source_id);

    # get the "main" GO/PO source id.
    $source_id = $self->get_source_id_for_source_name( $source_name, 'main' );

    $species_id       ||= $self->get_species_id_for_filename($file);
    my %species_id2name = $self->species_id2name();
    my $species         = $species_id2name{$species_id}->[0];

    $notify->("Parsing $source_name for $species\n");

    # These Ontology xrefs depend on gene ID xrefs (TAIR_LOCUS) having been
    # previously loaded - get the source.
    my $gene_source_id = $self->get_source_id_for_source_name('TAIR_LOCUS');

    #
    # Get mappimg from GO terms to descriptions from the obo file.
    #
    # The gene_association file in GAF format
    # http://www.geneontology.org/GO.format.gaf-2_0.shtml

    # The obo file (xml or plain)
    my %id_to_term;
    if ( my $file_desc = shift @$files ) {
        my $obo_io = $self->get_filehandle($file_desc);
        if ( !defined $obo_io ) {
            print STDERR "ERROR: Could not open $file_desc\n";
            return 1;    # 1 error
        }

        my $term = undef;
        my $desc = undef;

        if ( $file_desc =~ /xml/ ) {    # xml file format
            while ( $_ = $obo_io->getline() ) {
                if (/\<id\>([G|P]O:\d+)\<\/id\>/) {
                    $term = $1;
                }
                elsif (/\<name\>(.*)\<\/name\>/) {
                    if ( defined($term) ) {
                        $id_to_term{$term} = $1;
                    }
                    $term = undef;
                }
            }
        }
        else {
            while ( $_ = $obo_io->getline() ) {
                if (/^id:\s+([G|P]O:\d+)/) {
                    $term = $1;
                }
                elsif (/^name:\s+(.*)/) {
                    if ( defined($term) ) {
                        $id_to_term{$term} = $1;
                    }
                    $term = undef;
                }
            }
        }
        $obo_io->close();

        $notify->(
            sprintf "Mapped %s terms to descriptions\n", 
            scalar keys %id_to_term
        );
    }

    my %wrongtype;

    # Get a list of taxon IDs given that some species have multiple e.g. strains
    my %sp2tax  = $self->species_id2taxonomy();
    my @tax_ids = @{ $sp2tax{ $species_id } };

    # Loop for each taxon
    my $count         = 0;
    my $miss_parse    = 0;
    my $miss_get_xref = 0;

    for my $tax_id ( @tax_ids ) {
        $notify->("processing for taxon: $tax_id\n");
        my $taxon_line = "taxon:" . $tax_id;

        # Open the gene association GAF file
        my $gaf_io = $self->get_filehandle( $file );

        if ( !defined $gaf_io ) {
            print STDERR "ERROR: Could not open $file\n";
            return 1; # error
        }

        while ( my $line = $gaf_io->getline() ) {
            next if $line =~ /^!/;     # comment
            # next unless /$taxon_line/; # only process data for this taxon

            chomp;
            my (
                $DB,                 $DB_Object_ID,
                $DB_Object_Symbol,   $Qualifier,
                $GO_ID,              $DB_Reference,
                $Evidence_Code,      $With,
                $Aspect,             $DB_Object_Name,
                $DB_Object_Synonym,  $DB_Object_Type,
                $Taxon,              $Date,
                $Assigned_By,        $Annotation_Extension,
                $Gene_Product_Form_ID, @extra
            ) = split( /\t/, $line );

            # next unless $Taxon eq $taxon_line;   # Only process this species
            next if $Qualifier eq "NOT";         # Skip "NOT" terms entirely

            $DB_Object_Name =~ s/\'/\\\'/g;
            my $master = 0;

            # in po file stable_id=$array[9] # in GO asso file, stable_id=
            # ??[9] [10]

            # Find the stable_id in either DB_Object_Name or DB_Object_Synonym
            if ( my $stable_id = $DB_Object_Name || $DB_Object_Synonym ) {
                my $xref_id = $self->add_xref({
                    acc           => $stable_id, 
                    label         => '', 
                    desc          => $nasc_gene_id, 
                    source_id     => $nascg_source_id, 
                    species_id    => $species_id, 
                    info_type     => 'DIRECT',
                });

                $self->add_direct_xref( 
                    $tairg_xref_id, $gene_stable_id, 'Gene', 'DIRECT' 
                );
            $count++;
        } 

        $gaf_io->close();
    }

    $notify->(
        map { " - $_\n" }
        "$count $source_name dependent xrefs added",
        "$miss_parse lines did not contain a recognised stable_id",
        "$miss_get_xref lines contained a stable_id that did not " .
            "correspond with a known xref",
    );

    if ( defined $release_file ) {
        # Parse and set release information from $release_file.
        my $release_io = $self->get_filehandle($release_file);

        # Slurp mode.
        local $/;
        my $release = <$release_io>;
        $release_io->close();

        $release =~ tr/\n/ /;
        $release =~ s#.*The following table describes.*?of (POC.*?)<ul>.*#$1#;
        $release =~ s#<[^>]+>##g;

        $notify->("$source_name release: '$release'\n");
        $self->set_release( $source_id, $release );
    }

    return 0; # success
}

1;
