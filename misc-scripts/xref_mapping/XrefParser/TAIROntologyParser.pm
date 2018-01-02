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

package XrefParser::TAIROntologyParser;

=pod

=head1 NAME

XrefParser::TAIROntologyParser

=head1 DESCRIPTION

Parse ontology files.  File format (and example data):

   1: DB, database contributing the file (always "TAIR" for this file).
   2: DB_Object_ID  (TAIR's unique identifier for genes).
   3: DB_Object_Symbol, see below
   4: Qualifier (optional), one or more of 'NOT', 'contributes_to',
      'colocalizes_with' as qualifier(s) for a GO annotation, when needed,
      multiples separated by pipe (|)
   5: GO ID, unique numeric identifier for the GO term
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

E.g.:
   Field1: TAIR
   Field2: locus:2031476
   Field3: ENO1
   Field4: 
   Field5: GO:0000015
   Field6: TAIR:AnalysisReference:501748310
   Field7: IEA
   Field8: INTERPRO:IPR000941|INTERPRO:IPR020809|INTERPRO:IPR020810|INTERPRO:IPR020811
   Field9: C
  Field10: enolase 1
  Field11: AT1G74030|AT1G74030.1|F2P9.10|F2P9_10
  Field12: protein
  Field13: taxon:3702
  Field14: 20120418
  Field15: TAIR
  Field16: 
  Field17: TAIR:gene:2031475

   Field1: TAIR
   Field2: locus:2081840
   Field3: LTP12
   Field4: 
   Field5: PO:0001009
   Field6: TAIR:Publication:501710265|AGRICOLA_IND:IND23314018
   Field7: IDA
   Field8: 
   Field9: G
  Field10: AT3G51590
  Field11: AT3G51590|LTP12|lipid transfer protein 12|T18N14.1
  Field12: protein
  Field13: taxon:3702
  Field14: 20060215
  Field15: TAIR
  Field16: 
  Field17: TAIR:locus:2081840

=head1 AUTHOR

Ken Youens-Clark E<lt>kclark@cshl.eduE<gt>.

=cut

use strict;
use warnings;
use autodie;
use Data::Dumper;
use File::Basename 'basename';
use Readonly;
use List::MoreUtils 'uniq';
use XML::Simple 'XMLin';

use base qw( XrefParser::BaseParser );

Readonly my $DIRECT           => 'DIRECT';
Readonly my $GENE             => 'Gene';
Readonly my $TAIR_TRANSLATION => 'TAIR_TRANSLATION';
Readonly my @GAF_FIELDS       => qw(
    db                   
    db_object_id         
    db_object_symbol     
    qualifier            
    ont_id               
    db_reference         
    evidence_code        
    with                 
    aspect               
    db_object_name       
    db_object_synonym    
    db_object_type       
    taxon                
    date                 
    assigned_by          
    annotation_extension 
    gene_product_form_id 
);

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

    unless ( ref $files eq 'ARRAY' ) {
        die '"files" argument not an array?';
    }

    my ($ont_file)   = grep { /\.(xml|obo)/ } @$files;
    my ($assoc_file) = grep { $ont_file && !/$ont_file/ } @$files;

    if ( $assoc_file ) {
        $notify->(sprintf "%s parsing file '$assoc_file'\n", __PACKAGE__);
    }
    else {
        printf STDERR "%s called without a 'files' argument\n%s", 
            __PACKAGE__, Dumper($args);
        return 1; # error
    }

    # get the "main" GO/PO source id.
    my $source_name  = $self->get_source_name_for_source_id( $source_id );
    $source_id = $self->get_source_id_for_source_name( $source_name, 'main' );

    #
    # Ontologies are attached to the translations, so this gets all the 
    # translation xrefs as a hash with the translation IDs as the keys
    # and the xref IDs as the values.
    #
    my %master_xrefs = %{ 
        $self->get_valid_codes( $TAIR_TRANSLATION , $species_id ) 
    };

    $species_id       ||= $self->get_species_id_for_filename($assoc_file);
    my %species_id2name = $self->species_id2name();
    my $species         = $species_id2name{$species_id}->[0];

    $notify->("Parsing $source_name for $species\n");

    #
    # Get mappimg from GO terms to descriptions from the XML/OBO file
    #
    #
    my %id_to_term;
    if ( $ont_file ) {
        my $obo_io = $self->get_filehandle( $ont_file );
        if ( !defined $obo_io ) {
            print STDERR "ERROR: Could not open '$ont_file'\n";
            return 1;    # 1 error
        }

        my $term = undef;
        my $desc = undef;

        if ( $ont_file =~ /\.xml$/ ) { 
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
                    if ( defined $term ) {
                        $id_to_term{ $term } = $1;
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

    my $num_terms  = 0;
    my $count      = 0;
    my $miss_parse = 0;

    # Open the gene association GAF file
    my $gaf_io = $self->get_filehandle( $assoc_file );

    if ( !defined $gaf_io ) {
        print STDERR "ERROR: Could not open $assoc_file\n";
        return 1; # error
    }

    while ( my $line = $gaf_io->getline() ) {
        next if $line =~ /^!/;     # comment

        chomp $line;

        my @vals = split( /\t/, $line );
        my %rec  = map { $GAF_FIELDS[$_], $vals[$_] } 0..$#GAF_FIELDS;

        next if $rec{'db'}        ne 'TAIR'; # Only process TAIR annotations
        next if $rec{'qualifier'} eq 'NOT';  # Skip "NOT" terms entirely

        # Find the stable_id(s)
        my @stable_ids = uniq( 
            map  { defined $master_xrefs{ $_ } ? $_ : () }
            grep { /^AT\d+G/ }
            map  { s/\'/\\\'/g; $_ }
            map  { split /[|]/ }
            ( $rec{'db_object_name'}, $rec{'db_object_synonym'} )
        );

        if ( !@stable_ids ) {
            $miss_parse++;
            next;
        }

        my $ont_id = $rec{'ont_id'};
        my $desc   = $id_to_term{ $ont_id } || '';

        for my $stable_id ( @stable_ids ) {
            printf "%-70s\r", sprintf(
                "%10s: %s => %s", ++$num_terms, $ont_id, $stable_id
            );

            my @xref_ids = @{ $master_xrefs{ $stable_id } };

            for my $xref_id ( @xref_ids ) {
                $self->add_dependent_xref({
                    master_xref_id => $xref_id,
                    acc            => $ont_id,
                    label          => $ont_id,
                    desc           => $desc,
                    linkage        => $rec{'evidence_code'},
                    source_id      => $source_id,
                    species_id     => $species_id 
                });
            }
        }

        $count++;
    } 

    $gaf_io->close();

    $notify->(
        map { " - $_\n" }
        "$count $source_name dependent xrefs added",
        "$miss_parse lines did not contain a recognised stable_id",
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
