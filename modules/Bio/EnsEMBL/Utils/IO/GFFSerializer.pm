=pod

=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 NAME

GFFSerializer - Feature to GFF converter

=head1 SYNOPSIS

use Bio::EnsEMBL::Utils::IO::GFFSerializer;
use Bio::EnsEMBL::Utils::BiotypeMapper;

my $ontology_adaptor = $registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );
my $serializer = Bio::EnsEMBL::Utils::IO::GFFSerializer->new($ontology_adaptor,$output_fh);

my $variation_feature_adaptor = $registry->get_adaptor( $config{'species'}, 'variation', 'variationfeature' );
$serializer->print_metadata("Variation Features:");
my $iterator = $variation_feature_adaptor->fetch_Iterator_by_Slice($slice,undef,60000);
$serializer->print_feature_Iterator($iterator);

=head1 DESCRIPTION

Subclass of Serializer that can turn a feature into a line for the GFF3 format. Requires
a BiotypeMapper in order to translate biotypes to SO terms.

=cut

package Bio::EnsEMBL::Utils::IO::GFFSerializer;
use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Utils::BiotypeMapper;
use URI::Escape;
use Bio::EnsEMBL::Utils::IO::FeatureSerializer;

use base qw(Bio::EnsEMBL::Utils::IO::FeatureSerializer);

my %strand_conversion = ( '1' => '+', '0' => '?', '-1' => '-');

=head2 new

    Constructor
    Arg [1]    : Ontology Adaptor
    Arg [2]    : Optional File handle
    
    Returntype : Bio::EnsEMBL::Utils::IO::GFFSerializer

=cut

sub new {
    my $class = shift;
    my $self = {
        ontology_adaptor => shift,
        filehandle => shift,
    };
    bless $self, $class;
    if ( ref($self->{'ontology_adaptor'}) ne "Bio::EnsEMBL::DBSQL::OntologyTermAdaptor" ) {
        throw("GFF format requires an instance of Bio::EnsEMBL::DBSQL::OntologyTermAdaptor to function. See also Bio::EnsEMBL::Utils::BiotypeMapper");        
    }
    $self->{'mapper'} = new Bio::EnsEMBL::Utils::BiotypeMapper($self->{'ontology_adaptor'});
    
    if (!defined ($self->{'filehandle'})) {
        # no file handle, let the handle point to a copy of STDOUT instead
        open $self->{'filehandle'}, ">&STDOUT";
        $self->{'stdout'} = 1;
    }
    return $self;
}

=head2 print_feature

    Arg [1]    : Bio::EnsEMBL::Feature, subclass or related pseudo-feature
    Example    : $reporter->print_feature($feature,$slice_start_coordinate,"X")
    Description: Asks a feature for its summary, and generates a GFF3 
                 compliant entry to hand back again
                 Additional attributes are handed through to column 9 of the 
                 output using exact spelling and capitalisation of the 
                 feature-supplied hash.
    Returntype : none
=cut

sub print_feature {
    my $self = shift;
    my $feature = shift;
    my $biotype_mapper = $self->{'mapper'};

    my $text_buffer = "";
    if ($feature->can('summary_as_hash') ) {
        my %summary = %{$feature->summary_as_hash};
        my $row = "";
#    Column 1 - seqname, the name of the sequence/chromosome the feature is on. Landmark for start below
        if (!defined($summary{'seq_region_name'})) {$summary{'seq_region_name'} = "?";}
        $row .= $summary{'seq_region_name'}."\t";

#    Column 2 - source, complicated with Ensembl not being the originator of all data
        $row .= "EnsEMBL\t";

#   Column 3 - feature, the ontology term for the kind of feature this row is
        my $so_term = $biotype_mapper->translate_feature_to_SO_term($feature);
        $row .= $so_term."\t";

#    Column 4 - start, the start coordinate of the feature, here shifted to chromosomal coordinates
#    Start and end must be in ascending order for GFF. Circular genomes require the length of 
#   the circuit to be added on.
        if ($summary{'start'} > $summary{'end'}) {
            #assumes this is not a Compara circular sequence and can treat is as a Feature
            if ($feature->slice() && $feature->slice()->is_circular() ) {
                $summary{'end'} = $summary{'end'} + $feature->seq_region_length;
            }
            # non-circular, but end still before start
            else {$summary{'end'} = $summary{'start'};}
        }
        $row .= $summary{'start'} . "\t";

#    Column 5 - end, coordinates (absolute) for the end of this feature
        $row .= $summary{'end'} . "\t";

#    Column 6 - score, for variations only.
        if (exists($summary{'score'})) {
            $row .= $summary{'score'}."\t";
        }
        else {
            $row .= ".\t";
        }

#    Column 7 - strand, up or down
        if (exists($summary{'strand'})) {
            $row .= $strand_conversion{$summary{'strand'}}."\t";
        }
        else {
            $row .= ".\t";
        }

#   Column 8 - reading phase, necessary only for Exons
        if (exists($summary{'phase'})) {
          $row .= $summary{'phase'}."\t";
        }
        else {
          $row .= ".\t";
        }

#    Column 9 - the 'other' section for all GFF and GVF compliant attributes
#    We include Stable ID and biotype where possible to supplement the information in the other columns
        delete $summary{'seq_region_start'};
        delete $summary{'seq_region_name'};
        delete $summary{'start'};
        delete $summary{'end'};
        delete $summary{'strand'};
        delete $summary{'score'};
#   Slice the hash for specific keys in GFF-friendly order
        my @ordered_keys = qw(ID Name Alias Parent Target Gap Derives_from Note Dbxref Ontology_term Is_circular);
        my @ordered_values = @summary{@ordered_keys};
        while (my $key = shift @ordered_keys) {
            my $value = shift @ordered_values;
            if ($value) {
                $row .= $key."=".uri_escape($value,'\t\n\r;=%&,').";";
            }
            delete $summary{$key};
        }
#   Catch the remaining keys, containing whatever else the Feature provided
        my @keys = keys %summary;
        while(my $attribute = shift @keys) {
            if (ref $summary{$attribute} eq "ARRAY") {
                $row .= $attribute."=".join (',',map { uri_escape($_,'\t\n\r;=%&,') } grep { defined $_ } @{$summary{$attribute}});
            }
            else {
                if ($summary{$attribute}) { 
                  $row .= $attribute."=".uri_escape($summary{$attribute},'\t\n\r;=%&,'); 
                }
            }
            $row .= ';' if scalar(@keys) > 0;
        }
# trim off any trailing commas left by the ordered keys stage above:
        $text_buffer .= $row."\n";
    }
    else {
        warning("Feature failed to self-summarise");
    }
    #filehandle is inherited
    my $fh = $self->{'filehandle'};
    print $fh $text_buffer;
}

=head2 print_main_header

    Arg [1]    : Arrayref of slices going into the file.
    Description: Printing the header text or metadata required for GFF,
                 using a list of slices to be written
    Returntype : None
=cut

sub print_main_header {
    my $self = shift;
    my $arrayref_of_slices = shift;
    my $fh = $self->{'filehandle'};
    
    print $fh "##gff-version 3\n";
    foreach my $slice (@{$arrayref_of_slices}) {
        if (not defined($slice)) { warning("Slice not defined.\n"); return;}
        print $fh "##sequence-region   ",$slice->seq_region_name," ",$slice->start," ",$slice->end,"\n";
    }
}

sub print_metadata {
    my $self = shift;
    my $text = shift;
    my $fh = $self->{'filehandle'};
    print $fh "\n#".$text."\n";
}


1;
