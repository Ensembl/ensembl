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

=pod


=head1 NAME

Report Serializer - generating textual summary reports

=head1 SYNOPSIS

    use Bio::EnsEMBL::Registry;
    use Bio::EnsEMBL::Utils::IO::ReportSerializer;
    use IO::File;
    
    my $registry = 'Bio::EnsEMBL::Registry';
    $output_fh = IO::File->new($config{'output'},'w') or die;
    $serializer = new ReportSerializer($output_fh);
    my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );
    my $slice = $slice_adaptor->fetch_by_toplevel_location("6:1000000..1500000");
    
    $serializer->print_section_header($slice);
    $serializer->print_feature_list($slice->get_all_Genes);

=head1 DESCRIPTION

Subclass of Serializer that can turn a feature into a text block
Unsuited to very large slices, because it requires a select-all approach for features.

=cut

package Bio::EnsEMBL::Utils::IO::ReportSerializer;
use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception;
use URI::Escape;
use Bio::EnsEMBL::Utils::IO::FeatureSerializer;

use base qw(Bio::EnsEMBL::Utils::IO::FeatureSerializer);

my %strand_conversion = ( '1' => '+', '0' => '?', '-1' => '-');

my %feature_conversion = (     'Bio::EnsEMBL::Gene' => 'Gene',
                            'Bio::EnsEMBL::Transcript' => 'Transcript',
                            'Bio::EnsEMBL::Translation' => 'Translation',
                            'Bio::EnsEMBL::Variation::StructuralVariationFeature' => 'Structural Variation',
                            'Bio::EnsEMBL::Variation::VariationFeature' => 'Variation',
                            'Bio::EnsEMBL::Funcgen::RegulatoryFeature' => 'Regulatory Feature',
                            'Bio::EnsEMBL::Compara::ConstrainedElement' => 'Constrained Element',
                            'Feature' => 'Feature',
);

# Hash for selecting the correct attributes of unseen features for crude summary. This hash is 
# for fallback behaviour, slicing summary hashes for a limited set of values.
my %printables = ( 
                    'Bio::EnsEMBL::Gene' => ['id','biotype','start','end'],
                    'Bio::EnsEMBL::Transcript' => ['id','start','end'],
                    'Bio::EnsEMBL::Translation' => ['id'],
                    'Bio::EnsEMBL::Slice' => ['start','end','strand'],
                    'Bio::EnsEMBL::Variation::VariationFeature' => ['id','start','end','strand','seq_region_name'],
                    'Bio::EnsEMBL::Variation::StructuralVariationFeature' => ['id','start','end','strand','seq_region_name'],
                    'Bio::EnsEMBL::Funcgen::RegulatoryFeature' => ['id','start','end','strand'],
                    'Bio::EnsEMBL::Compara::ConstrainedElement' => ['id','start','end','strand','seq_region_name'],
                );

=head2 print_feature

    Arg [1]    : Bio::EnsEMBL::Feature, subclass or related pseudo-feature
    Example    : $reporter->print_feature($feature,$slice_start_coordinate,"X")
=cut

sub print_feature {
    my $self = shift;
    my $feature = shift;
    my $fh = $self->{'filehandle'};
    my $feature_type = ref($feature);
    
    if ($feature->can('summary_as_hash') ) {
        my %summary = %{$feature->summary_as_hash};
        my @values = @summary{ @{$printables{$feature_type}} };
        print $fh join(',',@values)."\n";
    }
    else {
        warning("Feature failed to self-summarise");
    }
}

=head2 print_feature_list

    Arg [1]    : Listref of Bio::EnsEMBL::Feature, subclass or related pseudo-feature
    Description: Relies on a list of similar features to print in a block together.
                 Overrides superclass method
                 Results are truncated after the first 100 features for brevity.
    Example    : $reporter->print_feature_list(\@features);
=cut

sub print_feature_list {
    my $self = shift;
    my $feature_list = shift;
    if (scalar(@$feature_list) > 0) {$self->{'achieved_something'} = 1;} #from superclass
    my $fh = $self->{'filehandle'};

    my $example_feature = $feature_list->[0];
    my $feature_type = ref($example_feature);
    my $feature_count = 0;
    unless (defined $feature_type) {$feature_type = "Feature"};
    if (scalar(@$feature_list) > 0) {
        print $fh "There are ",scalar(@$feature_list)," ",$feature_conversion{$feature_type},(scalar(@$feature_list) != 1) ? "s":""," in this region\n";
    }
    if (scalar(@$feature_list) > 100 ) { print $fh "Too many to display, results truncated to the first 100\n";}
    print $fh "\n";
    foreach my $feature (@$feature_list) {
        $feature_count++;
        my %attributes = %{$feature->summary_as_hash};
        
        if ($feature_count == 100) {last;}
        # Begin the feature-specific formatting code
        if ($feature_type eq "Bio::EnsEMBL::Gene") {
            print $fh "\tGene ".$feature_count.": ".$attributes{'external_name'}.",".$attributes{'id'}."\n";
            print $fh "\tBiotype: ".$attributes{'biotype'}."\n";
            print $fh "\tLocation: ".$attributes{'start'}."-".$attributes{'end'}." bp\n\n";
            
            print $fh "\tTranscripts and proteins\n";
            foreach my $transcript (@{$feature->get_all_Transcripts}) {
                my %tr_summary = %{$transcript->summary_as_hash};
                print $fh "\t\t ".$tr_summary{'id'};
                my $translation = $transcript->translation;
                if (defined $translation) {
                    my %pr_summary = %{$translation->summary_as_hash};
                    print $fh " - ".$pr_summary{'id'}."\n\n";
                }
                else {print $fh " - no protein\n\n";}
            }
            print $fh "\n";
        }
        elsif ($feature_type eq "Bio::EnsEMBL::Funcgen::RegulatoryFeature") {
            print $fh "\t".$attributes{'id'}."\n";
        }
        elsif ($feature_type eq "Bio::EnsEMBL::Compara::ConstrainedElement") {
            print $fh "\t".$attributes{'start'}."-".$attributes{'end'}."\n";
        } 
        elsif ( $feature_type eq "Bio::EnsEMBL::Variation::StructuralVariationFeature" 
            or $feature_type eq "Bio::EnsEMBL::Variation::VariationFeature") {
            print $fh "\tID: ".$attributes{'id'}."  Position: ".
                $attributes{'start'}."-".$attributes{'end'}." on strand ".$attributes{'strand'}." \n";
        }
        else {
            # slice favourite values out unformatted.
            my @values = @attributes{ @{$printables{$feature_type}} };
            print $fh $feature_type.join(',',@values)."\n";
            
        }
    }
}

# Just print individuals without awareness of list size and position.
sub print_feature_iterator {
    my $self = shift;
    my $feature_iterator = shift;
    while ($feature_iterator->has_next) {
        my $feature = $feature_iterator->next;
        $self->print_feature($feature);
    }
    $self->{'achieved_something'} = 1;
}

=head2 print_main_header

    Arg [1]    : Arrayref of slices going into the file.
    Description: Printing the header text for this report
                 Requires a slice list in order to report how many will be printed
    Returntype : None
=cut

sub print_main_header {
    my $self = shift;
    my $arrayref_of_slices = shift;
    my $fh = $self->{'filehandle'};

    my $regions = scalar @{$arrayref_of_slices};
    print $fh "Report for $regions region";
    if ($regions > 1) { print $fh "s";}
    print $fh "\n\n";
}


=head2 print_section_header 

    Arg [1]    : Bio::EnsEMBL::Slice
    Description: Prints a summary of the slice
                 Intended to be used prior to print_feature_list()
    Returntype : None

=cut

sub print_section_header {
    my $self = shift;
    my $slice = shift;
    my $fh = $self->{'filehandle'};

    print $fh "  Region: ",$slice->seq_region_name," ",$slice->start,"-",$slice->end," bp\n\n";

}

