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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 AUTHOR

Juguang Xiao <juguang@tll.org.sg>

=cut

=head1 NAME

Bio::EnsEMBL::Utils::Converter::bio_ens_featurePair

=head1 SYNOPISIS

=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Converter::bio_ens_featurePair;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::RepeatConsensus;
use Bio::EnsEMBL::ProteinFeature;
use Bio::EnsEMBL::Utils::Converter;
use Bio::EnsEMBL::Utils::Converter::bio_ens;
@ISA = qw(Bio::EnsEMBL::Utils::Converter::bio_ens);

sub _initialize {
    my ($self, @args) = @_;
    $self->SUPER::_initialize(@args);
    my ($translation_id) = $self->_rearrange([qw(TRANSLATION_ID)], @args);
    $self->translation_id($translation_id);
    
    # internal converter for seqFeature
    $self->{_bio_ens_seqFeature} = new Bio::EnsEMBL::Utils::Converter (
        -in => 'Bio::SeqFeature::Generic',
        -out => 'Bio::EnsEMBL::SeqFeature',
    );
}

sub _convert_single {
    my ($self, $pair) = @_;
    unless($pair && $pair->isa('Bio::SeqFeature::FeaturePair')){
        $self->throw('a Bio::SeqFeature::FeaturePair object needed');
    }
    
    if($self->out eq 'Bio::EnsEMBL::RepeatFeature'){
        return $self->_convert_single_to_repeatFeature($pair);
    }elsif($self->out eq 'Bio::EnsEMBL::FeaturePair'){
        return $self->_convert_single_to_featurePair($pair);
    }elsif($self->out eq 'Bio::EnsEMBL::ProteinFeature'){
        return $self->_convert_single_to_proteinFeature($pair);
    }else{
        my $output_module = $self->out;
        $self->throw("Cannot covert to [$output_module]");
    }
}

sub _convert_single_to_featurePair {
    my ($self, $pair) = @_;
    my $feature1 = $pair->feature1;
    my $feature2 = $pair->feature2;
    $self->{_bio_ens_seqFeature}->contig($self->contig);
    $self->{_bio_ens_seqFeature}->analysis($self->analysis);
    my $ens_f1 = $self->{_bio_ens_seqFeature}->_convert_single($feature1);
    my $ens_f2 = $self->{_bio_ens_seqFeature}->_convert_single($feature2);
    my $ens_fp = Bio::EnsEMBL::FeaturePair->new(
        -feature1 => $ens_f1,
        -feature2 => $ens_f2
    );
    return $ens_fp;
}

sub _convert_single_to_proteinFeature {
    my ($self, $pair) = @_;
    my $featurePair = $self->_convert_single_to_featurePair($pair);
    my $proteinFeature = Bio::EnsEMBL::ProteinFeature->new(
        -feature1 => $featurePair->feature1,
        -feature2 => $featurePair->feature2
    );
    $proteinFeature->seqname($self->translation_id);
    return $proteinFeature;
}

sub _convert_single_to_repeatFeature {
    my ($self, $pair) = @_;
    my $feature1 = $pair->feature1;
    my $feature2 = $pair->feature2;
    my $ens_repeatfeature = new Bio::EnsEMBL::RepeatFeature(
        -seqname => $feature1->seq_id,
        -start => $feature1->start,
        -end => $feature1->end,
        -strand => $feature1->strand,
        -source_tag => $feature1->source_tag,
    );
    
    my ($h_start, $h_end);
    if($feature1->strand == 1){
        $h_start = $feature2->start;
        $h_end = $feature2->end;
    }elsif($feature1->strand == -1){
        $h_start = $feature2->end;
        $h_end = $feature2->start;
    }else{
        $self->throw("strand cannot be outside of (1, -1)");
    }

    $ens_repeatfeature->hstart($h_start);
    $ens_repeatfeature->hend($h_end);
    my $repeat_name = $feature2->seq_id;
    my $repeat_class = $feature1->primary_tag;
    $repeat_class ||= $feature2->primary_tag;
    $repeat_class ||= "not sure";
    my $ens_repeat_consensus = 
        $self->_create_consensus($repeat_name, $repeat_class);
    $ens_repeatfeature->repeat_consensus($ens_repeat_consensus);
   
    my($contig) = ref ($self->contig) eq 'ARRAY' ? @{$self->contig} : $self->contig;

    $ens_repeatfeature->attach_seq($contig);
    $ens_repeatfeature->analysis($self->analysis);
    return $ens_repeatfeature;
}

sub _create_consensus{
    my ($self, $repeat_name, $repeat_class) = @_;
    my $consensus = new Bio::EnsEMBL::RepeatConsensus;
    $consensus->name($repeat_name);
    $consensus->repeat_class($repeat_class);
    return $consensus;
}

1;
