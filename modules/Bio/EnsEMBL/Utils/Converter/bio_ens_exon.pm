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

Bio::EnsEMBL::Utils::Converter::bio_ens_exon

=head1 SYNOPISIS

=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Converter::bio_ens_exon;

use strict;
use vars qw(@ISA %GTF_ENS_PHASE);
use Bio::EnsEMBL::Utils::Converter;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::Utils::Converter::bio_ens;
@ISA = qw(Bio::EnsEMBL::Utils::Converter::bio_ens);

BEGIN {
    %GTF_ENS_PHASE = (
        0  =>  0,
        1  =>  2,
        2  =>  1,
       '.' => -1
    );
}

sub _initialize {
    my ($self, @args) = @_;
    $self->SUPER::_initialize(@args);

    $self->{_bio_ens_seqFeature} = new Bio::EnsEMBL::Utils::Converter (
        -in => 'Bio::SeqFeature::Generic',
        -out => 'Bio::EnsEMBL::SeqFeature',
    );

    $self->{_bio_ens_featurePair} = new Bio::EnsEMBL::Utils::Converter (
        -in => 'Bio::SeqFeature::FeaturePair',
        -out => 'Bio::EnsEMBL::FeaturePair'
    );
}

sub _attach_supporting_feature {
    my ($self, $exon, $ens_exon) = @_;
    unless($exon->has_tag('supporting_feature')){
        return;
    }
    my ($sf) = $exon->each_tag_value('supporting_feature');
    unless(defined $sf){
        $self->warn("no supporting feature is attached in exon");
        return;
    }
#    $self->{_bio_ens_seqFeature}->contig($self->contig);
    $self->{_bio_ens_seqFeature}->analysis($self->analysis);
    my $ens_f1 = $self->{_bio_ens_seqFeature}->_convert_single($sf->feature1);
    my $ens_f2 = $self->{_bio_ens_seqFeature}->_convert_single($sf->feature2);
    $self->{_bio_ens_featurePair}->contig($self->contig);
    $self->{_bio_ens_featurePair}->analysis($self->analysis);
    my $ens_sf = $self->{_bio_ens_featurePair}->_convert_single($sf);
    my @align_feautre_args = (
        -feature1 => $ens_f1,
        -feature2 => $ens_f2,
        -features => [$ens_sf]
    );
    
    my $ens_supporting_feature = 
        Bio::EnsEMBL::DnaPepAlignFeature->new(@align_feautre_args);
    
    $ens_exon->add_supporting_features($ens_supporting_feature);
}

sub _convert_single {
    my ($self, $arg) = @_;
    unless($arg && $arg->isa('Bio::SeqFeature::Gene::Exon')){
        $self->throw("a Bio::SeqFeature::Gene::Exon object needed");
    }
    my $exon = $arg;
    my $ens_exon = Bio::EnsEMBL::Exon->new_fast(
        $self->contig, $exon->start, $exon->end, $exon->strand);

    my ($phase) = $exon->each_tag_value('phase');
    $ens_exon->phase($GTF_ENS_PHASE{$phase});
    my $ens_end_phase = 3 - ($exon->length - $phase) % 3;
    $ens_end_phase = 0 if $ens_end_phase == 3;
    $ens_exon->end_phase($ens_end_phase);
    if($self->contig->isa('Bio::EnsEMBL::RawContig')){
        $ens_exon->sticky_rank(1);
    }
    $self->_attach_supporting_feature($exon, $ens_exon);
    return $ens_exon;
}


1;
