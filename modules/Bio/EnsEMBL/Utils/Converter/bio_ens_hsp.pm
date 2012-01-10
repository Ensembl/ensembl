=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut


package Bio::EnsEMBL::Utils::Converter::bio_ens_hsp;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Utils::Converter::bio_ens;
use Bio::EnsEMBL::ProteinFeature;

@ISA = qw(Bio::EnsEMBL::Utils::Converter::bio_ens);

sub _initialize {
    my ($self, @args) = @_;
    $self->SUPER::_initialize(@args);

    # After super initialized, analysis and contig are ready.
    my $bio_ens_seqFeature_converter = new Bio::EnsEMBL::Utils::Converter(
        -in => 'Bio::SeqFeature::Generic',
        -out => 'Bio::EnsEMBL::SeqFeature',
        -analysis => $self->analysis,
        -contig => $self->contig
    );
    $self->_bio_ens_seqFeature_converter($bio_ens_seqFeature_converter);

}

sub _convert_single {
    my ($self, $hsp) = @_;

    unless(ref($hsp) && $hsp->isa('Bio::Search::HSP::GenericHSP')){
        $self->throw("a GenericHSP object needed");
    }
    
    my $in = $self->in;
    my $out = $self->out;

    if($out =~ /^Bio::EnsEMBL::ProteinFeature$/){
        return $self->_convert_single_to_proteinFeature($hsp);
    }elsif($out =~/^Bio::EnsEMBL::(DnaDna|DnaPep|PepDna)AlignFeature/){
        return $self->_convert_single_to_alignFeature($hsp);
    }else{
        $self->throw("[$in]->[$out], not implemented");
    }
}

sub _convert_single_to_featurePair {
    my ($self, $hsp) = @_;
    
    my $bio_ens_seqFeature_converter = $self->_bio_ens_seqFeature_converter;
    my $ens_feature1 = $bio_ens_seqFeature_converter->_convert_single(
        $hsp->feature1);
    my $ens_feature2 = $bio_ens_seqFeature_converter->_convert_single(
        $hsp->feature2);

    $ens_feature1->p_value($hsp->evalue);
    $ens_feature1->score($hsp->score);
    $ens_feature1->percent_id($hsp->percent_identity);
    $ens_feature2->p_value($hsp->evalue);
    $ens_feature2->score($hsp->score);
    $ens_feature2->percent_id($hsp->percent_identity);

    my $featurePair = Bio::EnsEMBL::FeaturePair->new(
        -feature1 => $ens_feature1,
        -feature2 => $ens_feature2
    );

    return $featurePair;
}

sub _convert_single_to_proteinFeature {
    my ($self, $hsp) = @_;
    
    my $ens_featurePair = $self->_convert_single_to_featurePair($hsp);
    my $ens_proteinFeature = Bio::EnsEMBL::ProteinFeature->new(
        -feature1 => $ens_featurePair->feature1,
        -feature2 => $ens_featurePair->feature2
    );
    $ens_proteinFeature->seqname($self->translation_id);
    return $ens_proteinFeature;
}

sub _convert_single_to_alignFeature {
    my ($self, $hsp) = @_;
    my $ens_featurePair = $self->_convert_single_to_featurePair($hsp);
    my $cigar_string = $hsp->cigar_string;
    my @args = (
        -feature1 => $ens_featurePair->feature1,
        -feature2 => $ens_featurePair->feature2,
        -cigar_string => $hsp->cigar_string
    );
    my $contig = $self->contig;
    # choose the AlignFeature based on the blast program
    my $program = $hsp->algorithm;

    $self->throw("HSP does not have algorithm value") unless(defined($program));
    my $align_feature;
    if($program =~ /blastn/i){
        $align_feature = new Bio::EnsEMBL::DnaDnaAlignFeature(@args);
#        $align_feature->attach_seq($contig);
    }elsif($program =~ /blastx/i){
        $align_feature = new Bio::EnsEMBL::DnaPepAlignFeature(@args);
#        $align_feature->attach_seq($contig);
    }else{
        $self->throw("\[$program\] is not supported yet");
    }
    return $align_feature;
}

sub _bio_ens_seqFeature_converter {
    my ($self) = shift ;
    return $self->{_bio_ens_seqFeature_converter} = shift if(@_);
    return $self->{_bio_ens_seqFeature_converter};
}
1;
