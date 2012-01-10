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

=head1 AUTHOR

Juguang Xiao <juguang@fugu-sg.org>

=cut

=head1 NAME

Bio::EnsEMBL::Utils::Converter::ens_bio_featurePair

=head1 SYNOPISIS

=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Converter::ens_bio_featurePair;

use strict;
use vars qw(@ISA);
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::FeaturePair;
use Bio::EnsEMBL::Utils::Converter::ens_bio;
@ISA = qw(Bio::EnsEMBL::Utils::Converter::ens_bio);

sub _convert_single {
    my ($self, @args) = @_;

}

# convert object from Bio::EnsEMBL::RepeatFeature
# to Bio::SeqFeature::FeaturePair

sub _convert_single_repeatFeature {
    my ($self, $ens_repeat) = @_;
    
    my $feature1 = new Bio::SeqFeature::Generic(
        -start => $ens_repeat->start,
        -end => $ens_repeat->end,
        -strand => $ens_repeat->strand,
        -source_tag => $ens_repeat->source_tag
        -primary_tag => $ens_repeat->repeat_class,
        -seq_id => $ens_repeat->seqname
    );
    
    my ($start2, $end2);
    if($ens_repeat->strand == 1){
        $start2 = $ens_repeat->hstart;
        $end2 = $ens_repeat->hend;
    }elsif($ens_repeat->strand == -1){
        $start2 = $ens_repeat->hend;
        $end2 = $ens_repeat->hstart;
    }else{
        $self->throw("strand cannot be out of range (1, -1)");
    }
        
    my $feature2 = new Bio::SeqFeature::Generic(
        -start => $start2,
        -end => $end2,
        -source_tag => $ens_repeat->source_tag,
        -primary_tag => $ens_repeat->repeat_class,
        -seq_id => $ens_repeat->repeat_name
    );
    
    my $output_module = $self->out;
    require "$output_module";
    return new Bio::SeqFeature::FeaturePair(
        -feature1 => $feature1,
        -feature2 => $feature2
    );
}

1;
