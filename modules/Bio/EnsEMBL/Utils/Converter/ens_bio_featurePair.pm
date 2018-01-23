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
    eval "require $output_module"; ## no critic
    return new Bio::SeqFeature::FeaturePair(
        -feature1 => $feature1,
        -feature2 => $feature2
    );
}

1;
