# Bio::EnsEMBL::Utils::Converter::ens_bio_featurePair
#
# Created and cared for by Juguang Xiao <juguang@fugu-sg.org>
# Created date: 5/3/2003
# 
# Copyright Juguang Xiao
# 
# You may distribute this module under the same terms as perl itself
#
# POD documentation
#

=head1 NAME

Bio::EnsEMBL::Utils::Converter::ens_bio_featurePair

=head1 SYNOPISIS



=head1 DESCRIPTION


=head1 FEEDBACK

=head2 Mailing Lists

=head2 Reporting Bugs


=head1 AUTHOR Juguang Xiao

Juguang Xiao <juguang@fugu-sg.org>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin ...

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
