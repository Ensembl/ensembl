# Bio::EnsEMBL::Utils::Converter::bio_ens_featurePair
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

Bio::EnsEMBL::Utils::Converter::bio_ens_featurePair

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

package Bio::EnsEMBL::Utils::Converter::bio_ens_featurePair;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::RepeatConsensus;

use Bio::EnsEMBL::Utils::Converter::bio_ens;
@ISA = qw(Bio::EnsEMBL::Utils::Converter::bio_ens);

sub _convert_single {
    my ($self, $pair) = @_;
    if($self->out eq 'Bio::EnsEMBL::RepeatFeature'){
        return $self->_convert_single_to_repeatFeature($pair);
    }else{
        my $output_module = $self->out;
        $self->throw("Cannot covert to [$output_module]");
    }
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
