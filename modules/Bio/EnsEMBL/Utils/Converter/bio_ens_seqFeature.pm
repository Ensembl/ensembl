# Bio::EnsEMBL::Utils::Converter::bio_ens_seqFeature
#
# Created and cared for by Juguang Xiao <juguang@tll.org.sg>
# Created date: 5/3/2003
# 
# Copyright Juguang Xiao
# 
# You may distribute this module under the same terms as perl itself
#
# POD documentation
#

=head1 NAME

Bio::EnsEMBL::Utils::Converter::bio_ens_seqFeature

=head1 SYNOPISIS

Please read Bio::EnsEMBL::Utils::Converter

=head1 DESCRIPTION


=head1 FEEDBACK

=head2 Mailing Lists

=head2 Reporting Bugs


=head1 AUTHOR Juguang Xiao

Juguang Xiao <juguang@tll.org.sg>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin ...

package Bio::EnsEMBL::Utils::Converter::bio_ens_seqFeature;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Utils::Converter::bio_ens;
@ISA = qw(Bio::EnsEMBL::Utils::Converter::bio_ens);

sub _convert_single {
    my ($self, $in) = @_;
    
    unless($in && defined($in) && ref($in) && $in->isa('Bio::SeqFeature::Generic')){
        $self->throw("a Bio::SeqFeature::Generic object needed");
    }
    
    my $seqFeature = $in;
    
    my $ens_seqFeature;
    my @args = (
        -start => $in->start,
        -end => $in->end,
        -strand => $in->strand,
        -score => $in->score
        -analysis => $self->analysis,
        -source_tag => $in->source_tag
        -seqname => $in->seq_id
    );

    my $output_module = $self->out;
    
    if($output_module eq 'Bio::EnsEMBL::SeqFeature'){
        
        $ens_seqFeature = new Bio::EnsEMBL::SeqFeature(@args);
    }elsif($self->out eq 'Bio::EnsEMBL::SimpleFeature'){
        $ens_seqFeature = new Bio::EnsEMBL::SimpleFeature(@args);
        # The field that there is in SimpleFeature, but not in SeqFeature.
        $ens_seqFeature->display_label('__NONE__');
    }elsif($self->out eq 'Bio::EnsEMBL::Exon'){
        $ens_seqFeature = Bio::EnsEMBL::Exon->new_fast(
            $self->contig, $seqFeature->start, $seqFeature->end, 
            $seqFeature->strand);
        
    }else{
        $self->throw("[$output_module] as -out, not supported");
    }
    
    $ens_seqFeature->attach_seq($self->contig);

    return $ens_seqFeature;
}

1;
