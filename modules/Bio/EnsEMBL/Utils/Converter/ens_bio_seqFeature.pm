# Bio::EnsEMBL::Utils::Converter::ens_bio_seqFeature
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

Bio::EnsEMBL::Utils::Converter::ens_bio_seqFeature

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

package Bio::EnsEMBL::Utils::Converter::ens_bio_seqFeature;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Utils::Converter::ens_bio;
@ISA = qw(Bio::EnsEMBL::Utils::Converter::ens_bio);

sub _convert_single {
    my ($self, $in) = @_;
    
    $self->throw("Input not defined") unless($in && defined($in));
    unless(ref($in) && $in->isa('Bio::EnsEMBL::SeqFeature')){
        $self->throw('A Bio::EnsEMBL::SeqFeature object needed');
    }

    my @args = (
        -start => $in->start,
        -end => $in->end,
        -strand => $in->strand,
        -score => $in->score,
        -source_tag => $in->source_tag,
        -seq_id => $in->seqname
    );
    
    my $seqFeature = new Bio::SeqFeature::Generic(@args);
    
    return $seqFeature;
}

1;
