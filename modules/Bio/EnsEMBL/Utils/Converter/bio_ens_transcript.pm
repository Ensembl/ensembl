# Bio::EnsEMBL::Utils::Converter::bio_ens_transcript
#
# Created and cared for by Juguang Xiao <juguang@tll.org.sg>
# Created date: 27/3/2003
# 
# Copyright Juguang Xiao
# 
# You may distribute this module under the same terms as perl itself
#
# POD documentation
#

=head1 NAME

Bio::EnsEMBL::Utils::Converter::bio_ens_transcript - the instance converter

=head1 SYNOPISIS



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

package Bio::EnsEMBL::Utils::Converter::bio_ens_transcript;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Utils::Converter::bio_ens;
@ISA = qw(Bio::EnsEMBL::Utils::Converter::bio_ens);

sub _convert_single {
    my ($self, $arg) = @_;
    unless($arg->isa('Bio::SeqFeature::Gene::Transcript')){
        $self->throw("A Bio::SeqFeature::Gene::Transcript object needed");
    }
    my $transcript = $arg;

    my @exons = $transcript->exons_ordered;
    $self->{_converter_for_exons}->contig($self->contig);
    $self->{_converter_for_exons}->analysis($self->analysis);
    
    my $ens_exons = $self->{_converter_for_exons}->convert(\@exons);
    
    my $ens_transcript = Bio::EnsEMBL::Transcript->new(@{$ens_exons});
    $ens_transcript->start($transcript->start);
    $ens_transcript->end($transcript->end);
#    $ens_transcript->strand($transcript->strand);
    return $ens_transcript;
}


sub _initialize {
    my ($self, @args) = @_;
    $self->SUPER::_initialize(@args);

    $self->{_converter_for_exons} = new Bio::EnsEMBL::Utils::Converter(
        -in => 'Bio::SeqFeature::Gene::Exon',
        -out => 'Bio::EnsEMBL::Exon'
    );

}

