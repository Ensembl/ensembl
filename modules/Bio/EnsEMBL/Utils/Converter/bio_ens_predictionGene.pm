# Bio::EnsEMBL::Utils::Converter::bio_ens_predictionGene
#
# Created and cared for by Juguang Xiao <juguang@tll.org.sg>
# Created date: 18/3/2003
# 
# Copyright Juguang Xiao
# 
# You may distribute this module under the same terms as perl itself
#
# POD documentation
#

=head1 NAME

Bio::EnsEMBL::Utils::Converter::bio_ens_predictionGene

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

package Bio::EnsEMBL::Utils::Converter::bio_ens_predictionGene;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Utils::Converter;
use Bio::EnsEMBL::PredictionTranscript;
use Bio::EnsEMBL::Utils::Converter::bio_ens;
@ISA = qw(Bio::EnsEMBL::Utils::Converter::bio_ens);

sub _initialize {
    my ($self, @args) = @_;

    $self->{_predictionExonConverter} = new Bio::EnsEMBL::Utils::Converter(
        -in => 'Bio::SeqFeature::Gene::Exon',
        -out => 'Bio::EnsEMBL::Exon',
    );
    $self->SUPER::_initialize(@args);
    
    $self->{_predictionExonConverter}->contig($self->contig);
    $self->{_predictionExonConverter}->analysis($self->analysis);
}

sub _convert_single {
    my ($self, $input) = @_;
    
    $self->throw("one argument needed") unless($input and defined($input));
    $self->throw("a Bio::Tools::Prediction::Gene object needed")
        unless(ref($input) && $input->isa('Bio::Tools::Prediction::Gene'));

    my $output = Bio::EnsEMBL::PredictionTranscript->new;
    $output->analysis($self->analysis);

    my @exons = sort {$a->start <=> $b->start} $input->exons;

    # Not sure on the correctivity of phase calculation.
    my $previous_end_phase = -1;
    foreach(@exons){
        my $length = $_->length;
        my $frame = $_->frame;
        my $phase = ($previous_end_phase+1) %3;
        my $end_phase = ($length-$frame)  %1;
        $previous_end_phase = $end_phase;
        $_->add_tag_value("phase", $phase);
        $_->add_tag_value("end_phase", $end_phase);
    }

    my @ens_exons = @{$self->{_predictionExonConverter}->convert(\@exons)};

    $output->add_Exon($_) foreach(@ens_exons);

    return $output;

}

=head2 contig
  Title   : contig
  Usage   : $self->contig
  Function: get and set for contig
  Return  : L<Bio::EnsEMBL::RawContig>
  Args    : L<Bio::EnsEMBL::RawContig>
=cut

sub contig {
    my ($self, $arg) = @_;
    if(defined($arg)){
        $self->throws("A Bio::EnsEMBL::RawContig object expected.") unless(defined $arg);
        $self->{_contig} = $arg;
        # assign it to the sub converter which converts exons
        $self->{_predictionExonConverter}->contig($arg);
    }
    return $self->{_contig};
}
1;
