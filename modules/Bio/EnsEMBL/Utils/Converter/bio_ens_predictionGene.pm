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

Juguang Xiao <juguang@tll.org.sg>

=cut

=head1 NAME

Bio::EnsEMBL::Utils::Converter::bio_ens_predictionGene

=head1 SYNOPISIS

=head1 DESCRIPTION

=head1 METHODS

=cut

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
