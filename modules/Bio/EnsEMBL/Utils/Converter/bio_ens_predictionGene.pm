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
