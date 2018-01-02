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

Bio::EnsEMBL::Utils::Converter::bio_ens_transcript - the instance converter

=head1 SYNOPISIS

=head1 DESCRIPTION

=head1 METHODS

=cut

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

