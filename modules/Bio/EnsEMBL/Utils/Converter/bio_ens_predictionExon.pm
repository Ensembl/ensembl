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

Bio::EnsEMBL::Utils::Converter::bio_ens_predictionExon

=head1 SYNOPISIS

  my $converter = new Bio::EnsEMBL::Utils::Converter(
    -in     => 'Bio::Tools::Prediction::Exon',
    -out    => 'Bio::EnsEMBL::Exon',
    -contig => $ens_contig
  );

=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Converter::bio_ens_predictionExon;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Utils::Converter;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Utils::Converter::bio_ens;
@ISA = qw(Bio::EnsEMBL::Utils::Converter::bio_ens);

sub _convert_single {
    my ($self, $input) = @_;
    

    $input || $self->throw("a input object needed");
    $self->throw("a Bio::Tools::Prediction::Exon object needed")
        unless($input->isa("Bio::Tools::Prediction::Exon"));

    my $output = Bio::EnsEMBL::Exon->new(
        -start => $input->start,
        -end => $input->end,
        -strand => $input->strand
    );

    $output->score($input->score);
    $output->p_value($input->significance);

    $output->phase($input->get_tag_values("phase")); # only first element is used
    $output->end_phase($input->get_tag_values("end_phase"));

    $output->contig($self->contig);

    return $output;
}


1;
