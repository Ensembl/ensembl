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

Bio::EnsEMBL::Utils::Converter::bio_ens_gene

=head1 SYNOPISIS

=head1 DESCRIPTION

This module is to convert from objects of
Bio::SeqFeature::Gene::GeneStructure to those of Bio::EnsEMBL::Gene

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::Converter::bio_ens_gene;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Utils::Converter::bio_ens;
@ISA = qw(Bio::EnsEMBL::Utils::Converter::bio_ens);

sub _initialize {
    my ($self, @args) = @_;
    $self->SUPER::_initialize(@args);
    
    $self->{_converter_for_transcripts} = new Bio::EnsEMBL::Utils::Converter(
        -in => 'Bio::SeqFeature::Gene::Transcript',
        -out => 'Bio::EnsEMBL::Transcript'
    );

}

sub _convert_single {
    my ($self, $input) = @_;

    unless($input->isa('Bio::SeqFeature::Gene::GeneStructure')){
        $self->throw("a Bio::SeqFeature::Gene::GeneStructure object needed");
    }
    my $gene = $input;
    my $ens_gene = Bio::EnsEMBL::Gene->new();
    $ens_gene->analysis($self->analysis);
    my @transcripts = $gene->transcripts;
    
    # contig is needed by exon and Supporting Feature; S.F. needs an analysis.
    $self->{_converter_for_transcripts}->contig($self->contig);
    $self->{_converter_for_transcripts}->analysis($self->analysis);
    
    my $ens_transcripts = $self->{_converter_for_transcripts}->convert(
        \@transcripts);
    
    foreach(@{$ens_transcripts}){
        $ens_gene->add_Transcript($_);
    }
    return $ens_gene;
}

1;
