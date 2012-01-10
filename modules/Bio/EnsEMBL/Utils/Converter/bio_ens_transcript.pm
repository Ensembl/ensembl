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

