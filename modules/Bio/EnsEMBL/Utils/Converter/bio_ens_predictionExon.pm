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
