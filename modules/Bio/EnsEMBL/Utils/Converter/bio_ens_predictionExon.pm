# Bio::EnsEMBL::Utils::Converter::bio_ens_predictionExon
#
# Created and cared for by Juguang Xiao <juguang@tll.org.sg>
# Created date: 19/3/2003
# 
# Copyright Juguang Xiao
# 
# You may distribute this module under the same terms as perl itself
#
# POD documentation
#

=head1 NAME

Bio::EnsEMBL::Utils::Converter::bio_ens_predictionExon

=head1 SYNOPISIS

my $converter = new Bio::EnsEMBL::Utils::Converter(
    -in => 'Bio::Tools::Prediction::Exon',
    -out => 'Bio::EnsEMBL::Exon'
    -contig => $ens_contig
);

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
