# Bio::EnsEMBL::Utils::Converter::bio_ens_gene
#
# Created and cared for by Juguang Xiao <juguang@tll.org.sg>
# Created date: 26/3/2003
# 
# Copyright Juguang Xiao
# 
# You may distribute this module under the same terms as perl itself
#
# POD documentation
#

=head1 NAME

Bio::EnsEMBL::Utils::Converter::bio_ens_gene

=head1 SYNOPISIS



=head1 DESCRIPTION

This module is to convert from objects of Bio::SeqFeature::Gene::GeneStructure
to those of Bio::EnsEMBL::Gene

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
