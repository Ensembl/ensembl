# Bio::EnsEMBL::Utils::Converter::bio_ens_analysis
#
# Created and cared for by Juguang Xiao <juguang@fugu-sg.org>
# Created date: 15/3/2003
# 
# Copyright Juguang Xiao
# 
# You may distribute this module under the same terms as perl itself
#
# POD documentation
#

=head1 NAME

Bio::EnsEMBL::Utils::Converter::bio_ens_analysis, a converter instance specific for Analysis.

=head1 SYNOPISIS

my $converter = Bio::EnsEMBL::Utils::Converter(
    -in => 'Bio::Pipeline::Analysis',
    -out => 'Bio::EnsEMBL::Analysis'
);
my $biopipe_analysis;
my ($ens_analysis) = @{$converter->convert([$biopipe_analysis])};

=head1 DESCRIPTION


=head1 FEEDBACK

=head1 AUTHOR Juguang Xiao

Juguang Xiao <juguang@fugu-sg.org>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin ...

package Bio::EnsEMBL::Utils::Converter::bio_ens_analysis;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Utils::Converter::bio_ens;
use Bio::EnsEMBL::Analysis;

@ISA = qw(Bio::EnsEMBL::Utils::Converter::bio_ens);

sub _convert_single {
    my ($self, $input) = @_;

    $self->throw("a Bio::Pipeline::Analysis object needed")
        unless(ref($input) && $input->isa('Bio::Pipeline::Analysis'));
    
    my $ens_analysis = Bio::EnsEMBL::Analysis->new(
        -logic_name => $input->logic_name,
        -db => $input->db,
        -db_version => $input->db_version,
        -db_file => $input->db_file,
        -program => $input->program,
        -program_version => $input->program_version,
        -program_file => $input->program_file,
        -parameters => $input->analysis_parameters,
        -module => $input->runnable,
        -gff_source => $input->gff_source,
        -gff_feature => $input->gff_feature,
        -id =>$input->dbID
    );

    return $ens_analysis;
}
;
