
package Bio::EnsEMBL::Utils::Converter::bio_ens_analysis;
use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Utils::Converter::bio_ens;
use Bio::EnsEMBL::Analysis;

@ISA = qw(Bio::EnsEMBL::Utils::Converter::bio_ens);

sub _convert_single {
    my ($self, $input) = @_;

    $self->throw("a Bio::Pipeline::Analysis object needed"){
        unless(ref($input) && $input->isa('Bio::Pipeline::Analysis');
    }

    my $ens_analysis = new Bio::EnsEMBL::Analysis(
        -logic_name => $input->logic_name
    );

    return $ens_analysis;
}
;
