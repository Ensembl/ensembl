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

Juguang Xiao <juguang@fugu-sg.org>

=cut

=head1 NAME

Bio::EnsEMBL::Utils::Converter::bio_ens_analysis, a converter instance
specific for Analysis.

=head1 SYNOPISIS

  my $converter = Bio::EnsEMBL::Utils::Converter(
    -in  => 'Bio::Pipeline::Analysis',
    -out => 'Bio::EnsEMBL::Analysis'
  );
  my $biopipe_analysis;
  my ($ens_analysis) = @{ $converter->convert( [$biopipe_analysis] ) };

=head1 DESCRIPTION

=head1 METHODS

=cut


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
