# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=pod

=encoding UTF-8

=head1 NAME

Bio::EnsEMBL - API to connect to and work with EnsEMBL genomic databases

=head1 SYNOPSIS

    use Bio::EnsEMBL::Registry;

    Bio::EnsEMBL::Registry->load_registry_from_db(
        -host => 'ensembldb.ensembl.org',
        -user => 'anonymous',
        -species => 'homo sapiens',
        -group => 'core'
    );
    my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor(
        'homo sapiens', 'Core', 'Slice'
    );
    my $slice = $slice_adaptor->fetch_by_gene_stable_id('ENSG00000101266');
    say $slice->display_id();

=head1 DESCRIPTION

L<Bio::EnsEMBL> is the namespace for the EnsEMBL Core API.
The Core API provides access to the EnsEMBL genomic databases.

Most people will want to use L<Bio::EnsEMBL::Registry> as an entry point.

=head1 SEE ALSO

L<https://www.ensembl.org/info/docs/api/index.html>

=head1 SUPPORT

Please email comments or questions to the public EnsEMBL developers list at
L<http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the EnsEMBL help desk at
L<http://www.ensembl.org/Help/Contact>.

=head1 COPYRIGHT AND LICENCE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Copyright [2016-2024] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0

=cut

use strict;
use warnings;
package Bio::EnsEMBL;
# ABSTRACT: API to connect to and work with EnsEMBL genomic databases

1;
