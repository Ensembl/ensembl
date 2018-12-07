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

=cut

=head1 NAME

Bio::EnsEMBL::ApiVersion

=head1 SYNOPSIS

  use Bio::EnsEMBL::ApiVersion;

  printf( "The API version used is %s\n", software_version() );

=head1 DESCRIPTION

The module exports the software_version() subroutine which returns the
release version of the Ensembl Core API.

=cut

package Bio::EnsEMBL::ApiVersion;

use strict;
use warnings;

use Exporter;

use base qw( Exporter );

our @EXPORT = qw( software_version );

my $API_VERSION = 96;

sub software_version { return $API_VERSION }

1;
