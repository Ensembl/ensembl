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

my $API_VERSION = 70;

sub software_version { return $API_VERSION }

1;
