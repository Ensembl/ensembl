=head1 LICENSE

  Copyright (c) 1999-2010 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Lite::DBAdaptor

=head1 SYNOPSIS

  $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -user   => 'anonymous',
    -dbname => 'homo_sapiens_lite_20_34c',
    -host   => 'ensembldb.ensembl.org',
    -driver => 'mysql'
  );

  $snp_adaptor = $db->get_SNPAdaptor();

  @snps = @{ $snp_adaptor->fetch_all_by_Slice($slice) }

=head1 DESCRIPTION

This is a database connection to the denormalised lite database.  It
allows for the rapid creation of drawable objects that are too slow
to retreive from normalised data sources. Formerly this included many
Ensembl objects such as genes, transcript, exons, etc. but is now
limited to just SNPs.

=head1 METHODS

=cut

package Bio::EnsEMBL::Lite::DBAdaptor;

use vars qw(@ISA);
use strict;

@ISA = qw(Bio::EnsEMBL::DBSQL::DBAdaptor);

sub get_available_adaptors{
  my %pairs = ("SNP", "Bio::EnsEMBL::Lite::SNPAdaptor");
  return (\%pairs);
}



1;
