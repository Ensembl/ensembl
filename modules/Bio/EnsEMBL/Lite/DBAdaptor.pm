
=Head1 NAME - Bio::EnsEMBL::Lite::DBAdaptor

=head1 SYNOPSIS

    $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -user   => 'anonymous',
        -dbname => 'homo_sapiens_lite_20_34c',
        -host   => 'ensembldb.ensembl.org',
        -driver => 'mysql');

    $snp_adaptor = $db->get_SNPAdaptor();

    @snps = @{$snp_adaptor->fetch_all_by_Slice($slice)}

=head1 DESCRIPTION

This is a database connection to the denormalised lite database.
It allows for the rapid creation of drawable objects that are too slow
to retreive from normalised data sources. Formerly this included many
Ensembl objects such as genes, transcript, exons, etc. but is now limited
to just SNPs.

=head1 CONTACT

This module is part of the Ensembl project: www.ensembl.org
Post questions to the Ensembl development list: <ensembl-dev@ebi.ac.uk>

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
