
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

use Bio::EnsEMBL::DBSQL::DBConnection;

use Bio::EnsEMBL::Utils::Exception qw(deprecate);

@ISA = qw(Bio::EnsEMBL::DBSQL::DBConnection);



=head2 get_SNPAdaptor

  Arg [1]    : none
  Example    : my $snp_adaptor = $db->get_SNPAdaptor();
  Description: Retrieves a SNPAdaptor that can be used to retrieve SNP
               information from the lite database.
  Returntype : Bio::EnsEMBL::Lite::SNPAdaptor
  Exceptions : none
  Caller     : webcode, snpview

=cut

sub get_SNPAdaptor {
  my $self = shift;

  return $self->_get_adaptor("Bio::EnsEMBL::Lite::SNPAdaptor");
}


=head2 get_GeneAdaptor

  Description: DEPRECATED, use core GeneAdaptor instead.

=cut

sub get_GeneAdaptor {
  my $self = shift;

  deprecate("Genes are no longer stored in the lite database.\n" .
            "Use the core GeneAdaptor instead.");

  my $core = $self->get_db_adaptor('core');

  return ($core) ? $core->get_GeneAdaptor() : undef;
}


=head2 get_RepeatFeatureAdaptor

  Description: DEPRECATED, use core RepeatFeatureAdaptor instead.

=cut

sub get_RepeatFeatureAdaptor {
  my $self = shift;

  deprecate("Repeats are no longer stored in the lite database.\n" .
            "Use the core RepeatAdaptor instead.");


  my $core = $self->get_db_adaptor('core');

  return ($core) ? $core->get_RepeatFeatureAdaptor : undef;
}


=head2 get_ChromosomeAdaptor

  Description: DEPRECATED, use core SliceAdaptor instead.

=cut

sub get_ChromosomeAdaptor {
  my $self = shift;

  deprecate("Chromosomes are no longer stored in the lite database.\n" .
            "Use the core SliceAdaptor instead.");

  my $core = $self->get_db_adaptor('core');

  return ($core) ? $core->get_SliceAdaptor() : undef;
}

=head2 get_DensityAdaptor

  Description: DEPRECATED, use core DenstiyFeatureAdaptor instead.

=cut

sub get_DensityAdaptor {
  my $self = shift;

  deprecate("DensityFeatures are no longer stored in the lite database.\n" .
            "Use the core DensityFeatureAdaptor instead.");

  my $core = $self->get_db_adaptor('core');

  return ($core) ? $core->get_DensityFeatureAdaptor() : undef;
}



1;
