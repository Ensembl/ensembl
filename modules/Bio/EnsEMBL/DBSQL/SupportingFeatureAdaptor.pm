# EnsEMBL Adaptor for retrieving SupportingFeatures
#
# Author: Graham McVicker
# 
# Date : 11-Oct-2002
#

=head1 NAME

Bio::EnsEMBL::DBSQL::SupportingFeatureAdaptor - Retrieves supporting features 
                                                from the database.

=head1 SYNOPSIS

$supporting_feature_adaptor = $database_adaptor->get_SupportingFeatureAdaptor;
@supporting_features = $supporting_feature_adaptor->fetch_by_Exon($exon);

=head1 CONTACT

Post questions to ensembl developer mailing list : <ensembl-dev@ebi.ac.uk>

=cut


use strict;

package Bio::EnsEMBL::DBSQL::SupportingFeatureAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use vars qw(@ISA);

#inherits from BaseAdaptor
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);



=head2 fetch_by_Exon

  Arg [1]    : Bio::EnsEMBL::Exon $exon 
               The exon to fetch supporting features for
  Example    : @supporting = $supporting_feature_adaptor->fetch_by_Exon($exon);
  Description: Retrieves supporting features (evidence) for a given exon. 
  Returntype : list of Bio::EnsEMBL::BaseAlignFeatures
  Exceptions : none
  Caller     : Bio::EnsEMBL::Exon

=cut

sub fetch_by_Exon {
  my ( $self, $exon )  = @_;

  # if exon is sticky, get supporting from components
  if( $exon->isa( 'Bio::EnsEMBL::StickyExon' )) {
    my @component_exons = $exon->each_component_exon();
    for my $component_exon ( @component_exons ) {
      $self->fetch_evidence_by_Exon( $component_exon );
    }
    return;
  }

  unless($exon->dbID) {
    $self->warn("exon has no dbID can't fetch evidence from db" .
		"no relationship exists\n");
    return 0;
  }


  my $sth = $self->prepare("SELECT feature_type, feature_id
                            FROM   supporting_feature
                            WHERE  exon_id = ?");

  $sth->execute($exon->dbID());
  
  my $prot_adp = $self->db->get_ProteinAlignFeatureAdaptor;
  my $dna_adp = $self->db->get_DnaAlignFeatureAdaptor;
  
  my @out;

  while(my ($type, $feature_id) = $sth->fetchrow){      
    if($type eq 'protein_align_feature'){
      push @out, $prot_adp->fetch_by_dbID($feature_id);
    }elsif($type eq 'dna_align_feature'){
      push @out, $dna_adp->fetch_by_dbID($feature_id);
    }
  }

  return @out;
}


1;
