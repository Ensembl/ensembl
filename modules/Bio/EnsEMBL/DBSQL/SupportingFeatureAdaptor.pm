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
@supporting_feats = @{$supporting_feat_adaptor->fetch_all_by_Exon($exon)};

=head1 CONTACT

Post questions to ensembl developer mailing list : <ensembl-dev@ebi.ac.uk>

=cut


use strict;

package Bio::EnsEMBL::DBSQL::SupportingFeatureAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use vars qw(@ISA);

#inherits from BaseAdaptor
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);



=head2 fetch_all_by_Exon

  Arg [1]    : Bio::EnsEMBL::Exon $exon 
               The exon to fetch supporting features for
  Example    : @sfs = @{$supporting_feat_adaptor->fetch_all_by_Exon($exon)};
  Description: Retrieves supporting features (evidence) for a given exon. 
  Returntype : list of Bio::EnsEMBL::BaseAlignFeatures in the same coordinate
               system as the $exon argument
  Exceptions : warning if $exon is not in the database (i.e. dbID not defined)
               throw if a retrieved supporting feature is of unknown type 
  Caller     : Bio::EnsEMBL::Exon

=cut

sub fetch_all_by_Exon {
  my ( $self, $exon )  = @_;

  my $out = [];

  #throw if this is a sticky exon
  if($exon->isa('Bio::EnsEMBL::StickyExon')) {
    $self->throw('Expected Exon but got StickyExon. ' .
		 'Call get_all_component_Exons first');
  }

  unless($exon->dbID) {
    $self->warn("exon has no dbID can't fetch evidence from db " .
		"no relationship exists\n");
    return [];
  }


  my $sth = $self->prepare("SELECT sf.feature_type, sf.feature_id
                            FROM   supporting_feature sf
                            WHERE  exon_id = ?");

  $sth->execute($exon->dbID());
  
  my $prot_adp = $self->db->get_ProteinAlignFeatureAdaptor;
  my $dna_adp = $self->db->get_DnaAlignFeatureAdaptor;
  
  my $feature;
  while(my ($type, $feature_id) = $sth->fetchrow){      
    if($type eq 'protein_align_feature'){
      $feature = $prot_adp->fetch_by_dbID($feature_id);
    }elsif($type eq 'dna_align_feature'){
      $feature = $dna_adp->fetch_by_dbID($feature_id);
    }else {
      $self->throw("Unknown feature type [$type]\n");
    }

    if($exon->contig()->isa("Bio::EnsEMBL::Slice")) {
      #tranform to slice coords
      $feature->transform($exon->contig());
    } else {
      #we might have to convert the features coordinate system
      next unless($feature->contig->dbID == $exon->contig->dbID);
    }
    push @$out, $feature;
  }
  return $out;
}



=head2 fetch_by_Exon

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use fetch_all_by_Exon instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_by_Exon {
  my ($self, @args) = @_;

  $self->warn("fetch_by_Exon has been renamed fetch_all_by_Exon\n" . caller);

  return $self->fetch_all_by_Exon(@args);
}




1;
