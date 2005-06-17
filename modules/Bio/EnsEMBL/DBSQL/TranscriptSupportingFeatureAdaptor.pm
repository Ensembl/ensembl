# EnsEMBL Adaptor for retrieving SupportingFeatures
#
# Author: Graham McVicker
# 
# Date : 11-Oct-2002
#

=head1 NAME

Bio::EnsEMBL::DBSQL::TranscriptSupportingFeatureAdaptor - Retrieves supporting features 
                                                from the database.

=head1 SYNOPSIS

$supporting_feature_adaptor = $database_adaptor->get_TranscriptSupportingFeatureAdaptor;
@supporting_feats = @{$supporting_feat_adaptor->fetch_all_by_Transcript($transcript)};

=head1 CONTACT

Post questions to ensembl developer mailing list : <ensembl-dev@ebi.ac.uk>

=cut


use strict;

package Bio::EnsEMBL::DBSQL::TranscriptSupportingFeatureAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);

#inherits from BaseAdaptor
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);



=head2 fetch_all_by_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript 
               The transcript to fetch supporting features for
  Example    : @sfs = @{$supporting_feat_adaptor->fetch_all_by_Transcript($transcript)};
  Description: Retrieves supporting features (evidence) for a given transcript. 
  Returntype : list of Bio::EnsEMBL::BaseAlignFeatures in the same coordinate
               system as the $transcript argument
  Exceptions : warning if $transcript is not in the database (i.e. dbID not defined)
               throw if a retrieved supporting feature is of unknown type 
  Caller     : Bio::EnsEMBL::Transcript

=cut

sub fetch_all_by_Transcript {
  my ( $self, $transcript )  = @_;

  my $out = [];

  unless($transcript->dbID) {
    warning("Cannot retrieve evidence for transcript without dbID");
    return [];
  }

  my $sth = $self->prepare("SELECT tsf.feature_type, tsf.feature_id
                            FROM   transcript_supporting_feature tsf
                            WHERE  transcript_id = ?");


  $sth->execute($transcript->dbID());

  my $prot_adp = $self->db->get_ProteinAlignFeatureAdaptor;
  my $dna_adp  = $self->db->get_DnaAlignFeatureAdaptor;

  my $feature;
  while(my ($type, $feature_id) = $sth->fetchrow){
    if($type eq 'protein_align_feature'){
      $feature = $prot_adp->fetch_by_dbID($feature_id);
    } elsif($type eq 'dna_align_feature'){
      $feature = $dna_adp->fetch_by_dbID($feature_id);
    } else {
      warning("Unknown feature type [$type]\n");
    }

    if(!$feature) {
      warning("Supporting feature $type $feature_id does not exist in DB");
    }

    my $new_feature = $feature->transfer($transcript->slice());

    push @$out, $new_feature if( $new_feature );
  }

  $sth->finish();

  return $out;
}

1;

