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
  Arg [2]    : (optional)boolean $sticky_component_flag
               Indicates $exon is a component exon to a sticky exon
  Example    : @supporting = $supporting_feature_adaptor->fetch_by_Exon($exon);
  Description: Retrieves supporting features (evidence) for a given exon. 
  Returntype : list of Bio::EnsEMBL::BaseAlignFeatures in the same coordinate
               system as the $exon argument
  Exceptions : warning if $exon is not in the database (i.e. dbID not defined)
               throw if a retrieved supporting feature is of unknown type 
  Caller     : Bio::EnsEMBL::Exon

=cut

sub fetch_by_Exon {
  my ( $self, $exon, $sticky_component_flag )  = @_;

  my @out;

  # if exon is sticky, get supporting from components
  if( $exon->isa( 'Bio::EnsEMBL::StickyExon' )) {
    foreach my $component_exon ( $exon->each_component_Exon() ) {
      push @out, $self->fetch_by_Exon( $component_exon, 1 );
    }
    return @out;
  }

  unless($exon->dbID) {
    $self->warn("exon has no dbID can't fetch evidence from db" .
		"no relationship exists\n");
    return 0;
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
    
    #we might have to convert the features coordinate system
    unless($feature->contig->name() eq $exon->contig->name()) {
      if($exon->contig()->isa("Bio::EnsEMBL::Slice")) {
	#tranform to slice coords
	$feature->transform($exon->contig());
      } else {
	#this feature is on the wrong contig
	#this is okay only if this is a sticky component exon
	next;
      }
    }
    push @out, $feature;
  }
  return @out;
}


1;
