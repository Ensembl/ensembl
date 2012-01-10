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

Bio::EnsEMBL::DBSQL::SupportingFeatureAdaptor - Retrieves supporting
features from the database.

=head1 SYNOPSIS

  my $sfa =
    $registry->get_adaptor( 'Human', 'Core', 'SupportingFeature' );

  my @supporting_feats = @{ $sfa->fetch_all_by_Exon($exon) };

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::SupportingFeatureAdaptor;

use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);

#inherits from BaseAdaptor
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);



=head2 fetch_all_by_Exon

  Arg [1]    : Bio::EnsEMBL::Exon $exon 
               The exon to fetch supporting features.
  Example    : @sfs =
                @{ $supporting_feat_adaptor->fetch_all_by_Exon($exon) };
  Description: Retrieves supporting features (evidence) for a given
               exon.
  Returntype : List of Bio::EnsEMBL::BaseAlignFeatures in the same
               coordinate system as the $exon argument
  Exceptions : Warning if $exon is not in the database (i.e. dbID
               not defined).
               Throw if a retrieved supporting feature is of unknown
               type.
  Caller     : Bio::EnsEMBL::Exon
  Status     : Stable

=cut

sub fetch_all_by_Exon {
  my ( $self, $exon )  = @_;

  my $out = [];

  unless($exon->dbID) {
    warning("Cannot retrieve evidence for exon without dbID");
    return [];
  }

  my $sth = $self->prepare("SELECT sf.feature_type, sf.feature_id
                            FROM   supporting_feature sf
                            WHERE  exon_id = ?");

  $sth->bind_param(1,$exon->dbID,SQL_INTEGER);
  $sth->execute();

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
    } else {
      my $new_feature = $feature->transfer($exon->slice());
      push @$out, $new_feature if( $new_feature );
    }
  }

  $sth->finish();

  return $out;
}

=head2 store

  Arg [1]    : Int $exonsID
               The dbID of an EnsEMBL exon to associate with
               supporting features.
  Arg [2]    : Ref to array of Bio::EnsEMBL::BaseAlignFeature
               (the support)
  Example    : $sfa->store($exon_id, \@features);
  Description: Stores a set of alignment features and associates an
               EnsEMBL exon with them
  Returntype : none
  Exceptions : thrown when invalid dbID is passed to this method
  Caller     : TranscriptAdaptor
  Status     : Stable

=cut

sub store {
  my ( $self, $exon_dbID, $aln_objs ) = @_;

  my $pep_check_sql = 
      "SELECT protein_align_feature_id " . 
      "FROM protein_align_feature " . 
      "WHERE seq_region_id = ? " . 
      "AND   seq_region_start = ? " . 
      "AND   seq_region_end   = ? " .
      "AND   seq_region_strand = ? " . 
      "AND   hit_name = ? " . 
      "AND   hit_start = ? " . 
      "AND   hit_end   = ? " . 
      "AND   analysis_id = ? " . 
      "AND   cigar_line = ? ";

  my $dna_check_sql = 
      "SELECT dna_align_feature_id " . 
      "FROM dna_align_feature " . 
      "WHERE seq_region_id = ? " . 
      "AND   seq_region_start = ? " . 
      "AND   seq_region_end   = ? " .
      "AND   seq_region_strand = ? " . 
      "AND   hit_name = ? " . 
      "AND   hit_start = ? " . 
      "AND   hit_end   = ? " . 
      "AND   analysis_id = ? " . 
      "AND   cigar_line = ? " . 
      "AND   hit_strand = ? ";

  my $assoc_check_sql = 
      "SELECT * " .  
      "FROM  supporting_feature " . 
      "WHERE exon_id = $exon_dbID " . 
      "AND   feature_type = ? " . 
      "AND   feature_id   = ? ";

  my $assoc_write_sql = "INSERT into supporting_feature " . 
      "(exon_id, feature_id, feature_type) " . 
      "values(?, ?, ?)";

  my $pep_check_sth = $self->prepare($pep_check_sql);
  my $dna_check_sth = $self->prepare($dna_check_sql);
  my $assoc_check_sth = $self->prepare($assoc_check_sql);
  my $sf_sth = $self->prepare($assoc_write_sql);

  my $dna_adaptor = $self->db->get_DnaAlignFeatureAdaptor();
  my $pep_adaptor = $self->db->get_ProteinAlignFeatureAdaptor();

  foreach my $f (@$aln_objs) {
    # check that the feature is in toplevel coords

    if($f->slice->start != 1 || $f->slice->strand != 1) {
    #move feature onto a slice of the entire seq_region
      my $tls = $self->db->get_sliceAdaptor->fetch_by_region($f->slice->coord_system->name(),
                                                             $f->slice->seq_region_name(),
                                                             undef, #start
                                                             undef, #end
                                                             undef, #strand
                                                             $f->slice->coord_system->version());
      $f = $f->transfer($tls);

      if(!$f) {
        throw('Could not transfer Feature to slice of ' .
              'entire seq_region prior to storing');
      }
    }

    if(!$f->isa("Bio::EnsEMBL::BaseAlignFeature")){
      throw("$f must be an align feature otherwise" .
            "it can't be stored");
    }
    
    my ($sf_dbID, $type, $adap, $check_sth);
    
    my @check_args = ($self->db->get_SliceAdaptor->get_seq_region_id($f->slice),
                      $f->start,
                      $f->end,
                      $f->strand,
                      $f->hseqname,
                      $f->hstart,
                      $f->hend,
                      $f->analysis->dbID,
                      $f->cigar_string);

    if($f->isa("Bio::EnsEMBL::DnaDnaAlignFeature")){
      $adap = $dna_adaptor;      
      $check_sth = $dna_check_sth;
      $type = 'dna_align_feature';
      push @check_args, $f->hstrand;
    } elsif($f->isa("Bio::EnsEMBL::DnaPepAlignFeature")){
      $adap = $pep_adaptor;
      $check_sth = $pep_check_sth;
      $type = 'protein_align_feature';
    } else {
      warning("Supporting feature of unknown type. Skipping : [$f]\n");
      next;
    }

    $check_sth->execute(@check_args);
    $sf_dbID = $check_sth->fetchrow_array;
    if (not $sf_dbID) {
      $adap->store($f);
      $sf_dbID = $f->dbID;
    }

    # now check association
    $assoc_check_sth->execute($type,
                              $sf_dbID);
    if (not $assoc_check_sth->fetchrow_array) {    
      $sf_sth->bind_param(1, $exon_dbID, SQL_INTEGER);
      $sf_sth->bind_param(2, $sf_dbID, SQL_INTEGER);
      $sf_sth->bind_param(3, $type, SQL_VARCHAR);
      $sf_sth->execute();
    }
  }

  $dna_check_sth->finish;
  $pep_check_sth->finish;
  $assoc_check_sth->finish;
  $sf_sth->finish;
  
}



1;

