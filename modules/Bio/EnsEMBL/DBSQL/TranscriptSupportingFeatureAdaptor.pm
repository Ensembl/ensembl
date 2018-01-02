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

Bio::EnsEMBL::DBSQL::TranscriptSupportingFeatureAdaptor - Retrieves
supporting features from the database.

=head1 SYNOPSIS

  $supporting_feature_adaptor =
    $database_adaptor->get_TranscriptSupportingFeatureAdaptor;

  @supporting_feats =
    @{ $supporting_feat_adaptor->fetch_all_by_Transcript($transcript) };

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::TranscriptSupportingFeatureAdaptor;

use strict;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);

#inherits from BaseAdaptor
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);



=head2 fetch_all_by_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript $transcript 
               The transcript to fetch supporting features for
  Arg [2]    : String (optional)
               Feature type to filter upon (either 'dna_align_feature' or 'protein_align_feature')
  Example    : @sfs = @{$supporting_feat_adaptor->fetch_all_by_Transcript($transcript)};
  			   @sfs = @{$supporting_feat_adaptor->fetch_all_by_Transcript($transcript, $feature_type)};	
  Description: Retrieves supporting features (evidence) for a given transcript. 
  Returntype : list of Bio::EnsEMBL::BaseAlignFeatures in the same coordinate
               system as the $transcript argument
  Exceptions : warning if $transcript is not in the database (i.e. dbID not defined)
               throw if a retrieved or requested supporting feature is of unknown type 
  Caller     : Bio::EnsEMBL::Transcript
  Status     : Stable

=cut

sub fetch_all_by_Transcript {
  my ( $self, $transcript, $feature_type )  = @_;

  my $out = [];
  my $out_feature_type = {};

  unless($transcript->dbID) {
    warning("Cannot retrieve evidence for transcript without dbID");
    return [];
  }
  
  if(defined $feature_type && $feature_type !~ /(dna)|(protein)_align_feature/) {
    throw("feature type must be dna_align_feature or protein_align_feature");
  }

  my $sth = $self->prepare("SELECT tsf.feature_type, tsf.feature_id
                            FROM   transcript_supporting_feature tsf
                            WHERE  transcript_id = ?");


  $sth->bind_param(1,$transcript->dbID,SQL_INTEGER);
  $sth->execute();

  my $prot_adp = $self->db->get_ProteinAlignFeatureAdaptor;
  my $dna_adp  = $self->db->get_DnaAlignFeatureAdaptor;

  my $feature;
  while(my ($type, $feature_id) = $sth->fetchrow){
   if ($type eq 'protein_align_feature') {
        $feature = $prot_adp-> fetch_by_dbID($feature_id);
    }
    elsif($type eq 'dna_align_feature') {
        $feature = $dna_adp-> fetch_by_dbID($feature_id);
    } else {
        warning("Unknown feature type [$type]\n");
    }
    
    if(!$feature) {
      warning("Supporting feature $type $feature_id does not exist in DB");
    } else {
      my $new_feature = $feature->transfer($transcript->slice());
 
      push @{$out_feature_type->{$type}}, $new_feature if ($new_feature);
    }
    }
	
	$sth->finish();
   
    if(defined $feature_type){
      return $out_feature_type->{$feature_type};
    }else{
      while(my ($feature_type, $new_features) = each(%$out_feature_type)){
      	push @$out, @{$new_features};
  	}
  	
    }
  
   return $out;
}



=head2 store
  Arg [2]    : Int $transID
               The dbID of an EnsEMBL transcript to associate with supporting
               features
  Arg [1]    : Ref to array of Bio::EnsEMBL::BaseAlignFeature (the support)
  Example    : $dbea->store($transcript_id, \@features);
  Description: Stores a set of alignment features and associates an EnsEMBL transcript
               with them
  Returntype : none
  Exceptions : thrown when invalid dbID is passed to this method
  Caller     : TranscriptAdaptor
  Status     : Stable

=cut

sub store {
  my ( $self, $tran_dbID, $aln_objs ) = @_;

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
      "AND   cigar_line = ? " .
      "AND   hcoverage = ? ";

  my $dna_check_sql = 
      "SELECT dna_align_feature_id " . 
      "FROM  dna_align_feature " . 
      "WHERE seq_region_id = ? " . 
      "AND   seq_region_start = ? " . 
      "AND   seq_region_end   = ? " .
      "AND   seq_region_strand = ? " . 
      "AND   hit_name = ? " . 
      "AND   hit_start = ? " . 
      "AND   hit_end   = ? " . 
      "AND   analysis_id = ? " . 
      "AND   cigar_line = ? " .
      "AND   hcoverage = ? " . 
      "AND   hit_strand = ? ";

  my $assoc_check_sql = 
      "SELECT * " .  
      "FROM  transcript_supporting_feature " . 
      "WHERE transcript_id = $tran_dbID " . 
      "AND   feature_type = ? " . 
      "AND   feature_id   = ? ";

  my $assoc_write_sql = "INSERT into transcript_supporting_feature " . 
      "(transcript_id, feature_id, feature_type) " . 
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
                      $f->cigar_string,
		      $f->hcoverage);
    
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
      $sf_sth->bind_param(1, $tran_dbID, SQL_INTEGER);
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

