#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor

=head1 SYNOPSIS

  use Bio::EnsEMBL::DBSQL::DBAdaptor;
  use Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor;

  $db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -user => 'root', 
                                            -db => 'pog' ,
                                            -host => 'caldy' ,
                                            -driver => 'mysql' );
  my $pfa = $db->get_ProteinFeatureAdaptor();

  my @prot_feats = @{$pfa->fetch_all_by_translation_id(1231)};

  my $prot_feat = $pfa->fetch_by_dbID(523);

  
=head1 DESCRIPTION


=head1 CONTACT

  Post questions to the EnsEMBl development mailing list : 
  ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut


package Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor;

use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::ProteinFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);


use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 fetch_all_by_translation_id

  Arg [1]    : int $transl
               the internal id of the translation corresponding to protein 
               whose features are desired 
  Example    : @prot_feats = $prot_feat_adaptor->fetch_by_translation_id(1234);
  Description: Gets all protein features present on a peptide using the
               translations internal identifier.  This method will return
               an unsorted list of all protein_feature types.  The feature
               types may be distinguished using the logic name attribute of
               the attached analysis objects.   
  Returntype : listref of Bio::EnsEMBL::ProteinFeatures
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub fetch_all_by_translation_id {
  my($self,$translation_id) = @_;

  if(!$translation_id) {
    throw("translation_id argument is required\n");
  }

  my @features;
  my $analysis_adaptor = $self->db()->get_AnalysisAdaptor();
    
  my $sth = $self->prepare
    ("SELECT protein_feature_id, p.seq_start, p.seq_end, p.analysis_id, " .
     "       p.score, p.perc_ident, p.evalue, p.hit_start, p.hit_end, " .
     "       p.hit_id, x.display_label, i.interpro_ac " .
     "FROM   protein_feature p " .
     "LEFT JOIN interpro AS i ON p.hit_id = i.id " .
     "LEFT JOIN xref AS x ON x.dbprimary_acc = i.interpro_ac " .
     "WHERE p.translation_id = ?");

  $sth->bind_param(1,$translation_id,SQL_INTEGER);
  $sth->execute();

  while( my $row = $sth->fetchrow_arrayref) {
    my ($dbID, $start, $end, $analysisid, $score, $perc_id, $evalue, $hstart,
        $hend,$hid,$desc, $interpro_ac) = @$row;

    my $analysis = $analysis_adaptor->fetch_by_dbID($analysisid);

    if(!$analysis) {
      warning("Analysis with dbID=$analysisid does not exist\n" .
	      "but is referenced by ProteinFeature $dbID");
    }

    my $feat = Bio::EnsEMBL::ProteinFeature->new
      (-DBID        => $dbID,
       -ADAPTOR     => $self,
       -SEQNAME     => $translation_id,
       -START       => $start,
       -END         => $end,
       -ANALYSIS    => $analysis,
       -PERCENT_ID  => $perc_id,
       -P_VALUE     => $evalue,
       -SCORE       => $score, 
       -HSTART      => $hstart,
       -HEND        => $hend,
       -HSEQNAME    => $hid,
       -IDESC       => $desc,
       -INTERPRO_AC => $interpro_ac);

    push(@features,$feat);
  }

  $sth->finish();

  return \@features;    
}


=head2 fetch_by_dbID

  Arg [1]    : int $protfeat_id
               the unique database identifier of the protein feature to obtain
  Example    : my $feature = $prot_feat_adaptor->fetch_by_dbID();
  Description: Obtains a protein feature object via its unique id
  Returntype : Bio::EnsEMBL::ProteinFeauture
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub fetch_by_dbID{
  my ($self,$protfeat_id) = @_;
    
  my $sth = $self->prepare(
	      "SELECT p.seq_start, p.seq_end, p.analysis_id, " .
        "       p.score, p.perc_ident, p.evalue, " .
        "       p.hit_start, p.hit_end, p.hit_id, " .
        "       x.display_label, i.interpro_ac " .
        "FROM   protein_feature p " .
        "LEFT JOIN interpro AS i ON p.hit_id = i.id " .
        "LEFT JOIN xref AS x ON x.dbprimary_acc = i.interpro_ac " .
        "WHERE  p.protein_feature_id = ?");
  
  $sth->bind_param(1,$protfeat_id,SQL_INTEGER);
  my $res = $sth->execute();
   
  if($sth->rows == 0) {
    $sth->finish();
    return undef;
  }

  my ($start, $end, $analysis_id, $score, $perc_ident, $pvalue, $hstart, 
      $hend, $hseqname, $idesc, $interpro_ac) = $sth->fetchrow_array();

  $sth->finish();

  my $analysis = $self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id);

  return Bio::EnsEMBL::ProteinFeature->new
    (-ADAPTOR    => $self,
     -DBID        => $protfeat_id,
     -START       => $start,
     -END         => $end,
     -HSTART      => $hstart,
     -HEND        => $hend,
     -HSEQNAME    => $hseqname,
     -ANALYSIS    => $analysis,
     -SCORE       => $score,
     -P_VALUE     => $pvalue,
     -PERCENT_ID  => $perc_ident,
     -IDESC       => $idesc,
     -INTERPRO_AC => $interpro_ac);
}




=head2 store

  Arg [1]    : Bio::EnsEMBL::ProteinFeature $feature
               The feature to be stored
  Arg [2]    : int $translation_id
            
  Example    : $protein_feature_adaptor->store($protein_feature);
  Description: Stores a protein feature in the database
  Returntype : int - the new internal identifier of the stored protein feature
  Exceptions : thrown if arg is not a Bio::EnsEMBL:
  Caller     : none
  Status     : Stable

=cut

sub store {
  my ($self,$feature, $translation_id) = @_;

  if(!ref($feature) || !$feature->isa('Bio::EnsEMBL::ProteinFeature')) {
    throw("ProteinFeature argument is required");
  }

  if(!$translation_id) {
    deprecate("Calling ProteinFeatureAdaptor without a translation_id is " .
              "deprecated.  Pass a translation_id argument rather than " .
              "setting the ProteinFeature seqname to be the translation " .
              "id");
    $translation_id = $feature->seqname();
  }

  my $db = $self->db();

  if($feature->is_stored($db)) {
    warning("ProteinFeature " . $feature->dbID() . " is already stored in " .
            "this database - not storing again");
  }

  my $analysis = $feature->analysis();
  if (!defined($analysis)) {
    throw("Feature doesn't have analysis. Can't write to database");
  }

  if(!$analysis->is_stored($db)) {
    $db->get_AnalysisAdaptor->store($analysis);
  }

  my $sth = 
    $self->prepare("INSERT INTO protein_feature " .
                   "        SET translation_id = ?, " .
                   "            seq_start      = ?, ".
                   "            seq_end        = ?, ".
                   "            analysis_id    = ?, ".
                   "            hit_start      = ?, ".
                   "            hit_end        = ?, ".
                   "            hit_id         = ?, ".
                   "            score          = ?, ".
                   "            perc_ident     = ?, ".
                   "            evalue         = ?");

  $sth->bind_param(1,$translation_id,SQL_INTEGER);
  $sth->bind_param(2,$feature->start,SQL_INTEGER);
  $sth->bind_param(3,$feature->end,SQL_INTEGER);
  $sth->bind_param(4,$analysis->dbID,SQL_INTEGER);
  $sth->bind_param(5,$feature->hstart,SQL_INTEGER);
  $sth->bind_param(6,$feature->hend,SQL_INTEGER);
  $sth->bind_param(7,$feature->hseqname,SQL_VARCHAR);
  $sth->bind_param(8,$feature->score,SQL_DOUBLE);
  $sth->bind_param(9,$feature->percent_id,SQL_FLOAT);
  $sth->bind_param(10,$feature->p_value,SQL_DOUBLE);

  $sth->execute();
  
  my $dbID = $sth->{'mysql_insertid'};

  $feature->adaptor($self);
  $feature->dbID($dbID);
  
  $sth->finish();

  return $dbID;          
}



sub fetch_by_translation_id {
  deprecate("Use fetch_all_by_translation_id instead.");
  fetch_all_by_translation_id(@_);
}

sub fetch_all_by_feature_and_dbID {
  my $self = shift;
  my $feature = shift;
  my $translation_id = shift;
  deprecate("Use fetch_all_by_translation_id instead.");

  print STDERR "translation_id = $translation_id feature = $feature\n";

  my $features = $self->fetch_all_by_translation_id($translation_id);

  my @out;
  foreach my $f (@$features) {
    my $logic_name = lc($f->analysis->logic_name());
    print STDERR "LOGIC_NAME = $logic_name | FEATURE = $feature\n";
    push(@out, $f) if($logic_name eq lc($feature));
  }

  return \@out;
}

1;

