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

Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor

=head1 SYNOPSIS

  use Bio::EnsEMBL::Registry;

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
  );

  $pfa = Bio::EnsEMBL::Registry->get_adaptor( "human", "core",
    "proteinfeature" );

  my @prot_feats = @{ $pfa->fetch_all_by_translation_id(1231) };

  my $prot_feat = $pfa->fetch_by_dbID(523);

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
  Example    : @prot_feats =
                  @{ $prot_feat_adaptor->fetch_by_translation_id(1234) };
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
  my ($self, $translation_id) = @_;

  my $constraint = "pf.translation_id = ?";
  $self->bind_param_generic_fetch($translation_id, SQL_INTEGER);
  my $features = $self->generic_fetch($constraint);

  return $features;
} ## end sub fetch_all_by_translation_id

=head2 fetch_all_by_logic_name

  Arg [1]    : string $logic_name
               the logic name of the type of features to obtain
  Example    : $fs = $a->fetch_all_by_logic_name('foobar');
  Description: Returns a listref of features created from the database.
               only features with an analysis of type $logic_name will
               be returned.  If the logic name is invalid (not in the
               analysis table), a reference to an empty list will be
               returned.
  Returntype : listref of Bio::EnsEMBL::ProteinFeatures
  Exceptions : thrown if no $logic_name
  Caller     : General
  Status     : Stable

=cut

sub fetch_all_by_logic_name {
  my ( $self, $logic_name ) = @_;

  if ( !defined($logic_name) ) {
    throw("Need a logic_name");
  }

  my $constraint = $self->_logic_name_to_constraint( '', $logic_name );

  if ( !defined($constraint) ) {
    warning("Invalid logic name: $logic_name");
    return [];
  }

  return $self->generic_fetch($constraint);
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

sub fetch_by_dbID {
  my ($self, $protfeat_id) = @_;

  my $sth = $self->prepare("SELECT p.translation_id, p.seq_start, p.seq_end, p.analysis_id, "
                            . "p.score, p.perc_ident, p.evalue, "
                            . "p.hit_start, p.hit_end, p.hit_name, p.hit_description, "
                            . "x.description, x.display_label, i.interpro_ac "
                            . "FROM   protein_feature p "
                            . "LEFT JOIN interpro AS i ON p.hit_name = i.id "
                            . "LEFT JOIN xref AS x ON x.dbprimary_acc = i.interpro_ac "
                            . "WHERE  p.protein_feature_id = ?");

  $sth->bind_param(1, $protfeat_id, SQL_INTEGER);
  my $res = $sth->execute();
   
  my ($translation_id, $start, $end, $analysis_id, $score, $perc_ident, $pvalue, $hstart,
      $hend, $hseqname, $hdesc, $idesc, $ilabel, $interpro_ac) = $sth->fetchrow_array();

  if($sth->rows == 0) {
    $sth->finish();
    return undef;
  }

  $sth->finish();

  my $analysis = $self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id);

  return
	Bio::EnsEMBL::ProteinFeature->new(-ADAPTOR     => $self,
									  -DBID        => $protfeat_id,
									  -START       => $start,
									  -END         => $end,
									  -HSTART      => $hstart,
									  -HEND        => $hend,
									  -HSEQNAME    => $hseqname,
									  -HDESCRIPTION => $hdesc,
									  -ANALYSIS    => $analysis,
									  -SCORE       => $score,
									  -P_VALUE     => $pvalue,
									  -PERCENT_ID  => $perc_ident,
									  -IDESC       => $idesc,
                                      -ILABEL      => $ilabel,
									  -INTERPRO_AC => $interpro_ac,
									  -TRANSLATION_ID  => $translation_id);
} ## end sub fetch_by_dbID

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
  my ($self, $feature, $translation_id) = @_;

  if (!ref($feature) || !$feature->isa('Bio::EnsEMBL::ProteinFeature')) {
	throw("ProteinFeature argument is required");
  }

  if (!$translation_id) {
	deprecate("Calling ProteinFeatureAdaptor without a translation_id is " . "deprecated.  Pass a translation_id argument rather than " . "setting the ProteinFeature seqname to be the translation " . "id");
	$translation_id = $feature->seqname();
  }

  my $db = $self->db();

  if ($feature->is_stored($db)) {
	warning("ProteinFeature " . $feature->dbID() . " is already stored in " . "this database - not storing again");
  }

  my $analysis = $feature->analysis();
  if (!defined($analysis)) {
	throw("Feature doesn't have analysis. Can't write to database");
  }
  if (!$analysis->is_stored($db)) {
	$db->get_AnalysisAdaptor->store($analysis);
  }

  my $insert_ignore = $self->insert_ignore_clause();

  my $sth = $self->prepare("
    ${insert_ignore} INTO protein_feature
                ( translation_id,
                  seq_start,
                  seq_end,
                  analysis_id,
                  hit_start,
                  hit_end,
                  hit_name,
                  hit_description,
                  score,
                  perc_ident,
                  evalue     )
         VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ? )
  ");

  $sth->bind_param(1,  $translation_id,        SQL_INTEGER);
  $sth->bind_param(2,  $feature->start,        SQL_INTEGER);
  $sth->bind_param(3,  $feature->end,          SQL_INTEGER);
  $sth->bind_param(4,  $analysis->dbID,        SQL_INTEGER);
  $sth->bind_param(5,  $feature->hstart,       SQL_INTEGER);
  $sth->bind_param(6,  $feature->hend,         SQL_INTEGER);
  $sth->bind_param(7,  $feature->hseqname,     SQL_VARCHAR);
  $sth->bind_param(8,  $feature->hdescription, SQL_LONGVARCHAR);
  $sth->bind_param(9,  $feature->score,        SQL_DOUBLE);
  $sth->bind_param(10, $feature->percent_id,   SQL_FLOAT);
  $sth->bind_param(11, $feature->p_value,      SQL_DOUBLE);

  $sth->execute();

  if (defined($sth->err) && $sth->err eq 0){ # is a warning if 0 and defined
      warning('SQL warning : ' . $sth->errstr ."\n");
  }
 
  my $dbID = $self->last_insert_id('protein_feature_id', undef, 'protein_feature');

  $feature->adaptor($self);
  $feature->dbID($dbID);

  $sth->finish();

  return $dbID;
} ## end sub store

sub _tables {
  my $self = shift;

  return (['protein_feature', 'pf'], ['interpro', 'ip'], ['xref', 'x']);
}

sub _left_join {
  return (['interpro', "pf.hit_name = ip.id"], ['xref', "x.dbprimary_acc = ip.interpro_ac"]);
}

sub _columns {
  my $self = shift;

  return qw( pf.protein_feature_id
             pf.translation_id pf.seq_start pf.seq_end
             pf.hit_start pf.hit_end pf.hit_name pf.hit_description
             pf.analysis_id pf.score pf.evalue pf.perc_ident
             x.description x.display_label ip.interpro_ac);
}


#  Arg [1]    : StatementHandle $sth
#  Example    : none
#  Description: PROTECTED implementation of abstract superclass method.
#               responsible for the creation of ProteinFeatures
#  Returntype : listref of Bio::EnsEMBL::ProteinFeatures
#  Exceptions : none
#  Caller     : internal
#  Status     : At Risk

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my($dbID, $translation_id, $start, $end,
     $hstart, $hend, $hid, $hdesc,
     $analysis_id, $score, $evalue, $perc_id, 
     $desc, $ilabel, $interpro_ac);

  $sth->bind_columns(\$dbID, \$translation_id, \$start, \$end, 
                     \$hstart, \$hend, \$hid, \$hdesc,
                     \$analysis_id, \$score, \$evalue, \$perc_id,
                     \$desc, \$ilabel, \$interpro_ac);

  my $analysis_adaptor = $self->db->get_AnalysisAdaptor();

  my @features;
  while($sth->fetch()) {
    my $analysis = $analysis_adaptor->fetch_by_dbID($analysis_id);

    push( 
      @features,
        my $feat = Bio::EnsEMBL::ProteinFeature->new(
           -DBID         => $dbID,
           -ADAPTOR      => $self,
           -SEQNAME      => $translation_id,
           -START        => $start,
           -END          => $end,
           -ANALYSIS     => $analysis,
           -PERCENT_ID   => $perc_id,
           -P_VALUE      => $evalue,
           -SCORE        => $score,
           -HSTART       => $hstart,
           -HEND         => $hend,
           -HSEQNAME     => $hid,
           -HDESCRIPTION => $hdesc,
           -IDESC        => $desc,
           -ILABEL       => $ilabel,
           -INTERPRO_AC  => $interpro_ac));

  }
  return \@features;
}

1;

