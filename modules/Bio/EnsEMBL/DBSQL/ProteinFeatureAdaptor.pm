#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::Feature_Obj
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::Feature_Obj - MySQL database adapter class for EnsEMBL 
Feature Objects

=head1 SYNOPSIS

  use Bio::EnsEMBL::DBSQL::DBAdaptor;
  use Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor;

  $db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -user => 'root', -dbname => 'pog' , -host => 'caldy' , -driver => 'mysql' );
  my $feature_obj=Bio::EnsEMBL::ProteinFeatureAdaptor->new($obj);

  
=head1 DESCRIPTION

This object deals with protein feature objects. It contains methods to fetch 
prtein features from the database, write protein feature into the database, 
and delete these protein features from the databases.

A protein feature is linked to a peptide. This linked is made through the 
translation id stored in the translation table. For example the method 
fetch_by_translationID will return all of the feature for this given peptide.

The Obj object represents a database that is implemented somehow 
(you shouldn\'t care much as long as you can get the object). 

=head1 CONTACT

Emmanuel Mongin: mongin@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal 
methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor;

use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::EnsEMBL::DBSQL::BaseAdaptor
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use Bio::EnsEMBL::ProteinFeature;
use Bio::EnsEMBL::SeqFeature;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);



=head2 fetch_all_by_translation_id

  Arg [1]    : int $transl
               the internal id of the translation corresponding to protein 
               whose features are desired 
  Example    : @prot_feats = $prot_feat_adaptor->fetch_by_translation_id(1234);
  Description: Get all of the protein feature objects of a peptide             
  Returntype : listref of Bio::EnsEMBL::ProteinFeatures
  Exceptions : none
  Caller     : ?

=cut

sub fetch_by_translation_id {
  my($self,$transl) = @_;

  my @features;
  my $analysis_adaptor = $self->db()->get_AnalysisAdaptor();
  
  
  my $sth = $self->prepare(
	      "SELECT p.seq_start, p.seq_end, p.analysis_id, 
                      p.score, p.perc_ident, p.evalue, 
                      p.hit_start, p.hit_end, p.hit_id, 
                      x.display_label, i.interpro_ac
               FROM protein_feature p,analysis a 
               LEFT JOIN interpro AS i ON p.hit_id = i.id
               LEFT JOIN xref AS x ON x.dbprimary_acc = i.interpro_ac
               WHERE p.translation_id = ? 
               AND p.analysis_id = a.analysis_id 
               AND a.gff_feature = ?");

  $sth->execute("$transl", 'domain');

  while( my $row = $sth->fetchrow_arrayref) {
    my ($start, $end, $analysisid, $score, $perc_id, $evalue, $hstart,
	$hend,$hid,$desc, $interpro_ac) = @$row;
    
    my $analysis = $analysis_adaptor->fetch_by_dbID($analysisid);
    
    my $feat = Bio::EnsEMBL::ProteinFeature->new();
    $feat->seqname($transl);
    $feat->start($start);
    $feat->end($end);
    $feat->analysis($analysis);
    $feat->percent_id($perc_id);
    $feat->p_value($evalue);
    
    $feat->hstart($hstart);
    $feat->hend($hend);
    $feat->hseqname($hid);
    
    $feat->idesc($desc);
    $feat->interpro_ac($interpro_ac);

    push(@features,$feat);
  }
  
  $sth = $self->prepare ("SELECT p.seq_start, p.seq_end, p.analysis_id, 
                                   p.score, p.perc_ident, p.evalue, 
                                   p.hit_start, p.hit_end, p.hit_id 
                            FROM   protein_feature p, analysis a 
                            WHERE  a.analysis_id = p.analysis_id 
                                   AND p.translation_id = ? 
                                   AND a.gff_feature != ?");

  $sth->execute("$transl", 'domain');

  while( my $row = $sth->fetchrow_arrayref) {
    my ($start,$end,$analysisid,$score,$perc_id,
	$evalue,$hstart,$hend,$hid) = @$row;
    
    my $analysis = $analysis_adaptor->fetch_by_dbID($analysisid);
    
    
    my $feat = Bio::EnsEMBL::ProteinFeature->new();
    $feat->seqname($transl);
    $feat->start($start);
    $feat->end($end);
    $feat->analysis($analysis);
    $feat->percent_id($perc_id);
    $feat->p_value($evalue);
    
    $feat->hstart($hstart);
    $feat->hend($hend);
    $feat->hseqname($hid);
    
    
    push(@features,$feat);	
  }
  
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

=cut

sub fetch_by_dbID{
   my ($self,$protfeat_id) = @_;
   
   my $features;
   my $sth = 
     $self->prepare("select * from protein_feature where id = '$protfeat_id'");
   my $res = $sth->execute;
   
   my $rowhash = $sth->fetchrow_hashref;
    
   if (!defined $rowhash->{'id'}) {
       $self->throw("This dbID: $protfeat_id, does not exist in the database");
   }

   my $feature = $self->_set_protein_feature($rowhash);
   
   return $feature;
}



=head2 fetch_by_feature_and_dbID

  Arg [1]    : string $feature
               analysis gff_source
  Arg [2]    : int $transl
               unique translation identifier
  Example    : none
  Description: Retrieves protein features of a protein using the dbID of the
               corresponding translation, and the analysis gff_source for 
               the feature
  Returntype : listref of Bio::EnsEMBL::ProteinFeatures
  Exceptions : none
  Caller     : ?

=cut

sub fetch_all_by_feature_and_dbID{
    my ($self,$feature,$transl) = @_;
    my @features;
    my $analysis_adaptor = $self->db()->get_AnalysisAdaptor();

    #The call has to be specific for the Interpro components because there 
    #is one join to make on the interpro table and then with the xref table
    if (($feature eq "PRINTS") || ($feature eq "Pfam") || 
	($feature eq "PROSITE") || ($feature eq "PROFILE")) {

      my $sth = $self->prepare("SELECT p.seq_start, p.seq_end, p.analysis_id, 
                                       p.score, p.perc_ident, p.evalue, 
                                       p.hit_start, p.hit_end, p.hit_id, 
                                       x.display_label 
                                FROM protein_feature p,analysis a 
                                      left join interpro as i on p.hit_id = i.id
                                      left join xref as x on x.dbprimary_acc = i.interpro_ac
                                WHERE p.translation_id = '$transl' 
                                      AND p.analysis_id = a.analysis_id 
                                      AND a.gff_source = '$feature'");
      
      $sth->execute();
            
      while( my $arrayref = $sth->fetchrow_arrayref) {
	my ($start,$end,$analysisid,$score,$perc_id,$evalue,
	    $hstart,$hend,$hid,$desc) = @{$arrayref};

	my $analysis = $analysis_adaptor->fetch_by_dbID($analysisid);

	my $feat1 = new Bio::EnsEMBL::SeqFeature ( -seqname => $transl,
						   -start => $start,
						   -end => $end,
						   -score => $score, 
						   -analysis => $analysis,
						   -percent_id => $perc_id,
						   -p_value => $evalue);
	    
	my $feat2 = new Bio::EnsEMBL::SeqFeature (-start => $hstart,
						  -end => $hend,
						  -analysis => $analysis,
						  -seqname => $hid);
	
	my $feature = new Bio::EnsEMBL::ProteinFeature(-feature1 => $feat1,
						       -feature2 => $feat2,);
	
	$feature->idesc($desc);
	
	if ($feature) {
	  push(@features,$feature);
	}
	
      }
    }

    #Superfamily has also a description attached to it but there is 
    #no join needed with the interpro table. but its also considered 
    #as a domain feature
    elsif ($feature eq "superfamily") {
      my $sth = $self->prepare ("SELECT p.seq_start, p.seq_end, p.analysis_id,
                                        p.score, p.perc_ident, p.evalue, 
                                        p.hit_start, p.hit_end, p.hit_id, 
                                        x.display_label, x.dbprimary_acc 
                                 FROM protein_feature AS p, analysis AS a, 
                                      xref AS x 
                                 WHERE a.gff_source = '$feature' AND 
                                       p.translation_id = '$transl' AND 
                                       a.analysis_id = p.analysis_id AND 
                                       x.dbprimary_acc = p.hit_id");
      $sth->execute();
      
      
      while( my $arrayref = $sth->fetchrow_arrayref) {
	    
	my ($start,$end,$analysisid,$score,$perc_id,$evalue,
	    $hstart,$hend,$hid,$desc,$interproac) = @{$arrayref};

	my $analysis = $analysis_adaptor->fetch_by_dbID($analysisid);

	my $feat1 = new Bio::EnsEMBL::SeqFeature ( -seqname => $transl,
						   -start => $start,
						   -end => $end,
						   -score => $score, 
						   -analysis => $analysis,
						   -percent_id => $perc_id,
						   -p_value => $evalue);
	    
	my $feat2 = new Bio::EnsEMBL::SeqFeature (-start => $hstart,
						  -end => $hend,
						  -analysis => $analysis,
						  -seqname => $hid);
	
	my $feature = new Bio::EnsEMBL::ProteinFeature(-feature1 => $feat1,
						       -feature2 => $feat2,);
	
	$feature->idesc($desc);
	$feature->interpro_ac($interproac);
	
	if ($feature) {
	  push(@features,$feature);
	}
	    
      }
    }
    
    #Get all of the other features...Coils, TMHMM, ...
    else {
      my $sth = $self->prepare ("SELECT p.seq_start, p.seq_end, p.analysis_id,
                                        p.score, p.perc_ident, p.evalue, 
                                        p.hit_start, p.hit_end, p.hit_id 
                                 FROM protein_feature p,analysis a 
                                 WHERE a.analysis_id = p.analysis_id 
                                       AND p.translation_id = '$transl' 
                                       AND a.gff_feature != 'domain' 
                                       AND a.gff_source = '$feature'");
      
      $sth->execute();
      my @a = $sth->fetchrow();

      $sth->execute();
      while( my $arrayref = $sth->fetchrow_arrayref) {
	my ($start, $end, $analysisid, $score, $perc_id, 
	    $evalue, $hstart,$hend,$hid) = @{$arrayref};

	my $analysis = $analysis_adaptor->fetch_by_dbID($analysisid);
	    
	my $feat1 = new Bio::EnsEMBL::SeqFeature ( -seqname => $transl,
						   -start => $start,
						   -end => $end,
						   -score => $score, 
						   -analysis => $analysis,
						   -percent_id => $perc_id,
						   -p_value => $evalue);
	
	my $feat2 = new Bio::EnsEMBL::SeqFeature (-start => $hstart,
						  -end => $hend,
						  -analysis => $analysis,
						  -seqname => $hid);
	
	my $feature = new Bio::EnsEMBL::ProteinFeature(-feature1 => $feat1,
						       -feature2 => $feat2);
	    
	if ($feature) {
	  push(@features,$feature);
	}
      }
    }	
    
    return \@features;
}



=head2 store

  Arg [1]    : Bio::EnsEMBL::ProteinFeature $feature
  Example    : $protein_feature_adaptor->store($protein_feature);
  Description: Stores a protein feature in the database
  Returntype : none
  Exceptions : thrown if arg is not a Bio::EnsEMBL::SeqFeatureI
  Caller     : none

=cut

sub store {
    my ($self,$feature) = @_;
    my $analysis;

   
    if( ! $feature->isa('Bio::EnsEMBL::FeaturePair') ) {
	$self->throw("A FeaturePair object required!");
    }
	
    
    eval {
	$feature->validate_prot_feature();
    };

    if (!defined($feature->analysis)) {
	$self->throw("Feature " . $feature->seqname . 
		     "doesn't have analysis. Can't write to database");
    } else {
	$analysis = $feature->analysis;
    }
    
    
    my $analysisid = $self->db->get_AnalysisAdaptor->store($analysis);
  

    my $homol = $feature->feature2;
      
    my $sth = 
      $self->prepare("INSERT INTO protein_feature(protein_feature_id, 
                                                  translation_id, seq_start,
                                                  seq_end, analysis_id,
                                                  hit_start, hit_end, 
                                                  hit_id, score,
                                                  perc_ident, evalue) ".
			       "VALUES ('NULL',"
			       ."'".$feature->seqname    ."',"
			       .$feature->start          .","
			       .$feature->end            .","
			       .$analysisid              .","
			       .$homol->start            .","
			       .$homol->end              .","
			       ."'".$homol->seqname      ."',"
			       .$feature->score         .","
			       .$feature->percent_id    .","
			       .$feature->p_value   .")");
    
    $sth->execute();

}

=head2 remove

  Arg [1]    : int $id 
               the database id of the protein feature to remove
  Example    : $protein_feature_adaptor->remove($prot_feat->dbID());
  Description: deletes a protein feature from the database
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub remove {
    my ($self,$id) = @_;
    my $sth = $self->prepare("delete from protein_feature where id = $id");

    $sth->execute;
}

=head2 _set_protein_feature

  Arg [1]    : DBI row hash $rowhash
  Example    : $prot_feat = $self->_set_protein_feature($row_hash);
  Description: PRIVATE creates a ProteinFeature object from a SQL hashref
  Returntype : Bio::EnsEMBL::ProteinFeature
  Exceptions : none
  Caller     : internal

=cut

sub _set_protein_feature {
   my ($self,$rowhash) = @_;

   my $analysis = 
     $self->db->get_AnalysisAdaptor->fetch_by_dbID($rowhash->{'analysis'});
   
   my $feat1 =
     new Bio::EnsEMBL::SeqFeature ( -seqname => $rowhash->{'translation'},
				    -start => $rowhash->{'seq_start'},
				    -end => $rowhash->{'seq_end'},
				    -score => $rowhash->{'score'}, 
				    -analysis => $analysis,
				    -percent_id => $rowhash->{'perc_id'},
				    -p_value => $rowhash->{'evalue'});
   
   my $feat2 = new Bio::EnsEMBL::SeqFeature (-start => $rowhash->{'hstart'},
					     -end => $rowhash->{'hend'},
					     -analysis => $analysis,
					     -seqname => $rowhash->{'hid'});
   
   my $feature = new Bio::EnsEMBL::ProteinFeature(-feature1 => $feat1,
						  -feature2 => $feat2,);
   return $feature;
}


1;







