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

Bio::EnsEMBL::DBSQL::Feature_Obj - MySQL database adapter class for EnsEMBL Feature Objects

=head1 SYNOPSIS

  use Bio::EnsEMBL::DBSQL::Obj;
  use Bio::EnsEMBL::DBSQL::Protein_Feature_Adaptor;

  $db = new Bio::EnsEMBL::DBSQL::Obj( -user => 'root', -db => 'pog' , -host => 'caldy' , -driver => 'mysql' );
  my $feature_obj=Bio::EnsEMBL::Protein_Feature_Adaptor->new($obj);

  
=head1 DESCRIPTION

This object deals with protein feature objects. It contains methods to fetch prtein features from the database, write protein feature into the database, and delete these protein features from the databases.
A protein feature is linked to a peptide. This linked is made through the translation id stored in the translation table. For example the method fetch_by_translationID will return all of the feature for this given peptide.

The Obj object represents a database that is implemented somehow (you shouldn\'t care much as long as you can get the object). 

=head1 CONTACT

Emmanuel Mongin: mongin@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are 
usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::DBSQL::Protein_Feature_Adaptor;

use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::EnsEMBL::DBSQL::BaseAdaptor
use Bio::EnsEMBL::DBSQL::BaseAdaptor;


use Bio::EnsEMBL::Protein_FeaturePair;
use Bio::EnsEMBL::SeqFeature;


@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

sub _feature_obj {
    my($self,$dbobj) = @_;
    if( !defined $self->{'feature_obj'}) {
	my $feat_obj = $self->db->get_Feature_Obj;
	$self->{'feature_obj'} = $feat_obj;
    }
    return $self->{'feature_obj'};
}

=head2 fetch_by_translationID

 Title   : fetch_by_translationID
 Usage   :@features = $prot_feat-> fetch_by_translationID($transl_id)
 Function:Get all of the protein feature objects of one peptide (this identifiant for this peptide is given by the transaltion id)
 Example :
 Returns : Protein feature objects
 Args    :


=cut



sub fetch_by_translationID {
    my($self,$transl) = @_;

    my @features;
    my $sth = $self->prepare ("select * from protein_feature where translation = '$transl'");
    my $res = $sth->execute;

    while (my $rowhash = $sth->fetchrow_hashref) {   
	my $feature = $self->_set_protein_feature($rowhash);
 
	push(@features,$feature);
	
    }
    return @features;    
}

=head2 fetch_by_dbID

 Title   : fetch_by_dbID
 Usage   :$feature = $prot_feat->fetch_by_dbID($id)
 Function:Get a protein feature object
 Example :
 Returns :Protein feature object 
 Args    :


=cut

sub fetch_by_dbID{
   my ($self,$protfeat_id) = @_;
   
   my $features;
   my $sth = $self->prepare ("select * from protein_feature where id = '$protfeat_id'");
   my $res = $sth->execute;
   
   my $rowhash = $sth->fetchrow_hashref;
    
   if (!defined $rowhash->{'id'}) {
       $self->throw("This dbID: $protfeat_id, does not exist in the database");
   }

   my $feature = $self->_set_protein_feature($rowhash);
   
   return $feature;
}



=head2 write_Protein_feature

 Title   : write_Protein_feature
 Usage   : $obj->write_Protein_feature($feature)
 Function: writes a protein feature object           
 Example :
 Returns : 
 Args    :


=cut

sub write_Protein_feature{
    my ($self,$feature) = @_;
    my $analysis;

   
    if( ! $feature->isa('Bio::EnsEMBL::SeqFeatureI') ) {
	$self->throw("Feature $feature is not a feature!");
    }
	
    
    eval {
	$feature->validate_prot_feature();
    };

    if (!defined($feature->analysis)) {
	$self->throw("Feature " . $feature->seqname . "doesn't have analysis. Can't write to database");
    } else {
	$analysis = $feature->analysis;
    }
    
    
    my $analysisid = $self->_feature_obj->write_Analysis($analysis);
  
    my $homol = $feature->feature2;
      
    my $sth = $self->prepare(  "insert into protein_feature(id,translation,seq_start,seq_end,analysis,hstart,hend,hid,score,perc_id,evalue) ".
			       "values ('NULL',"
			       ."'".$feature->seqname    ."',"
			       .$feature->start          .","
			       .$feature->end            .","
			       .$analysisid              .","
			       .$homol->start            .","
			       .$homol->end              .","
			       ."'".$homol->seqname      ."',"
			       .$feature->score         .","
			       .$feature->percent_id    .","
			       ."'".$feature->p_value   ."')");
    $sth->execute();

}

=head2 write_Protein_feature_by_translationID

 Title   : write_Protein_feature_by_translationID
 Usage   :$obj->write_Protein_feature_by_translation($pep,@features)
 Function: Write all of the protein features into the database of a particular peptide
 Example :
 Returns : nothing
 Args    :


=cut

sub write_Protein_feature_by_translationID {
    my ($self,$pep,@features) = @_;
    
    my $analysis;
   
#Check if the translation id exist in the database, throw an exeption if not.
    my $sth1 = $self->prepare("select id from translation where id = '$pep'");
    $sth1->execute;
   
    if ($sth1->rows == 0) {
	$self->throw("This translation id: $pep does not exist in the database");
    }

    FEATURE :
	foreach my $features(@features) {	
	   
	    if( ! $features->isa('Bio::EnsEMBL::SeqFeatureI') ) {
		$self->throw("Feature $features is not a feature!");
	    }
	    
	    eval {
		$features->validate_prot_feature();
	    };
	    
	   
	    
	    if ($@) {
		print STDERR "Feature for peptide ". $features->seqname." is not a protein feature, skipped\n";
		next FEATURE;
	    }
	    
	    
	    if (!defined($features->analysis)) {
		$self->throw("Feature " . $features->seqname . "doesn't have analysis. Can't write to database");
	    } else {
		$analysis = $features->analysis;
	    }
	    
	    my $analysisid = $self->_feature_obj->write_Analysis($analysis);
	    
	    if ( $features->isa('Bio::EnsEMBL::FeaturePair') ) {
		my $homol = $features->feature2;
		
		my $sth = $self->prepare(  "insert into protein_feature(id,translation,seq_start,seq_end,analysis,hstart,hend,hid,score,perc_id,evalue) ".
					   "values ('NULL',"
					   ."'".$pep                 ."',"
					   .$features->start          .","
					   .$features->end            .","
					   .$analysisid              .","
					   .$homol->start            .","
					   .$homol->end              .","
					   ."'".$homol->seqname      ."',"
					   .$features->score         .","
					   .$features->percent_id    .","
					   .$features->p_value       .")");
		$sth->execute();
	    }
	}
}
    


=head2 delete_by_translationID

 Title   : delete_by_translationID
 Usage   :
 Function: deletes all protein features for a particular peptide
 Example :
 Returns : 
 Args    :


=cut

sub delete_by_translationID {
    my ($self,$trans) = @_;
    my $sth = $self->prepare("delete from protein_feature where translation = '$trans'");

    $sth->execute;
}

=head2 delete_by_dbID

 Title   : delete_by_dbID
 Usage   :
 Function: deletes a protein feature
 Example :
 Returns : 
 Args    :


=cut

sub delete_by_dbID {
    my ($self,$id) = @_;
    my $sth = $self->prepare("delete from protein_feature where id = $id");

    $sth->execute;
}

=head2 get_interproac_by_signature_id

 Title   : get_interproac_by_signature_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_interproac_by_signature_id{
   my ($self,$sign_id) = @_;
   my $sth = $self->prepare("select interpro_ac from interpro where id = '$sign_id'");
   
   $sth->execute;
   
   my $interpro_ac = $sth->fetchrow;
   
   my $feat = new Bio::EnsEMBL::SeqFeature;
   
   $feat->external_db($interpro_ac);
   
   return $feat;

}


=head2 _set_protein_feature

 Title   : _set_protein_feature
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _set_protein_feature{
   my ($self,$rowhash) = @_;

my $analysis = $self->_feature_obj->get_Analysis($rowhash->{'analysis'});

#my $analysis = $self->_feature_obj->get_Analysis(7);

   
   my $feat1 = new Bio::EnsEMBL::SeqFeature ( -seqname => $rowhash->{'translation'},
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
   
   my $feature = new Bio::EnsEMBL::Protein_FeaturePair(-feature1 => $feat1,
					       -feature2 => $feat2,);
   return $feature;

}




















