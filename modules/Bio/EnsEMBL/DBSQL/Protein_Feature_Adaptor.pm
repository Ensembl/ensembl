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

This object deals with protein feature objects.

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

# Object preamble - inheriets from Bio::Root::Object
use Bio::EnsEMBL::DBSQL::BaseAdaptor;


use DBI;
use Bio::EnsEMBL::FeaturePair;
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

=head2 Fetch_protfeature_by_translation

 Title   : Fetch_protfeature_by_translation
 Usage   :@features = $prot_feat-> Fetch_protfeature_by_translation($transl_id)
 Function:Get all of the protein feature objects of one peptide
 Example :
 Returns : Protein feature objects
 Args    :


=cut



sub Fetch_protfeature_by_translation {
    my($self,$transl) = @_;

    my @features;
    my $sth = $self->prepare ("select * from protein_feature where translation = '$transl'");
    my $res = $sth->execute;

    while (my $rowhash = $sth->fetchrow_hashref) {
	
	    my $analysis = $self->_feature_obj->get_Analysis($rowhash->{'analysis'});
	
	    my $feat1 = new Bio::EnsEMBL::SeqFeature ( -start => $rowhash->{'seq_start'},
						       -end => $rowhash->{'seq_end'},
						       -score => $rowhash->{'score'}, 
						       -analysis => $analysis,
						       -seqname => $transl);

	    my $feat2 = new Bio::EnsEMBL::SeqFeature (-start => $rowhash->{'hstart'},
						      -end => $rowhash->{'hend'},
						      -analysis => $analysis,
						      -seqname => $rowhash->{'hid'});

	    my $feature = new Bio::EnsEMBL::FeaturePair(-feature1 => $feat1,
							-feature2 => $feat2,);
	    
	    push(@features,$feature);
	
    }
    return @features;    
}

=head2 Fetch_protfeature_by_id

 Title   : Fetch_protfeature_by_id
 Usage   :$feature = $prot_feat-> Fetch_protfeature_by_translation($feature_id)
 Function:Get a protein feature object
 Example :
 Returns :Protein feature object 
 Args    :


=cut

sub Fetch_protfeature_by_id{
   my ($self,$protfeat_id) = @_;
   
   my $features;
   my $sth = $self->prepare ("select * from protein_feature where id = '$protfeat_id'");
   my $res = $sth->execute || die "Can't execute";
   
   my $rowhash = $sth->fetchrow_hashref;
       
   my $analysis = $self->_feature_obj->get_Analysis($rowhash->{'analysis'});
   
   my $feat1 = new Bio::EnsEMBL::SeqFeature ( -start => $rowhash->{'seq_start'},
					      -end => $rowhash->{'seq_end'},
					      -score => $rowhash->{'score'}, 
					      -analysis => $analysis,
					      -seqname => $rowhash->{'translation'});
   
   my $feat2 = new Bio::EnsEMBL::SeqFeature (-start => $rowhash->{'hstart'},
					     -end => $rowhash->{'hend'},
					     -analysis => $analysis,
					     -seqname => $rowhash->{'hid'});
   
   my $feature = new Bio::EnsEMBL::FeaturePair(-feature1 => $feat1,
					       -feature2 => $feat2,);
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
    
    my $sth = $self->prepare(  "insert into protein_feature(id,seq_start,seq_end,score,analysis,translation,hstart,hend,hid) ".
			       "values ('NULL',"
			       .$feature->start          .","
			       .$feature->end            .","
			       .$feature->score          .","
			       .$analysisid                .","
			       ."'".$feature->seqname        ."',"
			       .$homol->start            .","
			       .$homol->end              .","
			       ."'".$homol->seqname      ."')");
		


    $sth->execute();
    
}

=head2 write_Protein_feature_by_translation

 Title   : write_Protein_feature_by_translation
 Usage   :$obj->write_Protein_feature_by_translation($pep,@features)
 Function: Write all of the protein features into the database of a particular peptide
 Example :
 Returns : nothing
 Args    :


=cut

sub write_Protein_feature_by_translation {
    my ($self,$pep,@features) = @_;
    
    my $analysis;
    
    FEATURE :
	foreach my $features(@features) {	
	
	    if( ! $features->isa('Bio::EnsEMBL::SeqFeatureI') ) {
		$self->throw("Feature $features is not a feature!");
	    }
	    
	    eval {
		$features->validate_prot_feature();
	    };
	    
	   
	    
	    if ($@) {
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
		
		my $sth = $self->prepare(  "insert into protein_feature(id,seq_start,seq_end,score,analysis,translation,hstart,hend,hid) ".
					   "values ('NULL',"
					   .$features->start          .","
					   .$features->end            .","
					   .$features->score          .","
					   .$analysisid              .","
					   ."'".$pep                 ."',"
					   .$homol->start            .","
					   .$homol->end              .","
					   ."'".$homol->seqname      ."')");
		$sth->execute();
	    }
	}
}
    


=head2 delete

 Title   : delete
 Usage   :
 Function: deletes all protein features for a particular peptide
 Example :
 Returns : 
 Args    :


=cut

sub delete {
    my ($self,$trans) = @_;
    my $sth = $self->prepare("delete from protein_feature where translation = '$trans'");

    my $res = $sth->execute;
}

=head2 delete_by_id

 Title   : delete_by_id
 Usage   :
 Function: deletes a protein feature
 Example :
 Returns : 
 Args    :


=cut

sub delete_by_id {
    my ($self,$id) = @_;
    my $sth = $self->prepare("delete from protein_feature where id = $id");

    my $res = $sth->execute;
}








