
#
# BioPerl module for Protein_Adaptor
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Protein_Adaptor - DESCRIPTION of Object

=head1 SYNOPSIS
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBSQL::Protein_Adaptor;

$db = new Bio::EnsEMBL::DBSQL::Obj( -user => 'root', -db => 'pog' , -host => 'caldy' , -driver => 'mysql' );
my $protein_adaptor=Bio::EnsEMBL::Protein_Adaptor->new($obj);


=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::Protein_Adaptor;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::Root::Object;
use Bio::EnsEMBL::Protein;
use Bio::EnsEMBL::DBSQL::Protein_Feature_Adaptor;
use Bio::EnsEMBL::ExternalData::Family::FamilyAdaptor;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

sub _gene_obj {
    my($self) = @_;
    if( !defined $self->{'gene_obj'}) {
	my $feat_obj = $self->db->get_Gene_Obj;
	$self->{'gene_obj'} = $feat_obj;
    }
    
    return $self->{'gene_obj'};
}

sub _protfeat_obj {
    my($self) = @_;
    if( !defined $self->{'protfeat_obj'}) {
	my $feat_obj = $self->db->get_Protfeat_Adaptor;
	$self->{'protfeat_obj'} = $feat_obj;
    }
    
    return $self->{'protfeat_obj'};
}

sub _familyAdaptor {
   my($self) = @_;
    if( !defined $self->{'familyAdaptor'}) {
	my $feat_obj = $self->db->get_FamilyAdaptor;
	$self->{'familyAdaptor'} = $feat_obj;
    }
    
    return $self->{'familyAdaptor'};
}
 
   

=head2 fetch_Protein_by_dbid

 Title   : fetch_Protein_by_dbid
 Usage   :$obj->fetch_Protein_by_dbid($id)
 Function:Built a whole Protein object for a given protein id
 Example :
 Returns : Protein Object
 Args    :Protein id (ENSPXXXX)


=cut

sub fetch_Protein_by_dbid{
   my ($self,$id) = @_;

   my $query = "select id from transcript where translation = '$id'";
   my $sth = $self->prepare($query);
   $sth ->execute();
   my $transid = $sth->fetchrow;
   
   my $query1 = "select g.created,g.modified from gene as g, transcript as t where t.gene = g.id and t.translation = '$id'";
   my $sth1 = $self->prepare($query1);
   $sth1 ->execute();
     
   my @date = $sth1->fetchrow;
   my $created = $date[0];
   my $modified = $date[1];


   my $transcript = $self->fetch_Transcript_by_dbid($transid);
   my @dblinks = $self->fetch_DBlinks_by_dbid($id);
   my @prot_feat = $self->fetch_Protein_features_by_dbid($id);
   #my $family = $self->fetch_Family_by_dbid($id);
   
   my $sequence = $transcript->translate->seq;
   
   my $moltype = "protein";
   
   #This has to be changed, the description may be take from the protein family description line
   my $desc = "Protein predicted by Ensembl";

   my $protein = Bio::EnsEMBL::Protein->new ( -seq =>$sequence,
					      -accession_number  => $id,
					      -display_id => $id,
					      -primary_id => $id,
					      -id => $id,
					      -desc => $desc,
					      -moltype => $moltype,
					      );

   my $ann  = Bio::Annotation->new;

   $protein ->annotation($ann);
   $protein->add_date($created);
   $protein->add_date($modified);

   foreach my $link (@dblinks) {
       $protein->annotation->add_DBLink($link);
   }
   my %seen;
#Get Interpro data and make with it a dblink object
   foreach my $feat (@prot_feat) {
       

       my $pfam = $feat->hseqname;
       my $query2 = "select interpro_ac from interpro where id = '$pfam'";
       my $sth2 = $self->prepare($query2);
       $sth2 ->execute();
       my $interpro = $sth2->fetchrow;
       
       

       if (! defined ($seen{$interpro})) {
	   my $dblink = Bio::Annotation::DBLink->new();
	   $dblink->database('InterPro');
	   $dblink->primary_id($interpro);
	   $protein->annotation->add_DBLink($dblink);
	   $seen{$interpro} = 1;
       }
   }


   foreach my $feat (@prot_feat) {
       
       $protein->add_Protein_feature($feat);
   }
   return $protein;
}

=head2 fetch_Transcript_by_dbid

 Title   : fetch_Transcript_by_dbid
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_Transcript_by_dbid{
   my ($self,$transcript_id) = @_;
   
   my $transcript = $self->_gene_obj->get_Transcript($transcript_id);
         
   return $transcript;

}


=head2 fetch_DBlinks_by_dbid

 Title   : fetch_DBlinks
 Usage   :$prot_adaptor->fetch_DBlinks_by_dbid($protein_id)
 Function:Get all of the DBlinks for one protein given a protein object and attach them to the object
 Example :
 Returns :an array of dblinks
 Args    :Protein id


=cut



sub fetch_DBlinks_by_dbid {
    my($self,$protein_id) = @_;
    my @links;

    my $query = "select dbl.external_db,dbl.external_id from transcriptdblink as dbl, transcript as t where dbl.transcript_id = t.id and t.translation = '$protein_id'";
    my $sth = $self->prepare($query);
    $sth ->execute();
    while( (my $hash = $sth->fetchrow_hashref()) ) {
	
	my $dblink = Bio::Annotation::DBLink->new();
	$dblink->database($hash->{'external_db'});
	$dblink->primary_id($hash->{'external_id'});
	

	push(@links,$dblink);
    }
    return @links;
}


=head2 fetch_Protein_features_by_dbid

 Title   : fetch_Protein_features_by_dbid
 Usage   :$prot_adaptor->fetch_Protein_features_by_dbid($protein_id)
 Function:Get all of the protein features for a given protein object and attach it to the object
 Example :
 Returns :nothing
 Args    :Protein Object


=cut

sub fetch_Protein_features_by_dbid{
   my ($self,$protein_id) = @_;

     
#This call a method contained in Protein_Feature_Adaptor, returns feature objects for the given protein
   my @features = $self->_protfeat_obj->fetch_by_translationID($protein_id);

   return @features;

}

=head2 fetch_Family_by_dbid

 Title   : fetch_Family_by_dbid
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_Family_by_dbid{
   my ($self,$protein_id) = @_;
   
   #This call a method contained in FamilyAdaptor, perhaps we should one day put all of these objects together; a big protein object

   my $family = $self->_familyAdaptor->get_Family_of_Ensembl_pep_id($protein_id);

   return $family;

}



=head2 fetch_by_DBlink

 Title   : fetch_by_DBlink
 Usage   :my @proteins = $obj->fetch_by_DBlink($dblinkid)
 Function:Get the proteins corresponding to the given DBlink. In most of the case only one protein will be returned
 Example :
 Returns : an array of protein objects being linked to this given DBlink
 Args    :Dblink id (external_id in transcriptdblink table)


=cut

sub fetch_by_DBlink{
   my ($self,$dblink) = @_;
   my @proteins;
   my $query = "select t.translation from transcript as t, transcriptdblink as tdb where tdb.external_id = '$dblink' and tdb.transcript_id = t.id";
    my $sth = $self->prepare($query);
    $sth ->execute();
    while( (my $pepid = $sth->fetchrow) ) {
	my $pep = $self->fetch_Protein($pepid);

	push(@proteins,$pep);
    }
   return @proteins;
}

=head2 fetch_by_feature

 Title   : fetch_by_feature
 Usage   :my @proteins = $obj->fetch_by_feature($feature_id)
 Function:This method should be in theory in the Protein_Feature_Adaptor object but has been placed here for convenience (the Protein_Feature_Adaptor object may be transfered to this Protein_Adaptor object...to be discussed). The function of this method is to retrieve all of the proteins which have a given feature (retrieved by feature id)
 Example :
 Returns : An Array of protein objects
 Args    : feature id (hid in the protein_feature table)


=cut

sub fetch_by_feature{
   my ($self,$feature) = @_;
   my @proteins;
   
   my $query = "select translation from protein_feature where hid = '$feature'";
   my $sth = $self->prepare($query);
   $sth ->execute();
   while( (my $pepid = $sth->fetchrow) ) {
       my $pep = $self->fetch_Protein;
       push(@proteins,$pepid);
   }
   return @proteins;
}

=head2 fetch_by_array_feature

 Title   : fetch_by_hid_array_feature
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_by_array_feature{
   my ($self,@feature) = @_;

   my $nb = scalar @feature;
   my %seen;
   my @result;

   if (@feature) {
       foreach my $dbl(@feature) {
	   my @protein_linked = $self->fetch_by_feature($dbl);
	   foreach my $prot (@protein_linked) {
	       if ($seen{$prot}) {
		   my $count = $seen{$prot};
		   $count = $count++;
		   $seen{$prot} = $count;
	       }
	       else {
		   $seen{$prot} = 1;
	       }
	   }
       }
   }
   foreach my $key (keys (%seen)) {
       if ($seen{$key} = $nb) {
	   push (@result, $key);
       }
   }
   
   return @result;
   
}







