
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

my $protein = $protein_adaptor->fetch_Protein_by_dbid;

If the SNPs are wanted:

my $snpdb = Bio::EnsEMBL::ExternalData::SNPSQL::DBAdapter->new(
							       -dbname=>$snpdbname,
							       -user=>$dbuser,
							       -host=>$host);    
 
$protein_adaptor->snp_obj($snpdb);


=head1 DESCRIPTION

This Object inherit from BaseAdaptor following the new adaptor rules. It also has pointers to 3 different objects: Obj.pm, Gene_Obj.pm and FamilyAdaptor.pm (which is not currently used). This pointers allow the object to use the methods contained in these objects. SNPs db adaptor may also be added in these pointers.
The main method may be fetch_Protein_by_dbid, which return a complete protein object. Different methods are going to be develloped in this object to allow complexe queries at the protein level in Ensembl.

=head1 CONTACT

mongin@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::Protein_Adaptor;
use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Protein;
use Bio::EnsEMBL::DBSQL::Protein_Feature_Adaptor;
use Bio::Species;
use Bio::EnsEMBL::DBSQL::DBEntryAdaptor;
use Bio::EnsEMBL::Utils::Eprof qw(eprof_start eprof_end);

# this is required at the right point if needed
#use Bio::EnsEMBL::ExternalData::GeneSNP;

#use Bio::EnsEMBL::ExternalData::Family::FamilyAdaptor;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

#The pointer are defined bellow.

sub _gene_obj {
    my($self) = @_;
    if( !defined $self->{'gene_obj'}) {
	my $feat_obj = $self->db->gene_Obj;
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
 
sub _dbEntryAdaptor {
     my($self) = @_;
     if( !defined $self->{'dbEntryAdaptor'}) {
	 my $dbentryadaptor = Bio::EnsEMBL::DBSQL::DBEntryAdaptor->new($self);
	 $self->{'dbEntryAdaptor'} = $dbentryadaptor;
    }
     return $self->{'dbEntryAdaptor'};
}
   
=head2 snp_obj

 Title   : snp_obj
 Usage   : $obj->snp_obj($newval)
 Function: 
 Returns : value of snp_obj
 Args    : newvalue (optional)


=cut

sub snp_obj{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'snp_obj'} = $value;
    }
    return $obj->{'snp_obj'};

}

=head2 fetch_Protein_by_transcriptId

 Title   : fetch_Protein_by_transcriptId
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_Protein_by_transcriptId{
   my ($self,$transid) = @_;
   my $query = "select translation from transcript where id = '$transid'";
   my $sth = $self->prepare($query);
   $sth->execute();
   my @row = $sth->fetchrow;
   return $self->fetch_Protein_by_dbid($row[0]);
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

#Get the transcript id from the translation id 
   my $query = "select id,gene from transcript where translation = '$id'";
   my $sth = $self->prepare($query);
   $sth ->execute();
   my @rowid = $sth->fetchrow;

   my $transid = $rowid[0];
   my $geneid = $rowid[1];

   if (!defined $transid) {
       $self->throw("$id does not have a transcript id");
   }

   if (!defined $geneid) {
       $self->throw("$id does not have a gene id");
   }

#Get the different dates (created and modified) for the corresponding gene   
   my $query1 = "select g.created,g.modified from gene as g, transcript as t where t.gene = g.id and t.translation = '$id'";
   my $sth1 = $self->prepare($query1);
   $sth1 ->execute();
     
   my @date = $sth1->fetchrow;
   my $created = $date[0];
   my $modified = $date[1];

#Add created and modified tag to the date
   ($created) = $created =~ /(\d+-\d+-\d+)/;
   $created = $created." (Created)";

   ($modified) = $modified =~ /(\d+-\d+-\d+)/;
   $modified = $modified." (Modified)";

#Get the transcript object (this will allow us to get the aa sequence of the protein
   my $transcript = $self->fetch_Transcript_by_dbid($transid);

#Get all of the Dblink for the given Peptide id
   my @dbentry = $self->fetchDBentry_by_dbID($id);

#Get all of the Protein Features for the given Protein
   my @prot_feat = $self->fetch_Protein_features_by_dbid($id);

#Get all of the family (at the Transcript level), not implemented yet
   #my $family = $self->fetch_Family_by_dbid($id);

#Get all SNPs ?????? method which would be nice to implement


#Get the aa sequence out of the transcript object   
   #my $sequence = $transcript->translate->seq;
   my $sequence = $transcript->translate->seq;

#Calculate the length of the Peptide   
   my $length = length($sequence);
   if ($length == 0) {
       $self->throw("Transcript".$transcript->id." does not have any amino acid sequence"); 
       #return 0;
   }
  
#Define the moltype
   my $moltype = "protein";

#Define the specie (here by default human, but will have to find something else when other organisms come into Ensembl) 
   my @class = ( "Eukaryota", "Metazoa", "Chordata", "Craniata", "Vertebrata", "Euteleostomi", "Mammalia", "Eutheria", "Primates", "Catarrhini", "Hominidae","Homo" ,"sapiens (human)");
   @class = reverse(@class);
   my $common;
   my $sub_species;
   my $org;

   my $species = Bio::Species->new();
   $species->classification( @class );
   $species->common_name( $common      ) if $common;
   $species->sub_species( $sub_species ) if $sub_species;
   $species->organelle  ( $org         ) if $org;

   #This has to be changed, the description may be take from the protein family description line
   my $desc = "Protein predicted by Ensembl";
  
#Create the Protein object
   my $protein = Bio::EnsEMBL::Protein->new ( -seq =>$sequence,
					      -accession_number  => $id,
					      -display_id => $id,
					      -primary_id => $id,
					      -id => $id,
					      -desc => $desc,

					      );

   $protein->transcriptac($transid);                                              
   $protein->geneac($geneid);    
   $protein->adaptor($self);	

#Add the species object to protein object
   $protein->species($species);

#Add the date of creation of the protein to the annotation object
   my $ann  = Bio::Annotation->new;
      
   $protein ->annotation($ann);
   $protein->add_date($created);
   $protein->add_date($modified);

#Give a gene name to the peptide, here , its Ensembl gene id   
#$ann->gene_name($geneid); 

#Add the DBlinks to the annotation object
   foreach my $link (@dbentry) {
       if ($link){
	   $protein->annotation->add_DBLink($link);
       }
   }
   my %seen1;
   my %seen2;

#Get Protein feature data and creat a DBlink object
   foreach my $feat (@prot_feat) {
       
#Get the the accession number of the feature matching to the given Protein
       my $featid = $feat->hseqname;

#This is supposed to be the rigth way 
       my $dbdesc = $feat->analysis->db;


#In the case of a protein feature being an Interpro signature (eg: Pfam, Prints,...) this  signature will also be added to the object as a dblink (this is a requirement to dump peptides in SP format and to use their automatic annotation)

#If protein feature is an Interpro signature and has not been already put into DBlink, add it.       
       if ((! defined ($seen1{$featid})) && (($dbdesc eq "Pfam") || ($dbdesc eq "PRINTS") || ($dbdesc eq "PROSITE"))) {
	   my $newdbentry = Bio::EnsEMBL::DBEntry->new();
	   $newdbentry->primary_id($featid);
	   $newdbentry->dbname($dbdesc);
	   
#To work with SP dump the signature id has to be given, because we don't store it, an X is given instead
	   $newdbentry->optional_id("X");
	   $protein->annotation->add_DBLink($newdbentry);
	   $seen1{$featid} = 1;
       }

#Get the Interpro AC for the given signature
	    my $query2 = "select interpro_ac from interpro where id = '$featid'";
	    my $sth2 = $self->prepare($query2);
	    $sth2 ->execute();
	    my $interpro = $sth2->fetchrow;
       
       if (!defined $interpro) {
	   #$self->throw("$featid does not have Interpro accession number");
       }

#If the Interpro accession number has not already been put into DBlink, add it 
       if (! defined ($seen2{$interpro}) && defined $interpro) {
	   my $dblink = Bio::EnsEMBL::DBEntry->new();
	   $dblink->dbname('InterPro');
	   $dblink->primary_id($interpro);
	   $protein->annotation->add_DBLink($dblink);
	   $seen2{$interpro} = 1;
	   
       }
   }
   


#Add the Ensembl gene id (ENSG) as a DBlink to the object
   my $dblink = Bio::EnsEMBL::DBEntry->new();
   $dblink->dbname('EnsEMBL');
   $dblink->primary_id($geneid);
   $protein->annotation->add_DBLink($dblink);

#Add each protein features to the protein object
   foreach my $feat (@prot_feat) {
       if ($feat) {
	   $protein->add_Protein_feature($feat);
       }
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

#Call the method get_Transcript on Gene_Obj.om   
   my $transcript = $self->_gene_obj->get_Transcript($transcript_id);
         
   return $transcript;

}

=head2 fetchDBentry_by_dbid

 Title   : fetchDBentry_by_dbid
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetchDBentry_by_dbID{
   my ($self,$protein_id) = @_;
   my @entries = $self->_dbEntryAdaptor->fetch_by_translation($protein_id);
   return @entries;
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

   if( !$family->isa('Bio::EnsEMBL::ExternalData::Family::Family') ) {
       $self->throw(" $family is not a family object");
   }

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
       my $pep = $self->fetch_Protein_by_dbid($pepid);
       push(@proteins,$pep);
   }
   return @proteins;
}

=head2 fetch_by_array_feature

 Title   : fetch_by_array_feature
 Usage   :my @proteins = $obj->fetch_by_array_feature(@domaines)
 Function:This method allow to query proteins matching different domains, this will only return the protein which have the domain queried
 Example :
 Returns : An array of protein objects
 Args    :Interpro signatures


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
       if ($seen{$key} == $nb) {
	   push (@result, $key);
       }
   }
   
   return @result;
   
}

=head2 get_Introns

 Title   : get_Introns
 Usage   : $protein_adaptor->get_Introns($proteinid)
 Function: Return an array of protein features introns
 Example :
 Returns : 
 Args    : Protein ac


=cut

sub get_Introns{
    my ($self,$protid) = @_;
    my @array_features;
    my $previous_ex_end=0;
    my $count = 1;

    my $query = "select id from transcript where translation = '$protid'";
    my $sth = $self->prepare($query);
    $sth ->execute();
    my @rowid = $sth->fetchrow;
    my $transid = $rowid[0];
    
    if (!defined $transid) {
	$self->throw("No transcript can be retrieved with this translation id: $protid")
	}
    
    my $transcript = $self->fetch_Transcript_by_dbid($transid);
    my ($starts,$ends) = $transcript->pep_coords;
    
    my @exons = $transcript->each_Exon();
    my $nbex = scalar(@exons);
    
    foreach my $ex(@exons) {
	my $length;
	my $exid = $ex->id();
	my ($ex_start, $ex_end) = $self->get_exon_global_coordinates($exid);
	if ($previous_ex_end != 0) {
	    my $intron_start = $previous_ex_end;
	    my $intron_end = $ex_start;

	    my $feat1 = new Bio::EnsEMBL::SeqFeature ( -seqname => $protid,
						       -start => $starts->[$count],
						       -end => $starts->[$count],
						       -score => 0, 
						       -percent_id => "NULL",
						       -p_value => "NULL");
	    
	    my $feat2 = new Bio::EnsEMBL::SeqFeature (-start => $intron_start,
						      -end => $intron_end,
						      -seqname => "Intron");
	    
	    my $feature = new Bio::EnsEMBL::Protein_FeaturePair(-feature1 => $feat1,
								-feature2 => $feat2,);
	    
	    push(@array_features,$feature);
	    
	    $count++;
	}
	$previous_ex_end = $ex_end;
	
    }
     return @array_features; 
}



=head2 get_exon_global_coordinates

 Title   : get_exon_global_coordinates
 Usage   : $protein_adaptor->get_exon_global_coordinates($exon_id)
 Function: Get the start and end position of the exon in globale coordinates
 Example :
 Returns : Start and end of the exon
 Args    : Exon id


=cut

sub get_exon_global_coordinates{
   my ($self,$exid) = @_;

#This sql satement calculate the start and end of the exon in global coordinates
   my $query ="SELECT
               IF(sgp.raw_ori=1,(e.seq_start+sgp.chr_start-sgp.raw_start),
                 (sgp.chr_start+sgp.raw_end-e.seq_end)) as start,
               IF(sgp.raw_ori=1,(e.seq_end+sgp.chr_start-sgp.raw_start),
                 (sgp.chr_start+sgp.raw_end-e.seq_start)) as end
               FROM       static_golden_path as sgp ,exon as e
               WHERE      sgp.raw_id=e.contig
               AND        e.id='$exid'";
              

   my $sth = $self->prepare($query);
   $sth->execute;
  
   my @row = $sth->fetchrow;

   

   my $start = $row[0];
   my $end = $row[1];

   return ($start,$end);
}


=head2 get_snps

 Title   : get_snps
 Usage   : $protein_adaptor->get_snps($protein)
 Function:Get all of the snps for a peptide (for now return them in the format of a Protein_FeatureObject)
 Example :
 Returns :An array of Protein_FeaturePair object containing the snps  
 Args    :protein object


=cut
    
sub get_snps {
    my ($self,$protein) = @_;
    my $db = $self->db;
    my $gene;
    my $transcript;
    my @ex_snps;
    my $count = 0;
    my @array_features;
    
#Get the transcript, gene accession numbers for the given peptide
    my $transid = $protein->transcriptac();
    my $geneid = $protein->geneac();
    my $protid = $protein->id();

#For now a static golden path is built (this is not the fastest solution but the easiest one). Should be changed to direct sql statements.
    $db->add_ExternalFeatureFactory($self->snp_obj()); 
    $db->static_golden_path_type('UCSC');

    my $virtual = $db->get_StaticGoldenPathAdaptor();
    
#Fetch a virtual contig for the given transcript. All of the external are retrieved for this virtual contig as well as all genes.   

    &eprof_start('snp_vc_build');

    my $vc = $virtual->fetch_VirtualContig_of_transcript($transid,0);

    &eprof_end('snp_vc_build');

    my @snips = $vc->get_all_ExternalFeatures();
    my ($gene) = $vc->get_Genes($geneid);


    
#Get which virtual transcript holds the transcript we want to study   


   my @transcripts = $gene->each_Transcript();

   foreach my $trans  (@transcripts) {
       if ($trans->id eq $transid) {
	   $transcript = $trans;
       }
   }


   my $genesnp;
   eval {
       require Bio::EnsEMBL::ExternalData::GeneSNP;
       $genesnp = new Bio::EnsEMBL::ExternalData::GeneSNP
	   (-gene => $gene,
	    -contig => $vc
	    );
   };
   if( $@ ) {
       $self->throw("Unable to build GeneSNP object. Probably don't have ensembl-external code in your PERL5LIB. Download ensembl-external and either install or move to PERL5LIB\n\nReal Exception $@");
   }
       

   $genesnp->transcript($transcript);
   
   my @seq_diff = $genesnp->snps2transcript(@snips);
       
    foreach my $diff (@seq_diff) {
	foreach my $var ($diff->each_Variant) {
	    if($var->isa('Bio::Variation::AAChange') ) {
		print STDERR $var->label, "\n";
		print STDERR $var->trivname, "\n";
		print STDERR $var->allele_ori->seq, "\n";
		print STDERR $var->allele_mut->seq, "\n";
		foreach my $all ($var->each_Allele) {
		    print STDERR $all->seq, "\n";
		}
	    }
	}
    }

    push(@array_features,@seq_diff);
   
    return @array_features;
}




