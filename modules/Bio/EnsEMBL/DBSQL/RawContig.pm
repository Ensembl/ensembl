#
# BioPerl module for Contig
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DB::RawContig - Handle onto a database stored raw contiguous DNA 

=head1 SYNOPSIS

    # get a contig object somehow, eg from an DB::Obj

    @genes = $contig->get_all_Genes();
    @sf    = $contig->get_all_RepeatFeatures();
    @sf    = $contig->get_all_SimilarityFeatures();

    $contig->id();
    $contig->length();
    $primary_seq = $contig->primary_seq();

=head1 DESCRIPTION

A RawContig is physical piece of DNA coming out a sequencing project,
ie a single product of an assembly process. A RawContig defines an atomic
coordinate system on which features and genes are placed (remember that
genes can cross between atomic coordinate systems).

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::RawContig;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::RootI;

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBSQL::Feature_Obj;
use Bio::EnsEMBL::DBSQL::Gene_Obj;
use Bio::EnsEMBL::DB::RawContigI;

use Bio::EnsEMBL::Repeat;
use Bio::EnsEMBL::ContigOverlap;
use Bio::EnsEMBL::FeatureFactory;
use Bio::EnsEMBL::Chromosome;
use Bio::EnsEMBL::DBSQL::DBPrimarySeq;
use Bio::PrimarySeq;

@ISA = qw(Bio::EnsEMBL::DB::RawContigI Bio::Root::RootI);

sub new {
    my( $pkg, @args ) = @_;
    
    my $self = bless {}, $pkg;

    my (
        $dbobj,
        $id,
        $perlonlysequences,
        $contig_overlap_source,
        $overlap_distance_cutoff,
        ) = $self->_rearrange([qw(
            DBOBJ
	    ID
	    PERLONLYSEQUENCES
            CONTIG_OVERLAP_SOURCE
            OVERLAP_DISTANCE_CUTOFF
	    )], @args);

    $id    || $self->throw("Cannot make contig db object without id");
    $dbobj || $self->throw("Cannot make contig db object without db object");
    $dbobj->isa('Bio::EnsEMBL::DBSQL::Obj') || $self->throw("Cannot make contig db object with a $dbobj object");

    $self->id($id);
    $self->dbobj($dbobj);
    $self->_got_overlaps(0);
    $self->fetch();
    $self->perl_only_sequences($perlonlysequences);
    $self->contig_overlap_source($contig_overlap_source);
    $self->overlap_distance_cutoff($overlap_distance_cutoff);

    return $self;
}

=head2 fetch

 Title   : fetch
 Usage   : $contig->fetch($contig_id)
 Function: fetches the data necessary to build a Rawcontig object
 Example : $contig->fetch(1)
 Returns : Bio::EnsEMBL::DBSQL::RawContig object
 Args    : $contig_id


=cut

sub fetch {
    my ($self) = @_;
 
    my $id=$self->id;

#    my $sth = $self->dbobj->prepare("select c.id,c.internal_id,cl.embl_version " . 
#                           "from dna as d,contig as c,clone as cl " .
#                           "where d.id = c.dna and c.id = '$id' and c.clone = cl.id");

    my $sth = $self->dbobj->prepare("
        SELECT contig.internal_id
          , contig.dna
          , clone.embl_version
          , contig.clone
        FROM dna
          , contig
          , clone
        WHERE contig.dna = dna.id
          AND contig.clone = clone.id
          AND contig.id = ?
        ");

    my $res = $sth->execute($id);

    if (my $row = $sth->fetchrow_arrayref) {  
        $self->internal_id($row->[0]);
        $self->dna_id($row->[1]);
        $self->seq_version($row->[2]);
	$self->cloneid    ($row->[3]);
    } else {
         $self->throw("Contig $id does not exist in the database or does not have DNA sequence");
    }

    return $self;
}


=head2 get_Genes_by_Type

 Title   : get_Genes_by_Type
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Genes_by_Type{
   my ($self,$type,$supporting) = @_;
   my $contig_id = $self->internal_id();   
    # prepare the SQL statement
unless ($type){$self->throw("I need a type argument e.g. ensembl")}; 

my $query="
        SELECT t.gene
        FROM transcript t,
             exon_transcript et,
             exon e,
             genetype gt
        WHERE e.contig = '$contig_id'
          AND et.exon = e.id
          AND t.id = et.transcript
          AND gt.gene_id=t.gene
          AND gt.type = '$type'
        ";


   my $sth = $self->dbobj->prepare($query);
   
   my $res = $sth->execute();
   my %got;
   my @gene_array;
   while (my $rowhash = $sth->fetchrow_hashref) { 
       if( ! exists $got{$rowhash->{'gene'}}) {  
	   push(@gene_array,$rowhash->{'gene'});
       }
           
       $got{$rowhash->{'gene'}} = 1;
   }       
   

#   print STDERR "Gene array is [@gene_array]\n";
   my $gene_obj = Bio::EnsEMBL::DBSQL::Gene_Obj->new($self->dbobj);             

   my @out;

   if (@gene_array) {
       @out = $gene_obj->get_array_supporting($supporting,@gene_array);
   }
   return @out;


}

=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Genes{
   my ($self, $supporting) = @_;
   my @out;
   my $contig_id = $self->internal_id();   
   my %got;
    # prepare the SQL statement
   my $sth = $self->dbobj->prepare("
        SELECT t.gene
        FROM transcript t,
             exon_transcript et,
             exon e
        WHERE e.contig = '$contig_id'
          AND et.exon = e.id
          AND t.id = et.transcript
        ");

    my $res = $sth->execute();
   
    while (my $rowhash = $sth->fetchrow_hashref) { 
            
        if( ! exists $got{$rowhash->{'gene'}}) {  
            
           my $gene_obj = Bio::EnsEMBL::DBSQL::Gene_Obj->new($self->dbobj);             
	   my $gene = $gene_obj->get($rowhash->{'gene'}, $supporting);
           if ($gene) {
	        push(@out, $gene);
           }
	   $got{$rowhash->{'gene'}} = 1;
        }       
    }
   
    if (@out) {
        return @out;
    }
    return;
}


=head2 has_genes

 Title   : has_genes
 Usage   :
 Function: returns 1 if there are genes, 0 otherwise.
 Example :
 Returns : 
 Args    :


=cut

sub has_genes{
   my ($self,@args) = @_;
   my $contig_id = $self->internal_id();   

   my $seen =0;
   my $sth = $self->dbobj->prepare("select id from exon where contig = '$contig_id' limit 1");
   $sth->execute();

   my $rowhash;
   while ( ($rowhash = $sth->fetchrow_hashref()) ) {
       $seen = 1;
       last;
   }
   return $seen;
}

=head2 primary_seq

 Title   : primary_seq
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub primary_seq{
   my ($self,@args) = @_;

   if( $self->perl_only_sequences == 1 ) {
       return $self->perl_primary_seq();
   }
   return $self->db_primary_seq();

}

=head2 db_primary_seq

 Title   : db_primary_seq
 Usage   : $dbseq = $contig->db_primary_seq();
 Function: Gets a Bio::EnsEMBL::DBSQL::DBPrimarySeq object out from the contig
 Example :
 Returns : Bio::EnsEMBL::DBSQL::DBPrimarySeq object
 Args    : 


=cut

sub db_primary_seq {
    my ($self) = @_;
    
    my $dbseq = Bio::EnsEMBL::DBSQL::DBPrimarySeq->new(
						       -dna => $self->dna_id,
						       -db_handle => $self->dbobj->_db_handle
						       );
    
    return $dbseq;
}

=head2 perl_primary_seq

 Title   : seq
 Usage   : $seq = $contig->perl_primary_seq();
 Function: Gets a Bio::PrimarySeqI object out from the contig
 Example :
 Returns : Bio::PrimarySeqI object
 Args    :


=cut

sub perl_primary_seq {
    my ($self) = @_;

    if ( $self->_seq_cache() ) {
        return $self->_seq_cache();
    }

    my $dna_id = $self->dna_id()
        or $self->throw("No dna_id in RawContig ". $self->id);
    my $sth = $self->dbobj->prepare(q{ SELECT sequence FROM dna WHERE id = ? });
    my $res = $sth->execute($dna_id);

    my($str) = $sth->fetchrow
        or $self->throw("No DNA sequence in RawContig " . $self->id . " for dna id " . $dna_id);

    # Shouldn't sequence integrity be checked on the way
    # into the datbase instead of here?
    $str =~ /[^ABCDGHKMNRSTVWY]/ && $self->warn("Got some non standard DNA characters here! Yuk!");
    $str =~ s/\s//g;
    $str =~ s/[^ABCDGHKMNRSTVWY]/N/g;

    my $ret = Bio::PrimarySeq->new( 
        -seq => $str, 
        -display_id => $self->id,           # eg: AC004092.00001
        -primary_id => $self->internal_id,  # eg: 874
        -moltype => 'dna',
        );
    $self->_seq_cache($ret);

    return $ret;
}

=head2 _seq_cache

 Title   : _seq_cache
 Usage   : $obj->_seq_cache($newval)
 Function: Used to cache the primary seq object to avoid 
           more than one trip to the database for the dna
 Returns : value of _seq_cache
 Args    : newvalue (optional)


=cut

sub _seq_cache{
   my $obj = shift;
   if( @_ ) {
       my $value = shift;
       $obj->{'_seq_cache'} = $value;
   }
   return $obj->{'_seq_cache'};

}


=head2 get_all_SeqFeatures

 Title   : get_all_SeqFeatures
 Usage   : foreach my $sf ( $contig->get_all_SeqFeatures
 Function: Gets all the sequence features on the whole contig
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SeqFeatures {
    my ($self) = @_;

    my @out;

    push(@out,$self->get_all_SimilarityFeatures);
    push(@out,$self->get_all_RepeatFeatures);
#   push(@out,$self->get_all_PredictionFeatures);

    return @out;
}


=head2 get_all_SimilarityFeatures_above_score

 Title   : get_all_SimilarityFeatures_above_score
 Usage   : foreach my $sf ( $contig->get_all_SimilarityFeatures_above_score(analysis_type, score) ) 
 Function:
 Example :
 Returns : 
 Args    : 


=cut

sub get_all_SimilarityFeatures_above_score{
    my ($self, $analysis_type, $score) = @_;

    $self->throw("Must supply analysis_type parameter") unless $analysis_type;
    $self->throw("Must supply score parameter") unless $score;
    
   my @array;

   my $id     = $self->internal_id();
   my $length = $self->length();

   my %analhash;

   #First of all, get all features that are part of a feature set with high enough score and have the right type

    my $statement = "SELECT feature.id, seq_start, seq_end, strand, feature.score, analysis, name, " .
		             "hstart, hend, hid, fset, rank, fset.score " .
		     "FROM   feature, fset_feature, fset, analysis " .
		     "WHERE  feature.contig ='$id' " .
		     "AND    fset_feature.feature = feature.id " .
		     "AND    fset.id = fset_feature.fset " .
                     "AND    feature.score > '$score' " .
                     "AND    feature.analysis = analysis.id " .
                     "AND    analysis.db = '$analysis_type' " .
                     "ORDER BY fset";
		     
   my $sth = $self->dbobj->prepare($statement);                                                                       
   $sth->execute();
   
   my ($fid,$start,$end,$strand,$f_score,$analysisid,$name,$hstart,$hend,$hid,$fset,$rank,$fset_score);
   my $seen = 0;
   
   # bind the columns
   $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$f_score,\$analysisid,\$name,\$hstart,\$hend,\$hid,\$fset,\$rank,\$fset_score);
   
   my $out;
   
   my $fset_id_str = "";

   while($sth->fetch) {

       my $analysis;

       if (!$analhash{$analysisid}) {
	   
	   my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->dbobj);
	   $analysis = $feature_obj->get_Analysis($analysisid);
	   $analhash{$analysisid} = $analysis;
       
       } else {
	   $analysis = $analhash{$analysisid};
       }
       
       if( !defined $name ) {
	   $name = 'no_source';
       }
       
       #Build fset feature object if new fset found
       if ($fset != $seen) {
#	   print(STDERR "Making new fset feature $fset\n");
	   $out =  new Bio::EnsEMBL::SeqFeature;
	   $out->id($fset);
	   $out->analysis($analysis);
	   $out->seqname ($self->id);
	   $out->score($fset_score);
	   $out->source_tag($name);
	   $out->primary_tag("FSET");

	   $seen = $fset;
	   push(@array,$out);
       }
       $fset_id_str = $fset_id_str . $fid . ",";       
       #Build Feature Object
       my $feature = new Bio::EnsEMBL::SeqFeature;
       $feature->seqname   ($self->id);
       $feature->start     ($start);
       $feature->end       ($end);
       $feature->strand    ($strand);
       $feature->source_tag($name);
       $feature->primary_tag('similarity');
       $feature->id         ($fid);
       
       if( defined $f_score ) {
	   $feature->score($f_score);
       }
       
       $feature->analysis($analysis);
       
       # Final check that everything is ok.
       $feature->validate();

       #Add this feature to the fset
       $out->add_sub_SeqFeature($feature,'EXPAND');

   }
   
   #Then get the rest of the features, i.e. featurepairs and single features that are not part of a fset
   $fset_id_str =~ s/\,$//;

   if ($fset_id_str) {
        $statement = "SELECT feature.id, seq_start, seq_end, strand, score, analysis, name, hstart, hend, hid " .
		     "FROM   feature, analysis " .
                     "WHERE  id not in (" . $fset_id_str . ") " .
                     "AND    feature.score > '$score' " . 
                     "AND    feature.analysis = analysis.id " .
                     "AND    analysis.db = '$analysis_type' " .
                     "AND    feature.contig = '$id' ";
                                     
       $sth = $self->dbobj->prepare($statement);
       
   } else {
        $statement = "SELECT feature.id, seq_start, seq_end, strand, score, analysis, name, hstart, hend, hid " .
		     "FROM   feature, analysis " .
                     "WHERE  feature.score > '$score' " . 
                     "AND    feature.analysis = analysis.id " .
                     "AND    analysis.db = '$analysis_type' " .
                     "AND    feature.contig = '$id' ";
                     
                     
                     
       $sth = $self->dbobj->prepare($statement);
   }

   $sth->execute();

   # bind the columns
   $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$f_score,\$analysisid,\$name,\$hstart,\$hend,\$hid);
   
   while($sth->fetch) {
       my $out;
       my $analysis;
              
       if (!$analhash{$analysisid}) {
	   
	   my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->dbobj);
	   $analysis = $feature_obj->get_Analysis($analysisid);
	   $analhash{$analysisid} = $analysis;
	   
       } else {
	   $analysis = $analhash{$analysisid};
       }
       
       if( !defined $name ) {
	   $name = 'no_source';
       }
       
       if( $hid ne '__NONE__' ) {
	   # is a paired feature
	   # build EnsEMBL features and make the FeaturePair
	 
	   $out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();


	   $out->set_all_fields($start,$end,$strand,$f_score,$name,'similarity',$self->id,
				$hstart,$hend,1,$f_score,$name,'similarity',$hid);

	   $out->analysis    ($analysis);
	   $out->id          ($hid);              # MC This is for Arek - but I don't
	                                          #    really know where this method has come from.
       } else {
	   $out = new Bio::EnsEMBL::SeqFeature;
	   $out->seqname   ($self->id);
	   $out->start     ($start);
	   $out->end       ($end);
	   $out->strand    ($strand);
	   $out->source_tag($name);
	   $out->primary_tag('similarity');
	   $out->id         ($fid);

	   if( defined $f_score ) {
	       $out->score($f_score);
	   }
	   $out->analysis($analysis);
       }
       # Final check that everything is ok.
       $out->validate();
       
      push(@array,$out);
      
   }
   
   return @array;
}



=head2 get_all_SimilarityFeatures

 Title   : get_all_SimilarityFeatures
 Usage   : foreach my $sf ( $contig->get_all_SimilarityFeatures($start,$end) ) 
 Function: Gets all the sequence similarity features on the whole contig
 Example :
 Returns : 
 Args    : 


=cut

sub get_all_SimilarityFeatures{
   my ($self) = @_;

   my @array;
   my @fps;

   my $id     = $self->internal_id();
   my $length = $self->length();

   my %analhash;

   #First of all, get all features that are part of a feature set

   #my $sth = $self->dbobj->prepare("select  p1.id, " .
   #                          "p1.seq_start, p1.seq_end, " . 
   #                           "p1.strand,p1.score,p1.analysis,p1.name,  " .
   #                          "p1.hstart,p1.hend,p1.hid,"  .
   #                          "p2.fset,p2.rank, " . 
   #                          "fs.score " .
   #                 "from    feature as p1,  " .
   #                 "        fset_feature as p2, " .
   #                 "        fset as fs " .
   #                 "where   p1.contig ='$id' " .
   #                 "and     p2.feature = p1.id " .
   #                 "and     fs.id = p2.fset " .
   #                 "order by p2.fset");

    my $statement = "SELECT feature.id, seq_start, seq_end, strand, feature.score, analysis, name, " .
		             "hstart, hend, hid, fset, rank, fset.score " .
		     "FROM   feature, fset_feature, fset " .
		     "WHERE  feature.contig =$id " .
		     "AND    fset_feature.feature = feature.id " .
		     "AND    fset.id = fset " .
                     "ORDER BY fset";
		     
                                    
   my $sth = $self->dbobj->prepare($statement);                                    
                                    
   $sth->execute();
   
   my ($fid,$start,$end,$strand,$f_score,$analysisid,$name,$hstart,$hend,$hid,$fset,$rank,$fset_score);
   my $seen = 0;
   
   # bind the columns
   $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$f_score,\$analysisid,\$name,\$hstart,\$hend,\$hid,\$fset,\$rank,\$fset_score);
   
   my $out;
   
   my $fset_id_str = "";

   while($sth->fetch) {

       my $analysis;

       if (!$analhash{$analysisid}) {
	   
	   my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->dbobj);
	   $analysis = $feature_obj->get_Analysis($analysisid);
	   $analhash{$analysisid} = $analysis;
       
       } else {
	   $analysis = $analhash{$analysisid};
       }
       
       if( !defined $name ) {
	   $name = 'no_source';
       }
       
       #Build fset feature object if new fset found
       if ($fset != $seen) {
#	   print(STDERR "Making new fset feature $fset\n");
	   $out =  new Bio::EnsEMBL::SeqFeature;
	   $out->id($fset);
	   $out->analysis($analysis);
	   $out->seqname ($self->id);
	   $out->raw_seqname ($self->id);
	   $out->score($fset_score);
	   $out->source_tag($name);
	   $out->primary_tag("FSET");

	   $seen = $fset;
	   push(@array,$out);
       }
       $fset_id_str = $fset_id_str . $fid . ",";       
       #Build Feature Object
       my $feature = new Bio::EnsEMBL::SeqFeature;
       $feature->seqname   ($self->id);
       $feature->raw_seqname ($self->id);
       $feature->start     ($start);
       $feature->end       ($end);
       $feature->strand    ($strand);
       $feature->source_tag($name);
       $feature->primary_tag('similarity');
       $feature->id         ($fid);
       
       if( defined $f_score ) {
	   $feature->score($f_score);
       }
       
       $feature->analysis($analysis);
       
       # Final check that everything is ok.
       $feature->validate();

       #Add this feature to the fset
       $out->add_sub_SeqFeature($feature,'EXPAND');

   }
   
   #Then get the rest of the features, i.e. featurepairs and single features that are not part of a fset
   $fset_id_str =~ s/\,$//;

   if ($fset_id_str) {
       $sth = $self->dbobj->prepare("select id,seq_start,seq_end,strand,score,analysis,name,hstart,hend,hid " .
				     "from feature where id not in (" . $fset_id_str . ") and contig = $id");
   } else {
       $sth = $self->dbobj->prepare("select id,seq_start,seq_end,strand,score,analysis,name,hstart,hend,hid ".
				     "from feature where contig = $id");
   }

   $sth->execute();

   # bind the columns
   $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$f_score,\$analysisid,\$name,\$hstart,\$hend,\$hid);
   
   while($sth->fetch) {
       my $out;
       my $analysis;
              
       if (!$analhash{$analysisid}) {
	   
	   my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->dbobj);
	   $analysis = $feature_obj->get_Analysis($analysisid);
	   $analhash{$analysisid} = $analysis;
	   
       } else {
	   $analysis = $analhash{$analysisid};
       }
       
       if( !defined $name ) {
	   $name = 'no_source';
       }
       
       if( $hid ne '__NONE__' ) {
	   # is a paired feature
	   # build EnsEMBL features and make the FeaturePair
	 
	   $out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();


	   $out->set_all_fields($start,$end,$strand,$f_score,$name,'similarity',$self->id,
				$hstart,$hend,1,$f_score,$name,'similarity',$hid);

	   $out->analysis    ($analysis);
	   $out->id          ($hid);              # MC This is for Arek - but I don't
	                                          #    really know where this method has come from.
       } else {
	   $out = new Bio::EnsEMBL::SeqFeature;
	   $out->seqname   ($self->id);
	   $out->raw_seqname   ($self->id);
	   $out->start     ($start);
	   $out->end       ($end);
	   $out->strand    ($strand);
	   $out->source_tag($name);
	   $out->primary_tag('similarity');
	   $out->id         ($fid);

	   if( defined $f_score ) {
	       $out->score($f_score);
	   }
	   $out->analysis($analysis);
       }
       # Final check that everything is ok.
       $out->validate();
       if( $out->can('attach_seq') ) {
	   $out->attach_seq($self->primary_seq);
       }

      push(@fps,$out);
      
   }
   
   my @extras = $self->get_extra_features;
   my @newfeatures = $self->filter_features(\@fps,\@extras);

#   print ("Size " . @array . " " . @newfeatures . "\n");
   push(@array,@newfeatures);

   return @array;
}

=head2 get_all_RepeatFeatures

 Title   : get_all_RepeatFeatures
 Usage   : foreach my $sf ( $contig->get_all_RepeatFeatures )
 Function: Gets all the repeat features on a contig.
 Example :
 Returns : 
 Args    :  


=cut

sub get_all_RepeatFeatures {
   my ($self) = @_;

   my @array;

   my $id     = $self->internal_id();
   my $length = $self->length();

   my %analhash;

   # make the SQL query
    my $statement = "select id,seq_start,seq_end,strand,score,analysis,hstart,hend,hid " . 
				    "from repeat_feature where contig = '$id'";
                                    

   my $sth = $self->dbobj->prepare($statement);

   $sth->execute();

   my ($fid,$start,$end,$strand,$score,$analysisid,$hstart,$hend,$hid);

   # bind the columns
   $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$score,\$analysisid,\$hstart,\$hend,\$hid);

   while( $sth->fetch ) {
       my $out;
       my $analysis;

       if (!$analhash{$analysisid}) {
	   
	   my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->dbobj);
	   $analysis = $feature_obj->get_Analysis($analysisid);

	   $analhash{$analysisid} = $analysis;

       } else {
	   $analysis = $analhash{$analysisid};
       }


       if( $hid ne '__NONE__' ) {
	   # is a paired feature
	   # build EnsEMBL features and make the FeaturePair

	   $out = Bio::EnsEMBL::FeatureFactory->new_repeat();
	   $out->set_all_fields($start,$end,$strand,$score,'repeatmasker','repeat',$self->id,
				$hstart,$hend,1,$score,'repeatmasker','repeat',$hid);

	   $out->analysis($analysis);

       } else {
	   $self->warn("Repeat feature does not have a hid. bad news....");
       }
       
       $out->validate();

      push(@array,$out);
  }
 
   return @array;
}

=head2 get_MarkerFeatures

  Title   : get_MarkerFeatures 
  Usage   : @fp = $contig->get_MarkerFeatures; 
  Function: Gets MarkerFeatures. MarkerFeatures can be asked for a Marker. 
            Its assumed, that when you can get MarkerFeatures, then you can 
            get the Map Code as well.
  Example : - 
  Returns : -
  Args : -

=cut


sub get_MarkerFeatures {
  my $self = shift;

  my $id = $self->internal_id;
  my @result = ();
  eval {
    require Bio::EnsEMBL::Map::MarkerFeature;

    # features for this contig with db=mapprimer
    my $sth = $self->dbobj->prepare
      ( "select f.seq_start, f.seq_end, f.score, f.strand, f.name, ".
        "f.hstart, f.hend, f.hid, f.analysis ".
        "from feature f, analysis a ".
        "where f.contig='$id' and ".
        "f.analysis = a.id and a.db='mapprimer'" );
    $sth->execute;
    
    my ($start, $end, $score, $strand, $hstart, 
        $name, $hend, $hid, $analysisid );
    my $analysis;
    my %analhash;

    $sth->bind_columns
      ( undef, \$start, \$end, \$score, \$strand, \$name, 
        \$hstart, \$hend, \$hid, \$analysisid );
        
    while( $sth->fetch ) {
      my ( $out, $seqf1, $seqf2 );
      
      if (!$analhash{$analysisid}) {
        $analysis = $self->dbobj->get_Analysis($analysisid);
        $analhash{$analysisid} = $analysis;
        
      } else {
        $analysis = $analhash{$analysisid};
      }
    
      $seqf1 = Bio::EnsEMBL::SeqFeature->new();
      $seqf2 = Bio::EnsEMBL::SeqFeature->new();
      $out = Bio::EnsEMBL::Map::MarkerFeature->new
	( -feature1 => $seqf1, -feature2 => $seqf2 );
      $out->set_all_fields
        ( $start,$end,$strand,$score,
          $name,'similarity',$self->id,
          $hstart,$hend,1,$score,$name,'similarity',$hid);
          $out->analysis($analysis);
      $out->mapdb( $self->dbobj->mapdb );
      $out->id ($hid);
      push( @result, $out );
    }
  };

  if( $@ ) {
    print STDERR ("Install the Ensembl-map package for this feature" );
  }
  return @result;
}


=head2 get_all_PredictionFeatures

 Title   : get_all_PredictionFeatures
 Usage   : foreach my $sf ( $contig->get_all_RepeatFeatures )
 Function: Gets all the repeat features on a contig.
 Example :
 Returns : 
 Args    : 


=cut

sub get_all_PredictionFeatures {
   my ($self) = @_;

   my @array;

   my $id     = $self->internal_id();
   my $length = $self->length();
   my $fsetid;
   my $previous;
   my %analhash;

   # make the SQL query
   my $query = "select f.id,f.seq_start,f.seq_end,f.strand,f.score,f.analysis,fset.id " . 
       "from feature f, fset fset,fset_feature ff where ff.feature = f.id and fset.id = ff.fset and contig = $id and name = 'genscan'";
    
   my $sth = $self->dbobj->prepare($query);
   
   $sth->execute();
   
   my ($fid,$start,$end,$strand,$score,$analysisid);
   
   # bind the columns
   $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$score,\$analysisid,\$fsetid);
   
   $previous = undef;
   my $current_fset;
   while( $sth->fetch ) {
       my $out;
       
       my $analysis;
	   
       if (!$analhash{$analysisid}) {

	   my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->dbobj);
	   $analysis = $feature_obj->get_Analysis($analysisid);

	   $analhash{$analysisid} = $analysis;
	   
       } else {
	   $analysis = $analhash{$analysisid};
       }


       if( $fsetid != $previous ) {
	   $current_fset = new Bio::EnsEMBL::SeqFeature;
	   $current_fset->source_tag('genscan');
	   $current_fset->primary_tag('prediction');
	   $current_fset->analysis($analysis);
	   $current_fset->seqname($self->id);
	   $current_fset->raw_seqname($self->id);
	   $current_fset->id($fsetid);
	   push(@array,$current_fset);
       }

       $out = new Bio::EnsEMBL::SeqFeature;
       
       $out->seqname   ($self->id);
       $out->raw_seqname($self->id);
       $out->start     ($start);
       $out->end       ($end);
       $out->strand    ($strand);
       $out->id($fsetid); # to make genscan peptide work
       $out->source_tag('genscan');
       $out->primary_tag('prediction');
       
       if( defined $score ) {
	   $out->score($score);
       }

       $out->analysis($analysis);

       # Final check that everything is ok.
       
       $out->validate();
       $current_fset->add_sub_SeqFeature($out,'EXPAND');
       $current_fset->strand($strand);
       $previous = $fsetid;
  }
 
   return @array;
}

=head2 get_all_ExternalFeatures

 Title   : get_all_ExternalFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    : 


=cut

sub get_all_ExternalFeatures{
   my ($self) = @_;
   
   my @out;
   my $acc;
   
   # this is not pretty.
   $acc = $self->id();
   $acc =~ s/\.\d+$//g;
   my $embl_offset = $self->embl_offset();

   foreach my $extf ( $self->dbobj->_each_ExternalFeatureFactory ) {
       push(@out,$extf->get_Ensembl_SeqFeatures_contig($self->id,$self->seq_version,1,$self->length));
       
       foreach my $sf ( $extf->get_Ensembl_SeqFeatures_clone($acc,$self->seq_version,$self->embl_offset,$self->embl_offset+$self->length()) ) {
	   my $start = $sf->start - $embl_offset;
	   my $end   = $sf->end   - $embl_offset;
	   $sf->start($start);
	   $sf->end($end);
	   push(@out,$sf);
       }
   }
   my $id = $self->id();
   foreach my $f ( @out ) {
       $f->seqname($id);
   }

   return @out;

}


=head2 length

 Title   : length
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub length{
   my ($self,@args) = @_;
   my $id= $self->internal_id();
    $self->throw("Internal ID not set") unless $id;
   if (! defined ($self->{_length})) {
       my $sth = $self->dbobj->prepare("select length from contig where internal_id = \"$id\" ");
       $sth->execute();
       
       my $rowhash = $sth->fetchrow_hashref();
       
       $self->{_length} = $rowhash->{'length'};
   }

   return $self->{_length};
}

sub cloneid {
    my ($self,$arg) = @_;
    
    if (defined($arg)) {
	$self->{_cloneid} = $arg;
    }

    return $self->{_cloneid};
}


=head2 seq_version

 Title   : seq_version
 Usage   : $obj->seq_version($newval)
 Function: 
 Example : 
 Returns : value of seq_version
 Args    : newvalue (optional)


=cut

sub seq_version{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'seq_version'} = $value;
    }
    return $obj->{'seq_version'};

}

=head2 embl_order

 Title   : order
 Usage   : $obj->embl_order
 Function: 
 Returns : 
 Args    : 


=cut

sub embl_order{
   my $self = shift;
   my $id = $self->id();
   my $sth = $self->dbobj->prepare("select corder from contig where id = \"$id\" ");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref();
   return $rowhash->{'corder'};
   
}




=head2 embl_offset

 Title   : embl_offset
 Usage   : 
 Returns : 
 Args    :


=cut

sub embl_offset{
   my $self = shift;
   my $id = $self->id();


   my $sth = $self->dbobj->prepare("select offset from contig where id = \"$id\" ");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref();
   return $rowhash->{'offset'};

}

=head2 embl_accession

 Title   : embl_accession
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub embl_accession{
   my $self = shift;

   return "AL000000";
}


=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Example : 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id {
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'id'} = $value;
    }
    return $self->{'id'};

}

=head2 internal_id

 Title   : internal_id
 Usage   : $obj->internal_id($newval)
 Function: 
 Example : 
 Returns : value of database internal id
 Args    : newvalue (optional)

=cut

sub internal_id {
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'internal_id'} = $value;
    }
    return $self->{'internal_id'};

}

=head2 dna_id

 Title   : dna_id
 Usage   : $obj->dna_id($newval)
 Function: Get or set the id for this contig in the dna table
 Example : 
 Returns : value of dna id
 Args    : newvalue (optional)


=cut

sub dna_id {
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'dna_id'} = $value;
    }
    return $self->{'dna_id'};

}


=head2 seq_date

 Title   : seq_date
 Usage   : $contig->seq_date()
 Function: Gives the unix time value of the dna table created datetime field, which indicates
           the original time of the dna sequence data
 Example : $contig->seq_date()
 Returns : unix time
 Args    : none


=cut

sub seq_date {
   my ($self) = @_; 

   my $id = $self->internal_id();
   my $sth = $self->dbobj->prepare("select UNIX_TIMESTAMP(d.created) from dna as d,contig as c where c.internal_id = $id and c.dna = d.id");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref(); 
   return $rowhash->{'UNIX_TIMESTAMP(d.created)'};
}


=head2 get_left_overlap

 Title   : get_left_overlap
 Usage   : $overlap_object = $contig->get_left_overlap();
 Function: Returns the overlap object of contig to the left.
           This could be undef, indicating no overlap
 Returns : A Bio::EnsEMBL::ContigOverlapHelper object
 Args    : None

=cut

sub get_left_overlap {
   my ($self,@args) = @_;

   if( $self->_got_overlaps == 0 ) {
       $self->_load_overlaps() ;
   }

   return $self->_left_overlap();
}


=head2 get_right_overlap

 Title   : get_right_overlap
 Usage   : $overlap_object = $contig->get_right_overlap();
 Function: Returns the overlap object of contig to the left.
           This could be undef, indicating no overlap
 Returns : A Bio::EnsEMBL::ContigOverlapHelper object
 Args    : None

=cut

sub get_right_overlap {
   my ($self,@args) = @_;

   if( $self->_got_overlaps == 0 ) {
       $self->_load_overlaps() ;
   }

   return $self->_right_overlap();
}



sub _db_obj {
   my ($self,@args) = @_;
   $self->warn("Someone is using a deprecated _db_obj call!");
   return $self->dbobj(@args);
}

=head2 dbobj

 Title   : dbobj
 Usage   :
 Function:
 Example :
 Returns : The Bio::EnsEMBL::DBSQL::ObjI object
 Args    :


=cut

sub dbobj {
   my ($self,$arg) = @_;

   if (defined($arg)) {
        $self->throw("[$arg] is not a Bio::EnsEMBL::DBSQL::Obj") unless $arg->isa("Bio::EnsEMBL::DBSQL::Obj");
        $self->{'_dbobj'} = $arg;
   }
   return $self->{'_dbobj'};
}
	
=head2 _got_overlaps

 Title   : _got_overlaps
 Usage   : $obj->_got_overlaps($newval)
 Function: 
 Returns : value of _got_overlaps
 Args    : newvalue (optional)


=cut

sub _got_overlaps {
    my($obj, $value) = @_;
    if (defined($value)) {
        $obj->{'_got_overlaps'} = $value;
    }
    return $obj->{'_got_overlaps'};
}

=head2 _load_overlaps

 Title   : _load_overlaps
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut


sub _load_overlaps {
    my( $self ) = @_;
    
    my $id = $self->id;
    my @over = $self->get_all_Overlaps;
    foreach my $lap (@over) {
        my( $end, $helper ) = $lap->make_ContigOverlapHelper($id);
        if ($end eq 'left') {
            $self->_left_overlap($helper);
        }
        elsif ($end eq 'right') {
            $self->_right_overlap($helper);
        }
        else {
            $self->throw("Weird, got: '$end', '$helper'");
        }
    }

    # Flag that we've visited the database to get overlaps
    $self->_got_overlaps(1);

    # sanity check ourselves
    if( $self->golden_start > $self->golden_end ) {
	$self->throw("This contig ".$self->id." has dodgy golden start/ends with start:".$self->golden_start." end:".$self->golden_end);
    }
}



{ # Begin bare block to keep @queries array private to get_all_Overlaps

    # Doing two queries seems to be quickest
    # Statements like:
    #   c.dna = o.dna_b_id OR c.dna = o.dna_a_id
    # make queries inordinately slow because MySQL
    # doesn't use indices on OR statements.
    my @queries = (
       q{SELECT cb.id
          , co.contig_a_position
          , co.contig_b_position
          , co.overlap_type
          , co.overlap_size
          , co.type
        FROM contigoverlap co
          , contig cb
        WHERE co.dna_b_id = cb.dna
          AND co.dna_a_id = ?
        },

       q{SELECT ca.id
          , co.contig_a_position
          , co.contig_b_position
          , co.overlap_type
          , co.overlap_size
          , co.type
        FROM contigoverlap co
          , contig ca
        WHERE co.dna_a_id = ca.dna
          AND co.dna_b_id = ?
        },
    );

    sub get_all_Overlaps {
        my ($self) = @_;

        my $id                  = $self->dna_id();
        my $version             = $self->seq_version();
        my $overlap_source_sub  = $self->contig_overlap_source();
        my $overlap_cutoff      = $self->overlap_distance_cutoff();
        
        my( @overlap );
        foreach my $i (0,1) {
            my $query_str = $queries[$i];

            my $sth = $self->dbobj->prepare($query_str);
            $sth->execute($id);

            while (my $row = $sth->fetchrow_arrayref) {
                
                my( $sister_id,
                    $pos_a,
                    $pos_b,
                    $type,
                    $distance,
                    $source,
                    ) = @$row;
                
                # Skip this overlap if it isn't from the right source
                next unless &$overlap_source_sub($source);
                
                # Skip overlaps with distances larger than the cutoff
                if ($overlap_cutoff > -1 and $distance > $overlap_cutoff) {
                    next;
                }
                
                # Make the other contig of the overlap
                my( $contig_a, $contig_b );
                if ($i == 0) {
                    $contig_a = $self;
                    $contig_b = $self->dbobj->get_Contig($sister_id);
                } else {
                    $contig_a = $self->dbobj->get_Contig($sister_id);
                    $contig_b = $self;
                }

                my $new_overlap = Bio::EnsEMBL::ContigOverlap->new(
                    '-contiga'      => $contig_a,
                    '-contigb'      => $contig_b,
                    '-positiona'    => $pos_a,
                    '-positionb'    => $pos_b,
                    '-overlap_type' => $type,
                    '-distance'     => $distance,
                    '-source'       => $source,
                    );
                push(@overlap, $new_overlap);
            }
        }
        if (@overlap > 2) {
            $self->throw("Got '". scalar(@overlap) ."' overlaps, which is too many for 1 contig!");
        } else {
            return @overlap;
        }
    }
} # End privacy block

=head2 _right_overlap

 Title   : _right_overlap
 Usage   : $obj->_right_overlap($newval)
 Function: 
 Example : 
 Returns : value of _right_overlap
 Args    : newvalue (optional)


=cut

sub _right_overlap {
   my ($obj,$value) = @_;


   if( defined $value) {
      $obj->{'_right_overlap'} = $value;
    }
    return $obj->{'_right_overlap'};

}

=head2 _left_overlap

 Title   : _left_overlap
 Usage   : $obj->_left_overlap($newval)
 Function: 
 Example : 
 Returns : value of _left_overlap
 Args    : newvalue (optional)


=cut

sub _left_overlap {
    my ($obj,$value) = @_;

    if (defined $value) {
        $obj->{'_left_overlap'} = $value;
    }
    return $obj->{'_left_overlap'};

}


sub feature_file {
    my ($self) = @_;
    
    return;

    if (!(defined($self->{_feature_file}))) {

	my $cloneid = $self->cloneid;
	
	if (defined($cloneid)) {
	    print STDERR "Fetching job for $cloneid\n";
	    my @job = $self->dbobj->get_JobsByInputId($cloneid);
	    print STDERR "Job is $job[0]\n";
	    
	    if (defined($job[0])) {
		$self->{_feature_file} = $job[0]->output_file;
	    }
	}
    }
    return $self->{_feature_file};

}

sub get_extra_features{
    my ($self) = @_;
    my @out;

    if (!defined($self->{_extras})) {
	$self->{_extras} = [];

#	print (STDERR "Output file " . $self->feature_file ."\n");
	if (defined($self->feature_file) && (-e $self->feature_file)) {
	    
	    my $object;
	    open (IN,"<" . $self->feature_file) || do {print STDERR ("Could not open output file\n")};
	    
	    while (<IN>) {
		$_ =~ s/\[//;
		$_ =~ s/\]//;
		$object .= $_;
	    }
	    close(IN);
	    
	    if (defined($object)) {
		my (@obj) = FreezeThaw::thaw($object);
		foreach my $array (@obj) {
#		    print STDERR "$array\n";
		    foreach my $f (@$array) {
			#    print STDERR "$f\n";
			if ($f->isa("Bio::EnsEMBL::FeaturePair")) {
			    $f->source_tag("vert_est2genome");
			    $f->primary_tag("similarity"); 
			    $f->analysis($self->analysis);
			    push(@out,$f);
#			print STDERR "Adding " . $f->hseqname . "\n";
			}
		    }
		}
	    }
	    
	}
	push(@{$self->{_extras}},@out);
    }
    return @{$self->{_extras}};
}

sub analysis {
    my ($self,$arg) = @_;

    if (!(defined($self->{_analysis}))) {
	my $ana = new Bio::EnsEMBL::Analysis(-db => 'vert',
					     -dbversion => 1,
					     -program => 'vert_est2genome',
					     -program_version => 1,
					     -gff_source => 'vert_est2genome',
					     -gff_features => 'similarity',
					     );
	$self->{_analysis} = $ana;
    }
    return $self->{_analysis};
}


sub filter_features {
    my ($self,$old,$new) = @_;

    my %newids;
    my %oldids;
    my %newfeatures;
    my %oldfeatures;
    my @out;

    foreach my $f (@$new) {
	$newids{$f->hseqname}++;
	push(@{$newfeatures{$f->hseqname}},$f);
    }

    foreach my $f (@$old) {
	$oldids{$f->hseqname}++;
	push(@{$oldfeatures{$f->hseqname}},$f);
	if (! (exists $newfeatures{$f->hseqname})) {
	    push(@out,$f);
	}
    }

    foreach my $newid (keys %newids) {
	if ($newids{$newid} > $oldids{$newid}) {
	    print(STDERR "Using new features for $newid\n");
	    push(@out,@{$newfeatures{$newid}});
	} else {
	    print(STDERR "Using old features for $newid\n");
	    push(@out,@{$oldfeatures{$newid}});
	}
    }
    return @out;
}

=head2 perl_only_sequences

 Title   : perl_only_sequences
 Usage   : $obj->perl_only_sequences($newval)
 Function: 
 Returns : value of perl_only_sequences
 Args    : newvalue (optional)


=cut

sub perl_only_sequences{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'perl_only_sequences'} = $value;
    }
    return $obj->{'perl_only_sequences'};

}


=head2 contig_overlap_source

 Title   : contig_overlap_source
 Usage   : my $source_sub = $contig->contig_overlap_source()
 Function: Gets or sets a subroutine which is used to
           decide which overlap sources are used to
           build virtual contigs.
 Returns : value of contig_overlap_source.
 Args    : ref to a subroutine


=cut

sub contig_overlap_source {
    my( $self, $sub ) = @_;
    
    if ($sub) {
        $self->throw("'$sub' is not a CODE reference")
            unless ref($sub) eq 'CODE';
        $self->{'_contig_overlap_source'} = $sub;
    }
    return $self->{'_contig_overlap_source'};
}


=head2 overlap_distance_cutoff

 Title   : overlap_distance_cutoff
 Usage   : my $cutoff = $contig->overlap_distance_cutoff()
 Function: Gets or sets an integer which is used when building
           VirtualContigs.  If the distance in a contig overlap
           is greater than the cutoff, then the overlap will
           not be returned.
 Returns : value of overlap_distance_cutoff
 Args    : positive integer


=cut


sub overlap_distance_cutoff {
    my( $self, $cutoff ) = @_;
    
    if (defined $cutoff) {
	if( $cutoff !~ /^\d+$/ && $cutoff != -1 ) {
	    $self->throw("'$cutoff' is not an positive integer");
	    }
        $self->{'_overlap_distance_cutoff'} = $cutoff;
    }
    return $self->{'_overlap_distance_cutoff'};
}


sub is_golden {
   my $self = shift;

   if( defined $self->get_left_overlap || defined $self->get_right_overlap ) {
       return 1;
   } 
   return 0;
}


=head2 set_attribute

 Title   : set_attribute
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub set_attribute{
   my ($self,$tag,$value) = @_;

   if( !$self->dbobj->extension_tables ) {
       # only warn
       $self->warn("attempting to set attribute with no extension tables. Skipping");
   }
   my $id = $self->internal_id;

   my $sth = $self->dbobj->prepare("insert into contigext (contig_id,tag,value) VALUES ($id,'$tag','$value')");
   $sth->execute();

}

=head2 get_attribute

 Title   : get_attribute
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_attribute{
   my ($self,$tag) = @_;

   if( !$self->dbobj->extension_tables ) {
       # only warn
       $self->warn("attempting to set attribute with no extension tables. Skipping");
   }
   if( !defined $tag ) {
       $self->throw("no tag passed to get attribute");
   }

   my $id = $self->internal_id;

   my $sth = $self->dbobj->prepare("select value from contigext where contig_id = $id and tag = '$tag'");
   $sth->execute();
   my ($value) = $sth->fetchrow_array();

   return $value;
}

#
# Static golden path tables
#


=head2 chromosome

 Title   : chromosome
 Usage   : $self->chromosome($newval)
 Function: 
 Returns : value of chromosome
 Args    :


=cut

sub chromosome{
    my $self = shift;

    if( defined $self->_chromosome) { return $self->_chromosome;}

    my $id  = $self->internal_id;
    my $type = $self->dbobj->static_golden_path_type();
    my $sth = $self->dbobj->prepare("select chr_name from static_golden_path where raw_id = $id and type = '$type'");
    $sth->execute;
    my ($value) = $sth->fetchrow_array();
    if( !defined $value) { return undef; }
    $self->_chromosome($value);
    return $value;

}

=head2 _chromosome

 Title   : chromosome
 Usage   : $self->_chromosome($newval)
 Function: 
 Returns : value of _chromosome
 Args    : newvalue (optional)


=cut

sub _chromosome{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'_chromosome'} = $value;
    }
    return $self->{'_chromosome'};

}

=head2 fpc_contig_name

 Title   : fpc_contig_name
 Usage   : $self->fpc_contig_name()
 Function: 
 Returns : value of fpc_contig
 Args    :


=cut

sub fpc_contig_name {
    my $self = shift;

    if( defined $self->_fpc_contig) { return $self->_fpc_contig;}

    my $id  = $self->internal_id;
    my $type = $self->dbobj->static_golden_path_type();
    my $sth = $self->dbobj->prepare("select fpcctg_name from static_golden_path where raw_id = $id and type = '$type'");
    $sth->execute;
    my ($value) = $sth->fetchrow_array();
    if( !defined $value) { return undef; }
    $self->_fpc_contig($value);
    return $value;


}

=head2 _fpc_contig

 Title   : fpc_contig
 Usage   : $self->_fpc_contig($newval)
 Function: 
 Returns : value of _fpc_contig
 Args    : newvalue (optional)


=cut

sub _fpc_contig{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'_fpc_contig'} = $value;
    }
    return $self->{'_fpc_contig'};

}

=head2 chr_start

 Title   : chr_start
 Usage   : $self->chr_start($newval)
 Function: 
 Returns : value of chr_start
 Args    :


=cut

sub chr_start{
    my $self = shift;

    if( defined $self->_chr_start) { return $self->_chr_start;}

    my $id  = $self->internal_id;
    my $type = $self->dbobj->static_golden_path_type();
    my $sth = $self->dbobj->prepare("select chr_start from static_golden_path where raw_id = $id and type = '$type'");
    $sth->execute;
    my ($value) = $sth->fetchrow_array();
    if( !defined $value) { return undef; }
    $self->_chromosome($value);
    return $value;


}

=head2 _chr_start

 Title   : chr_start
 Usage   : $self->_chr_start($newval)
 Function: 
 Returns : value of _chr_start
 Args    : newvalue (optional)


=cut

sub _chr_start{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'_chr_start'} = $value;
    }
    return $self->{'_chr_start'};

}

=head2 chr_end

 Title   : chr_end
 Usage   : $self->chr_end($newval)
 Function: 
 Returns : value of chr_end
 Args    :


=cut

sub chr_end{
    my $self = shift;

    if( defined $self->_chr_end) { return $self->_chr_end;}

    my $id  = $self->internal_id;
    my $type = $self->dbobj->static_golden_path_type();
    my $sth = $self->dbobj->prepare("select chr_end from static_golden_path where raw_id = $id and type = '$type'");
    $sth->execute;
    my ($value) = $sth->fetchrow_array();
    if( !defined $value) { return undef; }
    $self->_chr_end($value);
    return $value;

}

=head2 _chr_end

 Title   : chr_end
 Usage   : $self->_chr_end($newval)
 Function: 
 Returns : value of _chr_end
 Args    : newvalue (optional)


=cut

sub _chr_end{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'_chr_end'} = $value;
    }
    return $self->{'_chr_end'};

}


=head2 fpc_contig_start

 Title   : fpc_contig_start
 Usage   : $self->fpc_contig_start($newval)
 Function: 
 Returns : value of fpc_contig_start
 Args    :


=cut

sub fpc_contig_start{
    my $self = shift;

    if( defined $self->_fpc_contig_start) { return $self->_fpc_contig_start;}

    my $id  = $self->internal_id;
    my $type = $self->dbobj->static_golden_path_type();
    my $sth = $self->dbobj->prepare("select fpcctg_start from static_golden_path where raw_id = $id and type = '$type'");
    $sth->execute;
    my ($value) = $sth->fetchrow_array();
    if( !defined $value) { return undef; }
    $self->_chromosome($value);
    return $value;


}

=head2 _fpc_contig_start

 Title   : fpc_contig_start
 Usage   : $self->_fpc_contig_start($newval)
 Function: 
 Returns : value of _fpc_contig_start
 Args    : newvalue (optional)


=cut

sub _fpc_contig_start{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'_fpc_contig_start'} = $value;
    }
    return $self->{'_fpc_contig_start'};

}

=head2 fpc_contig_end

 Title   : fpc_contig_end
 Usage   : $self->fpc_contig_end($newval)
 Function: 
 Returns : value of fpc_contig_end
 Args    :


=cut

sub fpc_contig_end{
    my $self = shift;

    if( defined $self->_fpc_contig_end) { return $self->_fpc_contig_end;}

    my $id  = $self->internal_id;
    my $type = $self->dbobj->static_golden_path_type();
    my $sth = $self->dbobj->prepare("select fpcctg_end from static_golden_path where raw_id = $id and type = '$type'");
    $sth->execute;
    my ($value) = $sth->fetchrow_array();
    if( !defined $value) { return undef; }
    $self->_fpc_contig_end($value);
    return $value;

}

=head2 _fpc_contig_end

 Title   : fpc_contig_end
 Usage   : $self->_fpc_contig_end($newval)
 Function: 
 Returns : value of _fpc_contig_end
 Args    : newvalue (optional)


=cut

sub _fpc_contig_end{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'_fpc_contig_end'} = $value;
    }
    return $self->{'_fpc_contig_end'};

}



=head2 static_golden_start

 Title   : static_golden_start
 Usage   : $self->static_golden_start($newval)
 Function: 
 Returns : value of static_golden_start
 Args    :


=cut

sub static_golden_start{
    my $self = shift;

    if( defined $self->_static_golden_start) { return $self->_static_golden_start;}
    my $id  = $self->internal_id;
    my $type = $self->dbobj->static_golden_path_type();
    my $sth = $self->dbobj->prepare("select raw_start from static_golden_path where raw_id = $id and type = '$type'");
    $sth->execute;
    my ($value) = $sth->fetchrow_array();
    if( !defined $value) { return undef; }
    $self->_static_golden_start($value);
    return $value;


}

=head2 _static_golden_start

 Title   : static_golden_start
 Usage   : $self->_static_golden_start($newval)
 Function: 
 Returns : value of _static_golden_start
 Args    : newvalue (optional)


=cut

sub _static_golden_start{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'_static_golden_start'} = $value;
    }
    return $self->{'_static_golden_start'};

}

=head2 static_golden_end

 Title   : static_golden_end
 Usage   : $self->static_golden_end($newval)
 Function: 
 Returns : value of static_golden_end
 Args    :


=cut

sub static_golden_end{
    my $self = shift;

    if( defined $self->_static_golden_end) { return $self->_static_golden_end;}
    my $id  = $self->internal_id;
    my $type = $self->dbobj->static_golden_path_type();
    my $sth = $self->dbobj->prepare("select raw_end from static_golden_path where raw_id = $id and type = '$type'");
    $sth->execute;
    my ($value) = $sth->fetchrow_array();
    if( !defined $value) { return undef; }
    $self->_static_golden_end($value);
    return $value;


}

=head2 _static_golden_end

 Title   : static_golden_end
 Usage   : $self->_static_golden_end($newval)
 Function: 
 Returns : value of _static_golden_end
 Args    : newvalue (optional)


=cut

sub _static_golden_end{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'_static_golden_end'} = $value;
    }
    return $self->{'_static_golden_end'};

}

=head2 static_golden_ori

 Title   : static_golden_ori
 Usage   : $self->static_golden_ori($newval)
 Function: 
 Returns : value of static_golden_ori
 Args    :


=cut

sub static_golden_ori{
    my $self = shift;

    if( defined $self->_static_golden_ori) { return $self->_static_golden_ori;}
    my $id  = $self->internal_id;
    my $type = $self->dbobj->static_golden_path_type();
    my $sth = $self->dbobj->prepare("select raw_ori from static_golden_path where raw_id = $id and type = '$type'");
    $sth->execute;
    my ($value) = $sth->fetchrow_array();
    if( !defined $value) { return undef; }
    $self->_static_golden_ori($value);
    return $value;


}

=head2 _static_golden_ori

 Title   : static_golden_ori
 Usage   : $self->_static_golden_ori($newval)
 Function: 
 Returns : value of _static_golden_ori
 Args    : newvalue (optional)


=cut

sub _static_golden_ori{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'_static_golden_ori'} = $value;
    }
    return $self->{'_static_golden_ori'};

}

=head2 static_golden_type

 Title   : static_golden_type
 Usage   : $self->static_golden_type($newval)
 Function: 
 Returns : value of static_golden_type
 Args    :


=cut

sub static_golden_type{
    my $self = shift;

    return $self->dbobj->static_golden_path_type();

}


=head2 is_static_golden

 Title   : is_static_golden
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub is_static_golden{
   my ($self,@args) = @_;

   if( defined $self->fpc_contig_name ) {
       return 1;
   }

}

1;






