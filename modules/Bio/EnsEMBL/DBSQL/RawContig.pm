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
    $self->overlap_distance_cutoff($overlap_distance_cutoff);

    return $self;
}


sub direct_new {
    my( $pkg, @args ) = @_;
    my $self = bless {}, $pkg;

    my (
        $dbobj,
        $id,
        $perlonlysequences,
        $overlap_distance_cutoff,
	$internal_id,
	$dna_id,
	$seq_version,
	$cloneid,
	$chr_start,
	$chr_end
        ) = $self->_rearrange([qw(
				  DBOBJ
				  ID
				  PERLONLYSEQUENCES
				  OVERLAP_DISTANCE_CUTOFF
				  INTERNAL_ID
				  DNA_ID
				  SEQ_VERSION
				  CLONEID
				  CHR_START
				  CHR_END
	    )], @args);

    $id    || $self->throw("Cannot make contig db object without id");
    $dbobj || $self->throw("Cannot make contig db object without db object");
    $dbobj->isa('Bio::EnsEMBL::DBSQL::Obj') || $self->throw("Cannot make contig db object with a $dbobj object");

    if( !$internal_id || !$dna_id || !defined($seq_version) || !$cloneid || !defined $chr_start || !defined $chr_end) {
	$self->throw("you don't have all the data to make a direct new [$internal_id,$dna_id,$seq_version,$cloneid,$chr_start,$chr_end]!");
    }

    $self->id($id);
    $self->dbobj($dbobj);
    $self->_got_overlaps(0);
    $self->internal_id($internal_id);
    $self->dna_id($dna_id);
    $self->seq_version($seq_version);
    $self->cloneid    ($cloneid);
    $self->perl_only_sequences($perlonlysequences);
    $self->overlap_distance_cutoff($overlap_distance_cutoff);
    $self->_chr_start($chr_start);
    $self->_chr_end($chr_end);

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

    my $query = 
    "    SELECT contig.internal_id
          , contig.dna
          , clone.embl_version
          , clone.id
          , contig.offset
        FROM dna
          , contig
          , clone
        WHERE contig.dna = dna.id
          AND contig.clone = clone.internal_id
          AND contig.id = '$id'
        ";


    my $sth = $self->dbobj->prepare($query);    
    my $res = $sth->execute();

    if (my $row = $sth->fetchrow_arrayref) {  
        $self->internal_id($row->[0]);
        $self->dna_id($row->[1]);
        $self->seq_version($row->[2]);
	$self->cloneid    ($row->[3]);
	$self->embl_offset    ($row->[4]);
    } else {
         $self->throw("Contig $id does not exist in the database or does not have DNA sequence");
    }

    return $self;
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
 
my $query="
        SELECT t.gene
        FROM transcript t,
             exon_transcript et,
             exon e
        WHERE e.contig = '$contig_id'
          AND et.exon = e.id
          AND t.id = et.transcript
        ";


   return  $self->_gene_query($query,$supporting);
}

=head2 get_old_Genes

 Title   : get_old_Genes
 Usage   : my @mapped_Genes=$rc->get_old_Genes 
 Function: Used to get out old Genes (not modifying coordinates)
 Returns : an array of Bio::EnsEMBL::Exon objects
 Args    : none


=cut

sub get_old_Genes {
    my ($self) = @_;

    #This method requires a connection to a crossmatch database
    if (!$self->_crossdb) { $self->throw("You need a crossmatch database to call get_old_Genes!");}
    my $crossdb = $self->_crossdb;

    
    #The crossdb should be holding onto old and new dbs, we need the old one here...
    my $old_db;
    eval {
	$old_db=$self->_crossdb->old_dbobj;
    }; 
    if ($@) {
	$self->throw("The crossmatch database has to hold the old dna database to be able to call get_old_Genes! $@");
    }
    my $oldcontig;
    eval {
	$oldcontig = $old_db->get_Contig($self->id);
    };

    #If the clone does not exist, these are really new Genes
    if ($@) {
	#print STDERR "Contig ".$self->id." doesn't exist in old db, returning empty array...\n";
	return ();
    }
   
    my @genes=$oldcontig->get_all_Genes();
    my $size=scalar (@genes);
    #print STDERR "Returning $size old Genes as they are for contig ".$self->id."\n"; 
    return @genes;
}

=head2 get_all_Exons

 Title   : get_all_Exons
 Usage   :
 Function: returns all exons for this contig
 Example :
 Returns : 
 Args    :

=cut

sub get_all_Exons {

    my ($self)=@_;


    my $contig_id=$self->id;


    my $query="SELECT e.id, e.seq_start,e.seq_end,e.strand,e.phase,e.created,e.modified 
               FROM   exon e,contig c 
               WHERE  c.internal_id=e.contig and c.id ='$contig_id'";

    my $sth = $self->dbobj->prepare ($query);
    $sth->execute;

    my ($id,$start,$end,$strand,$phase,$created,$modified);
    $sth->bind_columns (undef,\$id,\$start,\$end,\$strand,\$phase,\$created,\$modified);
    
    my @exons;
    while ($sth->fetch){
	my $exon=Bio::EnsEMBL::Exon->new;
	
	$exon->id($id);
	$exon->start($start);
	$exon->end($end);
	$exon->strand($strand);
	$exon->seqname($self->id);
	$exon->contig_id($self->id);
	$exon->phase($phase);
	$exon->created($created);
	$exon->modified($modified);
	$exon->sticky_rank(1);

	push @exons,$exon;
    }
    return @exons;
}


=head2 get_old_Exons

 Title   : get_old_Exons
 Usage   : my @mapped_exons=$rc->get_old_Exons 
 Function: Used to get out exons in new coordinates
 Returns : an array of Bio::EnsEMBL::Exon objects
 Args    : none


=cut

sub get_old_Exons {
    my ($self,$logfile) = @_;

    #This method requires a connection to a crossmatch database
    if (!$self->_crossdb) { $self->throw("You need a crossmatch database to call get_old_exons!");}
    my $crossdb = $self->_crossdb;

    
    #The crossdb should be holding onto old and new dbs, we need the old one here...
    my $old_db;
    eval {
	$old_db=$self->_crossdb->old_dbobj;
    }; 
    if ($@) {
	$self->throw("The crossmatch database has to hold the old dna database to be able to call get_old_exons! $@");
    }
    my $oldclone;
    my $oldcontig;
    eval {
	$oldclone = $old_db->get_Clone($self->cloneid);
    };

    #If the clone does not exist, these are really new exons
    if ($@) {
	return ();
    }
   
    my $newclone= $self->dbobj->get_Clone($self->cloneid);
    #If the clones have the same version, the underlying dna hasn't changed,
    #therefore we just return the old exons...
    if ($oldclone->embl_version == $newclone->embl_version) {
	my $oldcontig;
	eval {
	    $oldcontig = $oldclone->get_Contig($self->id);
	};
	if ($@) {
	    print STDERR "Clones with id ".$oldclone->id." have the same version in old and new db, but contig ".$self->id." is not there! (CLONE VERSION BUG)\n";
	    return ();
	}
	my @exons=$oldcontig->get_all_Exons();
	my $size=scalar (@exons);
	#print STDERR "Returning $size old exons as they are for contig ".$self->id." on clone ".$oldclone->id."\n"; 
	return @exons; 
    }
    #We get out a SymmetricContigFeatureContainer from the crossdb and use it     #to retrieve feature pairs for this contig, then sort them
    my $sfpc = $crossdb->get_SymmetricContigFeatureContainer;
    my @fp=$sfpc->get_FeaturePair_list_by_rawcontig_id($self->id,$newclone->embl_version);
    my @sorted_fp= sort { $a->start <=> $b->start} @fp;
    
    my %validoldcontigs;
    my %fphash;
    my @old_exons;
    foreach my $fp ( @sorted_fp ) {
	my $contigid = $fp->hseqname;
	my $oldcontig=$old_db->get_Contig($contigid);
	push @old_exons, $oldcontig->get_all_Exons;
	$validoldcontigs{$contigid} = $fp->hseqname;
	if( !exists $fphash{$fp->hseqname} ) {
	    $fphash{$fp->hseqname} = [];
	}
	push(@{$fphash{$fp->hseqname}},$fp);
    }
    #We now need to get all the Genes for this clone on the old case
    # now perform the mapping

    my @mapped_exons;
    my @unmapped;
    EXON:foreach my $exon (@old_exons) {
	my $mapped=0;
	foreach my $fp ( @{$fphash{$validoldcontigs{$exon->seqname}}} ) {
	    if( $fp->hstart < $exon->start && $fp->hend > $exon->start ) {
		if( $fp->strand == $fp->hstrand ) {
		    # straightforward mapping
		    $exon->seqname($fp->seqname);
		    $exon->contig_id($fp->seqname);
		    $exon->start($fp->start + $exon->start - $fp->hstart);
		    $exon->end($fp->start + $exon->end - $fp->hstart);
		} else {
		    # Grrr strand hell.
		    my $oldstart = $exon->start;
		    my $oldend   = $exon->end;
		    $exon->seqname($fp->seqname);
		    $exon->contig_id($fp->seqname);
		    $exon->start($fp->hend - ($oldstart - $fp->hend));  
		    $exon->end  ($fp->hend - ($oldend   - $fp->hend));
		    $exon->strand( -1 * $exon->strand);
		}
		$mapped=1;
		push (@mapped_exons,$exon);
		next EXON;
	    }
	}
	if ($mapped == 0) {
	    print $logfile "LOST EXON: ".$exon->id." (In get_old_Exons)\n"; 
	}
    }
    return @mapped_exons;		
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
   my @out;
   my $contig_id = $self->internal_id();   
   my %got;
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


   return  $self->_gene_query($query,$supporting);
}


sub _gene_query{

 my ($self, $query,$supporting) = @_;

 my @out;
 my $contig_id = $self->internal_id();   
 my %got;
 # prepare the SQL statement
 my $sth = $self->dbobj->prepare($query);
 
 my $res = $sth->execute();
 
 while (my $rowhash = $sth->fetchrow_hashref) { 
     
     if( ! exists $got{$rowhash->{'gene'}}) {  
	 
	 my $gene_obj = Bio::EnsEMBL::DBSQL::Gene_Obj->new($self->dbobj);
         my $gene;
         #PL: this may be overcautious, e.g. when called by get_GeneByType()
         eval {
             $gene = $gene_obj->get($rowhash->{'gene'}, $supporting);
         };
         if ($@) {
             $self->warn("In RawContig, tried to get gene ".$rowhash->{'gene'}." but couldn't (data bug?)\n");
         }
         else {
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

    if( $self->use_rawcontig_acc() ) {
	my $fplist = Bio::EnsEMBL::Ext::RawContigAcc::FeaturePairList_by_Score($id,$analysis_type,$score);
	@array = $fplist->each_FeaturePair;
	return @array;
    }

   my %analhash;

   #First of all, get all features that are part of a feature set with high enough score and have the right type

    my $statement = "SELECT feature.id, seq_start, seq_end, strand, feature.score, analysis, name, " .
		             "hstart, hend, hid, evalue, perc_id, phase, end_phase, fset, rank, fset.score " .
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
   
   my ($fid,$start,$end,$strand,$f_score,$analysisid,$name,$hstart,$hend,$hid,$evalue,$perc_id,$phase,$end_phase,$fset,$rank,$fset_score);
   my $seen = 0;
   
   # bind the columns
  
    $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$f_score,\$analysisid,\$name,\$hstart,\$hend,\$hid,\$evalue,\$perc_id,\$phase,\$end_phase,\$fset,\$rank,\$fset_score);
   
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
       $feature->seqname    ($self->id);
       $feature->start      ($start);
       $feature->end        ($end);
       $feature->strand     ($strand);
       $feature->source_tag ($name);
       $feature->primary_tag('similarity');
       $feature->id         ($fid);
       $feature->p_value    ($evalue)       if (defined $evalue);
       $feature->percent_id ($perc_id)      if (defined $perc_id);
       $feature->phase      ($phase)        if (defined $phase);
       $feature->end_phase  ($end_phase)    if (defined $end_phase);
       
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
        $statement = "SELECT feature.id, seq_start, seq_end, strand, score, analysis, name, hstart, hend, hid, evalue, perc_id, phase, end_phase " .
		     "FROM   feature, analysis " .
                     "WHERE  feature.id not in (" . $fset_id_str . ") " .
                     "AND    feature.score > '$score' " . 
                     "AND    feature.analysis = analysis.id " .
                     "AND    analysis.db = '$analysis_type' " .
                     "AND    feature.contig = '$id' ";
                                     
       $sth = $self->dbobj->prepare($statement);
       
   } else {
        $statement = "SELECT feature.id, seq_start, seq_end, strand, score, analysis, name, hstart, hend, hid, evalue, perc_id, phase, end_phase " .
		     "FROM   feature, analysis " .
                     "WHERE  feature.score > '$score' " . 
                     "AND    feature.analysis = analysis.id " .
                     "AND    analysis.db = '$analysis_type' " .
                     "AND    feature.contig = '$id' ";
                     
                     
                     
       $sth = $self->dbobj->prepare($statement);
   }

   $sth->execute();

   # bind the columns
   $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$f_score,\$analysisid,\$name,\$hstart,\$hend,\$hid,\$evalue,\$perc_id,\$phase,\$end_phase);
   
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
	   $out->seqname    ($self->id);
	   $out->start      ($start);
	   $out->end        ($end);
	   $out->strand     ($strand);
	   $out->source_tag ($name);
	   $out->primary_tag('similarity');
	   $out->id         ($fid);
       $out->p_value    ($evalue)    if (defined $evalue);
       $out->percent_id ($perc_id)   if (defined $perc_id); 
       $out->phase      ($phase)     if (defined $phase);    
       $out->end_phase  ($end_phase) if (defined $end_phase);

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

sub get_all_SimilarityFeatures {
   my ($self) = @_;

   my @array;
   my @fps;

   my $id     = $self->internal_id();
   my $length = $self->length();

   my @genscan = $self->get_all_PredictionFeatures;

   push(@array,@genscan);
   my %analhash;

   #Then get the rest of the features, i.e. featurepairs and single features that are not part of a fset
   my ($fid,$start,$end,$strand,$f_score,$analysisid,$name,$hstart,$hend,$hid,$evalue,$perc_id,$phase,$end_phase);

   my $sth = $self->dbobj->prepare("select id,seq_start,seq_end,strand,score,analysis,name,hstart,hend,hid, evalue, perc_id, phase, end_phase ".
				"from feature where contig = $id");
   
   $sth->execute();

   # bind the columns
   $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$f_score,\$analysisid,\$name,\$hstart,\$hend,\$hid,\$evalue,\$perc_id,\$phase,\$end_phase);
   
   FEAT: while($sth->fetch) {
       my $out;
       my $analysis;
              
#       print STDERR  "\nID $fid, START $start, END $end, STRAND $strand, SCORE $f_score, EVAL $evalue, PHASE $phase, EPHASE $end_phase, ANAL $analysisid, FSET $fset\n";
       
       if (!$analhash{$analysisid}) {
	   
	   my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->dbobj);
	   eval {
	     $analysis = $feature_obj->get_Analysis($analysisid);
	     $analhash{$analysisid} = $analysis;
	   };
	   if ($@) {
	     print STDERR "Error fetching analysis $analysisid. Skipping [$@]\n";
	     next FEAT;
	   }
	   
       } else {
	   $analysis = $analhash{$analysisid};
       }
       
       if( !defined $name ) {
	 $name = 'no_source';
       } elsif ($name eq "genscan") {
	 next FEAT;
       }
       
       if( $hid ne '__NONE__' ) {
	   # is a paired feature
	   # build EnsEMBL features and make the FeaturePair
	 
	 $out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();


	   $out->set_all_fields($start,$end,$strand,$f_score,$name,'similarity',$self->id,
				            $hstart,$hend,1,$f_score,$name,'similarity',$hid, $evalue, $perc_id, $phase, $end_phase);

	   $out->analysis    ($analysis);
	   $out->id($fid);
	   # see comment below
	   #$out->id          ($hid);              # MC This is for Arek - but I don't
	                                          #    really know where this method has come from.
       } else {
	 $out = new Bio::EnsEMBL::SeqFeature;
	 $out->seqname    ($self->id);
	 $out->start      ($start);
	 $out->end        ($end);
	 $out->strand     ($strand);
	 $out->source_tag ($name);
	 $out->primary_tag('similarity');
	 $out->id         ($fid);
	 $out->p_value    ($evalue)    if (defined $evalue);
	 $out->percent_id ($perc_id)   if (defined $perc_id); 
	 $out->phase      ($phase)     if (defined $phase);    
	 $out->end_phase  ($end_phase) if (defined $end_phase); 
	   $out->raw_seqname   ($self->id);

	   if( defined $f_score ) {
	       $out->score($f_score);
	   }
	   $out->analysis($analysis);
       }
       # Final check that everything is ok.
       $out->validate() || $out->throw("Invalid data in $out");
       if( $out->can('attach_seq') ) {
	   $out->attach_seq($self->primary_seq);
       }

      push(@fps,$out);
      
   }
   
   push(@array,@fps);

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
	   $out->id($fid);
       } else {
	   $self->warn("Repeat feature does not have a hid. bad news....");
       }
       
       $out->validate();

      push(@array,$out);
  }
 
   return @array;
}                                       # get_all_RepeatFeatures

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
    my ($self)=@_;
  
    my $id = $self->internal_id;
    my @markers;

    eval {
        require Bio::EnsEMBL::Map::MarkerFeature;
        
        my $statement="SELECT f.seq_start, f.seq_end, f.score, f.strand, f.name, 
	                  f.hstart, f.hend, f.hid, f.analysis 
	           FROM   feature f, analysis a 
		   WHERE  f.contig='$id' 
                   AND    f.analysis = a.id and a.db='mapprimer'";
        
        @markers=$self->_create_MarkerFeatures($statement);	
	
    };

    if( $@ ) {
        print STDERR ("Problems retrieving map data. Most likely not connected to maps db\n$@\n" );
    }
    
    return @markers;
}                                       # get_MarkerFeatures


=head2 get_landmark_MarkerFeatures

  Title   : get_landmark_MarkerFeatures 
  Usage   : @fp = $contig->get_landmark_MarkerFeatures; 
  Function: Gets MarkerFeatures with identifiers like D8S509. 
            MarkerFeatures can be asked for a Marker. 
            Its assumed, that when you can get MarkerFeatures, then you can 
            get the Map Code as well.
  Example : - 
  Returns : -
  Args : -

=cut


sub get_landmark_MarkerFeatures {
      my ($self)=@_;
      
      
      my $dbname=$self->dbobj->dbname;
      my $mapsdbname=$self->dbobj->mapdbname;
      
      my $id = $self->internal_id;
      my @markers;
      
      
      eval {
          require Bio::EnsEMBL::Map::MarkerFeature;
          
          my $statement="SELECT f.seq_start, f.seq_end, f.score, f.strand, f.name, 
	                  f.hstart, f.hend, s.name, f.analysis 
                   FROM   $dbname.feature f, $dbname.analysis a, 
                          $mapsdbname.MarkerSynonym s,$mapsdbname.Marker m 
		   WHERE  f.contig='$id' 
                   AND    f.analysis = a.id 
                   AND    a.db='mapprimer'
                   AND    m.marker=s.marker 
                   AND    f.hid=m.marker 
                   AND    s.name regexp '^D[0-9,X,Y][0-9]?S'";
          
          @markers=$self->_create_MarkerFeatures($statement);	
          
      };
      
      if( $@ ) {
          print STDERR ("Problems retrieving map data. Most likely not connected to maps db\n$@\n" );
      }

      return @markers;
}                                       # get_landmark_MarkerFeatures



sub _create_MarkerFeatures
{
my ($self,$statement)=@_;

my @result; 
my $analysis;
my %analhash;
    

my $sth = $self->dbobj->prepare($statement);
$sth->execute;

my ($start, $end, $score, $strand, $hstart, 
    $name, $hend, $hid, $analysisid );


$sth->bind_columns
    ( undef, \$start, \$end, \$score, \$strand, \$name, 
      \$hstart, \$hend, \$hid, \$analysisid);


while( $sth->fetch ) {
    
    my ( $out, $seqf1, $seqf2 );
    
    if (!$analhash{$analysisid}) {
	
	my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->dbobj);
	$analysis = $feature_obj->get_Analysis($analysisid);
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

return @result;

}                                       # _create_MarkerFeatures

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
   my $query = "select f.id,f.seq_start,f.seq_end,f.strand,f.score,f.evalue,f.perc_id,f.phase,f.end_phase,f.analysis,f.hid ". 
       "from feature f where contig = $id and name = 'genscan'";

   my $sth = $self->dbobj->prepare($query);
   
   $sth->execute();
   
   my ($fid,$start,$end,$strand,$score,$evalue,$perc_id,$phase,$end_phase,$analysisid,$hid);
   
   # bind the columns
   $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$score,\$evalue,\$perc_id,\$phase,\$end_phase,\$analysisid,\$hid);
  
   $previous = -1;
   my $current_fset;
   my $count;

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


       if( $hid eq "Initial Exon" || $hid eq "Single Exon" || $previous eq "Single Exon" || $previous eq "Terminal Exon" || $previous == -1) {
	   $current_fset = new Bio::EnsEMBL::SeqFeature;
	   $current_fset->source_tag('genscan');
	   $current_fset->primary_tag('prediction');
	   $current_fset->analysis($analysis);
	   $current_fset->seqname($self->id);
	   $current_fset->id($count);
           $count++;
	   $current_fset->raw_seqname($self->id);
	   push(@array,$current_fset);
       }

       $out = new Bio::EnsEMBL::SeqFeature;
       
       $out->seqname   ($self->id);
       $out->raw_seqname($self->id);

       $out->start     ($start);
       $out->end       ($end);
       $out->strand    ($strand);
       $out->p_value   ($evalue)    if (defined $evalue);
       $out->percent_id($perc_id)   if (defined $perc_id); 
       $out->phase     ($phase)     if (defined $phase);    
       $out->end_phase ($end_phase) if (defined $end_phase);
        

	my $query="select fset from fset_feature where feature=$fid"; 
	my $sth = $self->dbobj->prepare($query);
   	$sth->execute();
	my $arr_ref=$sth->fetchrow_arrayref;

	$fsetid=$arr_ref->[0];

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
       $previous = $hid;
  }
 
   return @array;
}


=head2 get_genscan_peptides

 Title   : get_genscan_peptides
 Usage   : 
 Function: Returns genscan predictions as peptides
 Example :
 Returns : 
 Args    : 


=cut

#have written this to use the new phase and end_phase tag in SeqFeature
#Therefore this won't work with the older features, after all the old system was pretty ropey.
sub get_genscan_peptides {
    my ($self) = @_;
    my @transcripts;
    foreach my $fset ($self->get_all_PredictionFeatures) {
	my $trans = &Bio::EnsEMBL::DBSQL::Utils::fset2transcript($fset,$self);
	push(@transcripts,$trans);
    }

    return @transcripts;
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
   
   $acc = $self->cloneid();

   my $embl_offset = $self->embl_offset();

   foreach my $extf ( $self->dbobj->_each_ExternalFeatureFactory ) {

       if( $extf->can('get_Ensembl_SeqFeatures_contig') ) {

	   my @tmp = $extf->get_Ensembl_SeqFeatures_contig($self->internal_id,
							   $self->seq_version,
							   1,
							   $self->length,
							   $self->id);

	   push(@out,@tmp);
       }
       if( $extf->can('get_Ensembl_SeqFeatures_clone') ) {
       
	   foreach my $sf ( $extf->get_Ensembl_SeqFeatures_clone($acc,$self->seq_version,$self->embl_offset,$self->embl_offset+$self->length()) ) {

	       my $start = $sf->start - $embl_offset+1;
	       my $end   = $sf->end   - $embl_offset+1;
	       $sf->start($start);
	       $sf->end($end);
	       push(@out,$sf);
	   }
       }
   }
   my $id = $self->id();
   foreach my $f ( @out ) {
       $f->seqname($id);
   }

   return @out;

}                                       # get_all_ExternalFeatures


=head2 get_all_ExternalGenes

 Title   : get_all_ExternalGenes 
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_ExternalGenes {
   my ($self) = @_;
   my @out;
   my $acc;

   $acc = $self->cloneid();
   my $embl_offset = $self->embl_offset();

   foreach my $extf ( $self->dbobj->_each_ExternalFeatureFactory ) {
       if( $extf->can('get_Ensembl_Genes_clone') ) {
	   my @genes = $extf->get_Ensembl_Genes_clone($acc);
	   foreach my $gene (@genes){
	   foreach my $exon ( $gene->all_Exon_objects ) {
	       $exon->start($exon->start - $embl_offset+1);
	       $exon->end($exon->end - $embl_offset+1);
	   }
	   push(@out,@genes);
       }


       }
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

=head2 old_chromosome

 Title   : chromosome
 Usage   : $chr = $contig->chromosome( [$chromosome] )
 Function: get/set the chromosome of the contig.
 Example :
 Returns : the chromsome object
 Args    :

=cut

sub old_chromosome {
   my ($self,$chromosome ) = @_;
   my $id= $self->internal_id();

   if( defined( $chromosome )) {
       $self->{_chromosome} = $chromosome;
   } else {
       if (! defined ($self->{_chromosome})) {
	   my $sth = $self->dbobj->prepare("select chromosomeId from contig where internal_id = \"$id\" ");
	   $sth->execute();
	   
	   my $rowhash = $sth->fetchrow_hashref();
	   my $chrId = $rowhash->{'chromosomeId'};
	   $self->{_chromosome} = Bio::EnsEMBL::Chromosome->get_by_id
	       ( $chrId );
       }
   }
   return $self->{_chromosome};
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
   my ( $self, $arg )  = @_;
   my $id = $self->id();
   if( defined $arg ) {
     $self->{_embl_offset} = $arg;
     return;
   }

   if( defined $self->{_embl_offset} ) {
     return $self->{_embl_offset};
   } else {
     my $sth = $self->dbobj->prepare("select offset from contig where id = \"$id\" ");
     $sth->execute();
     my $rowhash = $sth->fetchrow_hashref();
     $self->{_embl_offset} = $rowhash->{'offset'};
     return $self->{_embl_offset};
   }
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

=head2 crossdb

 Title   : crossdb
 Usage   :
 Function:
 Example :
 Returns : The Bio::EnsEMBL::DBSQL::CrossMatchAdaptor object
 Args    :


=cut

sub _crossdb {
   my ($self,$arg) = @_;

   return $self->dbobj->_crossdb;
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

{ # (this brace is the beginning of a block that results in static compilation
  # of the SQL queries)

    # Certainly worth explaining here.
    #
    # The overlap type indicates which end on the two contigs this overlap is.
    # left means 5', right means 3'. There are four options. From these four
    # options we can figure out
    #    a) whether this overlap is on the 5' or the 3' of our contig
    #    b) which polarity the overlap on the next contig is 
    #
    # Polarity indicates whether the sequence is being read in the same 
    # direction as this contig. 
    # 
    # The sequence has to be appropiately versioned otherwise this gets complicated
    # in the update scheme.


    # Doing two queries seems to be quickest
    # Statements like:
    #   c.dna = o.dna_b_id OR c.dna = o.dna_a_id
    # make queries inordinately slow because MySQL
    # doesn't use indices on OR statements.
    my @queries = (
       q{SELECT c.id sister_id
          , o.contig_b_position sister_pos
          , o.contig_a_position self_pos
          , o.overlap_type
          , o.overlap_size
          , o.overlap_source
        FROM contigoverlap o
          , contig c
        WHERE c.dna = o.dna_b_id
          AND dna_a_id = ?},

       q{SELECT c.id sister_id
          , o.contig_a_position sister_pos
          , o.contig_b_position self_pos
          , o.overlap_type
          , o.overlap_size
          , o.overlap_source
        FROM contigoverlap o
          , contig c
        WHERE c.dna = o.dna_a_id
          AND dna_b_id = ?},
        );

    sub get_all_Overlaps {
        my ($self) = @_;

        my $id      = $self->dna_id();
        my $version = $self->seq_version();

        # Doing two queries seems to be quickest
        # Statements like:
        #   c.dna = o.dna_b_id OR c.dna = o.dna_a_id
        # seem to make queries inordinately slow.
        my $overlap_source_sub  = $self->dbobj->contig_overlap_source();
        my $overlap_cutoff      = $self->overlap_distance_cutoff();
        if( !defined $overlap_cutoff) {
	    $overlap_cutoff = 100000;
	}

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
	    #print(STDERR "Using new features for $newid\n");
	    push(@out,@{$newfeatures{$newid}});
	} else {
	    #print(STDERR "Using old features for $newid\n");
	    push(@out,@{$oldfeatures{$newid}});
	}
    }
    return @out;
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

sub fpc_contig_start {
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
 Returns : value of static_golden_start (in RawContig coordinates)
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
 Returns : value of static_golden_end (in RawContig coordinates)
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

=head2 use_rawcontig_acc

 Title   : use_rawcontig_acc
 Usage   : $obj->use_rawcontig_acc($newval)
 Function: 
 Returns : value of use_rawcontig_acc
 Args    : newvalue (optional)


=cut

sub use_rawcontig_acc{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'use_rawcontig_acc'} = $value;
    }
    return $obj->{'use_rawcontig_acc'};

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

=head2 SeqI implementing functions

=head2 species

 Title   : species
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub species{
   my ($self,@args) = @_;

   return undef;
}



1;
