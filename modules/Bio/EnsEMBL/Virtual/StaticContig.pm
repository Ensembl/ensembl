
#
# Ensembl module for Bio::EnsEMBL::Virtual::StaticContig
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Virtual::StaticContig - Virtual Contig specific to SQL databases with static_golden_path

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

This object inherits from Bio::EnsEMBL::Virtual::Contig, which in turn
inherits from Bio::EnsEMBL::DB::ContigI.  See
Bio::EnsEMBL::Virtual::Contig for examples of usage.

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Virtual::StaticContig;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI

use Bio::Root::RootI;
use Bio::EnsEMBL::Virtual::Contig;
use Bio::Annotation;
use Bio::Annotation::DBLink;
use Bio::EnsEMBL::VirtualGene;
use Bio::EnsEMBL::Utils::Eprof qw(eprof_start eprof_end);
use Bio::EnsEMBL::WebTranscript;

@ISA = qw(Bio::EnsEMBL::Virtual::Contig);


my $static_number = 0;

# overrides Bio::EnsEMBL::Virtual::Contig::new
sub new {
    my ($class,$global_start,$vc_start_position,$global_end,@contigs) = @_;
    
    my $self = {};
    bless $self,$class;
    $self->_make_datastructures(); # back to virtual contig


    # EMBL dumping support
    $self->{'date'} = [];
    $self->{_anal_hash}=();
    $self->annotation( Bio::Annotation->new());
    $self->{'additional_seqf'} = [];
   
    if( scalar(@contigs) == 0 ) {
	if( $global_end == -1 ) {
	    $self->throw("Cannot build a virtual contig from no raw contigs. Probably an error in the call to get raw contigs");
	} else {
	    # we have an all gap contig.
	    $self->_global_start($global_start);
	    $self->_global_end($global_end);
	    
	    $self->_vmap->length($global_end - $global_start+1);
	    $self->length($global_end - $global_start+1);
	}
    }


    # this loop is no longer easy because it has to deal with
    # right and left truncations now <sigh>

    foreach my $rc ( @contigs ) {

	my $rc_start;
	my $chr_start;
	my $chr_end;

	if( $rc->chr_start < $global_start ) {
	    if( $rc->static_golden_ori ==1 ) {
		# move start
		$rc_start = $rc->static_golden_start + ($global_start - $rc->chr_start);
	    } else {
		# don't need to move start, unless end - handled below
		$rc_start = $rc->static_golden_start;
	    }

	    $chr_start = $global_start;
	} else {
	    $rc_start = $rc->static_golden_start;
	    $chr_start = $rc->chr_start;
	}

	if( $global_end != -1 && $rc->chr_end > $global_end ) {
	    if( $rc->static_golden_ori == -1 ) {
		# need to move rstart
		$rc_start = $rc->static_golden_start + ($rc->chr_end - $global_end);
	    }
	    $chr_end = $global_end;
	} else {
	    $chr_end = $rc->chr_end;
	}


	$self->_vmap->create_MapContig($rc,
				       $chr_start - $global_start+$vc_start_position,
				       $chr_end   - $global_start+$vc_start_position,
				       $rc_start,
				       $rc->static_golden_ori);
    }

    $self->_global_start($global_start);

    # needs to handle overhangs...
    if( $global_end == -1 ) {
	@contigs = $self->_vmap->each_MapContig;
	my $last = pop @contigs;
	$self->_vmap->length($last->end);
	$self->length($last->end);
	$self->_global_end($last->end+$global_start);
    } else {
	$self->_vmap->length($global_end - $global_start+1);
	$self->length($global_end - $global_start+1);
	$self->_global_end($global_end);
    }
	
    
    $self->id("static".$static_number++);
    
    return $self;
}                                       # new

=head2 get_all_SimilarityFeatures

 Title   : get_all_SimilarityFeatures
 Usage   : foreach my $sf ( $contig->get_all_SimilarityFeatures ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SimilarityFeatures {
   my ($self) = @_;

   my $glob_start=$self->_global_start;
   my $glob_end=$self->_global_end;
   my $chr_name=$self->_chr_name;
   my $idlist  = $self->_raw_contig_id_list();
   
   unless ($idlist){
       return ();
   }
   
   my    $statement = "SELECT f.id, 
                        IF     (sgp.raw_ori=1,(f.seq_start+sgp.chr_start-sgp.raw_start-$glob_start),
                                 (sgp.chr_start+sgp.raw_end-f.seq_end-$glob_start)) as start,  
                        IF     (sgp.raw_ori=1,(f.seq_end+sgp.chr_start-sgp.raw_start-$glob_start),
                                 (sgp.chr_start+sgp.raw_end-f.seq_start-$glob_start)), 
                        IF     (sgp.raw_ori=1,f.strand,(-f.strand)),
                                f.score,f.analysis, f.name, f.hstart, f.hend, f.hid 
		        FROM   feature f, analysis a,static_golden_path sgp
                        WHERE  f.analysis = a.id 
                        AND    sgp.raw_id = f.contig
                        AND    f.contig in $idlist
                        AND    sgp.chr_end >= $glob_start 
		        AND    sgp.chr_start <=$glob_end 
		        AND    sgp.chr_name='$chr_name' 
                        AND    a.gff_feature = 'similarity'
                        ORDER  by start";
   

    my  $sth = $self->dbobj->prepare($statement);    
    $sth->execute(); 
    
    
    my ($fid,$start,$end,$strand,$f_score,$analysisid,$name,
	$hstart,$hend,$hid,$fset,$rank,$fset_score,$contig);
    
    $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$f_score,
                       \$analysisid,\$name,\$hstart,\$hend,\$hid);
    
    
    my @features;
    
    my $out;
    my %analhash;
    my $length=$self->length;
  FEATURE: 

    while($sth->fetch) {

	if (($end > $length) || ($start < 1)) {
	    next;
	}
	
	my @args=($fid,$start,$end,$strand,$f_score,$analysisid,$name,$hstart,$hend,$hid);

	my $analysis=$self->_get_analysis($analysisid);
	
	if( !defined $name ) {
	    $name = 'no_source';
	}
	
	$out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();   
	$out->set_all_fields($start,$end,$strand,$f_score,$name,'similarity',$contig,
			     $hstart,$hend,1,$f_score,$name,'similarity',$hid);

	$out->analysis($analysis);
        $out->id      ($hid);
	push(@features,$out);
    }
  
   return @features;

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
    my ($self, $analysis_type, $score, $bp) = @_;
    
    $self->throw("Must supply analysis_type parameter") unless $analysis_type;
    $self->throw("Must supply score parameter") unless $score;

    if( $self->_use_cext_get() ) {
	return $self->_cext_get_all_SimilarityFeatures_type($analysis_type);
    }
    
    my $glob_start = $self->_global_start;
    my $chr_name   = $self->_chr_name;
    my $idlist     = $self->_raw_contig_id_list;
    my $type       = $self->dbobj->static_golden_path_type;
    
    unless ($idlist){
	return ();
    }

    my    $statement = "SELECT f.id, 
                        IF     (sgp.raw_ori=1,(f.seq_start+sgp.chr_start-sgp.raw_start-$glob_start),
                                 (sgp.chr_start+sgp.raw_end-f.seq_end-$glob_start)) as start,  
                        IF     (sgp.raw_ori=1,(f.seq_end+sgp.chr_start-sgp.raw_start-$glob_start),
                                 (sgp.chr_start+sgp.raw_end-f.seq_start-$glob_start)), 
                        IF     (sgp.raw_ori=1,f.strand,(-f.strand)),
                                f.score,f.analysis, f.name, f.hstart, f.hend, f.hid 
		        FROM   feature f, analysis a,static_golden_path sgp
                        WHERE  f.score > $score
                        AND    f.analysis = a.id 
                        AND    sgp.raw_id = f.contig
                        AND    f.contig in $idlist
		        AND    a.db = '$analysis_type'  
                        AND    sgp.type = '$type'
		        AND    sgp.chr_name = '$chr_name' 
                        ORDER  by start";
    
    my  $sth = $self->dbobj->prepare($statement);    
    $sth->execute(); 
    
    
    my ($fid,$start,$end,$strand,$f_score,$analysisid,$name,
	$hstart,$hend,$hid,$fset,$rank,$fset_score,$contig);
    
    $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$f_score,
                       \$analysisid,\$name,\$hstart,\$hend,\$hid);
    
    
    my @features;
    
    my $out;
    my $length=$self->length;
  FEATURE: 

    while($sth->fetch) {

	if (($end > $length) || ($start < 1)) {
	    next;
	}
	
	my @args=($fid,$start,$end,$strand,$f_score,$analysisid,$name,$hstart,$hend,$hid);
	
	#exclude contained features
	if(  defined $bp && defined $out && $end < $out->end ) { next; }

	#Glob overlapping and close features
	if ( defined $bp && defined $out && $out->end+$bp >= $start ) {
		
	    # reset previous guys end to end
	    $out->end($end);
	    next;
	}

	
	my $analysis=$self->_get_analysis($analysisid);
	
	if( !defined $name ) {
	    $name = 'no_source';
	}
	
	$out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();   
	$out->set_all_fields($start,$end,$strand,$f_score,$name,'similarity',$contig,
			     $hstart,$hend,1,$f_score,$name,'similarity',$hid);

	$out->analysis($analysis);
        $out->id      ($hid);
	push(@features,$out);
    }
    return @features;
}


#
# This is for lucash's unigene only feature table...
#

sub get_all_unigene_features {
    my ($self, $bp) = @_;
    
    my $glob_start = $self->_global_start;
    my $chr_name   = $self->_chr_name;
    my $idlist     = $self->_raw_contig_id_list;
    my $type       = $self->dbobj->static_golden_path_type;
    
    unless ($idlist){
	return ();
    }

    my    $statement = "SELECT f.id, 
                        IF     (sgp.raw_ori=1,(f.seq_start+sgp.chr_start-sgp.raw_start-$glob_start),
                                 (sgp.chr_start+sgp.raw_end-f.seq_end-$glob_start)) as start,  
                        IF     (sgp.raw_ori=1,(f.seq_end+sgp.chr_start-sgp.raw_start-$glob_start),
                                 (sgp.chr_start+sgp.raw_end-f.seq_start-$glob_start)), 
                        IF     (sgp.raw_ori=1,f.strand,(-f.strand)),
                                f.score,f.analysis, f.name, f.hstart, f.hend, f.hid 
		        FROM   unigene_feature f, static_golden_path sgp
                        WHERE  sgp.raw_id = f.contig
                        AND    f.contig in $idlist
		        AND    sgp.chr_name = '$chr_name' 
                        ORDER  by start";
    
    my  $sth = $self->dbobj->prepare($statement);    
    $sth->execute(); 
    
    
    my ($fid,$start,$end,$strand,$f_score,$analysisid,$name,
	$hstart,$hend,$hid,$fset,$rank,$fset_score,$contig);
    
    $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$f_score,
                       \$analysisid,\$name,\$hstart,\$hend,\$hid);
    
    
    my @features;
    
    my $out;
    my $length=$self->length;
  FEATURE: 

    while($sth->fetch) {

	if (($end > $length) || ($start < 1)) {
	    next;
	}
	
	my @args=($fid,$start,$end,$strand,$f_score,$analysisid,$name,$hstart,$hend,$hid);
	
	#exclude contained features
	if(  defined $bp && defined $out && $end < $out->end ) { next; }

	#Glob overlapping and close features
	if ( defined $bp && defined $out && $out->end+$bp >= $start ) {
		
	    # reset previous guys end to end
	    $out->end($end);
	    next;
	}

	
	my $analysis=$self->_get_analysis($analysisid);
	
	if( !defined $name ) {
	    $name = 'no_source';
	}
	
	$out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();   
	$out->set_all_fields($start,$end,$strand,$f_score,$name,'similarity',$contig,
			     $hstart,$hend,1,$f_score,$name,'similarity',$hid);

	$out->analysis($analysis);
        $out->id      ($hid);
	push(@features,$out);
    }
    return @features;
}



=head2 get_all_SimilarityFeatures_above_pid

 Title   : get_all_SimilarityFeatures_above_pid
 Usage   : foreach my $sf ( $contig->get_all_SimilarityFeatures_above_pid(analysis_type, pid) ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SimilarityFeatures_above_pid{
    my ($self, $pid) = @_;
    
    $self->throw("Must supply pid parameter") unless $pid;
    
    my $glob_start=$self->_global_start;
    my $glob_end=$self->_global_end;
    my $chr_name=$self->_chr_name;
    
    my ($fid,$start,$end,$strand,$f_score,$analysisid,$name,
	$hstart,$hend,$hid,$fset,$rank,$fset_score,$contig,$chr_start,$chr_end,$raw_ori);
    

    # this needs to be rewritten properely EB
    
    my $type = $self->dbobj->static_golden_path_type;

    my $statement = "SELECT ".
      "f.id,
   if(s.raw_ori=1,(f.seq_start-s.raw_start+s.chr_start),(s.chr_start+s.raw_end-f.seq_start)),
   if(s.raw_ori=1,(f.seq_end  -s.raw_start+s.chr_start),(s.chr_start+s.raw_end-f.seq_end)),
   f.strand * s.raw_ori,
   f.score,f.analysis, f.name, f.hstart, f.hend, f.hid
               FROM feature f, static_golden_path s 
               WHERE f.perc_id > $pid 
               AND   s.raw_id  = f.contig
               AND NOT (s.chr_start > $glob_end) 
               AND NOT (s.chr_end < $glob_start) 
               AND   f.seq_start > s.raw_start 
               AND   f.seq_end   < s.raw_end
               AND   s.chr_name  = '$chr_name'
               AND   s.type = '$type'";


    #my    $statement = "SELECT f.id, f.seq_start+s.chr_start,f.seq_end+s.chr_start, 
    #                           f.strand, f.score,f.analysis, f.name, f.hstart, f.hend, f.hid,
    #                           s.chr_start,s.chr_end,s.raw_ori
#		        FROM   feature f, static_golden_path s
#                        WHERE  f.perc_id > $pid  
#                        AND    s.raw_id = f.contig
#                        AND    s.chr_end >= $glob_start 
#		        AND    s.chr_start <=$glob_end 
#		        AND    s.chr_name='$chr_name'";
    
    my  $sth = $self->dbobj->prepare($statement);    
    $sth->execute(); 

    $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$f_score,\$analysisid,
		       \$name,\$hstart,\$hend,\$hid);


    my @array;
    my %analhash;
    my $out;
    
  FEATURE: while($sth->fetch) {
      my $out;
      my $analysis;

      # clip and map to vc coordinates

      if ($start>=$glob_start && $end<=$glob_end){

    	$start=$start-$glob_start;
    	$end=$end-$glob_start;

	# create features

	  if (!$analhash{$analysisid}) 
	  {
	      my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->dbobj);
	      $analysis = $feature_obj->get_Analysis($analysisid);
	      $analhash{$analysisid} = $analysis;	   
	  } 
	  else {$analysis = $analhash{$analysisid};}
	  
	  if( !defined $name ) {
	      $name = 'no_source';
	  }
	  $out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();   
	  $out->set_all_fields($start,$end,$strand,$f_score,$name,'similarity',$contig,
				$hstart,$hend,1,$f_score,$name,'similarity',$hid);
	  $out->analysis    ($analysis);
	  $out->id          ($hid);              
	   
#	  $out = new Bio::EnsEMBL::SeqFeature;
	  $out->seqname   ($self->id);
	  $out->start     ($start);
	  $out->end       ($end);
	  $out->strand    ($strand);
	  $out->source_tag($name);
	  $out->primary_tag('similarity');
	  $out->id         ($hid);
	  
	  if( defined $f_score ) {
	      $out->score($f_score);
	  }
	  $out->analysis($analysis);
	  
	  $out->validate();

	  push(@array,$out);       
      }
  }

    return @array;
}



=head2 get_all_RepeatFeatures

 Title   : get_all_RepeatFeatures
 Usage   : foreach my $sf ( $contig->get_all_RepeatFeatures ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_RepeatFeatures {
    my ($self,$bp) = @_;
    
    my $glob_start=$self->_global_start;
    my $glob_end=$self->_global_end;
    my $chr_name=$self->_chr_name;
    my $length=$self->length;
    my $idlist  = $self->_raw_contig_id_list();
    
    $self->throw("No chromosome name") unless defined $chr_name;

    unless ($idlist){
	return ();
    }

    my $type = $self->dbobj->static_golden_path_type;

    my $statement = "SELECT rf.id,
                     IF     (sgp.raw_ori=1,(rf.seq_start+sgp.chr_start-sgp.raw_start-$glob_start),
                            (sgp.chr_start+sgp.raw_end-rf.seq_end-$glob_start)) as start,                                        
                     IF     (sgp.raw_ori=1,(rf.seq_end+sgp.chr_start-sgp.raw_start-$glob_start),
                            (sgp.chr_start+sgp.raw_end-rf.seq_start-$glob_start)), 
                     IF     (sgp.raw_ori=1,rf.strand,(-rf.strand)),                         
                            rf.score,rf.analysis,rf.hstart,rf.hend,rf.hid  
                     FROM   repeat_feature rf,static_golden_path sgp
                     WHERE  sgp.raw_id = rf.contig
                     AND    rf.contig in $idlist
                     AND    sgp.type = '$type'
		     AND    sgp.chr_name='$chr_name' 
                     ORDER  by start";
    
    my $sth = $self->dbobj->prepare($statement);
    $sth->execute();
    
    
    my ($fid,$start,$end,$strand,$score,$analysisid,$hstart,$hend,$hid);
    
    $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$score,\$analysisid,
		       \$hstart,\$hend,\$hid);
    
    
    my @features;
    my @distinct_features;  
    my %analhash;
    my $out;
    
    #Loop through repeat features, and glob them (i.e. join overlapping ones)
    while( $sth->fetch ) {
	
	if (($end > $length) || ($start < 1)) {
	    next;
	}
	
	#exclude contained features
	if(  defined $bp && defined $out && $end < $out->end ) { next; }
	
	#Glob overlapping and close features
	if ( defined $bp && defined $out && $out->end+$bp >= $start ) {
	    
	    # reset previous guys end to end
	    $out->end($end);
	    next;
	}

	my $analysis=$self->_get_analysis($analysisid);
	
	if($hid ne '__NONE__' ) {
	    # is a paired feature
	    # build EnsEMBL features and make the FeaturePair
    
	    $out = Bio::EnsEMBL::FeatureFactory->new_repeat();
	    $out->set_all_fields($start,$end,$strand,$score,'repeatmasker','repeat',$self->id,
				 $hstart,$hend,1,$score,'repeatmasker','repeat',$hid);
	    
	    $out->analysis($analysis);
	    
	} else {
	    #print(STDERR "Repeat feature does not have a hid. bad news....");
	}
	push(@features,$out);
    }
    return @features;
}



=head2 get_all_PredictionFeatures

 Title   : get_all_PredictionFeatures
 Usage   : foreach my $sf ( $contig->get_all_PredictionFeatures )
 Function: Gets all the repeat features on a contig.
 Example :
 Returns : 
 Args    : 


=cut

sub get_all_PredictionFeatures {
   my ($self) = @_;

   my @array;

#   my $id     = $self->internal_id();
   my $length = $self->length();

	
   if( exists $self->{'_genscan_cache'} ) {
       return @{$self->{'_genscan_cache'}};
   }


   my $glob_start=$self->_global_start;
   my $glob_end=$self->_global_end;
   my $chr_name=$self->_chr_name;
   my $idlist  = $self->_raw_contig_id_list();
   
   unless ($idlist){
       return ();
   }
   
   
   
   my $fsetid;
   my $previous;
   my $previous_contig;
   my %analhash;
   my $analysis_type='genscan';
   
    my $type = $self->dbobj->static_golden_path_type;
    
   my $query = "SELECT f.id, 
                        IF     (sgp.raw_ori=1,(f.seq_start+sgp.chr_start-sgp.raw_start-$glob_start),
                                 (sgp.chr_start+sgp.raw_end-f.seq_end-$glob_start)) as start,  
                        IF     (sgp.raw_ori=1,(f.seq_end+sgp.chr_start-sgp.raw_start-$glob_start),
                                 (sgp.chr_start+sgp.raw_end-f.seq_start-$glob_start)), 
                        IF     (sgp.raw_ori=1,f.strand,(-f.strand)),
                        f.score,f.evalue,f.perc_id,f.phase,f.end_phase,f.analysis,f.hid,f.contig 
                        FROM   feature f, analysis a,static_golden_path sgp 
                        WHERE    f.analysis = a.id 
                        AND    sgp.raw_id = f.contig
                        AND    f.contig in $idlist
		        AND    a.gff_source = '$analysis_type'  
                        AND    sgp.type = '$type'
		        AND    sgp.chr_name='$chr_name' 
                        ";
   
   my $sth = $self->dbobj->prepare($query);
   
   $sth->execute();
   
   my ($fid,$start,$end,$strand,$score,$evalue,$perc_id,$phase,$end_phase,$analysisid,$hid,$contig);
   
   # bind the columns
   $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$score,\$evalue,\$perc_id,\$phase,\$end_phase,\$analysisid,\$hid,\$contig);
   
   $previous = -1;
   my $current_fset;
   my $count;
   
   while( $sth->fetch ) {
       
       if (($end > $length) || ($start < 1)) {
	   next;
       }
       
       my $out;
       
       my $analysis;
       
       if (!$analhash{$analysisid}) {
	   
	   my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->dbobj);
	   $analysis = $feature_obj->get_Analysis($analysisid);
	   
	   $analhash{$analysisid} = $analysis;
	   
       } else {
	   $analysis = $analhash{$analysisid};
       }

       
       if( $hid eq "Initial Exon" || $hid eq "Single Exon" || $previous eq "Single Exon" || $previous eq "Terminal Exon" || $previous eq -1 || $previous_contig != $contig) {
	   $count++;
	   $current_fset = Bio::EnsEMBL::SeqFeature->new();
	   $current_fset->source_tag('genscan');
	   $current_fset->primary_tag('prediction');
	   $current_fset->analysis($analysis);
	   $current_fset->seqname($self->id);
	   $current_fset->id($count);
	   $current_fset->score(0.0);
	   
	   $current_fset->raw_seqname($self->id);
	   push(@array,$current_fset);
       }

       $out = Bio::EnsEMBL::SeqFeature->new;
       
       $out->seqname   ($self->id);
       $out->raw_seqname($self->id);

       $out->start     ($start);
       $out->end       ($end);
       $out->strand    ($strand);
       $out->score     ($score)     if (defined $score);
       $out->p_value   ($evalue)    if (defined $evalue);
       $out->percent_id($perc_id)   if (defined $perc_id); 
       $out->phase     ($phase)     if (defined $phase);    
       $out->end_phase ($end_phase) if (defined $end_phase);
        

	my $query="select fset from fset_feature where feature=$fid"; 
	my $sth = $self->dbobj->prepare($query);
   	$sth->execute();
	my $arr_ref=$sth->fetchrow_arrayref;

	my $fsetid=$arr_ref->[0];

      



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
       $previous_contig = $contig;
  }
 

   $self->{'_genscan_cache'} = \@array;

   return @array;
}


=head2 get_all_SimpleFeatures

    my @feat = $contig->get_all_SimpleFeatures

Returns a list of B<Bio::EnsEMBL::SeqFeature>
objects.  Features which overlap the ends of the
contig are truncated to the contig, not discarded
like

=head2 get_all_SimpleFeatures_by_feature_type

    my @cpg = $contig->get_all_SimpleFeatures_by_feature_type('cpg_island');

=cut

sub get_all_SimpleFeatures {
    my( $self ) = @_;
    
    return $self->_fetch_SimpleFeatures_SQL_clause('');
}

sub get_all_SimpleFeatures_by_feature_type {
    my( $self, $feature_type ) = @_;
    
    $self->throw("No feature type given") unless $feature_type;
    
    return $self->_fetch_SimpleFeatures_SQL_clause(qq{
        AND a.gff_feature = '$feature_type'
        });
}

sub get_all_SimpleFeatures_by_feature_type_above_score {
    my( $self, $feature_type, $score ) = @_;
    
    $self->throw("No feature type given") unless $feature_type;
    $self->throw("No score given") unless defined($score);
    
    return $self->_fetch_SimpleFeatures_SQL_clause(qq{
        AND a.type = '$feature_type'
        AND f.score >= $score
        });
}

sub get_all_SimpleFeatures_by_analysis_id {
    my( $self, $analysis_id ) = @_;
    
    $self->throw("No analysis ID given") unless $analysis_id;
    
    return $self->_fetch_SimpleFeatures_SQL_clause(qq{
        AND a.id = $analysis_id
        });
}

# Internal method used by the get_all_SimpleFeatures* methods
sub _fetch_SimpleFeatures_SQL_clause {
    my( $self, $sql_extra ) = @_;

    my $global_start    = $self->_global_start;
    my $chr_name        = $self->_chr_name;
    my $type            = $self->dbobj->static_golden_path_type;
    
    # Return if we don't have any raw contigs (all gap)
    my $idlist = $self->_raw_contig_id_list or return;

    # This is the generic SQL used by all the subroutines
    # which call this one.
    my $sql_begin = qq{
        SELECT f.id
          , IF (sgp.raw_ori = 1
              , (sgp.chr_start + f.seq_start - sgp.raw_start - $global_start)
              , (sgp.chr_start + sgp.raw_end - f.seq_end     - $global_start)) as start
          , IF (sgp.raw_ori = 1
              , (sgp.chr_start + f.seq_end   - sgp.raw_start - $global_start)
              , (sgp.chr_start + sgp.raw_end - f.seq_start   - $global_start))
          , sgp.raw_ori * f.strand
          , f.score
          , f.analysis
          , f.name
          , f.hid
        FROM feature f
          , analysis a
          , static_golden_path sgp
        WHERE f.analysis = a.id
          AND sgp.raw_id = f.contig
          AND f.contig IN $idlist
          AND sgp.type = '$type'
          AND sgp.chr_name = '$chr_name'
        };
    
    # All statements have this ORDER by clause on the end
    my $sql_order = qq{
        ORDER BY start
        };

    # Make the full statement and execute it
    my $sql = join(' ', $sql_begin, $sql_extra, $sql_order);
    my $sth = $self->dbobj->prepare($sql);    
    $sth->execute();     
    
    # Bind columns to variables for the fastest possible
    # retrieval of results.
    my ($fid,
        $start, $end, $strand,
        $f_score, $analysis_id, $name, $hid);
    $sth->bind_columns(undef,
        \$fid,
        \$start, \$end, \$strand,
        \$f_score, \$analysis_id, \$name, \$hid);
        
    my $length = $self->length;

    my( %anal, @features );
    while ($sth->fetch) {
        
        # Skip features outside our region
        if ($end < 1 or $start > $length) {
            next;
        }
        
        # Get an analysis object
        my( $analysis );
        unless ($analysis = $anal{$analysis_id}) {
            my $feature_obj = Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->dbobj);
	    $analysis = $feature_obj->get_Analysis($analysis_id);
	    $anal{$analysis_id} = $analysis;
        }
    
        # Truncate the feature coordinates to the region
        $start = 1       if $start < 1;
        $end   = $length if $end   > $length;
    
        # Make the feature
        my $feat = Bio::EnsEMBL::FeatureFactory->new_feature;
        $feat->id       ($fid);
        $feat->start    ($start);
        $feat->end      ($end);
        $feat->strand   ($strand);
        $feat->score    ($f_score);
        $feat->analysis ($analysis);
        
        push(@features, $feat);
    }
    return @features;
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
   my ($self,$glob) = @_;
   

   if (! defined $glob) {
       $glob=50;
   }

   if( $self->_external_feature_cache ) {
       return @{$self->{'_external_feature_cache_array'}};
   }
   
   &eprof_start("External-feature-get");
   
   my @web;
   my @das;
   my @std;

   my @features;
   my @contig_features;

   if( scalar($self->_vmap->get_all_RawContigs) == 0) {
       return();
   }
   ## Loop over the currently config'd EFFs and sort them onto lists that
   ## are 1. Lightweigth for the web, 2. normal or 3. those that don't know about
   ## Ensembl internal clone/contig IDs. 
   ## After sorting call each one to get a list of feature back in contig/clone coords.
   ## Note that they should always return lists (possible empty) or bad things happen.
      
   foreach my $extf ( $self->dbobj->_each_ExternalFeatureFactory ) {
       if( $extf->isa('Bio::EnsEMBL::DB::WebExternalFeatureFactoryI') ) {
	   		push(@web,$extf);
       } elsif( $extf->isa('Bio::EnsEMBL::ExternalData::DAS::DAS') ) {
	   		push(@das,$extf);
       } else {
	   		push(@std,$extf);
       }
   }

   if( scalar(@web) > 0 ) {
       my @clones;
       my %cloneh;
       my %featureh;
       foreach my $contig ( $self->_vmap->get_all_RawContigs) {
	   if( !defined $cloneh{$contig->cloneid} ) {
	       $cloneh{$contig->cloneid} = [];
	       $featureh{$contig->cloneid} = [];
	   }
	   my $string = $contig->cloneid.".".$contig->seq_version;
	  
	   push(@clones,$string);
	   push(@{$cloneh{$contig->cloneid}},$contig);
       }
       # get them out, push into array by clone
       foreach my $extf ( @web ) {
	   &eprof_start("external_get_web".$extf);
	   foreach my $feature ( $extf->get_Ensembl_SeqFeatures_clone_web($glob,@clones) ) {
	       my $clone = $feature->seqname;
	       $clone =~ s/\.\d+$//g;
	       push(@{$featureh{$clone}},$feature);
	   }
	   &eprof_end("external_get_web".$extf);
       }

       # loop over clone. Sort both feature and contig arrays.
       # then loop over features, changing identifiers and then push onto final array

       foreach my $clone ( keys %cloneh ) {
	   my @features = sort { $a->start <=> $b->start } @{$featureh{$clone}};
	   my @contigs  = sort { $a->embl_offset <=> $b->embl_offset } @{$cloneh{$clone}};
	   my $current_contig = shift @contigs;
          FEATURE :
	   foreach my $f ( @features ) {
	       while( $current_contig->length + $current_contig->embl_offset < $f->start ) {
		   $current_contig = shift @contigs;
                   if( !defined $current_contig ) { last FEATURE; }
	       }
	       if( $f->end < $current_contig->embl_offset ) {
		   next; # not on a contig on this golden path presumably
	       }
	       $f->start($f->start - $current_contig->embl_offset+1);
	       $f->end  ($f->end   - $current_contig->embl_offset+1);
	       $f->seqname($current_contig->id);
	       push(@contig_features,$f);
	   }
       }
   }
 
   if( scalar(@std) > 0 ) {
       foreach my $contig ( $self->_vmap->get_all_RawContigs) {       
	   foreach my $extf ( @std ) {
	       &eprof_start("external_get_std".$extf);

	       if( $extf->can('get_Ensembl_SeqFeatures_contig') ) {
		   foreach my $sf ($extf->get_Ensembl_SeqFeatures_contig($contig->internal_id,$contig->seq_version,1,$contig->length,$contig->id)) {
			$sf->seqname($contig->id);
			push(@contig_features,$sf);
		   }
	       }
	       if( $extf->can('get_Ensembl_SeqFeatures_clone') ) {
       
		   foreach my $sf ( $extf->get_Ensembl_SeqFeatures_clone($contig->cloneid,$contig->seq_version,$contig->embl_offset,$contig->embl_offset+$contig->length(),$contig->id) ) {
		       
		       my $start = $sf->start - $contig->embl_offset+1;
		       my $end   = $sf->end   - $contig->embl_offset+1;
		       $sf->start($start);
		       $sf->end($end);
		       $sf->seqname($contig->id);
		       push(@contig_features,$sf);
		   }
	       }

	       &eprof_end("external_get_std".$extf);
	   }
       }
   }
	

	## The DAS external feature factory is based on coordinates on contigs (at the moment)
	## The standara EFF system has been moved to use contig/clone internal IDs so we have to
	## make a special case for DAS (and possibly other) EFFs that know nothing about Ensembl
	## internal IDs	. There are probably more efficient ways to do this....what about DAS caching?
	if( scalar(@das) > 0 ) {
       foreach my $contig ( $self->_vmap->get_all_RawContigs) {       
	   foreach my $extf ( @das ) {
	       &eprof_start("external_get_std".$extf);

	       if( $extf->can('get_Ensembl_SeqFeatures_contig') ) {
		   foreach my $sf ($extf->get_Ensembl_SeqFeatures_contig($contig->id,$contig->seq_version,1,$contig->length,$contig->id)) {
			$sf->seqname($contig->id);
			push(@contig_features,$sf);
		   }
	       }
	       if( $extf->can('get_Ensembl_SeqFeatures_clone') ) {
       
		   foreach my $sf (
		  
$extf->get_Ensembl_SeqFeatures_clone($contig->cloneid,$contig->seq_version,$contig->embl_offset,$contig->embl_offset+$contig->length(),$contig->cloneid) ) {
		       
		       my $start = $sf->start - $contig->embl_offset+1;
		       my $end   = $sf->end   - $contig->embl_offset+1;
		       $sf->start($start);
		       $sf->end($end);
		       $sf->seqname($contig->id);
		       push(@contig_features,$sf);
		   }
	       }

	       &eprof_end("external_get_std".$extf);
	   }
       }
   }	    
   &eprof_end("External-feature-get");

   # ok. Now @contig_features are in contig coordinates. Map up.
   
   # this is the simplest way. We can do this more efficiently if need be

   
   &eprof_start("External-coordinate-lift");

   my @final;
   foreach my $f ( @contig_features ) {
       if( defined $self->_convert_seqfeature_to_vc_coords($f) ) {
	   push(@final,$f);
       }
   }


   &eprof_end("External-coordinate-lift");

   $self->{'_external_feature_cache_array'} = \@final;
   $self->_external_feature_cache(1);

   return @final;

}

=head2 _external_feature_cache

 Title   : _external_feature_cache
 Usage   : $obj->_external_feature_cache($newval)
 Function: 
 Returns : value of _external_feature_cache
 Args    : newvalue (optional)


=cut

sub _external_feature_cache{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_external_feature_cache'} = $value;
    }
    return $obj->{'_external_feature_cache'};

}


sub get_all_FPCClones {
    my $self = shift;

    my $glob_start = $self->_global_start;
    my $glob_end   = $self->_global_end;
    my $length     = $self->length;
    my $chr_name   = $self->_chr_name;

    
    my $mapdb;
    
    eval {
	$mapdb = $self->dbobj->mapdb();
    };
    if( $@ || !defined $mapdb ) {
	$self->warn("in get_all_FPCClones, unable to get mapdb. Returning empty list [$@]");
	return ();
    }

    my $fpcmap = $mapdb->get_Map('FPC');
    $chr_name =~ s/chr//g;

    my $chr = $fpcmap->get_ChromosomeMap($chr_name);

    my @fpcclones = $chr->get_Clones_by_start_end($glob_start,$glob_end);

    foreach my $clone ( @fpcclones ) {
        my $newstart = $clone->start - $glob_start;
       	$clone->start($newstart);
    }

    return @fpcclones;
}

   

sub get_landmark_MarkerFeatures_old {

my ($self) = @_;

my $glob_start = $self->_global_start;
my $glob_end   = $self->_global_end;
my $length     = $self->length;
my $chr_name   = $self->_chr_name;
my $dbname     = $self->dbobj->dbname;
my $mapsdbname = $self->dbobj->mapdbname;
my @markers;

my $idlist  = $self->_raw_contig_id_list();
unless ($idlist){
    return ();
}

my $glob = 100;

eval {
    require Bio::EnsEMBL::Map::MarkerFeature;

    my $type = $self->dbobj->static_golden_path_type;
    my $statement= "   SELECT 
                       IF     (sgp.raw_ori=1,(f.seq_start+sgp.chr_start-sgp.raw_start-$glob_start),
                              (sgp.chr_start+sgp.raw_end-f.seq_end-$glob_start)) as start,                                        
                       IF     (sgp.raw_ori=1,(f.seq_end+sgp.chr_start-sgp.raw_start-$glob_start),
                              (sgp.chr_start+sgp.raw_end-f.seq_start-$glob_start)), 
                              f.score, 
                       IF     (sgp.raw_ori=1,f.strand,(-f.strand)), 
                              f.name, f.hstart, f.hend, 
                              f.hid, f.analysis, c.name 
                       FROM   feature f,
                              contig_landmarkMarker c,
                              static_golden_path sgp 
                       WHERE    f.hid=c.marker  
                       AND    f.contig in $idlist 
                       AND    sgp.raw_id=f.contig 
                       AND    sgp.type = '$type'
                       AND    sgp.chr_name='$chr_name' 
                       GROUP BY f.hid ORDER BY start ";
    

	my $sth = $self->dbobj->prepare($statement);
	$sth->execute;
	
	my ($start, $end, $score, $strand, $hstart, 
	    $name, $hend, $hid, $analysisid,$synonym);
	
	my $analysis;
	my %analhash;
	
	$sth->bind_columns
	    ( undef, \$start, \$end, \$score, \$strand, \$name, 
	      \$hstart, \$hend, \$hid, \$analysisid,\$synonym);
	
    my $prev;
      while( $sth->fetch ) {

	    #clipping
	    if (($end > $length) || ($start < 1)) {
		next;
	    }
	    
	    if( defined $prev && $prev->end+$glob > $start && $synonym eq $prev->id) {
		next;
	    }

	    my @args=($start,$end,$score,$strand,$name,$hstart,$hend,$hid,
		      $analysisid,$synonym);
	    
	    
	    my $out=$self->_create_Marker_features(@args);
	    if (defined $out){
		push (@markers,$out);
		$prev = $out;
	    } 
	}
};


if($@){$self->warn("Problems retrieving map data\nMost likely not connected to maps db\n$@\n");}

return @markers;

}





sub get_landmark_MarkerFeatures{
    my ($self,$glob) = @_;

    my $chr_name   = $self->_chr_name;

    if( !defined $glob ) {
        $glob = 500000;
    }


    my $glob_start = $self->_global_start;
    my $glob_end   = $self->_global_end;
    my $length     = $self->length;


   my $statement= " SELECT  start,
			    end,
			    strand,
			    name 
		    FROM    contig_landmarkMarker 
		    WHERE   chr_name = '$chr_name'
                    AND     start >= $glob_start 
                    AND     end <= $glob_end 
		    ORDER BY start
		";
   
   $statement =~ s/\s+/ /g;
   
   my $sth = $self->dbobj->prepare($statement);
   $sth->execute;
   
   my ($start, $end, $strand, $name);
   
   my $analysis;
   my %analhash;
   
   $sth->bind_columns
       ( undef, \$start, \$end,  \$strand, \$name);
   
   my @out;
   my $prev;
   while( $sth->fetch ) {
	
       #############################
       # change to local coordinates
       #############################
       $start=$start-$glob_start;
       $end=$end-$glob_start;
       
       if( defined $prev && $prev->end + $glob > $start  && $prev->id eq $name ) {           
           next;
       }

       my $sf = Bio::EnsEMBL::SeqFeature->new();
       $sf->start($start);
       $sf->end($end);
       $sf->strand($strand);
       $sf->id($name);
       push(@out,$sf);
       $prev = $sf;
   } 

   return @out;
}





=head2 next_landmark_Marker

 Title   : next_landmark_Marker
 Usage   : $obj->next_landmark_Marker
 Function: retrieves next marker  
 Returns : marker feature
 Args    : golden path position, chromosome, Mb limit (optional)


=cut



sub next_landmark_Marker {

my ($self,$start,$chr_name,$Mb)=@_;

$self->throw("Must supply golden path position") unless $start;
$self->throw("Must supply chromosome") unless $chr_name;
if (!defined $Mb){$Mb=1000000;}
my $length = $self->length;
my $glob_start=$self->_global_start;
my $glob_end=$self->_global_end;
   $chr_name=$self->_chr_name;
my $dbname=$self->dbobj->dbname;
my $mapsdbname=$self->dbobj->mapdbname;
my $type = $self->dbobj->static_golden_path_type;

my @markers;

eval {
    require Bio::EnsEMBL::Map::MarkerFeature;

    my $end;
    my $limit;
    unless ($#markers>=0 || $end >255000000){

	$limit=$limit+$Mb;
	$end=$start+$limit;

	my $statement=   "SELECT    
                          IF        (sgp.raw_ori=1,(f.seq_start+sgp.chr_start-sgp.raw_start-$glob_start),
                                    (sgp.chr_start+sgp.raw_end-f.seq_end-$glob_start)) as start,                                        
                          IF        (sgp.raw_ori=1,(f.seq_end+sgp.chr_start-sgp.raw_start-$glob_start),
                                    (sgp.chr_start+sgp.raw_end-f.seq_start-$glob_start)),                                       
                                    f.score, 
                          IF        (sgp.raw_ori=1,f.strand,(-f.strand)),
                                    f.name, f.hstart, f.hend, 
                                    f.hid, f.analysis, s.name  
                          FROM      $dbname.feature f, $dbname.analysis a, 
                                    $mapsdbname.MarkerSynonym s, 
                                    $dbname.static_golden_path sgp
                          WHERE     sgp.raw_id=f.contig  
                          AND       f.hid=s.marker
                          AND       sgp.chr_name = '$chr_name' 
                          AND       sgp.type = '$type'
                          AND       f.analysis = a.id 
                          AND       a.db='mapprimer'
                          AND       sgp.chr_start > $start 
                          AND       sgp.chr_start < $end
                          AND       sgp.chr_start + f.seq_start - sgp.raw_start > $start  
                          AND       (s.name regexp '^D[0-9,X,Y][0-9]?S' OR s.name regexp '^AFM') 
                          ORDER BY  start limit 1";



	my $sth = $self->dbobj->prepare($statement);
	$sth->execute;

	my ($score, $strand, $hstart, $name, $hend, $hid, $analysisid,$synonym);

	my $analysis;
	my %analhash;

	$sth->bind_columns
	    ( undef, \$start, \$end, \$score, \$strand, \$name, 
	      \$hstart, \$hend, \$hid, \$analysisid,\$synonym);


	while( $sth->fetch ) {

	    my @args=($start,$end,$score,$strand,$name,$hstart,$hend,$hid,
		      $analysisid,$synonym);

	    my $out=$self->_create_Marker_features(@args);
	    if (defined $out){
		push (@markers,$out);
	    } 
	}
    }
};

if($@){$self->warn("Problems retrieving map data\nMost likely not connected to maps db\n$@\n");}

return $markers[0];
}


=head2 _previous_landmark_Marker

 Title   : previous_landmark_Marker
 Usage   : $obj->previous_landmark_Marker
 Function: retrieves previous marker  
 Returns : marker feature
 Args    : golden path position, chromosome, Mb limit (optional) 


=cut



sub previous_landmark_Marker
{

my ($self,$start,$chr_name,$Mb)=@_;

$self->throw("Must supply golden path position") unless $start;
$self->throw("Must supply chromosome") unless $chr_name;
if (!defined $Mb){$Mb=1000000;}


my $glob_start=$self->_global_start;
my $glob_end=$self->_global_end;
   $chr_name=$self->_chr_name;
my $dbname=$self->dbobj->dbname;
my $mapsdbname=$self->dbobj->mapdbname;
my $type = $self->dbobj->static_golden_path_type;

my @markers;


eval {
    require Bio::EnsEMBL::Map::MarkerFeature;

    my $end;
    my $limit;
    unless ($#markers>=0 || $end==1){
   
	$limit=$limit+$Mb; 
	$end=$start-$limit;
	if ($end<0){$end=1;}


	my $statement=   "SELECT    
                          IF        (sgp.raw_ori=1,(f.seq_start+sgp.chr_start-sgp.raw_start-$glob_start),
                                    (sgp.chr_start+sgp.raw_end-f.seq_end-$glob_start)) as start,                                        
                          IF        (sgp.raw_ori=1,(f.seq_end+sgp.chr_start-sgp.raw_start-$glob_start),
                                    (sgp.chr_start+sgp.raw_end-f.seq_start-$glob_start)),                                       
                                    f.score, 
                          IF        (sgp.raw_ori=1,f.strand,(-f.strand)), 
                                    f.name, f.hstart, f.hend, 
                                    f.hid, f.analysis, s.name  
                          FROM      $dbname.feature f, $dbname.analysis a, 
                                    $mapsdbname.MarkerSynonym s, 
                                    $dbname.static_golden_path sgp
                          WHERE     sgp.raw_id=f.contig  
                          AND       f.hid=s.marker
                          AND       sgp.chr_name='$chr_name'
                          AND       sgp.type = '$type'
                          AND       f.analysis = a.id 
                          AND       a.db='mapprimer'                       
                          AND       sgp.chr_start<$start 
                          AND       sgp.chr_start>=$end 
                          AND       sgp.chr_start+f.seq_start-sgp.raw_start<$start  
                          AND       (s.name regexp '^D[0-9,X,Y][0-9]?S' OR s.name regexp '^AFM')  
                          ORDER BY  start limit 1";



	my $sth = $self->dbobj->prepare($statement);
	$sth->execute;
	
	my ($score, $strand, $hstart, $name, $hend, $hid, $analysisid,$synonym);
	
	my $analysis;
	my %analhash;
	
	$sth->bind_columns
	    ( undef, \$start, \$end, \$score, \$strand, \$name, 
	      \$hstart, \$hend, \$hid, \$analysisid,\$synonym);
	
        
	while( $sth->fetch ) {
	    
	    my @args=($start,$end,$score,$strand,$name,$hstart,$hend,$hid,
		      $analysisid,$synonym);
	    
	    my $out=$self->_create_Marker_features(@args);
	    if (defined $out){
		push (@markers,$out);
	    } 
	}
    }
};

if($@){$self->warn("Problems retrieving map data\nMost likely not connected to maps db\n$@\n");}

return $markers[0];


}


=head2 get_all_Genes_exononly

 Title   : get_all_Genes_exononly
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Genes_exononly{
   my ($self) = @_;

   
   my $glob_start=$self->_global_start;
   my $glob_end=$self->_global_end;
   my $chr_name=$self->_chr_name;
   my $idlist  = $self->_raw_contig_id_list();
   my $type = $self->dbobj->static_golden_path_type;
   
   unless ($idlist){
       return ();
   }

   if( exists $self->{'_all_Genes_exononly'} ) {
       return @{$self->{'_all_Genes_exononly'}};
   }

   my $query = "
        SELECT e.id,e.sticky_rank,et.rank,et.transcript,t.gene, 
        IF     (sgp.raw_ori=1,(e.seq_start+sgp.chr_start-sgp.raw_start-$glob_start),
                 (sgp.chr_start+sgp.raw_end-e.seq_end-$glob_start)) as start,  
        IF     (sgp.raw_ori=1,(e.seq_end+sgp.chr_start-sgp.raw_start-$glob_start),
                 (sgp.chr_start+sgp.raw_end-e.seq_start-$glob_start)), 
        IF     (sgp.raw_ori=1,e.strand,(-e.strand))
        FROM   exon e,static_golden_path sgp,exon_transcript et,transcript t
        WHERE  t.id = et.transcript
        AND    et.exon = e.id 
        AND    sgp.raw_id = e.contig
	AND    sgp.chr_name = '$chr_name'
        AND    sgp.type = '$type'
        AND    e.contig in $idlist
        AND    sgp.chr_end >= $glob_start
	AND    sgp.chr_start <= $glob_end
        ORDER  BY t.gene,t.id,et.rank,e.sticky_rank";

   my $sth = $self->dbobj->prepare($query);
   $sth->execute();

   my ($exonid,$stickyrank,$rank,$transcriptid,$geneid,$start,$end,$strand);
   $sth->bind_columns(undef,\$exonid,\$stickyrank,\$rank,\$transcriptid,\$geneid,\$start,\$end,\$strand);
  
   my $current_transcript;
   my $current_gene;
   my $current_transcript_id='';
   my $current_gene_id='';
   my $previous_exon;

   my @out;
   my @trans;
   my $length = $glob_end - $glob_start;

   while( $sth->fetch ) {

       if (($end > $length) || ($start < 1)) {
	   next;
       }
       if( $geneid ne $current_gene_id ) {
	   # make a new gene
	   $current_gene = Bio::EnsEMBL::Gene->new;
	   $current_gene->id($geneid);
	   push(@out,$current_gene);
	   $current_gene_id = $geneid;
       }

       if( $transcriptid ne $current_transcript_id ) {
	   # make a new transcript
	   $current_transcript = Bio::EnsEMBL::WebTranscript->new();
	   $current_gene->add_Transcript($current_transcript);
	   push(@trans,$current_transcript);
	   if( $rank == 1 ) {
	       $current_transcript->is_start_exon_in_context('dummy',1);
	   } else {
	       $current_transcript->is_start_exon_in_context('dummy',0);
	   }

	   $current_transcript_id = $transcriptid;
	   $current_transcript->id($transcriptid);
       }

       if( $stickyrank > 1 ) {
	   if( !defined $previous_exon ) {
	       $self->warn("Really bad news - half-on-half off Sticky Exon. Faking it");
	   }
	   if( $previous_exon->end < $end ) {
	       $previous_exon->end($end);
	       next;
	   }

       }


       my $exon = Bio::EnsEMBL::Exon->new();
       $exon->start($start);
       $exon->end($end);
       $exon->strand($strand);
       $exon->id($exonid);
       $exon->seqname($self->id);
       $previous_exon = $exon;
       $current_transcript->add_Exon($exon);
       $current_transcript->end_exon_rank($rank);
       
   }

   #
   # We need to make another quick trip to the database for each
   # transcript to discover whether we have all of it or not
   #

   foreach my $trans ( @trans ) {
       my $sth2 = $self->dbobj->prepare("select max(rank) from exon_transcript where transcript = '".$trans->id."'");
       $sth2->execute;
       my ($rank) = $sth2->fetchrow_array();
       if( $rank == $trans->end_exon_rank) {
	   $trans->is_end_exon_in_context('dummy',1);
       } else {
	   $trans->is_end_exon_in_context('dummy',0);
       }
                                                              
   }


   #
   # This can obviously be optimised to a better single trip
   # to the database
   #

   my $gene_obj = $self->dbobj->gene_Obj;

   foreach my $g ( @out ) {
       $gene_obj->_get_dblinks($g);
       $gene_obj->_get_description($g);
   }

   $self->{'_all_Genes_exononly'} = \@out;

   return @out;

}





=head2 get_all_VirtualGenes_startend

 Title   : get_all_VirtualGenes_startend
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut




sub get_all_VirtualGenes_startend
{

    my $self = shift;
    
    my $gene;
    my @genes;
    
    if( $self->_cached_virtualgenes_startend == 1 ) {
	return @{$self->{'_virtualgenes_startend'}};
    }

    my $glob_start  = $self->_global_start;
    my $glob_end    = $self->_global_end;
    my $chr_name    = $self->_chr_name;
    my $idlist      = $self->_raw_contig_id_list();
    my $type = $self->dbobj->static_golden_path_type;

    unless ($idlist){
	return ();
    }
    
    $self->throw ("I need a chromosome name") unless defined $chr_name;
    $self->throw ("I need a chromosome end") unless defined $glob_end;
    $self->throw ("I need a chromosome start") unless defined $glob_start;
    
    &eprof_start("virtualgene-sql-get");

    my $query ="SELECT     STRAIGHT_JOIN t.gene,
                       MIN(IF(sgp.raw_ori=1,(e.seq_start+sgp.chr_start-sgp.raw_start-$glob_start),
                                  (sgp.chr_start+sgp.raw_end-e.seq_end-$glob_start))) as start,
                       MAX(IF(sgp.raw_ori=1,(e.seq_start+sgp.chr_start-sgp.raw_start-$glob_start),
                                  (sgp.chr_start+sgp.raw_end-e.seq_start-$glob_start))) as end 
            FROM       static_golden_path sgp ,exon e,exon_transcript et,transcript t 
            WHERE      sgp.raw_id=e.contig
            AND        e.contig in $idlist 
            AND        e.id=et.exon 
            AND        t.id=et.transcript 
            AND        sgp.chr_end >= $glob_start   
            AND        sgp.chr_start <=$glob_end 
            AND        sgp.chr_name='$chr_name'
            AND        sgp.type = '$type'
            GROUP BY   t.gene;";
    
    my $sth = $self->dbobj->prepare($query);
    $sth->execute;

    &eprof_end("virtualgene-sql-get");


    &eprof_start("virtualgene-build");

				# 
    my ($gene_id,$start,$end);	# 
    $sth->bind_columns(undef,\$gene_id,\$start,\$end);

    while ($sth->fetch){
	if( $end < 1 ) { 
	    # clip this gene to the left
	    next;
	}

        if( $start > $self->length ) {
	    # clip this gene to the right
	    next;
	}

 
	if (($end > $self->length)) {$end=$self->length;}
	if (($start < 1)) {$start=1;}

	$gene=Bio::EnsEMBL::Gene->new();
	$gene->id($gene_id);

	&eprof_start("virtualgene-externaldb");

	
#Get the DBlinks for the given gene
    my $entryAdaptor = $self->dbobj->get_DBEntryAdaptor();
    my @gene_xrefs = $entryAdaptor->fetch_by_gene($gene_id);

    foreach my $genelink (@gene_xrefs) {
        $gene->add_DBLink($genelink);
    }

	my $query1 = "select t.id from transcript t where t.gene = '$gene_id';";
	my $sth1 = $self->dbobj->prepare($query1);
	$sth1->execute;
	
	while (my $transid = $sth1->fetchrow) {
	    
	    $transid =~ s/T/P/;
	    
	    my @transcript_xrefs = $entryAdaptor->fetch_by_translation($transid);
	    
	    foreach my $translink(@transcript_xrefs) {
		
		$gene->add_DBLink($translink);
	    }
	}
#End for fetching the DBlinks
	

	&eprof_end("virtualgene-externaldb");

	my $genestr=1;
	my $vg = Bio::EnsEMBL::VirtualGene->new(-gene => $gene,
						-contig => $self, 
						-start => $start, 
						-end => $end, 
						-strand => $genestr
						);

	push @genes,$vg;
    }


    &eprof_end("virtualgene-build");

    $self->{'_virtualgenes_startend'} = \@genes;
    $self->_cached_virtualgenes_startend(1);

    return @genes;

}

=head2 _cached_virtualgenes_startend

 Title   : _cached_virtualgenes_startend
 Usage   : $obj->_cached_virtualgenes_startend($newval)
 Function: 
 Example : 
 Returns : value of _cached_virtualgenes_startend
 Args    : newvalue (optional)


=cut

sub _cached_virtualgenes_startend{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_cached_virtualgenes_startend'} = $value;
    }
    return $obj->{'_cached_virtualgenes_startend'};

}






=head2 get_all_VirtualTranscripts_startend

 Title   : get_all_VirtualTranscripts_startend
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut




sub get_all_VirtualTranscripts_startend {
    my ($self) = @_;

    my $transcript;
    my @transcripts;

    my $glob_start=$self->_global_start;
    my $glob_end=$self->_global_end;
    my $chr_name=$self->_chr_name;
    my $idlist  = $self->_raw_contig_id_list();
    my $type = $self->dbobj->static_golden_path_type;

    unless ($idlist){
        return ();
    }

    $self->throw ("I need a chromosome name") unless defined $chr_name;
    $self->throw ("I need a chromosome end") unless defined $glob_end;
    $self->throw ("I need a chromosome start") unless defined $glob_start;

    my $query ="SELECT     STRAIGHT_JOIN t.id,
                           MIN(IF(sgp.raw_ori=1,(e.seq_start+sgp.chr_start-sgp.raw_start-$glob_start),
                                      (sgp.chr_start+sgp.raw_end-e.seq_end-$glob_start))) as start,
                           MAX(IF(sgp.raw_ori=1,(e.seq_end+sgp.chr_start-sgp.raw_start-$glob_start),
                                      (sgp.chr_start+sgp.raw_end-e.seq_start-$glob_start))) as end 
                FROM       static_golden_path sgp ,exon e,exon_transcript et,transcript t 
                WHERE      sgp.raw_id=e.contig
                AND        e.contig in $idlist 
                AND        e.id=et.exon 
                AND        t.id=et.transcript 
                AND        sgp.chr_end >= $glob_start   
                AND        sgp.chr_start <=$glob_end 
                AND        sgp.chr_name='$chr_name'
                AND        sgp.type = '$type'
                GROUP BY   t.id;";


    my $sth = $self->dbobj->prepare($query);
    $sth->execute;
    
    my ($transcript_id,$start,$end);
    $sth->bind_columns(undef,\$transcript_id,\$start,\$end);

    while ($sth->fetch){

	if (($end > $self->length)) {$end=$self->length;}
	if (($start < 1)) {$start=1;}

	my $gene=Bio::EnsEMBL::Gene->new();
	$gene->id($transcript_id);

#	my $query = "select external_db,external_id from genedblink where gene_id = '$transcript_id'";
#	my $sth = $self->dbobj->prepare($query);
#	my $res = $sth ->execute();
#	while( (my $hash = $sth->fetchrow_hashref()) ) {
#	    my $dblink = Bio::Annotation::DBLink->new();
#	    $dblink->database($hash->{'external_db'});
#	    $dblink->primary_id($hash->{'external_id'});
#	    $transcript->add_DBLink($dblink);
#	}

	my $transcriptstr=1;
# I don't think we need virtual transcript object :-)
	my $tr = Bio::EnsEMBL::VirtualGene->new(-gene => $gene,
						-contig => $self, 
						-start => $start, 
						-end => $end, 
						-strand => $transcriptstr
						);
	push @transcripts,$tr;
    }

    return @transcripts;

}

=head2 get_Genes

 Title   : get_Genes
 Usage   : @genes = $vc->get_Genes(@geneids)
 Function: gets only the named genes in @geneids
           Used in geneview to speed up getting
           a gene in VC coords
 Example :
 Returns : 
 Args    :


=cut

sub get_Genes {
  my ($self,@gene_ids) = @_;

  my $gene_obj = $self->dbobj->gene_Obj();

  my @newgenes;
  my @genes = $gene_obj->get_array_supporting('without',@gene_ids);
  
  my %gene;
  
  foreach my $gene ( @genes ) {
    $gene{$gene->id()}= $gene;
  }
  
   # this delegates off to Virtual::Contig
   @newgenes=$self->_gene_query(%gene);

  $self->{'_static_vc_gene_get'} = \@newgenes;

  return @newgenes;
}
  

=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   : @genes = $vc->get_all_Genes()
 Function: accelerated get_all_Genes for statics
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Genes {
    my ($self) = @_;

    if( defined $self->{'_static_vc_gene_get'} ) {
        return @{$self->{'_static_vc_gene_get'}};
    }

    &eprof_start("total-static-gene-get");
    my $idlist  = $self->_raw_contig_id_list();
    
    if( $idlist !~ /\w/ ) { 
        return ();
    }

    my $query = "
        SELECT DISTINCT(t.gene)
        FROM exon e
          , exon_transcript et
          , transcript t
        WHERE e.contig IN $idlist
          AND e.id = et.exon
          AND et.transcript = t.id
        ";

   &eprof_start("gene-sql-get");

   my $sth = $self->dbobj->prepare($query);
   $sth->execute;
      
   &eprof_end("gene-sql-get");

   my ($gene_id,$start,$end);	# 
   $sth->bind_columns(undef,\$gene_id);

   my @gene_ids;

   while ($sth->fetch){
       push(@gene_ids,$gene_id);
   }

   if( scalar(@gene_ids) == 0 ) {
       &eprof_end("total-static-gene-get");
       return ();
   }

   &eprof_start("full-gene-get");

   my $gene_obj = $self->dbobj->gene_Obj();

   my @genes = $gene_obj->get_array_supporting('without',@gene_ids);

   my %gene;

   foreach my $gene ( @genes ) {
       $gene{$gene->id()}= $gene;
   }

   &eprof_end("full-gene-get");

   &eprof_start("gene-convert");
   
   # this delegates off to Virtual::Contig
   my @newgenes=$self->_gene_query(%gene);

   &eprof_end("gene-convert");

   $self->{'_static_vc_gene_get'} = \@newgenes;

   &eprof_end("total-static-gene-get");

   return @newgenes;
}

=head2 fetch_chromosome_length

 Title   : fetch_chromosome_length
 Usage   : $length_in_bp = $self->fetch_chromosome_length
 Function:
 Example :
 Returns : SV 
 Args    :


=cut

sub fetch_chromosome_length {
   my ($self,$chr) = @_;

   unless (defined $chr){
      $chr = $self->_chr_name();
   }
   my $kba = $self->dbobj->get_KaryotypeBandAdaptor();
   my $len = $kba->fetch_chromosome_length($chr);
   return($len);
}


=head2 fetch_karyotype_adaptor

 Title   : fetch_karyotype_adaptor
 Usage   : $band_obj = $self->fetch_karyotype_adaptor
 Function:
 Example :
 Returns : karyotype band adaptor 
 Args    :


=cut

sub fetch_karyotype_adaptor {
   my ($self,@args) = @_;

   return ($self->dbobj->get_KaryotypeBandAdaptor());

}


=head2 fetch_karyotype_band_byname

 Title   : fetch_karyotype_band_byname
 Usage   : $band_obj = $self->fetch_karyotype_band_byname
 Function:
 Example :
 Returns : karyotype band object 
 Args    :


=cut

sub fetch_karyotype_band_by_name {
   my ($self,$chr, $band) = @_;

   my $kadp = $self->dbobj->get_KaryotypeBandAdaptor();
   my $kband = $kadp->fetch_by_chromosome_name($chr, $band);

   return $kband; 
}



=head2 fetch_karyotype_band_startend

 Title   : fetch_karyotype_band_startend
 Usage   : $band_obj = $self->fetch_karyotype_band_startend
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_karyotype_band_start_end {
   my ($self,@args) = @_;

   my $kadp = $self->dbobj->get_KaryotypeBandAdaptor();
   my @bands = $kadp->fetch_by_chromosome_start_end($self->_chr_name,$self->_global_start,$self->_global_end);

   return @bands; 
}

=head2 fetch_karyotype_band

 Title   : fetch_karyotype_band
 Usage   : $band_obj = $self->fetch_karyotype_band
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_karyotype_band {
   my ($self,@args) = @_;

   my $kadp = $self->dbobj->get_KaryotypeBandAdaptor();
   my $band = $kadp->fetch_by_chromosome_position($self->_chr_name,$self->_global_start + int($self->length/2));

   return $band 
}
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




sub _get_analysis {
    my ($self,$analysisid)=@_;
    
    my $analysis;
    my $analhash=$self->{_anal_hash};
    if (!$analhash->{$analysisid}) {
	my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($self->dbobj);
	$analysis = $feature_obj->get_Analysis($analysisid);
	
	$self->{_anal_hash}->{$analysisid} = $analysis;
	
    } else {
	$analysis = $analhash->{$analysisid};
    }
    return $analysis;
}


# sub _gp_position has been removed since revision 1.18

sub _create_Marker_features {
    my ($self,@args)=@_;
    
    my ($start,$end,$score,$strand,$name,$hstart,$hend,$hid,$analysisid,$synonym)=@args;
    
    
    my ( $out, $seqf1, $seqf2 );
    my $analysis = $self->_get_analysis($analysisid);
    
    $seqf1 = Bio::EnsEMBL::SeqFeature->new();
    $seqf2 = Bio::EnsEMBL::SeqFeature->new();
    $out = Bio::EnsEMBL::Map::MarkerFeature->new
	( -feature1 => $seqf1, -feature2 => $seqf2 );
    
    $out->set_all_fields
	( $start,$end,$strand,$score,
	  $name,'similarity',$self->id,
	  $hstart,$hend,1,$score,$name,'similarity',$name);
    
    $out->analysis($analysis);
    $out->mapdb( $self->dbobj->mapdb );
    $out->id ($synonym);

    return $out;
}

=head2 _global_start

 Title   : _global_start
 Usage   : $obj->_global_start($newval)
 Function: 
 Returns : value of _global_start
 Args    : newvalue (optional)


=cut

sub _global_start{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_global_start'} = $value;
    }
    return $obj->{'_global_start'};

}


=head2 _global_end

 Title   : _global_end
 Usage   : $obj->_global_end($newval)
 Function: 
 Returns : value of _global_end
 Args    : newvalue (optional)


=cut

sub _global_end{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_global_end'} = $value;
    }
    return $obj->{'_global_end'};

}

=head2 _chr_name

 Title   : _chr_name
 Usage   : $obj->_chr_name($newval)
 Function: 
 Returns : value of _chr_name
 Args    : newvalue (optional)


=cut

sub _chr_name{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_chr_name'} = $value;
    }
    return $obj->{'_chr_name'};

}



=head2 EMBL Dumping support

These functions are just to support EMBL dumping

=cut

=head2 desc

 Title   : desc
 Usage   : $obj->desc($newval)
 Function: 
 Example : 
 Returns : value of desc
 Args    : newvalue (optional)


=cut

sub desc{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'desc'} = $value;
    }
    return $obj->{'desc'};

}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Example : 
 Returns : value of id
 Args    : newvalue (optional)

=cut

sub id{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'id'} = $value;
    }
    return $obj->{'id'};

}

=head2 htg_phase

 Title   : htg_phase
 Usage   : $obj->htg_phase($newval)
 Function: 
 Example : 
 Returns : value of htg_phase
 Args    : newvalue (optional)


=cut

sub htg_phase{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'htg_phase'} = $value;
    }
    return $obj->{'htg_phase'};

}

=head2 sv

 Title   : sv
 Usage   : $obj->sv($newval)
 Function: 
 Example : 
 Returns : value of sv
 Args    : newvalue (optional)


=cut

sub sv{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'sv'} = $value;
    }
    return $obj->{'sv'};

}


=head2 embl_id

 Title   : embl_id
 Usage   : $obj->embl_id($newval)
 Function: 
 Example : 
 Returns : value of embl_id
 Args    : newvalue (optional)


=cut

sub embl_id{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'embl_id'} = $value;
    }
    return $obj->{'embl_id'};

}

=head2 project_name

 Title   : project_name
 Usage   : $obj->project_name($newval)
 Function: 
 Example : 
 Returns : value of project_name
 Args    : newvalue (optional)


=cut

sub project_name{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'project_name'} = $value;
    }
    return $obj->{'project_name'};

}


=head2 add_SeqFeature

 Title   : add_SeqFeature
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_SeqFeature{
   my ($self,$sf) = @_;

   push(@{$self->{'additional_seqf'}},$sf);
}


=head2 top_SeqFeatures

 Title   : top_SeqFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub top_SeqFeatures{
   my ($self) = @_;
   my @sf;

   @sf = $self->SUPER::top_SeqFeatures();
   push(@sf,@{$self->{'additional_seqf'}});
   return @sf;
}


=head2 _raw_contig_id_list

 Title   : _raw_contig_id_list
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _raw_contig_id_list {
   my ($self,@args) = @_;
    
   my $string;

   if( defined $self->{'_raw_contig_id_list'} ) {
       return $self->{'_raw_contig_id_list'};
   }

   foreach my $c ( $self->_vmap->get_all_RawContigs) {
       $string .= $c->internal_id . ",";
   }

   chop $string;

   if ($string) { $string = "($string)";}

   $self->{'_raw_contig_id_list'} = $string;
   return $string;
		   
}


=head2 _use_cext_get

 Title   : _use_cext_get
 Usage   : $obj->_use_cext_get($newval)
 Function: 
 Example : 
 Returns : value of _use_cext_get
 Args    : newvalue (optional)


=cut

sub _use_cext_get{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_use_cext_get'} = $value;
    }
    return $obj->{'_use_cext_get'};

}

=head2 _cext_get_all_SimilarityFeatures_type

 Title   : _cext_get_all_SimilarityFeatures_type
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _cext_get_all_SimilarityFeatures_type{
   my ($self,$db) = @_;

   if( !exists $self->{'_cext_sim_cache'} ) {
       $self->_fill_cext_SimilarityFeature_cache();
   }

   return @{$self->{'_cext_sim_cache'}->{$db}};

}


=head2 _fill_cext_SimilarityFeature_cache

 Title   : _fill_cext_SimilarityFeature_cache
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _fill_cext_SimilarityFeature_cache{
   my ($self) = @_;

   my $host = $self->dbobj->host;
   my $user = $self->dbobj->username;
   my $dbname = $self->dbobj->dbname;
   my $pass = $self->dbobj->password;
   if( !defined $pass ) {
       $pass = '-';
   }

   &Bio::EnsEMBL::Ext::ContigAcc::prepare_Ensembl_cache($host,$user,$pass,$dbname);

   my $glob_start=$self->_global_start;
   my $glob_end=$self->_global_end;
   my $chr_name=$self->_chr_name;
   
   my $fpl = &Bio::EnsEMBL::Ext::ContigAcc::FeaturePairList_by_Score_VC($chr_name,$glob_start,$glob_end,'10');
   
   my $cache = {};
   $self->{'_cext_sim_cache'} = $cache;
   my $save;
   foreach my $f ( $fpl->each_FeaturePair ) {
       my $db = $f->analysis->db;
       if( !defined $cache->{$db} ) {
	   $cache->{$db} = [];
       }
       push(@{$cache->{$db}},$f);
       $save = $f;
   }

   $fpl = 0;
   &Bio::EnsEMBL::Ext::ContigAcc::release_Ensembl_cache();

}

=head2 get_all_coding_Snps

 Title   : get_all_coding_Snps
 Usage   : $vc->get_all_coding_Snps
 Function: get all coding SNPs for the given virual contig
 Example :
 Returns : 
 Args    : Nothing


=cut

sub get_all_coding_Snps{
   my ($self) = @_;
   my @snps;
   my %csnp;
   my %genestr;

   my $glob_start=$self->_global_start;
   my $glob_end=$self->_global_end;
   my $chr_name=$self->_chr_name;
   
   my $query = "select e.contig, e.gene_name, e.exon, s.refsnpid, s.snp_chrom_start from ensembl_lite100.gene_exon as e, ensembl_lite100.gene_snp as s where e.gene_chrom_start >= $glob_start and e.gene_chrom_end <= $glob_end and e.chr_name = '$chr_name' and e.gene = s.gene and s.snp_chrom_start>e.exon_chrom_start and s.snp_chrom_start<e.exon_chrom_end";

   my $sth = $self->dbobj->prepare($query);
   $sth->execute;

   while( my $rowhash = $sth->fetchrow_hashref) {
       my $contig_id = $rowhash->{'contig'};
       my $gene_id = $rowhash->{'gene_name'};
       my $ex_id = $rowhash->{'exon'};
       my $start = $rowhash->{'snp_chrom_start'};
       my $sn_id = $rowhash->{'refsnpid'};
       my $strand = $rowhash->{'exon_chrom_strand'};
       
       $csnp{$sn_id}->{'contig_id'} = $contig_id;
       $csnp{$sn_id}->{'gene_id'} = $gene_id;
       $csnp{$sn_id}->{'exon'} = $ex_id;
       $csnp{$sn_id}->{'snp_id'} = $sn_id;
       $csnp{$sn_id}->{'start'} = $start;
       $csnp{$sn_id}->{'strand'} = $strand;
   
   } 

   my $query2 = "select gene_name,exon, exon_chrom_start, exon_chrom_end, rank, exon_chrom_strand from ensembl_lite100.gene_exon as e where e.gene_chrom_start >= $glob_start and e.gene_chrom_end <= $glob_end";

   my $sth2 = $self->dbobj->prepare($query2);
   $sth2->execute;

   while( my $rowhash2 = $sth2->fetchrow_hashref) {
       my $geneid = $rowhash2->{'gene_name'};
       my $exonid = $rowhash2->{'exon'};
       my $exstart = $rowhash2->{'exon_chrom_start'};
       my $exend = $rowhash2->{'exon_chrom_end'};
       my $rank = $rowhash2->{'rank'};
       my $strand = $rowhash2->{'exon_chrom_strand'};

       $genestr{$exonid}->{'gene_id'} = $geneid;
       $genestr{$exonid}->{'start'} = $exstart;
        $genestr{$exonid}->{'end'} = $exend;
       $genestr{$exonid}->{'rank'} = $rank;
       $genestr{$exonid}->{'strand'} = $strand;
   }
      
   return (\%csnp,\%genestr);
}



1;

