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
use Bio::EnsEMBL::DBSQL::AnalysisAdaptor;
use Bio::EnsEMBL::DBSQL::AssemblyContigAdaptor;
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
    $self->{_anal_hash}={};
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
                                 (sgp.chr_start+sgp.raw_end-f.seq_start-$glob_start)) as end , 
                        IF     (sgp.raw_ori=1,f.strand,(-f.strand)) as strand,
                                f.score,f.analysis, f.name, f.hstart, f.hend, f.hid 
		        FROM   feature f, analysisprocess a,static_golden_path sgp
                        WHERE  f.analysis = a.analysisId 
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


sub get_all_SimilarityFeatures_by_analysis_id {
   my ( $self, $ana_id ) = @_;

   my $glob_start=$self->_global_start;
   my $glob_end  =$self->_global_end;
   my $chr_name  =$self->_chr_name;
   my $idlist    = $self->_raw_contig_id_list();
   my $type      = $self->dbobj->static_golden_path_type;

   
   unless ($idlist){
       return ();
   }
   
   my    $statement = "SELECT f.id, 
                        IF     (sgp.raw_ori=1,(f.seq_start+sgp.chr_start-sgp.raw_start-$glob_start),
                                 (sgp.chr_start+sgp.raw_end-f.seq_end-$glob_start))  as start,  
                        IF     (sgp.raw_ori=1,(f.seq_end+sgp.chr_start-sgp.raw_start-$glob_start),
                                 (sgp.chr_start+sgp.raw_end-f.seq_start-$glob_start))  as end , 
                        IF     (sgp.raw_ori=1,f.strand,(-f.strand)) as strand,
                                f.score,f.analysis, f.name, f.hstart, f.hend, f.hid 
		        FROM   feature f, analysisprocess a,static_golden_path sgp
                        WHERE  f.analysis = a.analysisId
                        AND    f.analysis = $ana_id
                        AND    sgp.raw_id = f.contig
                        AND    f.contig in $idlist
                        AND    sgp.chr_end >= $glob_start 
		        AND    sgp.chr_start <=$glob_end 
		        AND    sgp.chr_name ='$chr_name' 
                        AND    sgp.type = '$type'
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
  $self->throw("Must supply score parameter") unless defined($score);

  if( $self->_use_cext_get() ) {
    return $self->_cext_get_all_SimilarityFeatures_type($analysis_type);
  }
    
  my $glob_start = $self->_global_start;
  my $chr_name   = $self->_chr_name;
  my $idlist     = $self->_raw_contig_id_list;
  my $type       = $self->dbobj->static_golden_path_type;
  my $count = 0;
   
  unless ($idlist){
    print STDERR "There are no contigs to dump features on in get_all_SimilarityFeatures_above_score\n";
    return ();
  }

  #print STDERR "SIMI: got bump factor ",$bp,"\n";

  if( ! defined $self->{'_feature_cache'} || 
      $score < $self->{'_feature_cache_score'} ) {

    &eprof_start('similarity-query');

    #
    # EB. DO NOT remove the sort by hid as it is crucial for the 
    # glob factor below. Tony will have your testicles (unless you
    # are a lady in which case lets not talk about it)
    #

    # SMJS. Remove hit,start ordering if bp is undefined (no globbing should 
    # occur so hid and start ordering not needed). 
    # Note that the globbing code previously did combine features when $bp was 
    # undefined (see NOTES below), but I've fixed this.

    my $order_cols = defined($bp)  ? "hid, start" : "";

    # SMJS Added check for whether feature is at least partly within the
    # golden part of the contig - AND (f.seq_end <= sgp.raw_end) 
    # AND (f.seq_start >= sgp.raw_start). I haven't done this.

    # SMJS Strand used to be got with  IF (sgp.raw_ori=1,f.strand,(-f.strand)),
    
    my $statement = "SELECT f.id,
        IF (sgp.raw_ori=1,(f.seq_start+sgp.chr_start-sgp.raw_start-$glob_start),
            (sgp.chr_start+sgp.raw_end-f.seq_end-$glob_start)) as start,
        IF (sgp.raw_ori=1,(f.seq_end+sgp.chr_start-sgp.raw_start-$glob_start),
            (sgp.chr_start+sgp.raw_end-f.seq_start-$glob_start)), 
        sgp.raw_ori*f.strand,
        f.score,f.analysis, f.name, f.hstart, f.hend, f.hid, f.contig
        FROM   feature f, static_golden_path sgp
        WHERE  sgp.raw_id = f.contig
        AND    f.contig in $idlist
        AND    sgp.type = '$type'
        AND    f.score > $score
        AND    sgp.chr_name = '$chr_name'
        AND    f.seq_start >= sgp.raw_start
        AND    f.seq_end <= sgp.raw_end
        ";
        
    if ($order_cols) {
      $statement .= "ORDER  by $order_cols";
    }

   # PLEASE READ THE COMMENT ABOVE if you want to remove the sort by hid

    #open(T,">>/tmp/stat2.sql");
    #print T $statement,"\n";
    #close(T);
#    print STDERR $statement . "\n";

    my  $sth = $self->dbobj->prepare($statement);    
    $sth->execute(); 

    &eprof_end('similarity-query');
    &eprof_start('similarity-obj');

    my ($fid,$start,$end,$strand,$f_score,$analysisid,$name,
        $hstart,$hend,$hid,$fset,$rank,$fset_score,$contig);
    
    $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$f_score,
        	       \$analysisid,\$name,\$hstart,\$hend,\$hid,\$contig);
    
    
    $self->{'_feature_cache'} = {};
    $self->{'_feature_cache_score'} = $score;

    my $out;
    my $length=$self->length;
    
  FEATURE: 
    my $prev = undef;
      
    while($sth->fetch) {

      my @args=($fid,$start,$end,$strand,$f_score,$analysisid,$name,
                $hstart,$hend,$hid);

      if (($end > $length) || ($start < 1)) {
        next;
      }

# NOTE 1 Without the if (defined($bp)), even when $bp was undefined features 
# were being merged if $prev->end > $start. This is undesirable for gene 
# building as it will lead to incorrect hstart and hend values.

      if (defined($bp)) {

# NOTE 2 Without the ORDER BY hid, start above this will discard features which
# should be retained so protect with the if (defined($bp)) condition (because 
# when $bp is undefined no sorting on start is done [to speed up the query for 
# gene building]).

        if (defined $prev && $prev->hseqname eq $hid && 
            $prev->end + $bp > $start ) {
      
          $prev->end($end);
          next;
        }
      } 

      # &eprof_start('similarity-obj-creation');
      my $analysis=$self->_get_analysis($analysisid);
       
      #warn $analysis->logic_name;
        
      if( !defined $name ) {
        $name = 'no_source';
      }
          
      $out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();   
      $out->set_all_fields($start,$end,$strand,$f_score,$name,'similarity',
                           $contig, $hstart,$hend,1,$f_score,$name,
                           'similarity',$hid);

      $out->analysis($analysis);
      $out->id      ($hid);
        
      if( !defined $self->{'_feature_cache'}->{$analysis->db()} ) {
        $self->{'_feature_cache'}->{$analysis->db()} = [];
      }
      push( @{$self->{'_feature_cache'}->{$analysis->db()}}, $out );
        
      $prev = $out;
      $count++;
      # &eprof_end('similarity-obj-creation');
    }
      
    #print STDERR "FEATURE: got $count in entire call\n";
    &eprof_end('similarity-obj');
      
  }

  # now extract requested features

  my @features;
  foreach my $feature ( @{$self->{_feature_cache}->{$analysis_type}} ) {
    if( $feature->score() > $score ) {
      push( @features, $feature );
    } 
  }
 return @features;
}

=head2 get_all_SimilarityFeatures_by_strand

 Title   : get_all_SimilarityFeatures_by_strand
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SimilarityFeatures_by_strand{
   my ($self, $analysis_type, $score, $bp,$strand) = @_;

   &eprof_start("similarity-by-strand");
   &eprof_start("similarity-by-strand-call-to-score");
   my @features = $self->get_all_SimilarityFeatures_above_score($analysis_type, $score, $bp);
   &eprof_end("similarity-by-strand-call-to-score");
   my (@f);

   foreach my $f (@features) {
       if ($f->strand == 1 && $strand == 1) {
	   unshift(@f,$f);
       }
       elsif ($strand == -1 && $f->strand == -1) {
	   push(@f,$f);
       }
   }
   &eprof_end("similarity-by-strand");
   #print STDERR "Returning ",scalar(@f)," features for ",$analysis_type,"\n";
   return (@f);
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
                        AND    f.score > 300
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
    
  my $glob_start = $self->_global_start;
  my $glob_end   = $self->_global_end;
  my $chr_name   = $self->_chr_name;
    
  my ($fid,$start,$end,$strand,$f_score,$analysisid,$name,
      $hstart,$hend,$hid,$fset,$rank,$fset_score,$contig,$chr_start,$chr_end,
      $raw_ori);
    

  # this needs to be rewritten properely EB

  my $type = $self->dbobj->static_golden_path_type;

#    my $statement = "SELECT f.id,
#         if(s.raw_ori=1,(f.seq_start-s.raw_start+s.chr_start),(s.chr_start+s.raw_end-f.seq_start)),
#         if(s.raw_ori=1,(f.seq_end  -s.raw_start+s.chr_start),(s.chr_start+s.raw_end-f.seq_end)),
#         f.strand * s.raw_ori,
#         f.score,f.analysis, f.name, f.hstart, f.hend, f.hid
#      FROM feature f, static_golden_path s 
#      WHERE f.perc_id > $pid 
#      AND   s.raw_id  = f.contig
#      AND NOT (s.chr_start > $glob_end) 
#      AND NOT (s.chr_end < $glob_start) 
#      AND   f.seq_start > s.raw_start 
#      AND   f.seq_end   < s.raw_end
#      AND   s.chr_name  = '$chr_name'
#      AND   s.type = '$type'";

  my $idlist     = $self->_raw_contig_id_list;
  unless ($idlist){
    return ();
  } 

  my $statement = "SELECT f.id,
        IF (sgp.raw_ori=1,(f.seq_start+sgp.chr_start-sgp.raw_start),
            (sgp.chr_start+sgp.raw_end-f.seq_end)) as start,
        IF (sgp.raw_ori=1,(f.seq_end+sgp.chr_start-sgp.raw_start),
            (sgp.chr_start+sgp.raw_end-f.seq_start)), 
        f.strand * sgp.raw_ori,
        f.score,f.analysis, f.name, f.hstart, f.hend, f.hid, f.contig
        FROM   feature f, static_golden_path sgp
        WHERE  sgp.raw_id = f.contig
        AND    f.contig in $idlist
        AND    sgp.type = '$type'
        AND    sgp.chr_name = '$chr_name'
        AND    f.seq_start >= sgp.raw_start 
        AND    f.seq_end   <= sgp.raw_end
        AND    f.perc_id > $pid
        ORDER  by hid
       ";
   

  # print $statement . "\n";

    
  my  $sth = $self->dbobj->prepare($statement);    
  $sth->execute(); 

  $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$f_score,\$analysisid,
		     \$name,\$hstart,\$hend,\$hid,\$contig);


  my @array;
  my %analhash;
  my $out;
    

  FEATURE: while($sth->fetch) {
    my $out;
    my $analysis;

    # clip and map to vc coordinates

    if ($start >= $glob_start && $end <= $glob_end) {

      $start = $start-$glob_start;
      $end   = $end-$glob_start;

      # create features

      if (!$analhash{$analysisid}) {
        my $analysis_adp=new Bio::EnsEMBL::DBSQL::AnalysisAdaptor($self->dbobj);
        $analysis = $analysis_adp->fetch_by_dbID($analysisid);
        $analhash{$analysisid} = $analysis;           
      } 
      else {$analysis = $analhash{$analysisid};}
          
      if( !defined $name ) {
        $name = 'no_source';
      }

      # These conditions should not be needed, but incorrect data in the 
      # human_live database at one time lead to their inclusion
      if ($start > $end) {
        my $tmp = $start;
        $start = $end;
        $end = $tmp;
      }
      if ($hstart > $hend) {
        my $tmp = $hstart;
        $hstart = $hend;
        $hend = $tmp;
      }

      $out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();   
      $out->set_all_fields($start,$end,$strand,$f_score,$name,'similarity',
                           $contig,$hstart,$hend,1,$f_score,$name,'similarity',
                           $hid);

      $out->analysis    ($analysis);
      $out->id          ($hid);              
      $out->seqname     ($self->id);
      $out->start       ($start);
      $out->end         ($end);
      $out->strand      ($strand);
      $out->source_tag  ($name);
      $out->primary_tag ('similarity');
      $out->id          ($hid);
          
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

=head2 get_all_SNPFeatures

 Title   : get_all_SNPFeatures
 Usage   : foreach my $sf ( $contig->get_all_SNPFeatures ) 
 Function: Returns snp features from lite adaptor
 Example :
 Returns : 
 Args    :


=cut



sub get_all_SNPFeatures {
  my ($self,$bp) = @_;
  my @snps = $self->dbobj->get_LiteAdaptor->fetch_snp_features
    (
     $self->_chr_name, 
     $self->_global_start, 
     $self->_global_end,
     $bp
    );
  return @snps;
}

sub get_all_SNPFeatures_lite {
  my ($self,$bp) = @_;
  return $self->dbobj->get_LiteAdaptor->fetch_virtualsnps(
     $self->_chr_name, $self->_global_start, $self->_global_end, $bp );
}

sub get_all_virtualfeatures_lite {
  my ($self,$type,$score,$bp) = @_;
  return $self->dbobj->get_LiteAdaptor->fetch_virtualfeatures(
     $self->_chr_name, $self->_global_start, $self->_global_end, $type,$score,$bp );
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
   my $analysis_types="('genscan', 'fgenesh')";
   
   my $type = $self->dbobj->static_golden_path_type;
   my $query = "
    SELECT f.id
      , sgp.raw_start
      , sgp.raw_end
      , f.seq_start
      , f.seq_end
      , IF (sgp.raw_ori=1,(f.seq_start+sgp.chr_start-sgp.raw_start-$glob_start )
          , (sgp.chr_start+sgp.raw_end-f.seq_end-$glob_start) ) as start
      , IF (sgp.raw_ori=1,(f.seq_end+sgp.chr_start-sgp.raw_start-$glob_start)
          , (sgp.chr_start+sgp.raw_end-f.seq_start-$glob_start) )
      , IF (sgp.raw_ori=1,f.strand,(-f.strand))
      , f.score
      , f.evalue
      , f.perc_id
      , f.phase
      , f.end_phase
      , f.analysis
      , f.hid
      , f.contig,f.name
    FROM feature f
      , static_golden_path sgp
    WHERE sgp.raw_id = f.contig
      AND f.contig IN $idlist
      AND f.name IN $analysis_types
      AND sgp.type = '$type'
      AND sgp.chr_name='$chr_name'
    ORDER BY f.contig
      , f.strand * f.seq_start";

   my $sth = $self->dbobj->prepare($query);
  
   
   $sth->execute();
   
   
   my ($fid,$rawstart,$rawend,$seqstart,$seqend,$start,$end,$strand,$score,$evalue,$perc_id,$phase,$end_phase,$analysisid,$hid,$contig,$name);
   
   # bind the columns
   $sth->bind_columns(undef,\$fid,\$rawstart,\$rawend,\$seqstart,\$seqend,\$start,\$end,\$strand,\$score,\$evalue,\$perc_id,\$phase,\$end_phase,\$analysisid,\$hid,\$contig,\$name);
      
   $previous = -1;
   my $current_fset;
   my $fsetstart;
   my $count;
   my $prev;

   while( $sth->fetch ) {
       
       
       if (($end > $length) || ($start < 1)) {
	   next;
       }
       if ($seqstart < $rawstart || $seqend > $rawend) {
           next;
       }
       
       #MC This is one humdinger of a hack to get rid of duplicate genscans

       if (defined $prev && $start == $prev->start && $end == $prev->end) {
	 next;
       }
       
       
       my $out;
       
       my $analysis;
       
       if (!$analhash{$analysisid}) {
	 my $analysis_adp = new Bio::EnsEMBL::DBSQL::AnalysisAdaptor($self->dbobj);
	 $analysis = $analysis_adp->fetch_by_dbID($analysisid);
	 $analhash{$analysisid} = $analysis;	   
	   
       } else {
	   $analysis = $analhash{$analysisid};
       }
        
       # MC. Temporarily changed back the genscan fetching for a build.
#       if( $hid =~ /Initial/ || $hid =~ /Single/ || $previous =~ /Single/ || $previous =~ /Terminal/ || $previous eq -1 || $previous_contig != $contig) {
       if( $hid ne $previous|| $previous eq -1 || $previous_contig != $contig) {
#        if($phase != $previous || $previous_contig != $contig $previous == -1) {
	   $count++;
	   $current_fset = Bio::EnsEMBL::SeqFeature->new();
	   $current_fset->source_tag($name);
	   $current_fset->primary_tag('prediction');
	   $current_fset->analysis($analysis);
	   $current_fset->seqname($self->id);
	   $current_fset->id($count);
	   $current_fset->score(0.0);
	   $fsetstart = $seqstart;
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
        

       #MC Now fset ids have never worked well for genscans.  The new id
       #is contig internal_id . start in raw contig coords.  Ugly I know but 
       #it works.

       #my $query="select fset from fset_feature where feature=$fid"; 
       #my $sth = $self->dbobj->prepare($query);
       #$sth->execute();
       #my $arr_ref=$sth->fetchrow_arrayref;
       
       #my $fsetid=$arr_ref->[0];

       my $fsetid = $contig . "." .  $fsetstart;

       $out->id($fsetid); # to make genscan peptide work

       $out->source_tag($name);
       $out->primary_tag('prediction');
       
       if( defined $score ) {
	   $out->score($score);
       }

       $out->analysis($analysis);

       # Final check that everything is ok.
       
       $out->validate();

       $current_fset->add_sub_SeqFeature($out,'EXPAND');
       $current_fset->strand($strand);

       $previous        = $hid;
       $previous_contig = $contig;
       $prev            = $out;
  }
 

   $self->{'_genscan_cache'} = \@array;

   return @array;
}


sub get_all_PredictionFeatures_by_analysis_id {
    my ( $self, $ana_id ) = @_;

    $self->throw("No analysis_id given") unless $ana_id;

    my @array;
    #my $id     = $self->internal_id();
    my $length = $self->length();
   
   
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
    my $type = $self->dbobj->static_golden_path_type;
    my $transcript_pos = -1;
	
    # make the SQL query
    my $query =
                        "SELECT f.id
                          , sgp.raw_start
                          , sgp.raw_end
                          , f.seq_start
                          , f.seq_end
                          , IF (sgp.raw_ori=1,(f.seq_start+sgp.chr_start-sgp.raw_start-$glob_start)
                              , (sgp.chr_start+sgp.raw_end-f.seq_end-$glob_start))  as start
                          , IF (sgp.raw_ori=1,(f.seq_end+sgp.chr_start-sgp.raw_start-$glob_start)
                              , (sgp.chr_start+sgp.raw_end-f.seq_start-$glob_start))  as end
                          , IF (sgp.raw_ori=1,f.strand,(-f.strand))
                          , f.score
                          , f.evalue
                          , f.perc_id
                          , f.phase
                          , f.end_phase
                          , f.analysis
                          , f.hid
                          , f.contig
                          , f.name
                        FROM feature f
                          , static_golden_path sgp
                        WHERE sgp.raw_id = f.contig
                          AND f.contig IN $idlist
                          AND f.analysis = $ana_id
                          AND sgp.type = '$type'
                          AND sgp.chr_name='$chr_name'
                        ORDER BY f.contig,f.strand*f.seq_start
                        ";
   
   
   my $sth = $self->dbobj->prepare($query);
   
   $sth->execute();
  
   my ($fid, $rawstart, $rawend, $seqstart,
        $seqend,$start,$end,$strand,$score,
        $evalue,$perc_id,$phase,$end_phase,
        $analysisid,$hid,$contig,$name);
   
   
   # bind the columns
   $sth->bind_columns(undef,
   \$fid,\$rawstart,\$rawend,\$seqstart,
   \$seqend,\$start,\$end,\$strand,\$score,
   \$evalue,\$perc_id,\$phase,\$end_phase,
   \$analysisid,\$hid,\$contig,\$name);
   	
    my $current_fset;
    my $fsetstart;
    my $count = 1;
    my $prev;
    
    while ( $sth->fetch ) {
       if (($end > $length) || ($start < 1)) {
	   next;
       }
       if ($seqstart < $rawstart || $seqend > $rawend) {
           next;
       }
       
       #MC This is one humdinger of a hack to get rid of duplicate genscans

       if (defined $prev && $start == $prev->start && $end == $prev->end) {
	 next;
       }
        my $out;
	$prev     = $out;
        $count++;
        my $analysis;
        
       				
        if ( !$analhash{$analysisid} ) {
            my $analysis_adp =
              new Bio::EnsEMBL::DBSQL::AnalysisAdaptor( $self->dbobj );
            $analysis = $analysis_adp->fetch_by_dbID($analysisid);
            $analhash{$analysisid} = $analysis;

        }
        else {
            $analysis = $analhash{$analysisid};
        }
				
        # Oh boyoboy.  Yet another genscan hack to avoid duplicate genscans
        if ( defined($prev) && $start == $prev->start && $end == $prev->end ) {
            next;
        }
	        
        my ( $transcript_num, $exon_num ) = $hid =~/(\d+)\.(\d+)$/;
        if ($transcript_num =~/(\D)/){ warn 'Transcript Number contains non digit'};
		
	#print "Transcript :",$transcript_num,"\n";
		
	if ($transcript_num != $transcript_pos ) {			
            $current_fset = Bio::EnsEMBL::SeqFeature->new();  
            $current_fset->source_tag($name);
            $current_fset->primary_tag('prediction');
            $current_fset->analysis($analysis);
            $current_fset->seqname( $self->id );
            #$current_fset->id( $self->internal_id . "." . $fsetstart );
            $current_fset->score(0.0);
            $fsetstart = $seqstart;
            $current_fset->raw_seqname( $self->id );
            push ( @array, $current_fset );

            $transcript_pos = $transcript_num;			
        }

        $out = Bio::EnsEMBL::SeqFeature->new;

        $out->seqname( $self->id );
        $out->raw_seqname( $self->id );
        $out->start($start);
        $out->end($end);
        $out->strand($strand);
        $out->score($score)         if ( defined $score );
        $out->p_value($evalue)      if ( defined $evalue );
        $out->percent_id($perc_id)  if ( defined $perc_id );
        $out->phase($phase)         if ( defined $phase );
        $out->end_phase($end_phase) if ( defined $end_phase );
     
        $out->id($fsetid);    # to make genscan peptide work
        $out->source_tag($name);
        $out->primary_tag('prediction');

        if ( defined $score ) {
            $out->score($score);
        }

        $out->analysis($analysis);

        # Final check that everything is ok.
		
        $out->validate();
        $current_fset->add_sub_SeqFeature( $out, 'EXPAND' );
        $current_fset->strand($strand);
        $previous        = $hid;
        $previous_contig = $contig;
        $prev            = $out;
		
    }
    
    return @array;
}


=head2 get_all_SimpleFeatures

    my @feat = $contig->get_all_SimpleFeatures

Returns a list of B<Bio::EnsEMBL::SeqFeature>
objects.  Features which overlap the ends of the
contig are truncated to the contig, whereas the
other methods which return FeaturePairs discard
them.

=head2 get_all_SimpleFeatures_by_feature_type

    my @cpg = $contig->get_all_SimpleFeatures_by_feature_type(
        'cpg_island');

Same as B<get_all_SimpleFeatures>, but returns
only the features where the B<gff_feature> column
in the analsis table matches the string given in
the argument.

=head2 get_all_SimpleFeatures_by_feature_type_above_score

    my @cpg = $contig->get_all_SimpleFeatures_by_feature_type_above_score(
        'cpg_island', 400);

Extends B<get_all_SimpleFeatures_by_feature_type>,
taking an additional argument which is a number. 
Only features with a score higher than this
number are returned.

=head2 get_all_SimpleFeatures_by_analysis_id

    my @cpg = $contig->get_all_SimpleFeatures_by_analysis_id(21);

Same as B<get_all_SimpleFeatures>, but returns
only the features where the analysis_id matches
the integer argument.

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
        AND a.gff_feature = '$feature_type'
        AND f.score >= $score
        });
}

sub get_all_SimpleFeatures_by_analysis_id {
    my( $self, $analysis_id ) = @_;
    
    $self->throw("No analysis ID given") unless $analysis_id;
    
    return $self->_fetch_SimpleFeatures_SQL_clause(qq{
        AND a.analysisId = $analysis_id
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
              , (sgp.chr_start + sgp.raw_end - f.seq_start   - $global_start)) as end
          , sgp.raw_ori * f.strand
          , f.score
          , f.analysis
          , f.name
          , f.hid
          , f.perc_id
          , f.evalue
        FROM feature f
          , analysisprocess a
          , static_golden_path sgp
        WHERE f.analysis = a.analysisId
          AND sgp.raw_id = f.contig
          AND f.contig IN $idlist
          AND sgp.type = '$type'
          AND sgp.chr_name = '$chr_name'
        };
    
    # All statements have this ORDER by clause on the end
    my $sql_order = qq{
        ORDER BY seq_start
       };

    # Make the full statement and execute it
    my $sql = join(' ', $sql_begin, $sql_extra, $sql_order);
    my $sth = $self->dbobj->prepare($sql);    
    $sth->execute();     
    
    # Bind columns to variables for the fastest possible
    # retrieval of results.
    my ($fid,
        $start, $end, $strand,
        $f_score, $analysis_id, $name, $hid, $f_perc_id, $f_evalue);
    $sth->bind_columns(undef,
        \$fid,
        \$start, \$end, \$strand,
        \$f_score, \$analysis_id, \$name, \$hid, 
        \$f_perc_id, \$f_evalue);
        
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
	  my $analysis_adp = new Bio::EnsEMBL::DBSQL::AnalysisAdaptor($self->dbobj);
	  $analysis = $analysis_adp->fetch_by_dbID($analysis_id);
	    $anal{$analysis_id} = $analysis;
        }
    
        # Truncate the feature coordinates to the region
        $start = 1       if $start < 1;
        $end   = $length if $end   > $length;
    
        # Make the feature
        my $feat = Bio::EnsEMBL::FeatureFactory->new_feature;
        $feat->id        ($fid);
        $feat->start     ($start);
        $feat->end       ($end);
        $feat->strand    ($strand);
        $feat->score     ($f_score);
        $feat->analysis  ($analysis);
        $feat->source_tag($name);
        $feat->percent_id($f_perc_id) if ( defined $f_perc_id );
        $feat->p_value   ($f_evalue)  if ( defined $f_evalue ); 
        
        push(@features, $feat);
    }
    return @features;
}

=head2 get_all_DASFeatures

 Title   : get_all_DASFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut
### This is a hacky little function which gets the name of the contig
### from a finished clone - or returns undef otherwise.
### A finised clone is a clone with only one contig which starts at bp 1
### (So grab the last contig of the clone, and see if it starts at 1 <g>)
sub contig_from_clone {
    my ($self,$clone) = @_;
    my $sth = $self->dbobj->prepare(
            "select co.id from clone as cl, contig as co
              where co.clone = cl.internal_id and cl.id = '$clone'
              order by offset desc limit 1"
    );
    $sth->execute();
    my ($contig) = $sth->fetchrow_array();
    return $contig=~/^[^\.]+\.\d+\.1\./ ? $contig : undef;
}

sub get_all_DASFeatures{
   my ($self,@args) = @_;

   if( defined $self->{'_das_cached_features'} ) {
       return @{$self->{'_das_cached_features'}};
   }

   my @contig_features;
   my @chr_features;
   my @fpc_features;
   my @clone_features;
   my @genomic_features;

   my @rawcontigs = $self->_vmap->get_all_RawContigs();
   @rawcontigs = map { $_->id() } @rawcontigs;
   #foreach (@rawcontigs){
   #    print STDERR "Raw contig: ", $_, "\n";
   #}
    
    # need a call here to get a list of FPC contigs that overlap my VC
    # I also need to have their VC start end in the FPC coordinates.
    # and somehow pass all this stuff down to the DAS fetcher...eek!
   my @fpccontigs = (undef);
    
   #print STDERR "DAS ", $self->_chr_name,$self->_global_start,$self->_global_end, "\n";
   my @clones  = $self->get_all_Clones();
   #foreach (@clones){
   #    print STDERR "Clone: ", $_, "\n";
   #}
   my $chr_length = $self->fetch_chromosome_length();       
   foreach my $extf ( $self->dbobj->_each_DASFeatureFactory ) {
       
       if( $extf->can('get_Ensembl_SeqFeatures_DAS') ) {
	       foreach my $sf (
                $extf->get_Ensembl_SeqFeatures_DAS(
                    $self->_chr_name,$self->_global_start,$self->_global_end,
                    \@fpccontigs, \@clones,\@rawcontigs, $chr_length)
            ) {


# BAC.*_C are fly contigs....
# CRA_x are Celera mosquito contigs....

	           if( $sf->seqname() =~ /(\w+\.\d+\.\d+.\d+|BAC.*_C)|CRA_.*/ ) {
#                    warn ("Got a raw contig feature: ", $sf->seqname(), "\n");
 		            push(@contig_features,$sf);
               } elsif( $sf->seqname() =~ /chr[\d+|X|Y]/i) { 
#                    warn ("Got a chromosomal feature: ", $sf->seqname(), "\n");
 	                push(@chr_features, $sf);
               } elsif( $sf->seqname() =~ /^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|24|X|Y|2L|2R|3L|3R)$/o) {  # breaks on mouse!
#                    warn ("Got a chromosomal feature: ", $sf->seqname(), "\n");
 	                push(@chr_features, $sf);
               } elsif( $sf->seqname() =~ /ctg\d+|NT_\d+/i) { 
#                    warn ("Got a FPC contig feature: ", $sf->seqname(), "\n");
 	                push(@fpc_features, $sf);
               } elsif( $sf->seqname() =~ /^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|X)\.\d+\-\d+/i) { 
#                    warn ("Got a mouse clone feature: ", $sf->seqname(), "\n");
 	                push(@contig_features, $sf);
               } elsif( $sf->seqname() =~ /\w{1,2}\d+/i) { 
#                    print STDERR "CLONE >".$sf->seqname()."<\n";
                    if(my $contig_from_clone = $self->contig_from_clone($sf->seqname()) ) {
#                        print STDERR "CONTIG NAME FROM CLONE >$contig_from_clone<\n";
                        $sf->seqname($contig_from_clone);
 	                    push(@contig_features, $sf);
                    }
#                    warn ("Got a clone feature: ", $sf->seqname(), "\n");
               } elsif( $sf->das_type_id() eq '__ERROR__') { 
#                    Always push errors even if they aren't wholly within the VC
	                push(@genomic_features, $sf);
               } elsif( $sf->seqname() eq '') { 
                    #suspicious
	                warn ("Got a DAS feature with an empty seqname! (discarding it)\n");
	           } else {
		            warn ("Got a DAS feature with an unrecognized segment type: >", $sf->seqname(), "< >", $sf->das_type_id(), "<\n");
	           }
	       }
	   
       } else {
	        warn "StaticContig: Got a DAS feature factory that can't do get_Ensembl_SeqFeatures_DAS\n";
	        #$self->throw("Got a DAS feature factory that can't do get_Ensembl_SeqFeatures_DAS");
       }
   }
   
    my $chr_start = $self->_global_start();
    my $chr_end   = $self->_global_end();
   foreach my $sf ( @contig_features ) {
#            print STDERR "SEG ID: ",         $sf->seqname(), "\t";
#            print STDERR "ID: ",             $sf->das_id(), "\t";
#            print STDERR "DSN: ",            $sf->das_dsn(), "\t";
#            print STDERR "FEATURE START: ",  $sf->das_start(), "\t";
#            print STDERR "FEATURE END: ",    $sf->das_end(), "\t";
#            print STDERR "FEATURE STRAND: ", $sf->das_strand(), "\t";
#            print STDERR "FEATURE TYPE: ",   $sf->das_type_id(), "\n";
            my( $st, $en, $str) = $self->_convert_start_end_strand_vc( $sf->seqname(), $sf->das_start,$sf->das_end, $sf->das_strand );
            $sf->das_start(  $st );
            $sf->das_end(    $en );
            $sf->das_strand( $str );
#       if( defined $self->_convert_seqfeature_to_vc_coords($sf) ) {
        if($sf->das_start <= $chr_end-$chr_start && $sf->das_end >= 1) {
            push(@genomic_features, $sf);
        }
#       }
       #else {
       # print STDERR "Binning ", $sf->seqname(), "\n";
       #}
   }

    
   my $xx = 1;
   foreach my $sf ( @chr_features ) {
#       print STDERR "$xx BEFORE: ", $sf->seqname() , "\t";
#      print STDERR "$xx BEFORE: ", $sf->seqname() , "\t";
#           print STDERR "SC SEG ID: ",         $sf->seqname(), "\t";
#           print STDERR "SC DSN: ",            $sf->das_dsn(), "\t";
#           print STDERR "SC FEATURE START: ",  $sf->das_start(), "\t";
#           print STDERR "SC FEATURE END: ",    $sf->das_end(), "\t";
#           print STDERR "SC FEATURE STRAND: ", $sf->das_strand(), "\t";
#            print STDERR "SC FEATURE TYPE: ",   $sf->das_type_id(), "\n";
#       print STDERR "FEATURE START: ",  $sf->start() , "\t";
#       print STDERR "FEATURE END: ",    $sf->end() , "\t";
#       print STDERR "FEATURE STRAND: ", $sf->strand() , "\t";
#       print STDERR "FEATURE ID: ",   $sf->das_feature_id(), "\n";
       if( defined $self->_convert_chrfeature_to_vc_coords($sf) ) {
#       		print STDERR "$xx AFTER: ", $sf->seqname() , "\t";
#            print STDERR "FEATURE START: ",  $sf->das_start(), "\t";
#            print STDERR "FEATURE END: ",    $sf->das_end(), "\t";
#            print STDERR "FEATURE STRAND: ", $sf->das_strand(), "\t";
#       		print STDERR "FEATURE ID: ",   $sf->das_feature_id(), "\n";
#            
#			print STDERR "SEG ID: ",         $sf->seqname(), "\t";
#            print STDERR "DSN: ",            $sf->das_dsn(), "\t";
#            print STDERR "FEATURE TYPE: ",   $sf->das_type_id(), "\n";
	        push(@genomic_features, $sf);
       }
       $xx++;
   }
   $self->{'_das_cached_features'} = \@genomic_features;
   return @genomic_features;
}

=head2 _convert_chrfeature_to_vc_coords

 Title   : _convert_chrfeature_to_vc_coords
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _convert_chrfeature_to_vc_coords{
    my ($self, $f) = @_;

    my $chr_start = $self->_global_start();
    my $chr_end   = $self->_global_end();
    
    if($f->das_start > $chr_end || $f->das_end < $chr_start) {
#        print STDERR "DAS ERROR! Feature not on VC between $chr_start and $chr_end: START: ",
#                    $f->das_start,' END: ',$f->das_end,' STRAND: ',$f->das_strand, ' ID: ', $f->das_feature_id,
#                    "\n";
        return ();
    }

    $f->das_start( $f->das_start() - $chr_start +1 );
    $f->das_end(   $f->das_end() - $chr_start  + 1);
#    print STDERR "---> ",$f->das_start(), " , ", $f->das_end(),"\n";       
    return($f);
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
   my @rawcontigs = $self->_vmap->get_all_RawContigs();
   if( scalar(@rawcontigs) == 0) {
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
	        $self->throw("Should add DAS feature factories to add_DASFeatureFactory");
        } else {
	        push(@std,$extf);
        }
   }

   #Build needed arrays and hashes in one go
   my %int_ext;
   my @clones;
   my %cloneh;
   my %featureh;
   my @cintidlist;
   foreach my $contig (@rawcontigs) {
       $int_ext{$contig->internal_id}=$contig->id;
       push (@cintidlist,$contig->internal_id);
       if( !defined $cloneh{$contig->cloneid} ) {
	    $cloneh{$contig->cloneid} = [];
	    $featureh{$contig->cloneid} = [];
       }
       my $string = $contig->cloneid.".".$contig->seq_version;
       push(@clones,$string);
       push(@{$cloneh{$contig->cloneid}},$contig);
   }
   

   &eprof_start("External-feature-web-get");

   if( scalar(@web) > 0 ) {
      
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
    	   my @features = sort { $a->start <=> $b->start }             @{ $featureh{$clone} };
	       my @contigs  = sort { $a->embl_offset <=> $b->embl_offset } @{ $cloneh{$clone}   };
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

   &eprof_end("External-feature-web-get");

   &eprof_start("External-feature-std");

   #Standard EFFs now take a list of contig ids
   if( scalar(@std) > 0 ) {
       foreach my  $extf ( @std ) {
	   &eprof_start("external_get_contig_list".$extf);

	   if ( $extf->can('get_Ensembl_SeqFeatures_contig_list') ) {
               if (scalar(@cintidlist) == 1) {
                  my ($contig) = $self->_vmap->each_MapContig;
	          push(@contig_features,$extf->get_Ensembl_SeqFeatures_contig_list(\%int_ext,\@cintidlist,$contig->rawcontig_start, $contig->rawcontig_start + $self->_vmap->length));
                } else {
	          push(@contig_features,$extf->get_Ensembl_SeqFeatures_contig_list(\%int_ext,\@cintidlist));
                }
	   }
	   &eprof_end("external_get_contig_list".$extf);
       }
   }

   &eprof_end("External-feature-std");


   &eprof_end("External-feature-get");

   # ok. Now @contig_features are in contig coordinates. Map up.
   
   # this is the simplest way. We can do this more efficiently if need be
   
   
   &eprof_start("External-coordinate-lift");

   my @final;

### Hack by the web team - James Smith (js5) and Tony Cox (avc)
###
### The following code contains a hack specifically designed for external 
### datasources if there is a problem we return a single object with ->id 
### of "__ERROR__" and ->das_id of error code. This then gets passed 
### through even when the feature does not lie on the virtual contig.
###
### This allows the drawing code to display information about the track
### to indicate that it is empty because of no data or empty because of
### a failure to retrieve the information.

   #print STDERR "Got ",scalar(@contig_features),"before lift\n";
   foreach my $f ( @contig_features ) {
       if( defined $self->_convert_seqfeature_to_vc_coords($f) ) {
                push(@final, $f);
       } elsif( $f->id eq '__ERROR__') { #Always push errors even if they aren't wholly within the VC
            push(@final, $f);
       }
   }
   
   #print STDERR "Got ",scalar(@final),"after lift\n";

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
    my $flag = shift;
    
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

    my @fpcclones;
    if($flag eq 'FISH') {
        @fpcclones = $chr->get_Clones_by_start_end_including_FISH($glob_start,$glob_end);
    } else {
        @fpcclones = $chr->get_Clones_by_start_end($glob_start,$glob_end);
    }

    foreach my $clone ( @fpcclones ) {
        my $newstart = $clone->start - $glob_start;
       	$clone->start($newstart);
    }

    return @fpcclones;
}

sub get_cloneset_on_chromosome {
    my $self = shift;
    my $cloneset = shift;
    my $chr_name      = shift; 

    my $mapdb;
    eval {
        $mapdb = $self->dbobj->mapdb();
    };
    if( $@ || !defined $mapdb ) {
	    $self->warn("in get_all_clones_in_cloneset, unable to get mapdb. Returning empty list [$@]");
	    return ();
    }
    my $fpcmap = $mapdb->get_Map('FPC');
    $chr_name =~ s/chr//g;
    my $chr_name_X = $self->_chr_name;
    $chr_name_X =~ s/chr//g; 

    my $chr = $fpcmap->get_ChromosomeMap($chr_name||$chr_name_X);
    $cloneset = $chr->get_cloneset_on_chromosome($cloneset, $chr_name);
    return @$cloneset;
}

sub get_all_clones_in_cloneset {
    my $self = shift;
    my $cloneset = shift;
    
    my $mapdb;
    eval {
        $mapdb = $self->dbobj->mapdb();
    };
    if( $@ || !defined $mapdb ) {
	    $self->warn("in get_all_clones_in_cloneset, unable to get mapdb. Returning empty list [$@]");
	    return ();
    }
    my $glob_start = $self->_global_start;
    my $glob_end   = $self->_global_end;
    my $length     = $self->length;
    my $chr_name   = $self->_chr_name;
    my $fpcmap = $mapdb->get_Map('FPC');
    $chr_name =~ s/chr//g;

    my $chr = $fpcmap->get_ChromosomeMap($chr_name);

    my @cloneset = $chr->get_all_clones_in_cloneset($cloneset, $glob_start, $glob_end);
    return @cloneset;
}
sub get_landmark_MarkerFeatures_old {

my ($self) = @_;

$self->throw( "Method deprecated." );

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
	    
	    if( defined $prev && $prev->end+$glob > $start && uc($synonym) eq uc($prev->id)) {
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


   my $statement= " SELECT  chr_start,
			    chr_end,
			    chr_strand,
			    name 
		    FROM    landmark_marker 
		    WHERE   chr_name = '$chr_name'
                    AND     chr_start >= $glob_start 
                    AND     chr_end <= $glob_end 
		    ORDER BY chr_start
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
       
       if( defined $prev && $prev->end + $glob > $start  && uc($prev->id) eq uc($name) ) {           
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
                          FROM      $dbname.feature f, $dbname.analysisprocess a, 
                                    $mapsdbname.MarkerSynonym s, 
                                    $dbname.static_golden_path sgp
                          WHERE     sgp.raw_id=f.contig  
                          AND       f.hid=s.marker
                          AND       sgp.chr_name = '$chr_name' 
                          AND       sgp.type = '$type'
                          AND       f.analysis = a.analysisId 
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
                          FROM      $dbname.feature f, $dbname.analysisprocess a, 
                                    $mapsdbname.MarkerSynonym s, 
                                    $dbname.static_golden_path sgp
                          WHERE     sgp.raw_id=f.contig  
                          AND       f.hid=s.marker
                          AND       sgp.chr_name='$chr_name'
                          AND       sgp.type = '$type'
                          AND       f.analysis = a.analysisId 
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
 Function: Get all genes making sure there is no redundant exons
           Suitable if no translation has to be intended on transcript.
           Much faster than get_all_Genes.
 Example :
 Returns : Array of Bio::EnsEMBL::Gene
 Args    :

=cut

sub get_all_Genes_exononly{
  my ($self) = @_;

  my $glob_start = $self->_global_start;
  my $glob_end   = $self->_global_end;
  my $chr_name   = $self->_chr_name;
  my $idlist     = $self->_raw_contig_id_list();
  my $type       = $self->dbobj->static_golden_path_type;
   
  unless ($idlist){
    return ();
  }

  if( my $genes = $self->{'_all_Genes_exononly'} ) {
    return @$genes;
  }

  #  IF     (sgp.raw_ori=1,e.strand,(-e.strand))
  my $query = "
    SELECT e.exon_id,e.sticky_rank,e.phase,et.rank,et.transcript_id,t.gene_id, 
    IF     (sgp.raw_ori=1,(e.seq_start+sgp.chr_start-sgp.raw_start-$glob_start),
           (sgp.chr_start+sgp.raw_end-e.seq_end-$glob_start)) as start,  
    IF     (sgp.raw_ori=1,(e.seq_end+sgp.chr_start-sgp.raw_start-$glob_start),
           (sgp.chr_start+sgp.raw_end-e.seq_start-$glob_start)), 
    sgp.raw_ori*e.strand
    FROM   exon e,static_golden_path sgp,exon_transcript et,transcript t
    WHERE  t.transcript_id = et.transcript_id
    AND    et.exon_id = e.exon_id 
    AND    sgp.raw_id = e.contig_id
    AND    sgp.chr_name = '$chr_name'
    AND    sgp.type = '$type'
    AND    e.contig_id in $idlist
    AND    sgp.chr_end >= $glob_start
    AND    sgp.chr_start <= $glob_end
    ORDER  BY t.gene_id, t.transcript_id, et.rank, e.sticky_rank";

  my $sth = $self->dbobj->prepare($query);
  $sth->execute();

  my ($exonid,$stickyrank,$phase,$rank,$transcriptid,$geneid,
      $start,$end,$strand);
  $sth->bind_columns(undef,\$exonid,\$stickyrank,\$phase,\$rank,
                     \$transcriptid,\$geneid,\$start,\$end,\$strand);
  
  my $current_transcript;
  my $current_gene;
  my $current_transcript_id='';
  my $current_gene_id='';
  my $previous_exon;

  my @out;
  my @trans;
  my $length = $glob_end - $glob_start;
  my %exons;

  $query = "SELECT type FROM gene WHERE gene.gene_id = ?"; 
  my $sth_gene = $self->dbobj->prepare($query);

  while( $sth->fetch ) {
    next if (($end > $length) || ($start < 1));

    if( $geneid ne $current_gene_id ) {
      # make a new gene
      $current_gene = Bio::EnsEMBL::Gene->new;
      $current_gene->dbID($geneid);
      $current_gene->adaptor($self->dbobj->get_GeneAdaptor);
      $sth_gene->execute($geneid);

      my ($type) = $sth_gene->fetchrow_array;

      $current_gene->type($type);
      push(@out,$current_gene);
      $current_gene_id = $geneid;
    }

    if( $transcriptid ne $current_transcript_id ) {
      # make a new transcript
      $current_transcript = Bio::EnsEMBL::WebTranscript->new();
      $current_transcript->dbID($current_transcript_id);
      $current_transcript->adaptor($self->dbobj->get_TranscriptAdaptor);

      $current_gene->add_Transcript($current_transcript);
      push(@trans,$current_transcript);
      if( $rank == 1 ) {
        $current_transcript->is_start_exon_in_context('dummy',1);
      } else {
        $current_transcript->is_start_exon_in_context('dummy',0);
      }

      $current_transcript_id = $transcriptid;
      $current_transcript->dbID($transcriptid);
      $current_transcript->adaptor( $self->dbobj->get_TranscriptAdaptor );
    }
    if( $stickyrank > 1 && defined $previous_exon && 
        $previous_exon->dbID() == $exonid ) { 
      # This is the same EXON so amend its length
      $previous_exon->end($end);
      next;            ## else we can't do anything else so create the sticky 
                       ## exon anyway! (but this will possibly only return 
                       ## part of the sticky exon)
    } 
        
    unless( $exons{$exonid} ) { 
      my $exon = Bio::EnsEMBL::Exon->new();
      $exon->start($start);
      $exon->end($end);
      $exon->strand($strand);
      $exon->dbID($exonid);
      $exon->adaptor( $self->dbobj->get_ExonAdaptor() );
      $exon->seqname($self->id);
      $exon->phase($phase);
      $exon->attach_seq($self);
      $previous_exon = $exon;
      $exons{$exonid} = $exon;
    }
    $current_transcript->add_Exon($exons{$exonid});
    $current_transcript->end_exon_rank($rank);
  }

  #
  # We need to make another quick trip to the database for each
  # transcript to discover whether we have all of it or not
  #

  my $sth_trans = $self->dbobj->prepare("select max(rank) from exon_transcript where transcript_id = ?");

  foreach my $trans ( @trans ) {
    $sth_trans->execute($trans->dbID);
    my ($rank) = $sth_trans->fetchrow_array();
    if( $rank == $trans->end_exon_rank) {
      $trans->is_end_exon_in_context('dummy',1);
    } else {
      $trans->is_end_exon_in_context('dummy',0);
    }
  }

  $self->{'_all_Genes_exononly'} = \@out;

  return @out;
}





=head2 get_all_VirtualGenes_startend

 Title   : get_all_VirtualGenes_startend
 Usage   :
 Function: return VirtualGenes lying on this virtual contig
 Example :
 Returns : 
 Args    :

=cut

sub get_all_VirtualGenes_startend_lite {
	my  $self = shift;
	return $self->dbobj->get_LiteAdaptor->fetch_virtualgenes_start_end(
        $self->_chr_name, 
        $self->_global_start, 
        $self->_global_end
    ); 
}

sub get_all_VirtualTranscripts_startend_lite {
    my  $self = shift;
    return $self->dbobj->get_LiteAdaptor->fetch_virtualtranscripts_start_end(
        $self->_chr_name,
        $self->_global_start,
        $self->_global_end,
        @_
    );
}

sub get_all_VirtualTranscripts_startend_lite_coding {
    my  $self = shift;
    return $self->dbobj->get_LiteAdaptor->fetch_virtualtranscripts_coding_start_end(
        $self->_chr_name,
        $self->_global_start,
        $self->_global_end,
        @_
    );
}

sub get_all_VirtualGenscans_startend_lite {
    my  $self = shift;
    return $self->dbobj->get_LiteAdaptor->fetch_virtualgenscans_start_end(
        $self->_chr_name,
        $self->_global_start,
        $self->_global_end
    );
}

sub get_all_EMBLGenes_startend_lite {
	my  $self = shift;
	return $self->dbobj->get_LiteAdaptor->fetch_EMBLgenes_start_end(
        $self->_chr_name, 
        $self->_global_start, 
        $self->_global_end
    ); 
}

sub get_all_SangerGenes_startend_lite {
      my  $self = shift;
      return $self->dbobj->get_LiteAdaptor->fetch_SangerGenes_start_end(
        $self->_chr_name,
        $self->_global_start,
        $self->_global_end
    );
}

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

    my $query = "
        SELECT STRAIGHT_JOIN t.gene_id
          , MIN(IF(sgp.raw_ori=1,(e.seq_start+sgp.chr_start-sgp.raw_start-$glob_start)
                  , (sgp.chr_start+sgp.raw_end-e.seq_end-$glob_start))) as start
          , MAX(IF(sgp.raw_ori=1,(e.seq_end+sgp.chr_start-sgp.raw_start-$glob_start)
                  , (sgp.chr_start+sgp.raw_end-e.seq_start-$glob_start))) as end
          , IF (sgp.raw_ori=1,e.strand,(-e.strand))
          , gsi.stable_id
          , g.type
        FROM static_golden_path sgp ,exon e,exon_transcript et,transcript t
          , gene_stable_id gsi
          , gene g
        WHERE sgp.raw_id=e.contig_id
          AND e.contig_id IN $idlist
          AND e.exon_id=et.exon_id
          AND t.transcript_id=et.transcript_id
          AND sgp.chr_end >= $glob_start
          AND sgp.chr_start <=$glob_end
          AND sgp.chr_name='$chr_name'
          AND sgp.type = '$type'
          AND gsi.gene_id = t.gene_id
          AND t.gene_id = g.gene_id
        GROUP BY t.gene_id
        ";
# note: without the e.contig in $idlist, the query results and the query
# plan are identical. Scrap it ? PL

    my $sth = $self->dbobj->prepare($query);
    $sth->execute;

    &eprof_end("virtualgene-sql-get");


    &eprof_start("virtualgene-build");

				# 
    my ($gene_id,$start,$end,$strand, $stable_id, $type);	# 
    $sth->bind_columns(undef,\$gene_id,\$start,\$end,\$strand,\$stable_id,\$type);

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
	$gene->dbID($gene_id);
	$gene->stable_id( $stable_id );
        $gene->type($type);

	&eprof_start("virtualgene-externaldb");

	
#Get the DBlinks for the given gene
    my $entryAdaptor = $self->dbobj->get_DBEntryAdaptor();
    $entryAdaptor->fetch_by_gene($gene);


	&eprof_end("virtualgene-externaldb");

	my $vg = Bio::EnsEMBL::VirtualGene->new(-gene => $gene,
						-contig => $self, 
						-start => $start, 
						-end => $end, 
						-strand => $strand
						);

	push @genes,$vg;
    }


    &eprof_end("virtualgene-build");

    $self->{'_virtualgenes_startend'} = \@genes;
    $self->_cached_virtualgenes_startend(1);

    return @genes;

} # get_all_VirtualGenes_startend

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

    my $query ="SELECT     STRAIGHT_JOIN t.transcript_id,
                           MIN(IF(sgp.raw_ori=1,(e.seq_start+sgp.chr_start-sgp.raw_start-$glob_start),
                                      (sgp.chr_start+sgp.raw_end-e.seq_end-$glob_start))) as start,
                           MAX(IF(sgp.raw_ori=1,(e.seq_end+sgp.chr_start-sgp.raw_start-$glob_start),
                                      (sgp.chr_start+sgp.raw_end-e.seq_start-$glob_start))) as end 
                FROM       static_golden_path sgp ,exon e,exon_transcript et,transcript t   
                WHERE      sgp.raw_id=e.contig_id 
                AND        e.contig_id in $idlist 
                AND        e.exon_id=et.exon_id 
                AND        t.transcript_id=et.transcript_id  
                AND        sgp.chr_end >= $glob_start   
                AND        sgp.chr_start <=$glob_end 
                AND        sgp.chr_name='$chr_name'
                AND        sgp.type = '$type'
                GROUP BY   t.transcript_id;";


    my $sth = $self->dbobj->prepare($query);
    $sth->execute;
    
    my ($transcript_id,$start,$end);
    $sth->bind_columns(undef,\$transcript_id,\$start,\$end);

    while ($sth->fetch){

	if (($end > $self->length)) {$end=$self->length;}
	if (($start < 1)) {$start=1;}

	my $gene=Bio::EnsEMBL::Gene->new();
	$gene->dbID($transcript_id);

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
  my @genes = $gene_obj->get_array_supporting('without', @gene_ids);
  
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
        &eprof_end("total-static-gene-get");
        return ();
    }

    my $query = "
        SELECT DISTINCT(t.gene_id)
        FROM exon e
          , exon_transcript et
          , transcript t
        WHERE e.contig_id IN $idlist
          AND e.exon_id = et.exon_id
          AND et.transcript_id = t.transcript_id
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
    #print "@gene_ids\n";
   if( scalar(@gene_ids) == 0 ) {
       &eprof_end("total-static-gene-get");
       return ();
   }

   &eprof_start("full-gene-get");

   my $gadp = $self->dbobj->get_GeneAdaptor();

    my @genes;
   foreach my $geneid ( @gene_ids ) {
       push(@genes,$gadp->fetch_by_dbID($geneid));
   }

   my %gene;

    #my @temp_one = keys(%gene);
    
   foreach my $gene ( @genes ) {
       $gene{$gene->dbID()}= $gene;
   }

   &eprof_end("full-gene-get");

   &eprof_start("gene-convert");
    my @temp = keys(%gene);
    
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

    $chr ||= $self->_chr_name();
    
    my $cache_name = "_chr_length_$chr";
    $self->{ $cache_name } = $self->dbobj->get_KaryotypeBandAdaptor()->fetch_chromosome_length($chr)
        unless defined $self->{ $cache_name };
            
    return( $self->{ $cache_name } );
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
          if (!defined($self->{_analysis_adp})) {
	   $self->{_analysis_adp} = new Bio::EnsEMBL::DBSQL::AnalysisAdaptor($self->dbobj);
          }
	  $analysis = $self->{_analysis_adp}->fetch_by_dbID($analysisid);
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
   return $self->SUPER::top_SeqFeatures(), @{$self->{'additional_seqf'}};
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
    
   unless( defined $self->{'_raw_contig_id_list'} ) {
        my $string = join ',',
                     map { $_->internal_id } $self->_vmap->get_all_RawContigs;
        $self->{'_raw_contig_id_list'} = $string ? "($string)" : "";
   }
   return $self->{'_raw_contig_id_list'};
		   
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
    $pass = '-' unless defined $pass;

    &Bio::EnsEMBL::Ext::ContigAcc::prepare_Ensembl_cache($host,$user,$pass,$dbname);

    my $glob_start = $self->_global_start;
    my $glob_end   = $self->_global_end;
    my $chr_name   = $self->_chr_name;
   
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

=head2 get_all_Clones

 Title   : 
 Usage   : $sc->get_all_Clones()

 Function: Produces an array of accession IDs of the golden clones spanning 
           the static contig
 Example :
 Returns : An array of accession IDs
 Args    : 

 Originator: James Smith (js5)
=cut

sub get_all_Clones {
   	my $self = shift;

   	my $sth = $self->dbobj->prepare(
   		"select distinct cl.id
		   from static_golden_path as sgp, contig as c, clone as cl
		  where c.clone = cl.internal_id and c.internal_id = sgp.raw_id and
	                sgp.type = ? and sgp.chr_start < ? and
			sgp.chr_end > ? and sgp.chr_name = ?"
	);
   	$sth->execute(
   		$self->dbobj->static_golden_path_type,
		$self->_global_end,
		$self->_global_start,
		$self->_chr_name
	);
    my $res = $sth->fetchall_arrayref;
    return map {$_->[0]} @$res;
}

=head2 _sgp_select

 Title   : 
 Usage   : $sth = $self->prepare( "SELECT".
           _sgp_select( "e.seq_start", "e.seq_end", "e.strand" ).

 Function: Produces a string as part of sql-select statement
           which does the start end strand selection
 Example :
 Returns : 
 Args    : start end strand as column name strings


=cut

sub _sgp_select {
  my $self = shift;
  my ($start, $end, $strand ) = @_;

  my $glob_start = $self->_global_start();
  
  my $str = "IF (sgp.raw_ori=1,
                  ($start+sgp.chr_start-sgp.raw_start-$glob_start+1),
                  (sgp.chr_start+sgp.raw_end-$end-$glob_start+1)),  
             IF (sgp.raw_ori=1,
                  ($end+sgp.chr_start-sgp.raw_start-$glob_start+1),
                  (sgp.chr_start+sgp.raw_end-$start-$glob_start+1)), 
             IF (sgp.raw_ori=1,$strand,(-$strand)) ";
  return $str;
}

sub has_AssemblyContigs {
  my ($self) = @_;

  my $ad = new Bio::EnsEMBL::DBSQL::AssemblyContigAdaptor($self->dbobj);

  if ($ad->has_AssemblyContigs($self->dbobj->static_golden_path_type)) {
    return 1;
  } else {
    return 0;
  }
}

sub each_AssemblyContig {
  my ($self) = @_;

  my $chr      = $self->_chr_name;
  my $chrstart = $self->_global_start;
  my $chrend   = $self->_global_end;

  my $ad = new Bio::EnsEMBL::DBSQL::AssemblyContigAdaptor($self->dbobj);
  my @contigs = $ad->fetch_by_chr_start_end($chr,$chrstart,$chrend);

  return @contigs;

}

sub get_all_MapFrags {
    my $self = shift;
    my $mapset = shift;
    return $self->dbobj->get_MapFragAdaptor->fetch_mapset_chr_start_end( 
        $mapset, $self->_chr_name, $self->_global_start, $self->_global_end
    );
}    

sub has_MapSet {
    my( $self, $mapset_name ) = @_;
    return $self->dbobj->get_MapFragAdaptor->has_mapset( $mapset_name );
}

=head2 get_Haplotypes_start_end

 Title   : get_Haplotypes_start_end
 Usage   : foreach my $hap ( $contig->get_Haplotypes_start_end ) 
 Function: Returns "shallow" haplotype objects suitable for drawing
 Example :
 Returns : 
 Args    :


=cut

sub get_Haplotypes_start_end {
  my ($self,$adaptor) = @_;
  
  my @haps = $adaptor->fetch_Haplotype_by_chr_start_end(
                         $self->_chr_name, 
                         $self->_global_start, 
                         $self->_global_end,
                        );
                        
    # reset haplotype object  global coordinated to be vc-based...
    # to keep the drawing code happy
    
  foreach my $h (@haps){
      $h->start($h->start() - $self->_global_start + 1 );
      $h->end($h->end() - $self->_global_start + 1 );
  }
  
  return (@haps);
}



1;
