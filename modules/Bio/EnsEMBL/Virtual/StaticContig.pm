
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
   
# print STDERR "Static Contig created from ",scalar(@contigs), " contigs.\n";
    

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


=head2 get_all_SimilarityFeatures_above_score

 Title   : get_all_SimilarityFeatures_above_score
 Usage   : foreach my $sf ( $contig->get_all_SimilarityFeatures_above_score(analysis_type, score) ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SimilarityFeatures_above_score{
    my ($self, $analysis_type, $score,$bp) = @_;
    
    $self->throw("Must supply analysis_type parameter") unless $analysis_type;
    $self->throw("Must supply score parameter") unless $score;
    
    
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
                        WHERE  f.score > $score
                        AND    f.analysis = a.id 
                        AND    sgp.raw_id = f.contig
                        AND    f.contig in $idlist
		        AND    a.db = '$analysis_type'  
                        AND    sgp.chr_end >= $glob_start 
		        AND    sgp.chr_start <=$glob_end 
		        AND    sgp.chr_name='$chr_name' 
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

    unless ($idlist){
	return ();
    }

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
                     AND    sgp.chr_end >= $glob_start 
                     AND    sgp.chr_start <=$glob_end
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
	    print(STDERR "Repeat feature does not have a hid. bad news....");
	}
	push(@features,$out);
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

   
   my @web;
   my @std;

   my @features;
   my @contig_features;

   if( scalar($self->_vmap->get_all_RawContigs) == 0) {
       return();
   }

   foreach my $extf ( $self->dbobj->_each_ExternalFeatureFactory ) {
       if( $extf->isa('Bio::EnsEMBL::DB::WebExternalFeatureFactoryI') ) {
	   push(@web,$extf);
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
	   foreach my $feature ( $extf->get_Ensembl_SeqFeatures_clone_web($glob,@clones) ) {
	      
	       my $clone = $feature->seqname;
	       $clone =~ s/\.\d+$//g;
	       push(@{$featureh{$clone}},$feature);
	   }
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
	       if( $extf->can('get_Ensembl_SeqFeatures_contig') ) {
		   foreach my $sf ($extf->get_Ensembl_SeqFeatures_contig($contig->internal_id,$contig->seq_version,1,$contig->length,$contig->id)) {
			$sf->seqname($contig->id);
			push(@contig_features,$sf);
		   }
	       }
	       if( $extf->can('get_Ensembl_SeqFeatures_clone') ) {
       
		   foreach my $sf ( $extf->get_Ensembl_SeqFeatures_clone($contig->cloneid,$contig->seq_version,$contig->embl_offset,$contig->embl_offset+$contig->length()) ) {
		       
		       my $start = $sf->start - $contig->embl_offset+1;
		       my $end   = $sf->end   - $contig->embl_offset+1;
		       $sf->start($start);
		       $sf->end($end);
		       $sf->seqname($contig->id);
		       push(@contig_features,$sf);
		   }
	       }
	   }
       }
   }
		    

   # ok. Now @contig_features are in contig coordinates. Map up.
   
   # this is the simplest way. We can do this more efficiently if need be

   my @final;
   foreach my $f ( @contig_features ) {
       if( defined $self->_convert_seqfeature_to_vc_coords($f) ) {
	   push(@final,$f);
       }
   }

   return @final;

}




sub get_landmark_MarkerFeatures {

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

eval {
    require Bio::EnsEMBL::Map::MarkerFeature;


    my $statement= "   SELECT 
                       IF     (sgp.raw_ori=1,(f.seq_start+sgp.chr_start-sgp.raw_start-$glob_start),
                              (sgp.chr_start+sgp.raw_end-f.seq_end-$glob_start)),                                        
                       IF     (sgp.raw_ori=1,(f.seq_end+sgp.chr_start-sgp.raw_start-$glob_start),
                              (sgp.chr_start+sgp.raw_end-f.seq_start-$glob_start)), 
                              f.score, 
                       IF     (sgp.raw_ori=1,f.strand,(-f.strand)), 
                              f.name, f.hstart, f.hend, 
                              f.hid, f.analysis, c.name 
                       FROM   feature f,
                              contig_landmarkMarker c,
                              static_golden_path sgp 
                       WHERE  f.contig = c.contig
                       AND    f.hid=c.marker  
                       AND    f.contig in $idlist 
                       AND    sgp.raw_id=f.contig 
                       AND    sgp.chr_end >= $glob_start 
                       AND    sgp.chr_start <=$glob_end 
                       AND    sgp.chr_name='$chr_name'";
    
        $statement =~ s/\s+/ /g;
        #print STDERR "Doing Query ... $statement\n";

	my $sth = $self->dbobj->prepare($statement);
	$sth->execute;
	
	my ($start, $end, $score, $strand, $hstart, 
	    $name, $hend, $hid, $analysisid,$synonym);
	
	my $analysis;
	my %analhash;
	
	$sth->bind_columns
	    ( undef, \$start, \$end, \$score, \$strand, \$name, 
	      \$hstart, \$hend, \$hid, \$analysisid,\$synonym);
	
      while( $sth->fetch ) {

	    #clipping
	    if (($end > $length) || ($start < 1)) {
		next;
	    }
	    
	    my @args=($start,$end,$score,$strand,$name,$hstart,$hend,$hid,
		      $analysisid,$synonym);
	    
	    
	    my $out=$self->_create_Marker_features(@args);
	    if (defined $out){
		push (@markers,$out);
	    } 
	}
};


if($@){$self->warn("Problems retrieving map data\nMost likely not connected to maps db\n$@\n");}

return @markers;

}



=head2 next_landmark_Marker

 Title   : next_landmark_Marker
 Usage   : $obj->next_landmark_Marker
 Function: retrieves next marker  
 Returns : marker feature
 Args    : golden path position, chromosome, Mb limit (optional)


=cut



sub next_landmark_Marker
{

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
                          AND       sgp.chr_name='$chr_name' 
                          AND       f.analysis = a.id 
                          AND       a.db='mapprimer'
                          AND       sgp.chr_start>$start 
                          AND       sgp.chr_start <$end 
                          AND       sgp.chr_start+f.seq_start-sgp.raw_start>$start  
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
my ($self)=shift;

my $gene;
my @genes;

my $glob_start=$self->_global_start;
my $glob_end=$self->_global_end;
my $chr_name=$self->_chr_name;
my $idlist  = $self->_raw_contig_id_list();
    

unless ($idlist){
    return ();
}

$self->throw ("I need a chromosome name") unless defined $chr_name;
$self->throw ("I need a chromosome end") unless defined $glob_end;
$self->throw ("I need a chromosome start") unless defined $glob_start;

my $query ="SELECT     STRAIGHT_JOIN t.gene,
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
            GROUP BY   t.gene;";


    my $sth = $self->dbobj->prepare($query);
    $sth->execute;
    
    my ($gene_id,$start,$end);
    $sth->bind_columns(undef,\$gene_id,\$start,\$end);

    while ($sth->fetch){

	if (($end > $self->length)) {$end=$self->length;}
	if (($start < 1)) {$start=1;}

	$gene=Bio::EnsEMBL::Gene->new();
	$gene->id($gene_id);

	my $query = "select external_db,external_id from genedblink where gene_id = '$gene_id'";
	my $sth = $self->dbobj->prepare($query);
	my $res = $sth ->execute();
	while( (my $hash = $sth->fetchrow_hashref()) ) {
	    my $dblink = Bio::Annotation::DBLink->new();
	    $dblink->database($hash->{'external_db'});
	    $dblink->primary_id($hash->{'external_id'});
	    $gene->add_DBLink($dblink);
	}

	my $genestr=1;
	my $vg = Bio::EnsEMBL::VirtualGene->new(-gene => $gene,
						-contig => $self, 
						-start => $start, 
						-end => $end, 
						-strand => $genestr
						);
	push @genes,$vg;
    }

    return @genes;

}



=head2 fetch_karyotype_band

 Title   : fetch_karyotype_band
 Usage   : $label = $self->fetch_karyotype_band
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_karyotype_band {
   my ($self,@args) = @_;

   my $kadp = $self->dbobj->get_KaryotypeAdaptor();
   my $band = $kadp->get_band_label_by_position($self->_chr_name,$self->_global_start + ($self->length/2));

   return $band 
}




=head2 karyotype_band

 Title   : karyotype_band
 Usage   : $label = $self->karyotype_band
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub karyotype_band {
   my ($self,@args) = @_;

   my $kadp = $self->dbobj->get_KaryotypeAdaptor();
   my $band = $kadp->get_band_label_by_position($self->_chr_name,$self->_global_start + ($self->length/2));

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

=head2 species

 Title   : species
 Usage   : $obj->species($newval)
 Function: 
 Example : 
 Returns : value of species
 Args    : newvalue (optional)


=cut

sub species{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'species'} = $value;
    }
    return $obj->{'species'};

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

=head2 add_date

 Title   : add_date
 Usage   : $self->add_domment($ref)
 Function: adds a date
 Example :
 Returns : 
 Args    :


=cut

sub add_date {
   my ($self) = shift;
   foreach my $dt ( @_ ) {
       push(@{$self->{'date'}},$dt);
   }
}

=head2 each_date

 Title   : each_date
 Usage   : foreach $dt ( $self->each_date() )
 Function: gets an array of dates
 Example :
 Returns : 
 Args    :


=cut

sub each_date {
   my ($self) = @_;
   return @{$self->{'date'}}; 
}

=head2 annotation

 Title   : annotation
 Usage   : $obj->annotation($newval)
 Function: 
 Example : 
 Returns : value of annotation
 Args    : newvalue (optional)


=cut

sub annotation{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'annotation'} = $value;
    }
    return $obj->{'annotation'};

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





1;
