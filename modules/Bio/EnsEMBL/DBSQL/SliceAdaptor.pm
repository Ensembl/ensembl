
#
# Ensembl module for Bio::EnsEMBL::DBSQL::SliceAdaptor
#
# Cared for by Ewan Birney <ensembl-dev@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::SliceAdaptor - Adaptors for slices

=head1 SYNOPSIS
  



=head1 DESCRIPTION

Factory for getting out slices of assemblies. WebSlice is the highly
accelerated version for the web site.

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email ensembl-dev@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::Utils::Eprof qw(eprof_start eprof_end);
use vars qw(@ISA);
use strict;


# Object preamble - inherits from Bio::EnsEMBL::Root
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Slice;


@ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');


# new is inherieted from BaseAdaptor

=head2 new_slice

 Title   : new_slice
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub new_slice{
    my ($self,$chr,$start,$end,$strand,$type) = @_;


    my $slice = Bio::EnsEMBL::Slice->new( -chr_name  => $chr,
					  -chr_start => $start,
					  -chr_end   => $end,
					  -strand    => $strand,
					  -assembly_type      => $type,
					-adaptor => $self);

    return $slice;
}


=head2 new_web_slice

 Title   : new_web_slice
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub new_web_slice{
    my ($self,$chr,$start,$end,$strand,$type) = @_;
    
    die "Not implemented new slice yet";
    
}


sub fetch_all_repeat_features{
  my($self, $slice, $logic_name) = @_;

  if(!$slice){
    $self->throw("can't fetch all repeat features if con't have a slice to fetch them for\n");
  }

  my @repeats = $self->db->get_RepeatFeatureAdaptor->fetch_by_Slice($slice, $logic_name);

  return @repeats;

}

sub fetch_all_simple_features{
  my($self, $slice, $logic_name) = @_;

  if(!$slice){
    $self->throw("can't fetch all simple features if con't have a slice to fetch them for\n");
  }

  my @simple = $self->db->get_SimpleFeatureAdaptor->fetch_by_Slice($slice, $logic_name);

  return @simple;

}

sub fetch_all_prediction_transcripts{
  my($self, $slice, $logic_name) = @_;

  if(!$slice){
    $self->throw("can't fetch all simple features if con't have a slice to fetch them for\n");
  }

  my @prediction = $self->db->get_PredictionTranscriptAdaptor->fetch_by_Slice($slice, $logic_name);

  return @prediction;

}


sub fetch_all_similarity_features{
  my($self, $slice, $logic_name) = @_;

  if(!$slice){
    $self->throw("can't fetch all simple features if con't have a slice to fetch them for\n");
  }
  
  my @out;

  my @dnaalign = $self->db->get_DnaAlignFeatureAdaptor->fetch_by_Slice($slice, $logic_name);
  my @pepalign = $self->db->get_ProteinAlignFeatureAdaptor->fetch_by_Slice($slice, $logic_name);

  push(@out, @dnaalign);
  push(@out, @pepalign);

  return @out;
}


sub fetch_all_similarity_features_above_score{
  my($self, $slice, $score, $logic_name) = @_;

  if(!$slice){
    $self->throw("can't fetch all simple features if con't have a slice to fetch them for\n");
  }
  if(!$score){
    $self->throw("need score even if it 0\n");
  }
  my @out;

  my @dnaalign = $self->db->get_DnaAlignFeatureAdaptor->fetch_by_Slice_and_score($slice, $score, $logic_name);
  my @pepalign = $self->db->get_ProteinAlignFeatureAdaptor->fetch_by_Slice_and_score($slice, $score, $logic_name);

  push(@out, @dnaalign);
  push(@out, @pepalign);

  return @out;
}


sub fetch_all_similarity_features_above_pid{
  my($self, $slice, $pid, $logic_name) = @_;

  if(!$slice){
    $self->throw("can't fetch all simple features if con't have a slice to fetch them for\n");
  }
  if(!$pid){
    $self->throw("need percent_id even if it 0\n");
  }
  my @out;

  my @dnaalign = $self->db->get_DnaAlignFeatureAdaptor->fetch_by_Slice_and_pid($slice, $pid, $logic_name);
  my @pepalign = $self->db->get_ProteinAlignFeatureAdaptor->fetch_by_Slice_and_pid($slice, $pid, $logic_name);

  push(@out, @dnaalign);
  push(@out, @pepalign);

  return @out;
}



=head2 get_chr_start_end_of_contig

 Title   : get_chr_start_end_of_contig
 Usage   :
 Function: returns the chromosome name, absolute start and absolute end of the 
           specified contig
 Returns : returns chr,start,end
 Args    : contig id

=cut

sub get_chr_start_end_of_contig {
    my ($self,$contigid) = @_;

   if( !defined $contigid ) {
       $self->throw("Must have contig id to fetch Slice of contig");
   }
   
   my $type = $self->db->static_golden_path_type()
    or $self->throw("No assembly type defined");

   my $sth = $self->db->prepare("SELECT  c.name,
                        a.chr_start,
                        a.chr_end,
                        a.chromosome_id 
                    FROM assembly a, contig c 
                    WHERE c.name = '$contigid' 
                    AND c.contig_id = a.contig_id 
                    AND a.type = '$type'"
                    );
   $sth->execute();
   my ($contig,$start,$end,$chr_name) = $sth->fetchrow_array;

   if( !defined $contig ) {
     $self->throw("Contig $contigid is not on the golden path of type $type");
   }

   return ($chr_name,$start,$end);
}


=head2 fetch_Slice_by_chr_start_end

 Title   : fetch_Slice_by_chr_start_end
 Usage   :
 Function: create a Slice based on a segment of a chromosome and
           start/end
 Example :
 Returns : A Slice
 Args    : chromosome, start, end (in Chromosome coordinates)


=cut

sub fetch_Slice_by_chr_start_end {
    my ($self,$chr,$start,$end) = @_;

    if( !defined $end ) {   # Why defined?  Is '0' a valid end?
        $self->throw("must provide chr, start and end");
    }

    if( $start > $end ) {
        $self->throw("start must be less than end: parameters $chr:$start:$end");
    }

    my $slice;

    &eprof_start('Slice: staticcontig build');

    my $type = $self->db->static_golden_path_type();

    eval {
      $slice = Bio::EnsEMBL::Slice->new(
          -chr_name      => $chr,
          -chr_start     => $start,
          -chr_end       => $end,
          -assembly_type => $type,
          -adaptor       => $self->db->get_SliceAdaptor
      );
    } ;
    if( $@ ) {
      $self->throw("Unable to build a slice for $chr, $start,$end\n\nUnderlying exception $@\n");
    }
    &eprof_end('Slice: staticcontig build');

    return $slice;
}



=head2 fetch_Slice_by_contig

 Title   : fetch_Slice_by_contig
 Usage   : $slice = $slice_adaptor->fetch_Slice_by_contig('AC000012.00001',1000);
 Function: Creates a slice of the specified slice adaptor object.  If a context size is given, the slice is extended by that number of basepairs on either side of the contig.  Throws if the contig is not golden.
 Returns : Slice object 
 Args    : contig id, [context size in bp]


=cut

sub fetch_Slice_by_contig{
   my ($self,$contigid,$size) = @_;

   if( !defined $size ) {$size=0;}

   my ($chr_name,$start,$end) = $self->get_chr_start_end_of_contig($contigid); 

   return $self->fetch_Slice_by_chr_start_end(  $chr_name,
                            $start-$size,
                            $end+$size
                            );
  
}


=head2 fetch_Slice_by_clone

 Title   : fetch_Slice_by_clone
 Usage   : $slice = $slice_adaptor->fetch_Slice_by_clone('AC000012',1000);
 Function: Creates a Slice of the specified object.  If a context size is given, the Slice is extended by that number of basepairs on either side of the clone.  Throws if the clone is not golden.
 Returns : Slice object 
 Args    : clone id, [context size in bp]


=cut

sub fetch_Slice_by_clone{
   my ($self,$clone,$size) = @_;

   if( !defined $clone ) {
       $self->throw("Must have clone to fetch Slice of clone");
   }
   if( !defined $size ) {$size=0;}

   my $type = $self->db->static_golden_path_type()
    or $self->throw("No assembly type defined");

   my $sth = $self->db->prepare("SELECT  c.name,
                        a.chr_start,
                        a.chr_end,
                        a.chromosome_id 
                    FROM    assembly a, 
                        contig c, 
                        clone  cl
                    WHERE c.clone_id = cl.clone_id
                    AND cl.name = '$clone'  
                    AND c.contig_id = a.contig_id 
                    AND a.type = '$type' 
                    ORDER BY a.chr_start"
                    );
   $sth->execute();
 
   my ($contig,$start,$end,$chr_name); 
   my $counter; 
   my $first_start;
   while ( my @row=$sth->fetchrow_array){
       $counter++;
       ($contig,$start,$end,$chr_name)=@row;
       if ($counter==1){$first_start=$start;}      
   }

   if( !defined $contig ) {
       $self->throw("Clone is not on the golden path. Cannot build Slice");
   }
     
   my $slice = $self->fetch_Slice_by_chr_start_end(    $chr_name,
                            $first_start-$size,
                            $end+$size
                            );
   $slice->adaptor->db($self->db);
   return $slice;

}


=head2 get_Gene_chr_bp

 Title   : get_Gene_chr_bp
 Usage   : 
 Function: 
 Returns :  
 Args    :


=cut


sub get_Gene_chr_bp {
    my ($self,$geneid) =  @_;
   
   my $type = $self->db->static_golden_path_type()
    or $self->throw("No assembly type defined");

   my $sth = $self->db->prepare("SELECT  
   if(a.contig_ori=1,(e.contig_start-a.contig_start+a.chr_start),
                    (a.chr_start+a.contig_end-e.contig_end)),
   if(a.contig_ori=1,(e.contig_end-a.contig_start+a.chr_start),
                    (a.chr_start+a.contig_end-e.contig_start)),
     a.chromosome_id
  
                    FROM    exon e,
                        transcript tr,
                        exon_transcript et,
                        assembly a,
                        gene_stable_id gsi
                    WHERE e.exon_id=et.exon_id 
                    AND et.transcript_id =tr.transcript_id 
                    AND a.contig_id=e.contig_id 
                    AND a.type = '$type' 
                    AND tr.gene_id = gsi.gene_id
                    AND gsi.stable_id = '$geneid';" 
                    );
   $sth->execute();

   my ($start,$end,$chr);
   my @start;
   while ( my @row=$sth->fetchrow_array){
      ($start,$end,$chr)=@row;
       push @start,$start;
       push @start,$end;
   }   
   
   my @start_sorted=sort { $a <=> $b } @start;

   $start=shift @start_sorted;
   $end=pop @start_sorted;

   return ($chr,$start,$end); 
        
}


=head2 fetch_Slice_by_gene

 Title   : fetch_Slice_by_gene
 Usage   : $slice = $slice_adaptor->fetch_Slice_by_gene('ENSG00000012123',1000);
 Function: Creates a slice of the specified object.  If a context size is given, the slice is extended by that number of basepairs on either side of the gene.  Throws if the gene is not golden.
 Returns : Slice object 
 Args    : gene id, [context size in bp]


=cut

sub fetch_Slice_by_gene{
   my ($self,$geneid,$size) = @_;

   if( !defined $geneid ) {
       $self->throw("Must have gene id to fetch Slice of gene");
   }
   if( !defined $size ) {$size=0;}

   my ($chr_name,$start,$end) = $self->get_Gene_chr_bp($geneid);

   if( !defined $start ) {
       my $type = $self->adaptor->db->static_golden_path_type()
        or $self->throw("No assembly type defined");
       $self->throw("Gene is not on the golden path '$type'. Cannot build Slice.");
   }
     
   return $self->fetch_Slice_by_chr_start_end(  $chr_name,
                            $start-$size,
                            $end+$size
                            );
}


