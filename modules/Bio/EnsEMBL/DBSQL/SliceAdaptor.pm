
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
use vars qw(@ISA);
use strict;


# Object preamble - inherits from Bio::EnsEMBL::Root
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


@ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');


# new is inherited from BaseAdaptor


=head2 fetch_by_chr_start_end

 Title   : fetch_by_chr_start_end
 Usage   :
 Function: create a Slice based on a segment of a chromosome and
           start/end
 Example :
 Returns : A Slice
 Args    : chromosome, start, end (in Chromosome coordinates)


=cut

sub fetch_by_chr_start_end {
    my ($self,$chr,$start,$end) = @_;

    if( !defined $end ) {   # Why defined?  Is '0' a valid end?
        $self->throw("must provide chr, start and end");
    }

    if( $start > $end ) {
      $self->throw("start must be less than end: parameters $chr:$start:$end");
    }
    
    my $slice;
    my $type = $self->db->assembly_type();

    $slice = Bio::EnsEMBL::Slice->new(
          -chr_name      => $chr,
          -chr_start     => $start,
          -chr_end       => $end,
          -assembly_type => $type,
          -adaptor       => $self->db->get_SliceAdaptor()
	 );

    return $slice;
}



=head2 fetch_by_contig_name

 Title   : fetch_by_contig_name
 Usage   : $slice = $slice_adaptor->fetch_by_contig_name('AC000012.00001',1000);
 Function: Creates a slice of the specified slice adaptor object.  If a context size is given, the slice is extended by that number of basepairs on either side of the contig.  Throws if the contig is not golden.
 Returns : Slice object 
 Args    : contig id, [context size in bp]


=cut

sub fetch_by_contig_name {
   my ($self,$contigid,$size) = @_;

   if( !defined $size ) {$size=0;}

   my ($chr_name,$start,$end) = $self->_get_chr_start_end_of_contig($contigid);

   $start -= $size;
   $end += $size;

   if($start < 1) {
     $start  = 1;
   }

   return $self->fetch_by_chr_start_end($chr_name, $start, $end);
 }


=head2 fetch_by_fpc_name

 Title   : fetch_by_fpc_name
 Usage   :
 Function: create a Slice representing a complete FPC contig
 Example :
 Returns : 
 Args    : the FPC contig id.


=cut

sub fetch_by_fpc_name {
    my ($self,$fpc_name) = @_;

    my $type = $self->db->assembly_type();

    my $sth = $self->db->prepare("
        SELECT chromosome_id, superctg_ori, MIN(chr_start), MAX(chr_end)
        FROM assembly
        WHERE superctg_name = '$fpc_name'
        AND type = '$type'
        GROUP by superctg_name
        ");

    $sth->execute;

    my ($chr, $strand, $slice_start, $slice_end) = $sth->fetchrow_array;

    my $slice;

    $slice = new Bio::EnsEMBL::Slice($chr,$slice_start,$slice_end,
                                     $strand,$type);

    return $slice;
}



=head2 fetch_by_clone_accession

 Title   : fetch_by_clone_accession
 Usage   : $slice = $slice_adaptor->fetch_by_clone_accession('AC000012',1000);
 Function: Creates a Slice of the specified object.  If a context size is given, the Slice is extended by that number of basepairs on either side of the clone.  Throws if the clone is not golden.
 Returns : Slice object 
 Args    : clone id, [context size in bp]


=cut

sub fetch_by_clone_accession{
   my ($self,$clone,$size) = @_;

   if( !defined $clone ) {
     $self->throw("Must have clone to fetch Slice of clone");
   }
   if( !defined $size ) {$size=0;}

   my $type = $self->db->assembly_type()
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
     
   $first_start -= $size;
   $end += $size;

   if($first_start < 1) {
     $first_start = 1;
   }

   my $slice = $self->fetch_by_chr_start_end($chr_name, $first_start, $end);
   return $slice;
}



=head2 fetch_by_transcript_stable_id

 Title   : fetch_by_transcript_stable_id
 Usage   : $slice = $slice_adaptor->fetch_by_transcript_stable_id(
                                       'ENST00000302930',1000);
 Function: Creates a slice of the specified object.  If a context
           size is given, the slice is extended by that number of
	   basepairs on either side of the transcript.  Throws if
	   the transcript is not golden.
 Returns : Slice object 
 Args    : transcript stable ID, [context size in bp]


=cut

sub fetch_by_transcript_stable_id{
  my ($self,$transcriptid,$size) = @_;

  # Just get the dbID, then fetch slice by that
  my $ta = $self->db->get_TranscriptAdaptor;
  my $transcript_obj = $ta->fetch_by_stable_id($transcriptid);
  my $dbID = $transcript_obj->dbID;
  
  return $self->fetch_by_transcript_id($dbID, $size);
}


=head2 fetch_by_transcript_id

 Title   : fetch_by_transcript_id
 Usage   : $slice = $slice_adaptor->fetch_by_transcript_id(24,1000);
 Function: Creates a slice of the specified object.  If a context
           size is given, the slice is extended by that number of
	   basepairs on either side of the transcript.  Throws if
	   the transcript is not golden.
 Returns : Slice object 
 Args    : transcript dbID, [context size in bp]

=cut

sub fetch_by_transcript_id {
  my ($self,$transcriptid,$size) = @_;
     if( !defined $transcriptid ) {
       $self->throw("Must have transcriptid id to fetch Slice of transcript");
   }
   if( !defined $size ) {$size=0;}
   my $emptyslice = Bio::EnsEMBL::Slice->new( '-empty'   => 1,
                                              '-adaptor' => $self,
					      '-ASSEMBLY_TYPE' =>
					      $self->db->assembly_type);
   my $ta = $self->db->get_TranscriptAdaptor;
   my $transcript_obj = $ta->fetch_by_dbID($transcriptid);

   my %exon_transforms;
   for my $exon ( $transcript_obj->get_all_Exons() ) {
     my $newExon = $exon->transform( $emptyslice );
     $exon_transforms{ $exon } = $newExon;
   }
   $transcript_obj->transform( \%exon_transforms );


  my $start = $transcript_obj->start() - $size;
  my $end = $transcript_obj->end() + $size;

  if($start < 1) {
    $start = 1;
  }

  my $slice = $self->fetch_by_chr_start_end($emptyslice->chr_name,
					    $start, $end);
  return $slice;
}


=head2 fetch_by_gene_stable_id

 Title   : fetch_by_gene_stable_id
 Usage   : $slice = $slice_adaptor->fetch_by_gene_stable_id('ENSG00000012123',1000);
 Function: Creates a slice of the specified object.  If a context size is given, the slice is extended by that number of basepairs on either side of the gene.  Throws if the gene is not golden.
 Returns : Slice object 
 Args    : gene id, [context size in bp]


=cut

sub fetch_by_gene_stable_id{
   my ($self,$geneid,$size) = @_;

   if( !defined $geneid ) {
       $self->throw("Must have gene id to fetch Slice of gene");
   }
   if( !defined $size ) {$size=0;}

   my ($chr_name,$start,$end) = $self->_get_chr_start_end_of_gene($geneid);

   if( !defined $start ) {
     my $type = $self->db->assembly_type()
       or $self->throw("No assembly type defined");
     $self->throw("Gene is not on the golden path '$type'. Cannot build Slice.");
   }
     
   $start -= $size;
   $end += $size;
   
   if($start < 1) {
     $start = 1;
   }

   return $self->fetch_by_chr_start_end($chr_name, $start, $end);
}


=head2 fetch_by_chr_name

 Title   : fetch_by_chr_name
 Usage   : $slice = $slice_adaptor->fetch_by_chr_name('20');
 Function: Creates a slice of an entire chromosome 
 Returns : Slice object 
 Args    : chromosome name


=cut

sub fetch_by_chr_name{
   my ($self,$chr_name) = @_;

   unless( $chr_name ) {
       $self->throw("Chromosome name argument required");
   }

   my $chr_start = 1;
   
   #set the end of the slice to the end of the chromosome
   my $ca = $self->db()->get_ChromosomeAdaptor();
   my $chromosome = $ca->fetch_by_chrname($chr_name);
   my $chr_end = $chromosome->length();

   return $self->fetch_by_chr_start_end($chr_name, $chr_start, $chr_end);
}



=head2 _get_chr_start_end_of_contig

 Title   : _get_chr_start_end_of_contig
 Usage   :
 Function: returns the chromosome name, absolute start and absolute end of the 
           specified contig
 Returns : returns chr,start,end
 Args    : contig id

=cut

sub _get_chr_start_end_of_contig {
    my ($self,$contigid) = @_;

   if( !defined $contigid ) {
       $self->throw("Must have contig id to fetch Slice of contig");
   }
   
   my $type = $self->db->assembly_type()
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

=head2 _get_chr_start_end_of_gene

 Title   : get_Gene_chr_bp
 Usage   : 
 Function: 
 Returns :  
 Args    :


=cut


sub _get_chr_start_end_of_gene {
  my ($self,$geneid) =  @_;
  
  my $type = $self->db->assembly_type()
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
