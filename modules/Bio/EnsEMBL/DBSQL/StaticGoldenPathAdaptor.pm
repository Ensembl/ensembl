#
# Ensembl module for Bio::EnsEMBL::DBSQL::StaticGoldenPathAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::StaticGoldenPathAdaptor - Database adaptor for static golden path

=head1 SYNOPSIS

    # get a static golden path adaptor from the obj

    $adaptor = $db->get_StaticGoldenPathAdaptor();

    # these return sorted lists:

    @rawcontigs = $adaptor->fetch_RawContigs_by_fpc_name('ctg123');

    @rawcontigs = $adaptor->fetch_RawContigs_by_chr('chr2');

    #Create Virtual Contigs for fpc contigs or chromosomes

    $vc = $adaptor->fetch_VirtualContig_by_fpc_name('ctg123');

    $vc = $adaptor->fetch_VirtualContig_by_chr_name('chr2');


=head1 DESCRIPTION

Database adaptor for static golden path.  Affords access methods for retrieving virtual contigs.

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::StaticGoldenPathAdaptor;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root

use Bio::EnsEMBL::Root;
#use Bio::EnsEMBL::Virtual::StaticContig;
use Bio::EnsEMBL::RawContig;
use Bio::EnsEMBL::Utils::Eprof qw(eprof_start eprof_end);
use Bio::EnsEMBL::Slice;

@ISA = qw(Bio::EnsEMBL::Root);

# new() is written here 

sub new {
  my($class,@args) = @_;
  
  my $self = {};
  bless $self,$class;
  
  my ($db) = $self->_rearrange([qw(DB)],@args);

  if( !defined $db) {
      $self->throw("got no db. Aaaaah!");
  }

  $self->db($db);

# set stuff in self from @args
  return $self;
}


sub get_Gene_chr_MB {
    my ($self,$gene) = @_;

    my ($chr,$bp) = $self->get_Gene_chr_bp($gene);
    my $mbase = $bp/1000000;
    my $round = sprintf("%.1f",$mbase);   

    return ($chr,$round);
}



sub get_Gene_chr_bp {
    my ($self,$geneid) =  @_;
   
   my $type = $self->db->static_golden_path_type()
    or $self->throw("No assembly type defined");

   my $sth = $self->db->prepare("SELECT  
   if(a.contig_ori=1,(e.seq_start-a.contig_start+a.chr_start),
                    (a.chr_start+a.contig_end-e.seq_end)),
   if(a.contig_ori=1,(e.seq_end-a.contig_start+a.chr_start),
                    (a.chr_start+a.contig_end-e.seq_start)),
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
       $self->throw("Must have contig id to fetch VirtualContig of contig");
   }
   
   my $type = $self->db->static_golden_path_type()
    or $self->throw("No assembly type defined");

   my $sth = $self->db->prepare("SELECT  c.id,
                        a.chr_start,
                        a.chr_end,
                        a.chromosome_id 
                    FROM assembly a, contig c 
                    WHERE c.id = '$contigid' 
                    AND c.internal_id = a.contig_id 
                    AND a.type = '$type'"
                    );
   $sth->execute();
   my ($contig,$start,$end,$chr_name) = $sth->fetchrow_array;

   if( !defined $contig ) {
     $self->throw("Contig $contigid is not on the golden path of type $type");
   }

   return ($chr_name,$start,$end);
}


=head2 fetch_RawContigs_by_fpc_name

 Title   : fetch_RawContigs_by_fpc_name
 Usage   :
 Function: find all contigs belonging to the given FPC and lying on the
           Golden Path
 Example :
 Returns : returns an list of all rawContigs 
 Args    : the FPC id.


=cut

sub fetch_RawContigs_by_fpc_name {
   my ($self,$fpc) = @_;
   
   my $type = $self->db->static_golden_path_type()
    or $self->throw("No assembly type defined");

   # very annoying. DB obj wont make contigs by internalid. doh!
   my $sth = $self->db->prepare("SELECT  c.id 
                    FROM    assembly a, contig c 
                    WHERE c.internal_id = a.contig_id 
                    AND a.superctg_name = '$fpc' 
                    AND  a.type = '$type' 
                    ORDER BY a.superctg_start"
                    );
   $sth->execute;
   my @out;
   my $cid;

   while( ( my $cid = $sth->fetchrow_arrayref) ) {
       my $rc = $self->db->get_Contig($cid->[0]);
       push(@out,$rc);
   }
   if ($sth->rows == 0) {
       $self->throw("Could not find rawcontigs for fpc contig $fpc!");
   }
   return @out;
}

=head2 convert_chromosome_to_fpc
 
  Title   : convert_chromosome_to_fpc
  Usage   : ($fpcname,$start,$end) = $stadp->convert_chromosome_to_fpc('chr1',10000,10020)
  Function:
  Returns : 
  Args    :
 

=cut
 
sub convert_chromosome_to_fpc {
    my ($self,$chr,$start,$end) = @_;
 
    my $type = $self->db->static_golden_path_type()
     or $self->throw("No assembly type defined");
 
    my $sth = $self->db->prepare("
        SELECT superctg_name, chr_start
        FROM assembly
        WHERE chromosome_id = '$chr'
          AND superctg_start = 1
          AND chr_start <= $start
          AND type = '$type'
        ORDER BY chr_start DESC
        ");
    $sth->execute;
    my ($fpc,$startpos) = $sth->fetchrow_array;
 
    return ($fpc,$start-$startpos,$end-$startpos);
}

=head2 convert_fpc_to_chromosome

  Title   : convert_fpc_to_chromosome
  Usage   : ($chrname,$start,$end) = $stadp->convert_fpc_to_chromosome('ctg1234',10000,10020)
  Function:
  Returns : 
  Args    :
 

=cut
 
sub convert_fpc_to_chromosome {
    my ($self,$fpc,$start,$end) = @_;
 
    my $type = $self->db->static_golden_path_type()
     or $self->throw("No assembly type defined");
 
    my $sth = $self->db->prepare("
        SELECT chromosome_id, chr_start
        FROM assembly
        WHERE superctg_name = '$fpc'
            AND type = '$type'
        ORDER BY superctg_start LIMIT 1
        ");
    $sth->execute;
    my ($chr,$startpos) = $sth->fetchrow_array;
 
    if( !defined $chr ) {
        $self->throw("Couldn't find fpc contig $fpc in the database with $type golden path");
    }
    return ($chr,$start+$startpos,$end+$startpos) ;
}


=head2 convert_rawcontig_to_fpc

 Title   : convert_rawcontig_to_fpc
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub convert_rawcontig_to_fpc{
    my ($self,$rc,$start,$end,$strand) = @_;


    my $type = $self->db->static_golden_path_type()
     or $self->throw("No assembly type defined");
 
    my $sth = $self->db->prepare("
        SELECT a.superctg_name
          , a.superctg_start
          , a.contig_start
          , a.contig_ori
          , a.contig_end
        FROM assembly a
          , contig c
        WHERE c.internal_id = a.contig_id
          AND a.type = '$type'
          AND c.id = '$rc'
        ");
   $sth->execute;
   my ($fpc,$fpcstart,$rawstart,$rawori,$rawend) = $sth->fetchrow_array;
   
   if( $rawori == 1 ) {
       return ($fpc,$fpcstart+$start-$rawstart,$fpcstart+$end-$rawend,$strand);
   } else {
       return ($fpc,$fpcstart+($rawend - $end),$fpcstart+($rawend - $start),$strand*-1);
   }

}


=head2 fetch_RawContigs_by_chr_name

 Title   : fetch_RawContigs_by_chr_name
 Usage   :
 Function: get all the RawContigs on given chromosome belonging to the
           Golden Path
 Example :
 Returns : a list of RawContigs
 Args    : the chromosome name


=cut

sub fetch_RawContigs_by_chr_name{
   my ($self,$chr) = @_;

   my $type = $self->db->static_golden_path_type()
    or $self->throw("No assembly type defined");

   my $sth = $self->db->prepare("
        SELECT c.id
          , c.internal_id
          , c.dna
          , cl.id
          , cl.embl_version
          , a.chr_start
          , a.chr_end
          , a.contig_start
          , a.contig_end
          , a.contig_ori
          , c.offset
          , c.length
        FROM assembly a
          , contig c
          , clone cl
        WHERE cl.internal_id = c.clone
          AND c.internal_id = a.contig_id
          AND a.chromosome_id = '$chr'
          AND a.type = '$type'
        ");
#### fix this; should also return raw start/end, ori, offset and length
### see how done later on in this file
   $sth->execute;

   my @out;
   my $cid;
   while( ( my $array = $sth->fetchrow_arrayref) ) {

       my ($id,$internalid,$dna,$clone,$seq_version,$chr_start,$chr_end,$raw_start,$raw_end,$raw_ori,$offset,$contig_length) = @{$array};
       my $rc = Bio::EnsEMBL::DBSQL::RawContig->direct_new
       ( 
         -db    => $self->db,
         -id    => $id,
         -perlonlysequences => $self->db->perl_only_sequences,
         -internal_id => $internalid,
         -dna_id => $dna,
         -seq_version => $seq_version,
         -cloneid     => $clone,
             -chr_start   => $chr_start,
             -chr_end     => $chr_end,
             -raw_start   => $raw_start,
             -raw_end     => $raw_end,
             -raw_ori     => $raw_ori,
             -offset      => $offset,
             -contig_length => $contig_length
         );
       push(@out,$rc);
   }

   return @out;
}



=head2 fetch_RawContigs_by_chr_start_end

 Title   : fetch_RawContigs_by_chr_start_end
 Usage   :
 Function: return all RawContigs on given chromosome between start and
           end, on current Golden Path
 Example :
 Returns : list of RawContigs
 Args    : chromosome, start, end (in chromosome coordinates)


=cut

sub fetch_RawContigs_by_chr_start_end {
   my ($self,$chr,$start,$end) = @_;

   my $type = $self->db->static_golden_path_type()
    or $self->throw("No assembly type defined");
   
   # go for new go-faster method 
   # PL: is the query below correct? The 'NOT
   # <' looks odd, the join is different from e.g. the branch version
   # (e.g., look at 1.3.2.41). We'll keep it for now :-)
   #
   # JGRG: The NOT clauses are slightly counterintuitive,
   # but fetch every contig that overlaps start and end.
   # We search for every contig that doesn't begin after
   # our end (doesn't overlap) and doesn't end before
   # our start (doesn't overlap), and therefore get every
   # contig that overlaps.  This takes care of the condition
   # where our start and end lie within a contig.


   &eprof_start('VC: fetch_rc_get');
   my $sth = $self->db->prepare("
        SELECT c.id
          , c.internal_id
          , c.dna
          , cl.id
          , cl.embl_version
          , a.chr_start
          , a.chr_end
          , a.contig_start
          , a.contig_end
          , a.contig_ori
          , c.offset
          , c.length
        FROM assembly a
          , contig c
          , clone cl
        WHERE cl.internal_id = c.clone
          AND c.internal_id = a.contig_id
          AND a.chromosome_id = '$chr'
          AND a.type = '$type'
          AND NOT (a.chr_start > $end) 
          AND NOT (a.chr_end < $start) 
        ");
   $sth->execute;
   &eprof_end('VC: fetch_rc_get');

   my @out;
   my $cid;

   &eprof_start('VC: rc_build');

   while( ( my $array = $sth->fetchrow_arrayref) ) {

       my ($id,$internalid,$dna,$clone,$seq_version,$chr_start,$chr_end,$raw_start,$raw_end,$raw_ori,$offset,$contig_length) = @{$array};
       my $rc = Bio::EnsEMBL::DBSQL::RawContig->direct_new
       ( 
         -db    => $self->db,
         -id    => $id,
         -perlonlysequences => $self->db->perl_only_sequences,
         -internal_id => $internalid,
         -dna_id => $dna,
         -seq_version => $seq_version,
         -cloneid     => $clone,
             -chr_start   => $chr_start,
             -chr_end     => $chr_end,
             -raw_start   => $raw_start,
             -raw_end     => $raw_end,
             -raw_ori     => $raw_ori,
             -offset      => $offset,
             -contig_length => $contig_length
         );
       push(@out,$rc);
   }

   &eprof_end('VC: rc_build');

   return @out;
   
}


=head2 fetch_VirtualContig_by_chr_start_end

 Title   : fetch_VirtualContig_by_chr_start_end
 Usage   :
 Function: create a Virtual Contig based on a segment of a chromosome and
           start/end
 Example :
 Returns : A VirtualContig
 Args    : chromosome, start, end (in Chromosome coordinates)


=cut

sub fetch_VirtualContig_by_chr_start_end {
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


=head2 fetch_VirtualContig_of_clone

 Title   : fetch_VirtualContig_of_clone
 Usage   : $vc = $stadp->fetch_VirtualContig_of_clone('AC000012',1000);
 Function: Creates a virtual contig of the specified object.  If a context size is given, the vc is extended by that number of basepairs on either side of the clone.  Throws if the clone is not golden.
 Returns : Virtual Contig object 
 Args    : clone id, [context size in bp]


=cut

sub fetch_VirtualContig_of_clone{
   my ($self,$clone,$size) = @_;

   if( !defined $clone ) {
       $self->throw("Must have clone to fetch VirtualContig of clone");
   }
   if( !defined $size ) {$size=0;}

   my $type = $self->db->static_golden_path_type()
    or $self->throw("No assembly type defined");

   my $sth = $self->db->prepare("SELECT  c.id,
                        a.chr_start,
                        a.chr_end,
                        a.chromosome_id 
                    FROM    assembly a, 
                        contig c, 
                        clone  cl
                    WHERE c.clone = cl.internal_id
                    AND cl.id = '$clone'  
                    AND c.internal_id = ass.contig_id 
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
       $self->throw("Clone is not on the golden path. Cannot build VC");
   }
     
   my $vc = $self->fetch_VirtualContig_by_chr_start_end(    $chr_name,
                            $first_start-$size,
                            $end+$size
                            );
   $vc->db($self->db);
   return $vc;

}



=head2 fetch_VirtualContig_of_contig

 Title   : fetch_VirtualContig_of_contig
 Usage   : $vc = $stadp->fetch_VirtualContig_of_contig('AC000012.00001',1000);
 Function: Creates a virtual contig of the specified object.  If a context size is given, the vc is extended by that number of basepairs on either side of the contig.  Throws if the contig is not golden.
 Returns : Virtual Contig object 
 Args    : contig id, [context size in bp]


=cut

sub fetch_VirtualContig_of_contig{
   my ($self,$contigid,$size) = @_;

   if( !defined $size ) {$size=0;}

   my ($chr_name,$start,$end) = $self->get_chr_start_end_of_contig($contigid); 

   return $self->fetch_VirtualContig_by_chr_start_end(  $chr_name,
                            $start-$size,
                            $end+$size
                            );
  
}




=head2 fetch_VirtualContig_of_gene

 Title   : fetch_VirtualContig_of_gene
 Usage   : $vc = $stadp->fetch_VirtualContig_of_gene('ENSG00000012123',1000);
 Function: Creates a virtual contig of the specified object.  If a context size is given, the vc is extended by that number of basepairs on either side of the gene.  Throws if the gene is not golden.
 Returns : Virtual Contig object 
 Args    : gene id, [context size in bp]


=cut

sub fetch_VirtualContig_of_gene{
   my ($self,$geneid,$size) = @_;

   if( !defined $geneid ) {
       $self->throw("Must have gene id to fetch VirtualContig of gene");
   }
   if( !defined $size ) {$size=0;}

   my ($chr_name,$start,$end) = $self->get_Gene_chr_bp($geneid);

   if( !defined $start ) {
       my $type = $self->db->static_golden_path_type()
        or $self->throw("No assembly type defined");
       $self->throw("Gene is not on the golden path '$type'. Cannot build VC");
   }
     
   return $self->fetch_VirtualContig_by_chr_start_end(  $chr_name,
                            $start-$size,
                            $end+$size
                            );
}




=head2 fetch_VirtualContig_of_exon

 Title   : fetch_VirtualContig_of_exon
 Usage   : $vc = $stadp->fetch_VirtualContig_of_exon('ENSE00000648605',1000);
 Function: Creates a virtual contig of the specified object.  If a context size is given, the vc is extended by that number of basepairs on either side of the gene.  Throws if the object is not golden.
 Returns : Virtual Contig object 
 Args    : exon id, [context size in bp]

=cut

sub fetch_VirtualContig_of_exon{
   my ($self,$exonid,$size) = @_;


   $self->warn("Use of StaticGoldenPathAdaptor.fetch_VirtualContig_of_exon is deprecated.");

   return undef;
}


=head2 fetch_VirtualContig_of_feature

 Title   : fetch_VirtualContig_of_exon
 Usage   : $vc = $stadp->fetch_VirtualContig_of_feature('AC000001.1.1.2000',20,200,100);
 Function: Creates a virtual contig of the arbitrary feature.  If a context size is given, the vc is extended by that number of basepairs on either side of the gene.  Throws if the object is not golden.
 Returns : Virtual Contig object 
 Args    : contig id, seq_start, seq_end, [context size in bp]

=cut

sub fetch_VirtualContig_of_feature {
   my ($self,$contigid,$seq_start,$seq_end,$size) = @_;

   $self->warn("Use of StaticGoldenPathAdaptor.fetch_VirtualContig_of_feature is deprecated.");

   return undef;
}


=head2 get_location_of_feature

 Title   : get_location_of_feature
 Usage   : $vc = $stadp->get_location_of_feature('AC000001.1.1.2000',20,200);
 Function: Gets golden path co-ordinates of an arbitrary feature. Throws if contig is not golden... 
 Returns : array consisting of ( chromosome, orientation, start, end )
 Args    : contig id, seq_start, seq_end

=cut

sub get_location_of_feature {
   my ($self,$contigid,$seq_start,$seq_end) = @_;

   unless( defined $contigid && defined $seq_start && defined $seq_end ) {
       $self->throw("Must have feature details to get location of it");
   }

   my $type = $self->db->static_golden_path_type()
    or $self->throw("No assembly type defined");

   my $sth = $self->db->prepare("SELECT  
			if(a.contig_ori=1,
			    ($seq_start-a.contig_start+a.chr_start),
			    (a.chr_start+a.contig_end-$seq_end)),
   
			if(a.contig_ori=1,
			    ($seq_end-a.contig_start+a.chr_start),
			    (a.chr_start+a.contig_end-$seq_start)),

			    a.contig_ori,
			    a.chromosome_id
                    FROM    contig c, 
			    assembly a 
                    WHERE   c.id = '$contigid' 
			    AND c.internal_id = a.contig_id
                            AND a.type = '$type' 
                    ");
   $sth->execute();

   my ($start,$end,$raw_ori,$chr_name)=$sth->fetchrow_array;
   
   if( !defined $start ) {
       $self->throw("Contig $contigid is not on the current $type golden path.");
   }
     
   return ( $chr_name, $raw_ori, $start, $end );
}


=head2 fetch_VirtualContig_of_transcript

 Title   : fetch_VirtualContig_of_transcript
 Usage   : $vc = $stadp->fetch_VirtualContig_of_transcript('ENST00000012123',1000);
 Function: Creates a virtual contig of the specified object.  If a context size is given, the vc is extended by that number of basepairs on either side of the gene.  Throws if not golden.
 Returns : Virtual Contig object 
 Args    : transcript id, [context size in bp]


=cut

sub fetch_VirtualContig_of_transcript{
   my ($self,$transcriptid,$size) = @_;

   if( !defined $transcriptid ) {
       $self->throw("Must have gene id to fetch VirtualContig of transcript");
   }
   if( !defined $size ) {$size=0;}


   my $type = $self->db->static_golden_path_type()
    or $self->throw("No assembly type defined");

   my $sth = $self->db->prepare("SELECT  
   if(a.contig_ori=1,(e.seq_start-a.contig_start+a.chr_start),
                    (a.chr_start+a.contig_end-e.seq_end)),
   if(a.contig_ori=1,(e.seq_end-a.contig_start+a.chr_start),
                    (a.chr_start+a.contig_end-e.seq_start)),
     a.chromosome_id
  
                    FROM    exon e,
                        exon_transcript et,
                        assembly a,
                        transcript_stable_id tsi
                    WHERE tsi.stable_id = '$transcriptid'  
                    AND et.transcript_id = tsi.transcript_id
                    AND e.exon_id=et.exon_id 
                    AND a.contig_id=e.contig_id 
                    AND a.type = '$type' 
                    ");
   $sth->execute();

   my ($start,$end,$chr_name);
   my @start;
   while ( my @row=$sth->fetchrow_array){
      ($start,$end,$chr_name)=@row;
       push @start,$start;
       push @start,$end;
   }   
   
   my @start_sorted=sort { $a <=> $b } @start;

   $start=shift @start_sorted;
   $end=pop @start_sorted;

   if( !defined $start ) {
       $self->throw("Transcript is not on the golden path. Cannot build VC");
   }
     
   return $self->fetch_VirtualContig_by_chr_start_end(  $chr_name,
                            $start-$size,
                            $end+$size
                            );
   
}



=head2 fetch_VirtualContig_by_clone

 Title   : fetch_VirtualContig_by_clone
 Usage   : $vc = $stadp->fetch_VirtualContig_by_clone('AC000012',40000);
 Function: create a VirtualContig based on clone, and of a
           given length. The VC is centered around the start of the clone.
 Example :
 Returns : 
 Args    : clone name, size


=cut

sub fetch_VirtualContig_by_clone {
   my ($self,$clone,$size) = @_;

   if( !defined $size ) {
       $self->throw("Must have clone and size to fetch VirtualContig by clone");
   }

   my $type = $self->db->static_golden_path_type()
    or $self->throw("No assembly type defined");


   my $sth = $self->db->prepare("SELECT  c.id,
                        a.chr_start,
                        a.chromosome_id 
                    FROM assembly a,contig c,clone cl 
                    WHERE c.clone = cl.internal_id
                    AND cl.id = '$clone' 
                    AND c.internal_id = st.raw_id 
                    AND a.type = '$type' 
                    ORDER BY a.chr_start"
                    );
   $sth->execute();
   my ($contig,$start,$chr_name) = $sth->fetchrow_array;

   if( !defined $contig ) {
       $self->throw("Clone is not on the golden path. Cannot build VC");
   }


   my $halfsize = int($size/2);

   return $self->fetch_VirtualContig_by_chr_start_end(  $chr_name,
                            $start-$halfsize,
                            $start+$size-$halfsize
                            );
}




=head2 fetch_VirtualContig_by_contig

 Title   : fetch_VirtualContig_by_contig
 Usage   : $vc = $stadp->fetch_VirtualContig_by_clone('AC000012.00001',40000);
 Function: create a VirtualContig based on a RawContig, and of a
           given length. The VC is centered around the start of the clone.
 Example :
 Returns : 
 Args    : contigid (display_id, not internal one).

=cut

sub fetch_VirtualContig_by_contig {
   my ($self,$contigid,$size) = @_;

   if( !defined $size ) {
       $self->throw("Must have contig id and size to fetch VirtualContig by contig");
   }

   my $type = $self->db->static_golden_path_type()
    or $self->throw("No assembly type defined");

   my $sth = $self->db->prepare("SELECT  c.id,
                        a.chr_start,
                        a.chromosome_id 
                    FROM assembly a,contig c 
                    WHERE c.id = '$contigid' 
                    AND c.internal_id = a.contig_id 
                    AND a.type = '$type'"
                    );
   $sth->execute();
   my ($contig,$start,$chr_name) = $sth->fetchrow_array;

   if( !defined $contig ) {
     $self->throw("Contig $contigid is not on the golden path of type $type");
   }

   my $halfsize = int($size/2);
       return $self->fetch_VirtualContig_by_chr_start_end(  $chr_name,
                                $start-$halfsize,
                            $start+$size-$halfsize
                            );
}





=head2 fetch_VirtualContig_by_gene

 Title   : fetch_VirtualContig_by_gene
 Usage   : $vc = $stadp->fetch_VirtualContig_by_gene('ENSG00000012123',40000);
 Function: Creates a virtual contig of the specified size, centred around the given gene.
 Returns : Virtual Contig object 
 Args    : ensemblgene id, VC size in bp


=cut

sub fetch_VirtualContig_by_gene{
   my ($self,$geneid,$size) = @_;

   if( !defined $geneid ) {
       $self->throw("Must have gene id to fetch VirtualContig of gene");
   }
   if( !defined $size ) {$size=0;}

   my ($chr_name,$start) = $self->get_Gene_chr_bp($geneid);

   if( !defined $start ) {
       $self->throw("Gene is not on the golden path. Cannot build VC");
   }
     
   my $halfsize = int($size/2);

   return $self->fetch_VirtualContig_by_chr_start_end(  $chr_name,
                            $start-$halfsize,
                            $start+$size-$halfsize
                            );
}


=head2 fetch_VirtualContig_by_fpc_name

 Title   : fetch_VirtualContig_by_fpc_name
 Usage   :
 Function: create a VirtualContig representing a complete FPC contig
 Example :
 Returns : 
 Args    : the FPC contig id.


=cut

sub fetch_VirtualContig_by_fpc_name{
    my ($self,$fpc_name) = @_;

    my $type = $self->db->static_golden_path_type();

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

    &eprof_start('Slice: staticcontig build');

    eval {
      $slice = Bio::EnsEMBL::Slice->new($chr,$slice_start,$slice_end,$strand,$type);
    } ;
    if( $@ ) {
      $self->throw("Unable to build a slice using its fpc_name for for $chr, $slice_start,$slice_end\n\nUnderlying exception $@\n");
    }
    &eprof_end('Slice: staticcontig build');

    return $slice;
}


=head2 fetch_VirtualContig_list_sized

 Title   : fetch_VirtualContig_list_sized
 Usage   : @vclist = $stadaptor->fetch_VirtualContig_list_sized('ctg123',2000000,50000,4000000,100)
 Function: returns a list of virtual contigs from a FPC contig, split at gaps. The
           splitting happens as a greedy process:
              read as many contigs in until the first lenght threshold hits
              after this, split at the first gap length given
              If no gaps of this length are around, when the next length threshold is hit
              split at that gap.
 Returns : A list of VirtualContigs
 Args    : name,first lenght threshold, first gap size, second length threshold, second gap size


=cut

sub fetch_VirtualContig_list_sized {
   my ($self,$name,$length1,$gap1,$length2,$gap2) = @_;

   $self->warn("Use of StaticGoldenPathAdaptor.fetch_VirtualContig_list_sized is deprecated.");

   return undef;
}



=head2 fetch_VirtualContig_by_chr_name

 Title   : fetch_VirtualContig_by_chr_name
 Usage   :
 Function: create a VirtualContig representing the complete given chromosome
 Example :
 Returns : 
 Args    : chromosome name


=cut

sub fetch_VirtualContig_by_chr_name{
   my ($self,$name) = @_;

   my $vc = Bio::EnsEMBL::Virtual::StaticContig->new(1,1,-1,
                    $self->fetch_RawContigs_by_chr_name($name));
  
   $vc->db($self->db);
   $vc->_chr_name($name);
   return $vc; 
}


=head2 get_all_fpc_ids

 Title   : get_all_fpc_ids
 Usage   :
 Function:
 Returns : 
 Args    :


=cut


sub get_all_fpc_ids {
   my ($self,@args) = @_;

   my $type = $self->db->static_golden_path_type()
    or $self->throw("No assembly type defined");
   my $sth = $self->db->prepare("SELECT DISTINCT(superctg_name) 
                    FROM assembly 
                    WHERE type = '$type'"
                );
   $sth->execute();
   my @out;
   my $cid;
   while (my $rowhash = $sth->fetchrow_hashref){
       push (@out,$rowhash->{'fpcctg_name'});
   }
   if ($sth->rows == 0) {
       $self->throw("Could not find any fpc contigs in golden path $type!");
   }
   return @out;
}

sub get_chromosome_length {
    my ($self,$chrname) = @_;

    $self->throw("No chromosome name entered") unless defined($chrname);

    my $type = $self->db->static_golden_path_type()
     or $self->throw("No assembly type defined");
    
    my $sth = $self->db->prepare("
        SELECT MAX(chr_end)
        FROM assembly
        WHERE chromosome_id = '$chrname'
          AND type = '$type'
        ");

    $sth->execute;

    my ($len) = $sth->fetchrow;

    return $len;
}


=head2 db

 Title   : db
 Usage   : $obj->db($newval)
 Function: 
 Example : 
 Returns : value of db (i.e., the database handle)
 Args    : newvalue (optional)


=cut

sub db{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'db'} = $value;
    }
    return $obj->{'db'};

}


# sneaky

sub is_golden_static_contig {
    my ($self,$cid,$pos) = @_;
    my $type = $self->db->static_golden_path_type()
     or $self->throw("No assembly type defined");

    my $query = "
     SELECT c.id, a.contig_start, a.contig_end 
     FROM contig c, assembly a 
     WHERE c.id = '$cid' 
     AND a.contig_id = c.internal_id
     AND a.type = '$type'";

    my $sth = $self->db->prepare($query);
    $sth->execute;
    my $row = $sth->fetchrow_hashref;
    if ($row){
    if (defined($pos)) {
        if ($pos >= $row->{'raw_start'} && $pos <= $row->{'raw_end'}) {
        return 1;
        }
    } 
    else {
        return 1; 
    } 
    }

    return 0;
    
}

sub is_golden_static_clone {
    my ($self,$clone) = @_;
    my $type = $self->db->static_golden_path_type()
     or $self->throw("No assembly type defined");

    my $query = "   SELECT co.id 
            FROM contig co, clone cl, assembly a 
            WHERE cl.id = '$clone' 
            AND co.clone = cl.internal_id 
            AND a.contig_id = co.internal_id
            AND a.type = '$type'
        ";
       
    my $sth = $self->db->prepare($query);

    $sth->execute;

    return scalar($sth->fetchrow_array);
}
