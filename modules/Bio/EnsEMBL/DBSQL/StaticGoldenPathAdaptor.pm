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

    $adaptor = $dbobj->get_StaticGoldenPathAdaptor();

    # these return sorted lists:

    @rawcontigs = $adaptor->fetch_RawContigs_by_fpc_name('ctg123');

    @rawcontigs = $adaptor->fetch_RawContigs_by_chr('chr2');

    #Create Virtual Contigs for fpc contigs or chromosomes

    $vc = $adaptor->VirtualContig_by_fpc_name('ctg123');

    $vc = $adaptor->VirtualContig_by_chr('chr2');

    # can throw an exception: Not on Same Chromosome
    @rawcontigs = $adaptor->fetch_RawContigs_between_RawContigs($start_rc,$end_rc);


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

# Object preamble - inherits from Bio::Root::RootI

use Bio::Root::RootI;
use Bio::EnsEMBL::Virtual::StaticContig;

@ISA = qw(Bio::Root::RootI);

# new() is written here 

sub new {
  my($class,@args) = @_;
  
  my $self = {};
  bless $self,$class;
  
  my ($dbobj) = $self->_rearrange([qw( DBOBJ)],@args);

  if( !defined $dbobj) {
      $self->throw("got no dbobj. Aaaaah!");
  }

  $self->dbobj($dbobj);

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
    my ($self,$gene) =  @_;

    my $query = "
     SELECT STRAIGHT_JOIN p.chr_name, p.chr_start 
     FROM transcript tr, translation tl, exon e, static_golden_path p 
     WHERE tr.gene = '$gene' 
       AND tl.id = tr.translation 
       AND tl.start_exon = e.id 
       AND e.contig = p.raw_id";

    my $sth = $self->dbobj->prepare($query);

    $sth->execute();

    my ($chr,$mbase) = $sth->fetchrow_array;


    return ($chr,$mbase); 
        
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
   
   my $type = $self->dbobj->static_golden_path_type();

   # very annoying. DB obj wont make contigs by internalid. doh!
   my $sth = $self->dbobj->prepare("SELECT  c.id 
				    FROM    static_golden_path st,
					    contig c 
				    WHERE c.internal_id = st.raw_id 
				    AND st.fpcctg_name = '$fpc' 
				    AND  st.type = '$type' 
				    ORDER BY st.fpcctg_start"
				    );
   $sth->execute;
   my @out;
   my $cid;

   while( ( my $cid = $sth->fetchrow_arrayref) ) {
       my $rc = $self->dbobj->get_Contig($cid->[0]);
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
 
sub convert_chromosome_to_fpc{
    my ($self,$chr,$start,$end) = @_;
 
    my $type = $self->dbobj->static_golden_path_type();
 
    my $sth = $self->dbobj->prepare("SELECT fpcctg_name,
					    chr_start 
				    FROM static_golden_path 
				    WHERE chr_name = '$chr' 
				    AND	fpcctg_start = 1 
				    AND	chr_start <= $start 
				    ORDER BY chr_start DESC"
				    );
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
 
    my $type = $self->dbobj->static_golden_path_type();
 
    my $sth = $self->dbobj->prepare("SELECT chr_name,
					    chr_start 
				    FROM static_golden_path 
				    WHERE fpcctg_name = '$fpc' 
                                    ORDER by fpcctg_start LIMIT 1");
#				    AND fpcctg_start = 1" 
#				    );
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


    my $type = $self->dbobj->static_golden_path_type();
 
    my $sth = $self->dbobj->prepare("SELECT st.fpcctg_name,
					    st.fpcctg_start,
                                            st.raw_start,
                                            st.raw_ori,
                                            st.raw_end
				    FROM static_golden_path st,contig c 
				    WHERE c.id = '$rc' AND c.internal_id = st.raw_id"
				    );
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

   my $type = $self->dbobj->static_golden_path_type();
   
   # very annoying. DB obj wont make contigs by internalid. doh!
   my $sth = $self->dbobj->prepare("SELECT  c.id 
				    FROM    static_golden_path st,
					    contig c 
				    WHERE c.internal_id = st.raw_id 
				    AND st.chr_name = '$chr' 
				    AND  st.type = '$type' 
				    ORDER BY st.fpcctg_start"
				    );
   $sth->execute;
   my @out;
   my $cid;
   while( ( my $cid = $sth->fetchrow_arrayref) ) {
       my $rc = $self->dbobj->get_Contig($cid->[0]);
       push(@out,$rc);
   }
   if ($sth->rows == 0) {
       $self->throw("Could not find rawcontigs for chromosome $chr!");
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


   my $type = $self->dbobj->static_golden_path_type();
   
   $self->throw("I need a golden path type") unless ($type);
   
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
   my $sth = $self->dbobj->prepare("
        SELECT c.id
          , c.internal_id
          , c.dna
          , cl.id
          , cl.embl_version
        FROM static_golden_path st
          , contig c
          , clone cl
        WHERE cl.internal_id = c.clone
          AND c.internal_id = st.raw_id
          AND st.chr_name = '$chr'
          AND st.type = '$type'
          AND NOT (st.chr_start > $end) 
          AND NOT (st.chr_end < $start) 
        ");
   $sth->execute;

   my @out;
   my $cid;
   while( ( my $array = $sth->fetchrow_arrayref) ) {

       my ($id,$internalid,$dna,$clone,$seq_version) = @{$array};
       my $rc = Bio::EnsEMBL::DBSQL::RawContig->direct_new
	   ( 
	     -dbobj => $self->dbobj,
	     -id    => $id,
	     -perlonlysequences => $self->dbobj->perl_only_sequences,
	     -contig_overlap_source      => $self->dbobj->contig_overlap_source(),
	     -overlap_distance_cutoff    => $self->dbobj->overlap_distance_cutoff(),
	     -internal_id => $internalid,
	     -dna_id => $dna,
	     -seq_version => $seq_version,
	     -cloneid => $clone
	     );
       push(@out,$rc);
   }

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


    my @rc = $self->fetch_RawContigs_by_chr_start_end($chr,$start,$end);
	# Ewan's note - We can have zero rawcontigs - when there is a gap
        #or $self->throw("Got zero rawcontigs");


    my $vc;

    eval {
      $vc = Bio::EnsEMBL::Virtual::StaticContig->new($start,1,$end,@rc);
    } ;
    if( $@ ) {
      $self->throw("Unable to build a virtual contig at $chr, $start,$end\n\nUnderlying exception $@\n");
    }

    $vc->_chr_name($chr);
    $vc->dbobj($self->dbobj);
    return $vc;
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

   my $type = $self->dbobj->static_golden_path_type();

   my $sth = $self->dbobj->prepare("SELECT  c.id,
   					    st.chr_start,
					    st.chr_end,
					    st.chr_name 
				    FROM    static_golden_path st, 
					    contig c 
				    WHERE c.clone = '$clone' 
                                    AND c.internal_id = st.raw_id 
				    AND st.type = '$type' 
                                    ORDER BY st.fpcctg_start"
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
     
   my $vc = $self->fetch_VirtualContig_by_chr_start_end(	$chr_name,
   							$first_start-$size,
							$end+$size
							);
   $vc->dbobj($self->dbobj);
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

   if( !defined $contigid ) {
       $self->throw("Must have contig id to fetch VirtualContig of contig");
   }
   
   if( !defined $size ) {$size=0;}

   my $type = $self->dbobj->static_golden_path_type();

   my $sth = $self->dbobj->prepare("SELECT  c.id,
   					    st.chr_start,
					    st.chr_end,
					    st.chr_name 
                                    FROM static_golden_path st,contig c 
				    WHERE c.id = '$contigid' 
                                    AND c.internal_id = st.raw_id 
				    AND st.type = '$type'"
		   		    );
   $sth->execute();
   my ($contig,$start,$end,$chr_name) = $sth->fetchrow_array;

   if( !defined $contig ) {
     $self->throw("Contig $contigid is not on the golden path of type $type");
   }
   
   return $self->fetch_VirtualContig_by_chr_start_end(	$chr_name,
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


   my $type = $self->dbobj->static_golden_path_type();

   my $sth = $self->dbobj->prepare("SELECT  
   if(sgp.raw_ori=1,(e.seq_start-sgp.raw_start+sgp.chr_start),
                    (sgp.chr_start+sgp.raw_end-e.seq_end)),
   if(sgp.raw_ori=1,(e.seq_end-sgp.raw_start+sgp.chr_start),
                    (sgp.chr_start+sgp.raw_end-e.seq_start)),
     sgp.chr_name
  
				    FROM    exon e,
					    transcript tr,
					    exon_transcript et,
					    static_golden_path sgp 
				    WHERE e.id=et.exon 
				    AND et.transcript=tr.id 
				    AND sgp.raw_id=e.contig 
				    AND sgp.type = '$type' 
				    AND tr.gene = '$geneid';" 
		   		    );
   $sth->execute();

   my ($start,$end,$chr_name);
   my @start;
   while ( my @row=$sth->fetchrow_array){
      ($start,$end,$chr_name)=@row;
       print STDERR "Got $start-$end \n";
       push @start,$start;
       push @start,$end;
   }   
   
   my @start_sorted=sort { $a <=> $b } @start;

   $start=shift @start_sorted;
   $end=pop @start_sorted;

   if( !defined $start ) {
       $self->throw("Gene is not on the golden path. Cannot build VC");
   }
     
   return $self->fetch_VirtualContig_by_chr_start_end(	$chr_name,
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

   my $type = $self->dbobj->static_golden_path_type();


   my $sth = $self->dbobj->prepare("SELECT  c.id,
   					    st.chr_start,
					    st.chr_name 
				    FROM static_golden_path st,contig c,clone cl 
				    WHERE c.clone = cl.internal_id
                                    AND cl.id = '$clone' 
				    AND c.internal_id = st.raw_id 
				    AND st.type = '$type' 
				    ORDER BY st.fpcctg_start"
				    );
   $sth->execute();
   my ($contig,$start,$chr_name) = $sth->fetchrow_array;

   if( !defined $contig ) {
       $self->throw("Clone is not on the golden path. Cannot build VC");
   }


   my $halfsize = int($size/2);

   return $self->fetch_VirtualContig_by_chr_start_end(	$chr_name,
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

   my $type = $self->dbobj->static_golden_path_type();

   my $sth = $self->dbobj->prepare("SELECT  c.id,
   					    st.chr_start,
					    st.chr_name 
				    FROM static_golden_path st,contig c 
				    WHERE c.id = '$contigid' 
				    AND c.internal_id = st.raw_id 
				    AND st.type = '$type'"
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


   my $type = $self->dbobj->static_golden_path_type();

   my $sth = $self->dbobj->prepare("SELECT  STRAIGHT_JOIN (e.seq_start+sgp.chr_start),
					    sgp.chr_name 
				    FROM    transcript tr, 
					    exon_transcript et,
                                            exon e,
					    static_golden_path sgp 
				    WHERE e.id=et.exon 
                                    AND et.transcript=tr.id 
				    AND sgp.raw_id=e.contig 
				    AND tr.gene = '$geneid';" 
		   		    );
   $sth->execute();


   my ($start,$chr_name); 
   my @start;
   while ( my @row=$sth->fetchrow_array){
       ($start,$chr_name)=@row;

       push @start,$start;
   }
   
   my @start_sorted=sort @start;

   $start=pop @start_sorted;

   if( !defined $start ) {
       $self->throw("Gene is not on the golden path. Cannot build VC");
   }
     
   my $halfsize = int($size/2);

   return $self->fetch_VirtualContig_by_chr_start_end(	$chr_name,
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
   my ($self,$name) = @_;
   
   my @fpc = $self->fetch_RawContigs_by_fpc_name($name);
   my $start = $fpc[0];
   my $vc = Bio::EnsEMBL::Virtual::StaticContig->new(	$start->chr_start,
							1,
							-1,
							@fpc
						    );
 
   $vc->dbobj($self->dbobj);
   $vc->id($name);
   return $vc;
}

# depracated
=head2 fetch_VirtualContig_by_fpc_name_slice

 Title   : fetch_VirtualContig_by_fpc_name_slice
 Usage   : do not use; depracated. Use a construct with
           fetch_VirtualContig_list_sized() instead.

 Function: bit bizarre: start and end (in fpc coords) indicate which
           RawContigs to use, then construct a VC consisting of the _full_
           extent of these RCs. (As a result, its length is not simply
           end-start+1)
 Example :
 Returns : a Virtual contig, consisting of all the overlap of th
 Args    : fpc contig id, and start end in fpc coordinates?

=cut

sub fetch_VirtualContig_by_fpc_name_slice {
   my ($self,$name,$start,$end) = @_;
   
   $self->warn("Usage of StaticGoldenPathAdaptor.fetch_VirtualContig_by_fpc_name_slice is depracated. Use a construct with fetch_VirtualContig_list_sized() instead");
   
   if( !defined $end ) {
       $self->throw("must have start end to fetch by slice");
   }

   my @fpc = $self->fetch_RawContigs_by_fpc_name($name);
   my @finalfpc;

   foreach my $fpc ( @fpc ) {
       if( $fpc->fpc_contig_start >= $start && $fpc->fpc_contig_end <= $end ) {
	   push(@finalfpc,$fpc);
       }
   }
   if( scalar @finalfpc == 0 ) {
       $self->throw("No complete raw contigs between $start and $end");
   }

   $start = $finalfpc[0];
   my $vc = Bio::EnsEMBL::Virtual::StaticContig->new(	$start->chr_start,
						    $start->fpc_contig_start,
						    -1,
						    @finalfpc
						    );
   $vc->id("$name-$start-$end");
   $vc->dbobj($self->dbobj);
   return $vc;
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

   if( !defined $gap2 ) {
       $self->throw("Must fetch Virtual Contigs in sized lists");
   }
   my @fpc = $self->fetch_RawContigs_by_fpc_name($name);
   my $chr;
   if ($#fpc >= 0) { $chr = $fpc[0]->chromosome; }

   my @finalfpc;
   my @vclist;

   my $current_start = 1;
   my $prev = shift @fpc;
   push(@finalfpc,$prev);
   foreach my $fpc ( @fpc ) {
       $fpc->dbobj($self->dbobj);

       if( ( ($fpc->fpc_contig_end - $current_start+1) > $length1 && ($fpc->fpc_contig_start - $prev->fpc_contig_end -1) >= $gap1) ||
	   ( ($fpc->fpc_contig_end -$current_start+1) > $length2 && ($fpc->fpc_contig_start - $prev->fpc_contig_end -1) >= $gap2) ) {
	   # build new vc and reset stuff

	   my $start = $finalfpc[0];

	   my $vc = Bio::EnsEMBL::Virtual::StaticContig->new($start->chr_start,$start->fpc_contig_start,-1,@finalfpc);
	   $vc->id($name);
           $vc->dbobj($self->dbobj);
	   $vc->_chr_name($chr);

	   push(@vclist,$vc);
	   
	   $prev = $fpc;
	   $current_start = $prev->fpc_contig_start;
	   @finalfpc = ();
	   push(@finalfpc,$prev);
 	   $prev->dbobj($self->dbobj);
       } else {
	   push(@finalfpc,$fpc);
	   $prev = $fpc;
	   $fpc->dbobj($self->dbobj);
       }
   }
   # last contig

   my $start = $finalfpc[0];
   my $vc = Bio::EnsEMBL::Virtual::StaticContig->new($start->chr_start,$start->fpc_contig_start,-1,@finalfpc);
   $vc->dbobj($self->dbobj);
   $vc->_chr_name($chr);

   push(@vclist,$vc);

   return @vclist;
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
  
   $vc->dbobj($self->dbobj);
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

   my $type = $self->dbobj->static_golden_path_type();
   my $sth = $self->dbobj->prepare("SELECT DISTINCT(fpcctg_name) 
				    FROM static_golden_path 
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

    my $sth = $self->dbobj->prepare("SELECT max(chr_end)
                                     FROM   static_golden_path
                                     WHERE  chr_name = \'$chrname\'");

    $sth->execute;

    my $rowhash = $sth->fetchrow_hashref;

    if (defined($rowhash)) {
	my $len = $rowhash->{'max(chr_end)'};
	return $len;
    }
}


=head2 dbobj

 Title   : dbobj
 Usage   : $obj->dbobj($newval)
 Function: 
 Example : 
 Returns : value of dbobj (i.e., the database handle)
 Args    : newvalue (optional)


=cut

sub dbobj{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'dbobj'} = $value;
    }
    return $obj->{'dbobj'};

}


# sneaky

sub is_golden_static_contig {
    my ($self,$cid,$pos) = @_;

    my $query = "
     SELECT c.id,p.raw_start,p.raw_end 
     FROM contig c, static_golden_path p 
     WHERE c.id = '$cid' AND p.raw_id = c.internal_id";

    my $sth = $self->dbobj->prepare($query);
    $sth->execute;
#    my @contigs;
     foreach my $row ($sth->fetchrow_hashref) {
         my $contig = $row->{'id'};
         my $start  = $row->{'raw_start'};
         my $end    = $row->{'raw_end'};
         if (defined($pos)) {
             if ($pos >= $start && $pos <= $end) {
                 return 1;
                # push(@contigs,$contig);
             }
         } else {
             return 1; 
             # push(@contigs,$contig);
         } 
     }
     return 0;
#     return scalar(@contigs);
}

sub is_golden_static_clone {
    my ($self,$clone) = @_;

    my $query = "
    SELECT co.id 
    FROM contig co, clone cl, static_golden_path p 
    WHERE cl.id = '$clone' 
      AND co.clone = cl.internal_id 
      AND p.raw_id = co.internal_id";
   
    my $sth = $self->dbobj->prepare($query);

    $sth->execute;

    return scalar($sth->fetchrow_array);
}
