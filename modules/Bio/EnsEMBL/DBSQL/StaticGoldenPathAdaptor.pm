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

Describe the object here

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

Describe contact details here

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
   my $sth = $self->dbobj->prepare("select c.id from static_golden_path st,contig c where c.internal_id = st.raw_id AND st.fpcctg_name = '$fpc' AND  st.type = '$type' ORDER BY st.fpcctg_start");
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
   my $sth = $self->dbobj->prepare("select c.id from static_golden_path st,contig c where c.internal_id = st.raw_id AND st.chr_name = '$chr' AND  st.type = '$type' ORDER BY st.fpcctg_start");
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
   
   # very annoying. DB obj wont make contigs by internalid. doh!
   my $sth = $self->dbobj->prepare("select c.id from static_golden_path st,contig c where c.internal_id = st.raw_id AND st.chr_name = '$chr' AND  st.type = '$type' AND st.chr_end > $start AND st.chr_start < $end ORDER BY st.fpcctg_start");

   $sth->execute;
   my @out;
   my $cid;
   while( ( my $cid = $sth->fetchrow_arrayref) ) {
       my $rc = $self->dbobj->get_Contig($cid->[0]);
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

   if( !defined $end ) {
       $self->throw("must provide chr, start and end");
   }

   if( $start > $end ) {
       $self->throw("start must be less than end");
   }

   
   my @rc = $self->fetch_RawContigs_by_chr_start_end($chr,$start,$end);

   
   my $vc = Bio::EnsEMBL::Virtual::StaticContig->new($start,1,$end,@rc);

   $vc->_chr_name($chr);
   return $vc;
}

=head2 fetch_VirtualContig_by_clone

 Title   : fetch_VirtualContig_by_clone
 Usage   : $vc = $stadp->fetch_VirtualContig_by_clone('AC000012',40000);
 Function: create a VirtualContig based on clone, and of a
           given length. The thing is centered around the start of the clone.
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

   my $sth = $self->dbobj->prepare("select c.id,st.chr_start,st.chr_name from static_golden_path st,contig c,clone cl where cl.id = '$clone' AND cl.internal_id = c.clone AND c.internal_id = st.raw_id AND st.type = '$type' ORDER BY st.fpcctg_start");
   $sth->execute();
   my ($contig,$start,$chr_name) = $sth->fetchrow_array;

   if( !defined $contig ) {
       $self->throw("Clone is not on the golden path. Cannot build VC");
   }


   my $halfsize = int($size/2);
   if( $start > $size/2 ) {       
       return $self->fetch_VirtualContig_by_chr_start_end($chr_name,$start-$halfsize,$start+$size-$halfsize);
   } else {
       return $self->fetch_VirtualContig_by_chr_start_end($chr_name,1,$size);
   }
}

=head2 fetch_VirtualContig_by_contig

 Title   : fetch_VirtualContig_by_contig
 Usage   : $vc = $stadp->fetch_VirtualContig_by_clone('AC000012.00001',40000);
 Function: as fetch_VirtualContig_by_clone, but based on a RawContig. 
 Example :
 Returns : 
 Args    : contigid (display_id, not internal one).

=cut

sub fetch_VirtualContig_by_contig {
   my ($self,$contigid,$size) = @_;

   if( !defined $size ) {
       $self->throw("Must have clone and size to fetch VirtualContig by clone");
   }

   my $type = $self->dbobj->static_golden_path_type();


   # PL: could use shortcut, since Virtual<->Raw is a one-to-one thing?
   my $sth = $self->dbobj->prepare("select c.id,st.chr_start,st.chr_name from static_golden_path st,contig c where c.id = '$contigid' AND c.internal_id = st.raw_id AND st.type = '$type'");
   $sth->execute();
   my ($contig,$start,$chr_name) = $sth->fetchrow_array;

   my $halfsize = int($size/2);
   if( $start > $size/2 ) {       
       return $self->fetch_VirtualContig_by_chr_start_end($chr_name,$start-$halfsize,$start+$size-$halfsize);
   } else {
       return $self->fetch_VirtualContig_by_chr_start_end($chr_name,1,$size);
   }
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
   my $vc = Bio::EnsEMBL::Virtual::StaticContig->new($start->chr_start,1,-1,@fpc);
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
   my $vc = Bio::EnsEMBL::Virtual::StaticContig->new($start->chr_start,$start->fpc_contig_start,-1,@finalfpc);
   $vc->id("$name-$start-$end");
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
 Example :
 Returns : A list of VirtualContigs
 Args    : name,first lenght threshold, first gap size, second length threshold, second gap size


=cut

sub fetch_VirtualContig_list_sized {
   my ($self,$name,$length1,$gap1,$length2,$gap2) = @_;

   if( !defined $gap2 ) {
       $self->throw("Must fetch Virtual Contigs in sized lists");
   }
   my @fpc = $self->fetch_RawContigs_by_fpc_name($name);

   my @finalfpc;
   my @vclist;

   my $current_start = 1;
   my $prev = shift @fpc;
   push(@finalfpc,$prev);
   foreach my $fpc ( @fpc ) {
       #print STDERR "Looking at ",$fpc->id," ",$fpc->fpc_contig_start," ",$fpc->fpc_contig_end," ",($fpc->fpc_contig_start - $prev->fpc_contig_end -1),"\n";
       if( ( ($fpc->fpc_contig_end - $current_start+1) > $length1 && ($fpc->fpc_contig_start - $prev->fpc_contig_end -1) >= $gap1) ||
	   ( ($fpc->fpc_contig_end -$current_start+1) > $length2 && ($fpc->fpc_contig_start - $prev->fpc_contig_end -1) >= $gap2) ) {
	   # build new vc and reset stuff

	   my $start = $finalfpc[0];
	   #print STDERR "Building with ",$start->id," ",$start->fpc_contig_start,"\n";

	   my $vc = Bio::EnsEMBL::Virtual::StaticContig->new($start->chr_start,$start->fpc_contig_start,-1,@finalfpc);
	   $vc->id($name);
	   push(@vclist,$vc);
	   
	   $prev = $fpc;
	   $current_start = $prev->fpc_contig_start;
	   @finalfpc = ();
	   push(@finalfpc,$prev);
       } else {
	   push(@finalfpc,$fpc);
	   $prev = $fpc;
       }
   }
   # last contig

   my $start = $finalfpc[0];
   my $vc = Bio::EnsEMBL::Virtual::StaticContig->new($start->chr_start,$start->fpc_contig_start,-1,@finalfpc);
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

   return Bio::EnsEMBL::Virtual::StaticContig->new(1,1,-1,$self->fetch_RawContigs_by_chr_name($name));


}



=head2 get_all_fpc_ids

 Title   : get_all_fpc_ids
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_fpc_ids {
   my ($self,@args) = @_;

   my $type = $self->dbobj->static_golden_path_type();
   my $sth = $self->dbobj->prepare("select distinct(fpcctg_name) from static_golden_path where type = '$type'");
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
