#
# EnsEMBL module for Bio::EnsEMBL::DSQL::SequenceAdaptor
#
# Cared for by Arne Stabenau <stabenau@ebi.ac.uk>
#
# Copyright EMBL/EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::SequenceAdaptor - produce sequence strings from locations

=head1 SYNOPSIS

$seq_adptr = $database_adaptor->get_SequenceAdaptor();
$dna = $seq_adptr->fetch_by_RawContig_start_end_strand($contig, 1, 1000, -1);

=head1 DESCRIPTION

An adaptor for the retrieval of sequences of DNA from the database

=head1 CONTACT

Arne Stabenau - stabenau@ebi.ac.uk
Elia Stupka - elia@fugu-sg.org

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...
package Bio::EnsEMBL::DBSQL::SequenceAdaptor;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object
use Bio::EnsEMBL::DBSQL::BaseAdaptor;


@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 fetch_by_RawContig_start_end_strand

  Arg [1]    : Bio::EnsEMBL::RawContig $contig
  Arg [2]    : int $start
  Arg [3]    : int $end
               a -1 means until the end
  Arg [4]    : int $strand
               -1, 1 are possible values
  Arg [5]    : (optional) force read from compressed DNA table
  Example    : $dna = $seq_adp->fetch_by_RawContig_start_end_strand($contig, 1,
								    1000, -1);
  Description: retrieves the dna string from the database from the 
               given RawContig. 
  Returntype : string 
  Exceptions : thrown if start < 1
  Caller     : Bio::EnsEMBL::RawContig::seq(), RawContig::subseq()

=cut

sub fetch_by_RawContig_start_end_strand {
  my ( $self, $contig, $start, $end, $strand, $compressed ) = @_;
  my $sth;
  
  $compressed ||= 0;

  if( $start < 1 ) {
    $self->throw( "Wrong parameters" );
  }

  unless($contig && ref $contig && $contig->isa("Bio::EnsEMBL::RawContig")) {
    $self->throw("contig arg must be contig not a [$contig]");
  } 

  {

    my $query;

    if($compressed==1 || $self->is_compressed){

      # query for compressed mode, uses hex and a different range
      # (compression implemented by th, 14/7/03)
      # needs to getback dna_id as will have to check for Ns in dnan table
      my $start4=int(($start-1)/4)+1;
      
      if( $end == -1 ) { 
	$query = "SELECT c.length, HEX( SUBSTRING( d.sequence, $start4 ) ), 
                         d.n_line
                    FROM dnac d, contig c 
                   WHERE d.dna_id = c.dna_id 
                     AND c.contig_id = ?";
      
      } else {
	my $end4=int(($end-1)/4)+1;
	my $length4 = $end4 - $start4 + 1;
	if( $length4 < 1 ) {
	  $self->throw( "Wrong parameters" );
	}
      
	$query = "SELECT c.length, HEX( SUBSTRING( d.sequence, $start4, $length4 ) ), 
                         d.n_line
                    FROM dnac d, contig c 
                   WHERE d.dna_id = c.dna_id 
                     AND c.contig_id = ?";    
      }

    }else{
      
      if( $end == -1 ) { 

	$query = "SELECT c.length, SUBSTRING( d.sequence, $start )
                    FROM dna d, contig c 
                   WHERE d.dna_id = c.dna_id 
                     AND c.contig_id = ?";
      
      } else {
	my $length = $end - $start + 1;
	if( $length < 1 ) {
	  $self->throw( "Wrong parameters" );
	}
    
	$query = "SELECT c.length, SUBSTRING( d.sequence, $start, $length )
                    FROM dna d, contig c 
                   WHERE d.dna_id = c.dna_id 
                     AND c.contig_id = ?";    
      }
      
    }

    #use the dna db if defined
    if( defined $self->db()->dnadb() ) {
      $sth = $self->db()->dnadb()->prepare( $query );
    } else {
      $sth = $self->prepare( $query );
    }

    $sth->execute($contig->dbID());

    if( my $aref = $sth->fetchrow_arrayref() ) {
      my ( $length, $seq, $n_line ) = @$aref;
      my $lenx=length($seq);
      if($compressed || $self->is_compressed()){
	$seq=$self->_dna_uncompress($length,$n_line,$start,$end,\$seq);
      }
      $seq =~ s/\s//g;
      $seq = uc($seq);
      if( $strand == -1 ) {
	return $self->_reverse_comp( $seq );
      } else {
	return uc($seq);
      }
    } else {

      my $flag;
      if($self->is_compressed==0){
	# found nothing...perhaps in this entry is compressed - check
	# and if present, set is_compressed
	eval{
	  $query = "SELECT c.length, d.dna_id
                      FROM dnac d, contig c 
                     WHERE d.dna_id = c.dna_id 
                       AND c.contig_id = ?";
	  $sth = $self->prepare( $query );
	  $sth->execute($contig->dbID());
	  if( my $aref = $sth->fetchrow_arrayref() ) {
	    $flag=1;
	  }
	};
	if($flag==1){
	  $self->is_compressed(1);
	  print "Switched to compressed DNA reading mode\n";
	  redo;
	}
      }else{
	# found nothing...perhaps in this entry is uncompressed - check
	# and if present, set is_compressed
	eval{
	  $query = "SELECT c.length, d.dna_id
                      FROM dna d, contig c 
                     WHERE d.dna_id = c.dna_id 
                       AND c.contig_id = ?";
	  $sth = $self->prepare( $query );
	  $sth->execute($contig->dbID());
	  if( my $aref = $sth->fetchrow_arrayref() ) {
	    $flag=1;
	  }
	};
	if($flag==1){
	  $self->is_compressed(0);
	  print "Switched to uncompressed DNA reading mode\n";
	  redo;
	}
      }
      return undef;
    }
  }
}




=head2 fetch_by_Slice_start_end_strand

  Arg  [1]   : Bio::EnsEMBL::Slice slice
               The slice from which you want the sequence
  Arg  [2]   : int startBasePair 
               count from 1
  Arg  [3]   : int endBasePair 
               count from 1, -1 is last one
  Arg  [4]   : int strand 
               1, -1
  Example    : $dna = $seq_adptr->fetch_by_Slice_start_end_strand($slice, 1, 
                                                                  1000, -1);
  Description: retrieves from db the sequence for this slice
               uses AssemblyMapper to find the assembly
  Returntype : string 
  Exceptions : endBasePair should be less or equal to length of slice 
  Caller     : Bio::EnsEMBL::Slice::seq(), Slice::subseq() 

=cut

sub fetch_by_Slice_start_end_strand {
   my ( $self, $slice, $start, $end, $strand ) = @_;

   my $seq= "";

   if( !$slice ){
     $self->throw("need a slice to work\n");
   }

   unless 
     ($slice->isa("Bio::EnsEMBL::Slice")) {
       $self->throw("$slice isn't a slice");
     }
   
   if( $end == -1 ) {
     $end = $slice->chr_end() - $slice->chr_start() + 1;
   }

   # need to check the strand'edness of the slice as this
   # affects the direction in which the dna seq is grabbed
   if ( $slice->strand == 1 ) {
     $seq = $self->fetch_by_assembly_location
       (
	$slice->chr_start()+$start-1,
	$slice->chr_start()+$end-1,
	$strand,
	$slice->chr_name(),
	$slice->assembly_type() 
       );
   }
   elsif ( $slice->strand == -1 ) {
     $seq = $self->fetch_by_assembly_location
       (
	$slice->chr_end()-$end+1,
	$slice->chr_end()-$start+1,
	$strand * -1, #have to make strand relative to slice's strand
	$slice->chr_name(),
	$slice->assembly_type() 
       );
   }
   else {
     $self->throw("Incorrect strand set on slice $slice");
   }
   return $seq;
}



=head2 fetch_by_assembly_location

  Arg   [1]  : int $chrStart
  Arg   [2]  : int $chrEnd
  Arg   [3]  : int $strand
  Arg   [4]  : txt $chrName
  Arg   [5]  : txt $assemblyType
  Example    : $dna = $fetch_by_assembly_location( 1, 100, -1, 'X', NCBI30 );
  Description: retrieve specified sequence from db. Using AssemblyMapper. Gaps
               are filled with N 
  Returntype : string 
  Exceptions : Wrong parameters give undef as result 
  Caller     : general, fetch_by_Slice_start_end_strand 

=cut

sub fetch_by_assembly_location {
   my ( $self, $chrStart, $chrEnd, 
        $strand, $chrName, $assemblyType ) = @_;

   my $mapper = 
     $self->db->get_AssemblyMapperAdaptor->fetch_by_type($assemblyType);
   # $mapper->register_region($chrName,$chrStart,$chrEnd);
   
   my @coord_list = $mapper->map_coordinates_to_rawcontig
     ( $chrName, $chrStart, $chrEnd, $strand );
   
   # for each of the pieces get sequence
   my $seq = "";
   for my $segment ( @coord_list ) {
     if( $segment->isa( "Bio::EnsEMBL::Mapper::Coordinate" )) {

       my $contig = 
	 $self->db->get_RawContigAdaptor()->fetch_by_dbID($segment->id());

       my $contig_seq = $self->fetch_by_RawContig_start_end_strand
	 ( $contig,
	   $segment->start(),
	   $segment->end(),
	   $segment->strand() );

       $seq .= $contig_seq;
	   
     } else {
       # its a gap
       my $length = $segment->end() - $segment->start() + 1;
       $seq .= "N" x $length;
     }
   }
   
   return $seq;
}



=head2 store

  Arg [1]    : string $sequence the dna sequence to be stored in the database
  Arg [2]    : string $date create date to be associated with the dna sequence
               to be stored.
  Example    : $dbID = $seq_adaptor->store('ACTGGGTACCAAACAAACACAACA', $date); 
  Description: stores a dna sequence in the databases dna table and returns the
               database identifier for the new record.
  Returntype : int
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::RawContigAdaptor::store

=cut

sub store {
  my ($self, $sequence, $date) = @_;
  
  $sequence =~ tr/atgcn/ATGCN/;

  # save compressed
  if($self->is_compressed){
      return $self->store_compressed($sequence,$date);
  }
  
  my $statement = $self->prepare("
        INSERT INTO dna(sequence,created) 
        VALUES(?, FROM_UNIXTIME(?))
        "); 
  
  my $rv = $statement->execute($sequence, $date);  
  $self->throw("Failed to insert dna $sequence") unless $rv;
  
  $statement->finish;

  $statement = $self->prepare("SELECT last_insert_id()");
  $statement->execute();
  my ($id) = $statement->fetchrow();
  $statement->finish;
  
  return $id;
}


=head2 store_compressed

  Arg [1]    : string $sequence the dna sequence to be stored in the database
  Arg [2]    : string $date create date to be associated with the dna sequence
               to be stored.
  Arg [3]    : dbid (optional)
  Returntype : int
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::RawContigAdaptor::store_compressed, 
               Bio::EnsEMBL::DBSQL::RawContigAdaptor::store

=cut

sub store_compressed {
  my ($self, $sequence, $date, $dbid) = @_;

  # look for N's, save list of start and ends and convert to A's
  my $n_line;
  my $len=length($sequence);
  #print "DEBUG $sequence ($len)\n";
  {
      if($sequence=~/^([ACGT]*)([^ACGT]+)(.*)$/){
	  my $l1=length($1);
	  my $l2=length($2);
	  my $start=$l1+1;
	  my $end=$l1+$l2;
	  $sequence=$1.('A' x $l2).$3;
	  # might be more than one letter
	  my $insert=$2;
	  #print "DEBUG: $insert ($start,$end)\n";
	  {
	    if($insert=~/^(\w)(\1*)(.*)$/){
	      my $insert2=$1.$2;
	      my $end2=$start+length($insert2)-1;
	      $n_line.=" " if $n_line;
	      my $modifier='';
	      if($1 ne 'N'){$modifier=":$1";}
	      $n_line.="$start-$end2$modifier";
	      #print "DEBUG: write $start-$end2$modifier\n";
	      $insert=$3;
	      $start=$end2+1;
	      redo if $insert;
	    }
	  }
	  redo;
      }
  }
  #print "DEBUG $sequence\n";

  # convert ACGT -> 00011011
  my $sequence_comp=$self->_dna_compress($sequence);
  
  my $id;
  if($dbid){
      my $statement = $self->prepare("
            INSERT INTO dnac(dna_id,sequence,created,n_line) 
            VALUES(?, ?, FROM_UNIXTIME(?), ?)
            "); 
  
      my $rv = $statement->execute($dbid,$sequence_comp, $date, $n_line);  
      $self->throw("Failed to insert dna $sequence") unless $rv;  
      $statement->finish;
      $id=$dbid;
  }else{
      my $statement = $self->prepare("
            INSERT INTO dnac(sequence,created,flag_dnan) 
            VALUES(?, FROM_UNIXTIME(?), ?)
            "); 
  
      my $rv = $statement->execute($sequence_comp, $date, $n_line);  
      $self->throw("Failed to insert dna $sequence") unless $rv;  
      $statement->finish;

      $statement = $self->prepare("SELECT last_insert_id()");
      $statement->execute();
      ($id) = $statement->fetchrow();
      $statement->finish;
  }
  return $id;
}


=head2 _dna_compress

  Arg  1    : txt $dna_sequence
  Function  : build hex representation of string (ACGT -> 1B)
  Returntype: txt
  Exceptions: none
  Caller    : private to this module

=cut

{

  my %table=(
	     'AA'=>'0',
	     'AC'=>'1',
	     'AG'=>'2',
	     'AT'=>'3',
	     'CA'=>'4',
	     'CC'=>'5',
	     'CG'=>'6',
	     'CT'=>'7',
	     'GA'=>'8',
	     'GC'=>'9',
	     'GG'=>'A',
	     'GT'=>'B',
	     'TA'=>'C',
	     'TC'=>'D',
	     'TG'=>'E',
	     'TT'=>'F',
	     );

  sub _dna_compress {
    my $self = shift;
    my $seq = shift;
    # length must be multiple of 4, so pad if necesary
    my $len=length($seq);
    $seq.='A' x ($len-(int($len/4))*4);
    #print "DEBUG $seq\n";
    # process in blocks of 2 to do conversion
    $seq=~s/(\G..)/$table{$1}/g;
    return pack("H*",$seq);
  }
}


=head2 _dna_uncompress

  Arg [1]    : int $length
  Arg [2]    : int $n_line
  Arg [3]    : int $start
  Arg [4]    : int $end
  Arg [5]    : reference to $seq (hex string)
  Function   : converts hex string from DB into ACGT with appropriate truncation
  Returntype : string
  Exceptions : none
  Caller     : private to this module

=cut

{

  my %table=(
	     '0'=>'AA',
	     '1'=>'AC',
	     '2'=>'AG',
	     '3'=>'AT',
	     '4'=>'CA',
	     '5'=>'CC',
	     '6'=>'CG',
	     '7'=>'CT',
	     '8'=>'GA',
	     '9'=>'GC',
	     'A'=>'GG',
	     'B'=>'GT',
	     'C'=>'TA',
	     'D'=>'TC',
	     'E'=>'TG',
	     'F'=>'TT',
	     );

  sub _dna_uncompress {
    my($self,$length,$n_line,$start,$end,$rseq)=@_;

    # convert sequence back to ACGT (already in Hex)
    my $seq2=join('',map{$table{$_}}split(//,$$rseq));

    # calculate extent of $seq2 and truncate
    my $start4=(int(($start-1)/4)+1)*4-3;
    my $off=$start-$start4;
    my $len;
    if($end==-1){
      $end=$length;
    }
    $len=$end-$start+1;
    #print "truncate: $start,$start4,$end,$len\n";
    my $seq=substr($seq2,$off,$len);
    
    # mask with N's
    foreach my $range (split(/ /,$n_line)){
      my($st,$ed)=split(/\-/,$range);
      my $char='N';
      if($ed=~/(\d+):(\w)/){
	$ed=$1;
	$char=$2;
      }
      #print "before: $st-$ed $start-$end\n";
      # check in range
      next if($ed<$start || $st>$end);
      $ed=$end if $ed>$end;
      $st=$start if $st<$start;
      #print "after: $st-$ed\n";
      my $len=$ed-$st+1;
      $st-=$start;
      substr($seq,$st,$len)=$char x $len;
    }

    return $seq;
  }
}


=head2 is_compressed

  Arg [1]    : (optional) flag
               to set or unset
  Example    : $adaptor->is_compressed(1)
  Description: Getter/Setter for whether compressed DNA should be read/written
  Returntype : true/false
  Exceptions : none
  Caller     : Adaptors inherited from BaseAdaptor

=cut

sub is_compressed { 
    my( $obj, $flag ) = @_;
    if (defined $flag) {
        $obj->{'_is_truncated'} = $flag ? 1 : 0;
    }
    return $obj->{'_is_truncated'};
}


=head2 _reverse_comp

  Arg  1    : txt $dna_sequence
  Function  : build reverse complement string
  Returntype: txt
  Exceptions: none
  Caller    : private to this module

=cut

sub _reverse_comp {
  my $self = shift;
  my $seq = shift;

  $seq = reverse( $seq );
  $seq =~
    tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
  return $seq;
}


1;
