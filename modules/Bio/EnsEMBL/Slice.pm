#
# Ensembl module for Bio::EnsEMBL::Assembly::Slice
#
# Cared for by Ewan Birney <ensembl-dev@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Slice - Arbitary Slice of a genome

=head1 SYNOPSIS


   foreach $gene ( $slice->get_all_Genes ) {
      # do something with a gene
   }
       

=head1 DESCRIPTION



=head1 AUTHOR - Ewan Birney

=head1 CONTACT

This modules is part of the Ensembl project http://www.ensembl.org

Questions can be posted to the ensembl-dev mailing list:
ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Slice;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::PrimarySeqI;
use Bio::EnsEMBL::Tile;

@ISA = qw(Bio::EnsEMBL::Root Bio::PrimarySeqI);



sub new {
  my($class,@args) = @_;

  my $self = {};
  bless $self,$class;

  my ($chr,$start,$end,$strand,$type,$adaptor, $dbID, $empty) = 
    $self->_rearrange([qw(CHR_NAME 
			  CHR_START 
			  CHR_END 
			  STRAND 
			  ASSEMBLY_TYPE 
			  ADAPTOR 
			  DBID
                          EMPTY)],
		      @args);

  if( ! defined $empty ) {
    if( !defined $chr || !defined $start || !defined $end || !defined $type ) {
      print STDERR "Chr: " . $chr . "\t" . "Start: " . $start . "\t" . 
	"End: " . $end . "\t" . "Type: " . $type . "\n";
      $self->throw("Do not have all the parameters for slice");
    }
    $self->chr_name($chr);
    $self->chr_start($start);
    $self->chr_end($end);

    if(!defined $strand) {
      $strand = 1; #default slice strand is 1
    } else {
      unless($strand == 1 || $strand == -1) {
	$self->throw("Slice strand must be either -1 or 1 not [$strand].");
      }
    }
    $self->strand($strand);
  } else {
    $self->strand( 1 );
    $self->chr_start( 1 );
    
    # empty Slices are used to do mapping to chromosomal coords.
    # After the mapping, the Slice contains chr_name and is reference 
    # point for the mapped object
  }

  $self->assembly_type($type);
  $self->adaptor($adaptor);
  $self->dbID( $dbID );
  # set stuff in self from @args

  if( defined $adaptor && !defined $type ) {
    $self->assembly_type
      ( $adaptor->db()->get_MetaContainer()->get_default_assembly());
  }
  return $self;
}



=head2 adaptor

  Arg [1]    : (optional) Bio::EnsEMBL::DBSQL::SliceAdaptor $adaptor
  Example    : $adaptor = $slice->adaptor();
  Description: Getter/Setter for the slice object adaptor used
               by this slice for database interaction.
  Returntype : Bio::EnsEMBL::DBSQL::SliceAdaptor
  Exceptions : none
  Caller     : general

=cut

sub adaptor{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'adaptor'} = $value;
    }
    return $self->{'adaptor'};
}



=head2 dbID

  Arg [1]    : (optioanl) int $value
  Example    : none
  Description: Getter/Setter for the unique database identifier for this 
               slice. This is not currently useful since slices are 
               abstractions and not actually stored in a database.  This 
               function is present to mirror RawContigs dbID method and
               because it could in theory be used one day.
  Returntype : int
  Exceptions : none
  Caller     : none

=cut

sub dbID {
   my ( $self, $value ) = @_;
   if( defined $value ) {
     $self->{'dbID'} = $value;
   }
   return $self->{'dbID'};
}



=head2 name

  Arg [1]    : none
  Example    : $name = $slice->name();
  Description: Returns the name of this slice. The name is formatted as a 
               the following string: "$chr_name:$chr_start-$chr_end". 
               (e.g. 'X.10000-20000')
               This essentially allows slices to be easily compared and 
               can also act as a hash value. This is similar to the name 
               method in RawContig so for exons which can have either type 
               of sequence attached it provides a more common interface.
  Returntype : string
  Exceptions : none
  Caller     : 

=cut

sub name {
  my $self = shift;

  my $string = join '', $self->chr_name, '.', 
                        $self->chr_start, '-', $self->chr_end();

  if($self->strand == -1) {
    return "reverse($string)";
  }

  return $string;
}





=head2 id

  Arg [1]    : none 
  Example    : none
  Description: Here to mirror same method in RawContig.  Simply returns 
               the same thing as $slice->name() and generally name should be
               used instead.
  Returntype : string
  Exceptions : none
  Caller     : none

=cut

sub id {
   my $self = shift;

   return $self->name() || $self->dbID();
}



=head2 length

  Arg [1]    : none
  Example    : $length = $slice->length();
  Description: Returns the length of this slice
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub length {
  my ($self) = @_;

  return $self->chr_end() - $self->chr_start() + 1;
}



=head2 invert

  Arg [1]    : none
  Example    : $slice->invert;
  Description: Creates a copy of this slice on the opposite strand and 
               returns it.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : general

=cut

sub invert {
  my $self = shift;
  
  my %s = %$self;
  my $slice = bless \%s, ref $self;

  $slice->strand($self->strand * -1);

  return $slice;
}


=head2 seq

  Args      : none
  Function  : Returns the entire sequence string for this Slice.  This method
              will return the reverse complement of the sequence of another 
              slice on the same region but on the opposite strand.
              Note that the slice needs the adaptor to be set to obtain the 
              sequence.
  Returntype: txt
  Exceptions: none
  Caller    : general

=cut

sub seq {
  my $self = shift;
  my $seqAdaptor = $self->adaptor->db->get_SequenceAdaptor();
  return $seqAdaptor->fetch_by_Slice_start_end_strand( $self, 1, -1, 1 );
}



=head2 subseq

  Arg  1    : int $startBasePair
              relative to start of slice, which is 1.
  Arg  2    : int $endBasePair
              relative to start of slice.
  Arg  3    : (optional) int $strand
              The strand of the slice to obtain sequence from. Default value is
              1.
  Function  : returns string of dna sequence
  Returntype: txt
  Exceptions: end should be at least as big as start
              strand must be set
  Caller    : general

=cut

sub subseq {
  my ( $self, $start, $end, $strand ) = @_;

  if ( $end < $start ) {
    $self->throw("End coord is less then start coord");
  }

  $strand = 1 unless(defined $strand);

  if ( $strand != -1 && $strand != 1 ) {
    $self->throw("Invalid strand [$strand] in call to Slice::subseq.");
  }

  my $seqAdaptor = $self->adaptor->db->get_SequenceAdaptor();
  my $seq = $seqAdaptor->fetch_by_Slice_start_end_strand( $self, $start, 
                                                          $end, $strand );

  return $seq;
}



=head2 get_all_PredictionTranscripts

  Arg [1]    : (optional) string $logic_name
               The name of the analysis used to generate the prediction
               transcripts obtained.
  Example    : @transcripts = $slice->get_all_PredictionTranscripts();
  Description: Retrieves the list of prediction transcripts which overlap
               this slice.
  Returntype : listref of Bio::EnsEMBL::PredictionTranscript
  Exceptions : none
  Caller     : none

=cut

sub get_all_PredictionTranscripts {
   my ($self,$logic_name) = @_;

   my $pta = $self->adaptor()->db()->get_PredictionTranscriptAdaptor();

   return $pta->fetch_all_by_Slice($self, $logic_name);
}



=head2 get_all_DnaAlignFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the dna align features
               to obtain.
  Arg [2]    : (optional) float $score
               The mimimum score of the features to retrieve
  Example    : @dna_align_feats = $slice->get_all_DnaAlignFeatures()
  Description: Retrieves the DnaDnaAlignFeatures which overlap this slice
  Returntype : listref of Bio::EnsEMBL::DnaDnaAlignFeatures
  Exceptions : none
  Caller     :general

=cut

sub get_all_DnaAlignFeatures {
   my ($self, $logic_name, $score) = @_;

   my $dafa = $self->adaptor->db->get_DnaAlignFeatureAdaptor();

   return $dafa->fetch_all_by_Slice_and_score($self,$score, $logic_name);
}



=head2 get_all_ProteinAlignFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the protein align features
               to obtain.
  Arg [2]    : (optional) float $score
               The mimimum score of the features to retrieve
  Example    : @pep_align_feats = $slice->get_all_ProteinAlignFeatures()
  Description: Retrieves the PepDnaAlignFeatures which overlap this slice.
  Returntype : listref of Bio::EnsEMBL::PepDnaAlignFeatures
  Exceptions : none
  Caller     : general

=cut

sub get_all_ProteinAlignFeatures {
  my ($self, $logic_name, $score) = @_;

  my $pafa = $self->adaptor()->db()->get_ProteinAlignFeatureAdaptor();

  return $pafa->fetch_all_by_Slice_and_score($self, $score, $logic_name);
}



=head2 get_all_SimilarityFeatures

  Arg [1]    : (optional) string $logic_name
               the name of the analysis performed on the features to retrieve
  Arg [2]    : (optional) float $score
               the lower bound of the score of the features to be retrieved
  Example    : @feats = $slice->get_all_SimilarityFeatures()
  Description: Retrieves all dna_align_features and protein_align_features
               with analysis named $logic_name and with score above $score.
               It is probably faster to use get_all_ProteinAlignFeatures or
               get_all_DnaAlignFeatures if a sepcific feature type is desired.
  Returntype : listref of Bio::EnsEMBL::BaseAlignFeatures
  Exceptions : none
  Caller     : general

=cut

sub get_all_SimilarityFeatures {
  my ($self, $logic_name, $score) = @_;

  my @out = ();

  push @out, @{$self->get_all_ProteinAlignFeatures($logic_name, $score) };
  push @out, @{$self->get_all_DnaAlignFeatures($logic_name, $score) };

  return \@out;
}


=head2 get_all_SimpleFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the simple features
               to obtain.
  Arg [2]    : (optional) float $score
               The mimimum score of the features to retrieve
  Example    : @simple_feats = $slice->get_all_SimpleFeatures()
  Description: Retrieves the SimpleFeatures which overlap this slice.
  Returntype : listref of Bio::EnsEMBL::DnaDnaAlignFeature
  Exceptions : none
  Caller     : general

=cut

sub get_all_SimpleFeatures {
  my ($self, $logic_name, $score) = @_;

  my $sfa = $self->adaptor()->db()->get_SimpleFeatureAdaptor();

  return $sfa->fetch_all_by_Slice_and_score($self, $score, $logic_name);
}



=head2 get_all_RepeatFeatures

  Arg [1]    : (optional) string $logic_name
               The name of the analysis performed on the repeat features
               to obtain.
  Example    : @repeat_feats = $slice->get_all_RepeatFeatures()
  Description: Retrieves the RepeatFeatures which overlap this slice.
  Returntype : listref of Bio::EnsEMBL::RepeatFeatures
  Exceptions : none
  Caller     : general

=cut

sub get_all_RepeatFeatures {
   my ($self, $logic_name) = @_;

   my $rpfa = $self->adaptor()->db()->get_RepeatFeatureAdaptor();

   return $rpfa->fetch_all_by_Slice($self, $logic_name);
}

=head2 get_all_SNPs

  Args      : none
  Function  : returns all SNPs on this slice
  Returntype: listref of Bio::EnsEMBL::External::Variation
  Exceptions: none
  Caller    : GlyphSet_feature inherited objects

=cut

sub get_all_SNPs {
  my $self = shift;

  my $snpa = $self->adaptor()->db()->get_SNPAdaptor();

  return $snpa->fetch_all_by_Slice($self);
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
   my ($self, $empty_flag) = @_;

   #caching is performed on a per slice basis in the GeneAdaptor
   return $self->adaptor->db->get_GeneAdaptor()->fetch_all_by_Slice($self, $empty_flag);
}



=head2 get_Genes_by_source

 Title   : get_Genes_by_source
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub get_Genes_by_source {
  my $self = shift; 
  warn( "\$Slice->get_Genes_by_source(\$source): You don't want to do it like that!.\nYou want to do it like this \$Slice->get_all_Genes_by_source(\$source)");
  return $self->get_all_Genes_by_source( @_ );
}

sub get_all_Genes_by_source{
   my ($self, $source, $empty_flag) = @_;
   my @out = grep { $_->source eq $source } @{$self->get_all_Genes($empty_flag)};
   return \@out;
}



=head2 get_Genes_by_type

 Title   : get_Genes_by_type
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Genes_by_type {
  
  my $self = shift; 
  warn( "\$Slice->get_Genes_by_type(\$type): You don't want to do it like that!.\nYou want to do it like this \$Slice->get_all_Genes_by_type(\$type)");
  return $self->get_all_Genes_by_type( @_ );
}


sub get_all_Genes_by_type{
  my ($self, $type, $empty_flag) = @_;
  
  my @out = grep { $_->type eq $type } @{ $self->get_all_Genes($empty_flag)};
  
  return \@out;
}



=head2 chr_name

 Title   : chr_name
 Usage   : $obj->chr_name($newval)
 Function: 
 Example : 
 Returns : value of chr_name
 Args    : newvalue (optional)


=cut

sub chr_name{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'chr_name'} = $value;
    }
    return $self->{'chr_name'};

}

=head2 chr_start

 Title   : chr_start
 Usage   : $obj->chr_start($newval)
 Function: 
 Example : 
 Returns : value of chr_start
 Args    : newvalue (optional)


=cut

sub chr_start{
  my ($self,$value) = @_;
  if( defined $value) {
    $self->{'chr_start'} = $value;
  }
  return $self->{'chr_start'};
}



=head2 chr_end

 Title   : chr_end
 Usage   : $obj->chr_end($newval)
 Function: 
 Example : 
 Returns : value of chr_end
 Args    : newvalue (optional)

=cut

sub chr_end{
  my ($self,$value) = @_;
  if( defined $value) {
    $self->{'chr_end'} = $value;
  }
  return $self->{'chr_end'};
}



=head2 strand

 Title   : strand
 Usage   : $obj->strand($newval)
 Function: 
 Example : 
 Returns : value of strand
 Args    : newvalue (optional)

=cut

sub strand{
   my ($self,$value) = @_;

   if( defined $value) {
      $self->{'strand'} = $value;
    }
    return $self->{'strand'};

}



=head2 assembly_type

 Title   : assembly_type
 Usage   : $obj->assembly_type($newval)
 Function: 
 Example : 
 Returns : value of assembly_type
 Args    : newvalue (optional)


=cut

sub assembly_type{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'assembly_type'} = $value;
    }
    return $self->{'assembly_type'};

}


sub get_all_KaryotypeBands {
  my ($self) = @_;
  
  my $kadp = $self->adaptor->db->get_KaryotypeBandAdaptor();
  return $kadp->fetch_all_by_Slice($self);
}


sub get_Chromosome {
  my $self = shift @_;

  my $ca =  $self->adaptor->db->get_ChromosomeAdaptor();

  return $ca->fetch_by_chr_name($self->chr_name());
}


=head2 get_repeatmasked_seq

  Arg [1]    : string $logic_name (optional)
  Arg [2]    : int $soft_masking_enable (optional)
  Example    : $slice->get_repeatmasked_seq 
               or $slice->get_repeatmasked_seq('RepeatMask',1)
  Description: Returns Bio::PrimarySeq containing the masked (repeat replaced 
               by N) 
               or soft-masked (when Arg[2]=1, repeat in lower case while non
               repeat in upper case) sequence corresponding to the Slice 
               object.
               Will only work with database connection to get repeat features.
  Returntype : Bio::PrimarySeq
  Exceptions : none
  Caller     : general

=cut

sub get_repeatmasked_seq {
    my ($self,$logic_name,$soft_mask) = @_;

    if(!$logic_name){
      $logic_name = 'RepeatMask';
    }

    unless (defined $soft_mask) {
      $soft_mask = 0;
    }

    #$self->warn("Slice: get_repeatmasked_seq\n");      

    my $repeats = $self->get_all_RepeatFeatures($logic_name);
    my $dna = $self->seq();
    my $masked_dna = $self->_mask_features($dna,$repeats,$soft_mask);
    my $masked_seq = Bio::PrimarySeq->new('-seq'        => $masked_dna,
					  '-display_id' => $self->id,
					  '-primary_id' => $self->id,
					  '-moltype'    => 'dna'
					 );
    return $masked_seq;
}

=head2 _mask_features

  Arg [1]    : string $dna_string
  Arg [2]    : array_ref $repeats
               reference to a list Bio::EnsEMBL::RepeatFeature
               give the list of coordinates to replace with N or with 
               lower case
  Arg [3]    : int $soft_masking_enable (optional)
  Example    : 
  Description: replaces string positions described in the RepeatFeatures
               with Ns (default setting), or with the lower case equivalent 
               (soft masking)
  Returntype : string 
  Exceptions : none
  Caller     : get_repeatmasked_seq

=cut

sub _mask_features {
  my ($self,$dnastr,$repeats,$soft_mask) = @_;
  
  # explicit CORE::length call, to avoid any confusion with the Slice 
  # length method
  my $dnalen = CORE::length($dnastr);

 REP:foreach my $f (@{$repeats}) {
    my $start  = $f->start;
    my $end    = $f->end;
    my $length = ($end - $start) + 1;
    
    # check if we get repeat completely outside of expected slice range
    if ($end < 1 || $start > $dnalen) {
      warn "Repeat completely outside slice coordinates! " .
	"That should not happen! repeat_start $start or repeat_end $end not" .
	"within [1-$dnalen] slice range coordinates\n";
      next REP;
    }
    
    # repeat partly outside slice range, so correct
    # the repeat start and length to the slice size if needed
    if ($start < 1) { 
      $start = 1;
      $length = ($end - $start) + 1;
    }
    
    # repeat partly outside slice range, so correct
    # the repeat end and length to the slice size if needed
    if ($end > $dnalen) {
      $end = $dnalen;
      $length = ($end - $start) + 1;
    }

    $start--;
    
    my $padstr;
    
    if ($soft_mask) {
      $padstr = lc substr ($dnastr,$start,$length);
    } else {
      $padstr = 'N' x $length;
    }
    substr ($dnastr,$start,$length) = $padstr;

  }
  return $dnastr;
} 


sub get_all_MapFrags {
    my $self = shift;
    my $mapset = shift;

    my $mfa = $self->adaptor()->db()->get_MapFragAdaptor();

    return $mfa->fetch_all_by_mapset_chr_start_end($mapset, 
					       $self->chr_name,
					       $self->chr_start, 
					       $self->chr_end);
}    



sub has_MapSet {
  my( $self, $mapset_name ) = @_;
    
  my $mfa = $self->adaptor()->db()->get_MapFragAdaptor();

  return $mfa->has_mapset($mapset_name);
}



=head2 get_tiling_path

  Arg [1]    : none
  Example    : @tiles = @{$slice->get_tiling_path()};
  Description: Retrieve a listref of Bio::EnsEMBL::Tile objects representing
               the tiling path used to construct the contiguous slice sequence.
  Returntype : list reference of Bio::EnsEMBL::Tile objects
  Exceptions : none
  Caller     : general

=cut

sub get_tiling_path {
  my ($self) = @_;

  my $mapper = $self->adaptor()->db->get_AssemblyMapperAdaptor()->
    fetch_by_type($self->assembly_type());

  # Get the ids of the raw_contigs in this region specified in chrmsml coords 

  my @mapped = $mapper->map_coordinates_to_rawcontig
    (
     $self->chr_name(),
     $self->chr_start(),
     $self->chr_end(),
     $self->strand()
    );

  # Extract the IDS of the Coordinates, ommitting Gaps
  my @raw_contig_ids = ();
  foreach my $map_item (@mapped) {
    if($map_item->isa("Bio::EnsEMBL::Mapper::Coordinate" )) {
       push @raw_contig_ids, $map_item->id();
     } 
  }

  #Fetch filled raw contigs (non lazy-loaded) containing filled clone objects
  my $rca = $self->adaptor->db->get_RawContigAdaptor();
  my $raw_contigs = $rca->fetch_filled_by_dbIDs(@raw_contig_ids);

  my @tiling_path = ();
  my $current_start = 1;

  my($length, $slice_start, $slice_end, 
     $contig, $contig_start, $contig_end, $contig_ori);  

  foreach my $coord ( @mapped ) {
    $contig_start = $coord->start();
    $contig_end   = $coord->end();
    $length       = $contig_end - $contig_start + 1; 

    if ( $coord->isa("Bio::EnsEMBL::Mapper::Coordinate" ) ) {
      # create a tile for each coordinate
      $contig_ori  =  $coord->strand();
      $slice_start = $current_start;
      $slice_end   = $current_start + $length - 1;
      $contig      = $raw_contigs->{ $coord->id() };
      
      push @tiling_path, Bio::EnsEMBL::Tile->new_fast($self,
						      $slice_start,
						      $slice_end,
						      $contig,
						      $contig_start,
						      $contig_end,
						      $contig_ori);
						
 
      
      $current_start += $length;
    } else {
      # this is a gap, just add the length and discard it
      $current_start += $length;
    }
  }
  return \@tiling_path;
}
  


=head2 get_all_landmark_MarkerFeatures

  Arg [1]    : none
  Example    : my @marker_feats = @{$slice->get_all_landmark_MarkerFeatures};
  Description: Retrieves all landmark markers which lie on this slice. This
               query uses the lite database, and an empty list will be returned
               if the lite database is not available.
  Returntype : reference to a list of Bio::EnsEMBL::MarkerFeatures
  Exceptions : warning if lite database is not available
  Caller     : contigview?, general

=cut

sub get_all_landmark_MarkerFeatures {
  my $self = shift;

  my $lma = $self->adaptor()->db()->get_LandmarkMarkerAdaptor();
  if( ! defined $lma ) {
   $self->warn("Lite database must be available to retrieve landmark markers");
   return [];
  } else {
    return $lma->fetch_all_by_Slice( $self );
  }

}


### This is a hacky little function which gets the name of the contig
### from a finished clone - or returns undef otherwise.
### A finised clone is a clone with only one contig which starts at bp 1
### (So grab the last contig of the clone, and see if it starts at 1 <g>)
sub contig_from_clone {
    my ($self,$clone) = @_;
    my $sth = $self->adaptor->db->prepare(
            "select co.name, co.embl_offset from clone as cl, contig as co
              where co.clone_id = cl.clone_id and cl.name = '$clone'
              order by embl_offset desc limit 1"
    );
    $sth->execute();
    my ($contig,$offset) = $sth->fetchrow_array();
    return $offset==1 ? $contig : undef;
}

sub get_all_DASFeatures{
   my ($self,@args) = @_;

   if( defined $self->{'_das_cached_features'} ) {
       return $self->{'_das_cached_features'};
   }

   my %genomic_features;

   my $mapper = $self->adaptor()->db->get_AssemblyMapperAdaptor()->fetch_by_type($self->assembly_type());

   my @raw_contig_ids = $mapper->list_contig_ids( $self->chr_name,
                                                  $self->chr_start,
                                                  $self->chr_end );


   # need a call here to get a list of FPC contigs that overlap my VC
   # I also need to have their VC start end in the FPC coordinates.
   # and somehow pass all this stuff down to the DAS fetcher...eek!
   my @fpccontigs = (undef);
   
   my $rca = $self->adaptor()->db()->get_RawContigAdaptor();
   my $raw_Contig_Hash = $rca->fetch_filled_by_dbIDs( @raw_contig_ids );

   # provide mapping from contig names to internal ids
   my %contig_name_hash =
     map { ( $_->name(), $_) } values %$raw_Contig_Hash;
   my @raw_contig_names = keys %contig_name_hash;

   # retrieve all embl clone accessions
   my %clone_hash  =
     map {( $_->clone->embl_id(), 1 ) } values %$raw_Contig_Hash;
   my @clones = keys %clone_hash;


   my $chr_length = $self->get_Chromosome()->length();
   my $chr_start  = $self->chr_start();
   my $offset     = 1 - $chr_start;
   my $chr_end    = $self->chr_end();
   my $chr_strand = $self->strand();
   my $length     = $chr_end + $offset;

foreach my $extf ( $self->adaptor()->db()->_each_DASFeatureFactory ) {
       if( $extf->can('get_Ensembl_SeqFeatures_DAS') ) {
           my @contig_features;
           my @chr_features;
           my @fpc_features;
           my @clone_features;
           my $dsn = $extf->_dsn();

   #  warn("Fetching features from $dsn!!!" );
          foreach my $sf (
              @{ $extf->get_Ensembl_SeqFeatures_DAS(
                    $self->chr_name,$self->chr_start,$self->chr_end,
                    \@fpccontigs, \@clones,\@raw_contig_names, $chr_length) }
           ) {
# BAC.*_C are fly contigs....
# CRA_x are Celera mosquito contigs....
#	  warn( join ' - ', $sf->das_dsn, $sf->das_id, $sf->seqname, $sf->das_start, $sf->das_end, $sf->das_strand );
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
                     push(@{$genomic_features{$dsn}}, $sf);
               } elsif( $sf->seqname() eq '') {
                    #suspicious
                    warn ("Got a DAS feature with an empty seqname! (discarding it)\n");
               } else {
                   warn ("Got a DAS feature with an unrecognized segment type: >", $sf->seqname(), "< >", $sf->das_type_id(), "<\n");
               }
           }
           foreach my $sf ( @contig_features ) {
#warn( sprintf "SEG ID: %s\tID: %s\t DSN %s\t (%s-%s) STRAND: %s\tTYPE: %s\n", $sf->seqname, $sf->das_id, $sf->das_dsn, $sf->das_start, $sf->das_end, $sf->das_strand, $sf->das_type_id );
     # map to a chromosomal coordinate
               my( $coord ) = $mapper->map_coordinates_to_assembly( $contig_name_hash{ $sf->seqname() }->dbID(), $sf->das_start, $sf->das_end, $sf->das_strand );
#warn( join ' - ',  $sf->das_id, $contig_name_hash{ $sf->seqname() }->dbID(), $sf->das_start, $sf->das_end, $sf->das_strand );
     # if its not mappable than ignore the feature
               next if( $coord->isa( "Bio::EnsEMBL::Mapper::Gap" ));

               $sf->das_move( $coord->{'start'}+$offset, $coord->{'end'}+$offset, $coord->{'strand'} );
#warn( sprintf STDERR "SEG ID: %s\tID: %s\t DSN %s\t (%s-%s) STRAND: %s\tTYPE: %s\n", $sf->seqname, $sf->das_id, $sf->das_dsn, $sf->das_start, $sf->das_end, $sf->das_strand, $sf->das_type_id );
               next if $sf->das_start > $length || $sf->das_end < 1;
               push @{$genomic_features{$dsn}}, $sf
           }
           foreach my $sf ( @chr_features ) {
               $sf->das_shift( $offset );
               next if $sf->das_start > $length || $sf->das_end < 1; # skip anything thats fallen off the end
               push @{$genomic_features{$dsn}}, $sf
           }
       } else {
           warn "Slice: Got a DAS feature factory that can't do get_Ensembl_SeqFeatures_DAS\n";
       }
  
   }

   return $self->{'_das_cached_features'} = \%genomic_features;
}


=head2 get_all_ExternalFeatures

  Arg [1]    : (optional) string $track_name
               If specified only features from ExternalFeatureAdaptors with 
               the track name $track_name are retrieved.  
               If not set, all features from every ExternalFeatureAdaptor are 
               retrieved.
  Example    : @x_features = @{$slice->get_all_ExternalFeatures}
  Description: Retrieves features on this slice from external feature adaptors 
  Returntype : listref of Bio::SeqFeatureI implementing objects in slice 
               coordinates 
  Exceptions : none
  Caller     : general

=cut

sub get_all_ExternalFeatures {
   my ($self, $track_name) = @_;

   my $features = [];

   my $xfa_hash = $self->adaptor->db->get_ExternalFeatureAdaptors;
   my @xf_adaptors = ();

   if($track_name) {
     #use a specific adaptor
     push @xf_adaptors, $xfa_hash->{$track_name};
   } else {
     #use all of the adaptors
     push @xf_adaptors, values %$xfa_hash;
   }


   foreach my $xfa (@xf_adaptors) {
     push @$features, @{$xfa->fetch_all_by_Slice($self)};
   }

   return $features;
}



=head2 Methods included only for BioPerl compliance
=cut
###############################################################################

=head2 display_id

  Arg [1]    : none
  Example    : none
  Description: Only for BioPerl compliance.
  Returntype : string
  Exceptions : none
  Caller     : none

=cut

sub display_id{
  my $self = shift;

  return $self->id();
}

=head2 desc

  Arg [1]    : none
  Example    : none
  Description: Only for BioPerl compliance
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub desc{
  my $self = shift;
  return "Slice, no description";
}

=head2 moltype

  Arg [1]    : none
  Example    : none
  Description: Only for BioPerl compliance
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub moltype {
  my $self = shift;
  return 'DNA';
}

=head2 accession_number

  Arg [1]    : none
  Example    : none
  Description: Only for BioPerl compliance
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub accession_number {
  my $self = shift;
  return $self->dbID();
}


=head2 DEPRECATED methods
=cut
###############################################################################



=head2 primary_seq

 Title   : primary_seq
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

=head2 primary_seq

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use Bio::EnsEMBL:: instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub primary_seq{
   my ($self,@args) = @_;
   $self->warn("Call to deprecated method Bio::EnsEMBL::Slice::primary_seq" .
               "Slice is now a PrimarySeq and can be used directly as such\n");

   return $self;
}

=head2 get_all_Genes_exononly

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use get_all_Genes instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_all_Genes_exononly{
   my ($self) = @_;

   my ($p,$f,$l) = caller;
   $self->warn("$f:$l get_all_Genes_exononly has been deprecated. get_all_Genes called");

   return $self->get_all_Genes();
}


=head2 get_all_SangerGenes_startend_lite

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use get_Genes_by_source instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_all_SangerGenes_startend_lite {
  my $self = shift;

  $self->warn("Slice->get_all_SangerGenes_startend_lite deprecated" . 
	      " use get_allGenes() instead\n");
  
  return $self->get_Genes_by_source('sanger');
}


=head2 get_all_VirtualGenes_startend_lite

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use get_Genes_by_source instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut
  
sub get_all_VirtualGenes_startend_lite {
  my $self = shift;

  $self->warn("Slice->get_all_VirtualGenes_startend_lite deprecated" .
	      " use get_all_Genes() instead\n");

  return $self->get_all_Genes();
}


=head2 get_all_EMBLGenes_startend_lite

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use get_Genes_by_source instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_all_EMBLGenes_startend_lite {
  my $self = shift;

  $self->warn("Slice->get_all_EMBLGenes_startend_lite deprecated" .
	      " use get_Genes_by_source() instead\n");

  return $self->get_Genes_by_source('embl');
}



=head2 fetch_chromosome_length

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use get_Chromosome()->length() instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub fetch_chromosome_length {
  my ($self) = @_;

  $self->warn( "Call to deprecated method fetch_chromosome_length\n" .
	       "use \$slice->get_Chromosome()->length(); instead.\n" .
	       $self->stack_trace_dump());

  return $self->get_Chromosome()->length();
}



=head2 fetch_karyotype_band_start_end

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use get_KaryotypeBands instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut
        
sub fetch_karyotype_band_start_end {
   my ($self,@args) = @_;

   $self->warn( "Call to deprecated method fetch_karyotype_band_start_end\n" .
		"use \$slice->get_KaryotypeBands(); instead.\n" .
	       $self->stack_trace_dump());

   return $self->get_KaryotypeBands();
}



=head2 convert_Gene_to_raw_contig

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use Bio::EnsEMBL::Gene::transform instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub convert_Gene_to_raw_contig{
  my ($self,$gene) = @_;

  $self->warn("Call to deprecated method convert_Gene_to_raw_contig" . 
	      "use the gene's transform method directly\n");

  if(!$gene->isa("Bio::EnsEMBL::Gene")){
    $self->throw("trying to use the wrong method can called convert gene to RawContig coords on ".$gene."\n");
  }
     
  $gene->transform;
  
  return $gene;
}




=head2 get_landmark_MarkerFeatures

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use get_all_landmarkMarkerFeatures instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_landmark_MarkerFeatures {
  my ($self, @args) = @_;

  $self->warn("get_landmark_MarkerFeatures has been renamed get_all_landmarkMarkerFeatures\n" . caller);

  return $self->get_all_landmark_MarkerFeatures(@args);
}



1;
