
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

Describe the object here

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email ensembl-dev@ebi.ac.uk

Describe contact details here

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

@ISA = qw(Bio::EnsEMBL::Root Bio::PrimarySeqI);


# new() is written here 

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
    $self->id("$chr.$start-$end");
    
    #set strand to a default of 1 if it is not set
    if ( defined $strand) {
      $self->strand($strand);
    } else {
      $self->strand('1');
    }
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


=head2 First pass implementation

The first pass implementation tries to mimic precisely the important
parts of the old VirtualContig system for the pipeline. These
are the methods to implement

=cut



=head2 get_all_PredictionTranscripts

 Title   : get_all_PredictionTranscripts
 Usage   : $obj->get_all_PredictionTranscripts
 Function: Use to derive a list of prediction transcripts specific to the analysis type specified by the logic name.
 Example : my @pred_rm_trans = $obj->get_all_PredictionTranscripts('RepeatMasker');
 Returns : a list of Bio::EnsEMBL::PredictionTranscript objects
 Args    : a logic name - the name of the analysis that created or returned the prediction feature.


=cut

sub get_all_PredictionTranscripts{
   my ($self,$logic_name) = @_;

   my $pta = $self->adaptor()->db()->get_PredictionTranscriptAdaptor();

   return $pta->fetch_by_Slice($self, $logic_name);
}



  

=head2 get_all_SNPs

  Args      : none
  Function  : returns all SNPs on this slice
  Returntype: @Bio::EnsEMBL::External::Variation
  Exceptions: none
  Caller    : GlyphSet_feature inherited objects

=cut

sub get_all_SNPs {
  my $self = shift;

  my $snpa = $self->adaptor()->db()->get_SNPAdaptor();

  return $snpa->fetch_by_Slice($self);
}



=head2 get_all_DnaAlignFeatures_above_score

  Args      : $logic_name, $score
  Function  : returns all DnaAlignFeatures of type logic_name and above score
  Returntype: @Bio::EnsEMBL::DnaDnaAlignFeature
  Exceptions: none
  Caller    : GlyphSet_feature inherited objects

=cut

sub get_all_DnaAlignFeatures_above_score{
   my ($self,$logic_name, $score) = @_;

   if( !defined $score ) {
     $self->throw("No defined score.");
   }

   my $dafa = $self->adaptor->db->get_DnaAlignFeatureAdaptor();

   return $dafa->fetch_by_Slice_and_score($self,$score, $logic_name);
}


=head2 get_all_ProteinAlignFeatures_above_score

  Args      : $logic_name, $score
  Function  : getss all ProteinAlignFeatures of type logic_name and above score
  Returntype: @Bio::EnsEMBL::DnaPepAlignFeature
  Exceptions: none
  Caller    : GlyphSet_feature inherited objects
  
=cut

sub get_all_ProteinAlignFeatures_above_score {
  my ($self, $logic_name, $score) = @_;

  my $pafa = $self->adaptor()->db()->get_ProteinAlignFeatureAdaptor();

  return $pafa->fetch_by_Slice_and_score($self, $score, $logic_name);
}


=head2 get_all_SimpleFeatures_above_score

  Args      : $logic_name, $score
  Function  : retrieves all SimpleFeatures on this Slice of type logic_name 
              and above score
  Returntype: @Bio::EnsEMBL::SimpleFeature
  Exceptions: none
  Caller    : GlyphSet_feature inherited objects
  
=cut

sub get_all_SimpleFeatures_above_score {
  my ($self, $logic_name, $score) = @_;

  my $sfa = $self->adaptor()->db()->get_SimpleFeatureAdaptor();

  return $sfa->fetch_by_Slice_and_score($self, $score, $logic_name);
}


=head2 seq

  Args      : none
  Function  : returns the entire sequence string for this Slice
              needs the adaptor to be set.
  Returntype: txt
  Exceptions: none
  Caller    : general

=cut

sub seq {
  my $self = shift;
  my $seqAdaptor = $self->adaptor->db->get_SequenceAdaptor();
  my $seq = $seqAdaptor->fetch_by_Slice_start_end_strand( $self, 1, -1, 1 );

  return $seq;
}


=head2 subseq

  Arg  1    : int $startBasePair
              relative to start of slice, which is 1.
  Arg  2    : int $endBasePair
              relative to start of slice.
  Arg  3    : int $strand
  Function  : returns string of dna sequence
  Returntype: txt
  Exceptions: end should be at least as big as start
              strand must be set
  Caller    : general

=cut

sub subseq {
  my ( $self, $start, $end, $strand ) = @_;

  if ( $end < $start ) {
    $self->throw("End coord is less then start coord to call on Slice subseq.");
  }

  if ( !defined $strand || ( $strand != -1 && $strand != 1 )) {
#    $self->throw("Incorrect strand information set to call on Slice subseq.");
    $strand = 1;
  }

  my $seqAdaptor = $self->adaptor->db->get_SequenceAdaptor();
  my $seq = $seqAdaptor->fetch_by_Slice_start_end_strand( $self, $start, $end, $strand );

  return $seq;
}



=head2 get_all_SimilarityFeatures_above_pid

 Title   : get_all_SimilarityFeatures_above_pid
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SimilarityFeatures_above_pid{
   my ($self,@args) = @_;

   $self->throw("Ewan has not implemented this function! Complain!!!!");
}




=head2 get_all_RepeatFeatures

 Title   : get_all_RepeatFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_RepeatFeatures{
   my ($self, $logic_name) = @_;

   my $rpfa = $self->adaptor()->db()->get_RepeatFeatureAdaptor();

   return $rpfa->fetch_by_Slice($self, $logic_name);
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
   my ($self,@args) = @_;

   #if the genes have not been requested before, cache them now
   unless(defined $self->{_gene_cache}) {
     my $gene_adaptor = $self->adaptor->db->get_GeneAdaptor();
     my @genes = $gene_adaptor->fetch_by_Slice($self);
     foreach my $gene (@genes) {
       $self->{_gene_cache}->{$gene->dbID()} = $gene;
     }
   }
     
   return values(%{$self->{_gene_cache}});
}


=head2 get_Genes_by_source

 Title   : get_Genes_by_source
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Genes_by_source{
   my ($self,$source) = @_;
   my @genes = $self->get_all_Genes();
   
   my @out = ();

   foreach my $gene (@genes) {
     if($gene->source() eq $source) {
       push @out, $gene;
     }
   }

   return @out;
}

=head2 get_Genes_by_type

 Title   : get_Genes_by_type
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Genes_by_type{
   my ($self,$type) = @_;
   
   # Possibly this can be improved by selecting genes a query,
   # we expect that most times there will not be many genes in a region
   # however
   my @genes = $self->get_all_Genes();
   
   my @out = ();

   print STDERR "Slice : Getting genes of type $type \n";

   foreach my $gene (@genes) {
     print STDERR "Slice : Got gene of type " . $gene->type() . "\n";
     if($gene->type() eq $type) {
       push @out, $gene;
     }
   }

   return @out;
}


=head2 invert

 Title   : invert
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub invert{
   my ($self,@args) = @_;

   $self->throw("Ewan has not implemented this function! Complain!!!!");

}


=head2 primary_seq

 Title   : primary_seq
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub primary_seq{
   my ($self,@args) = @_;

    return $self;
}


=head2 convert_Gene_to_raw_contig

 Title   : convert_Gene_to_raw_contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub convert_Gene_to_raw_contig{
  my ($self,$gene) = @_;

  if(!$gene->isa("Bio::EnsEMBL::Gene")){
    $self->throw("trying to use the wrong method can called convert gene to RawContig coords on ".$gene."\n");
  }
     
  $gene->transform;
  
  return $gene;
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


=head2 adaptor

 Title   : adaptor
 Usage   : $obj->adaptor($newval)
 Function: 
 Example : 
 Returns : value of adaptor
 Args    : newvalue (optional)


=cut

sub adaptor{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'adaptor'} = $value;
    }
    return $self->{'adaptor'};

}


=head2 dbID

  Arg [1]   : int databaseInternalId
              A slice might exist in the database and will than have this
              internal id.
  Function  : attribute function
  Returntype: int
  Exceptions: none
  Caller    : DBSQL::SliceAdaptor

=cut

sub dbID {
   my ( $self, $value ) = @_;
   if( defined $value ) {
     $self->{'dbID'} = $value;
   }
   return $self->{'dbID'};
}

sub id {
   my ( $self, $value ) = @_;
   if( defined $value ) {
     $self->{'id'} = $value;
   }
   return $self->{'id'};
}

sub display_id{
  my ( $self, $value ) = @_;
  if( defined $value ) {
    $self->{'display_id'} = $value;
  }
  return $self->{'display_id'};
}

sub desc{
  my ( $self, $value ) = @_;
  if( defined $value ) {
    $self->{'desc'} = $value;
  }
  return $self->{'desc'};
}

sub moltype {
    my ($obj) = @_;
    return 'dna';
}

sub accession_number {
    my( $obj, $acc ) = @_;

    if (defined $acc) {
        $obj->{'accession_number'} = $acc;
    } else {
        $acc = $obj->{'accession_number'};
        $acc = 'unknown' unless defined $acc;
    }
    return $acc;
}



sub get_KaryotypeBands() {
  my ($self) = @_;
  
  my $kadp = $self->adaptor->db->get_KaryotypeBandAdaptor();
  my @bands = $kadp->fetch_by_chromosome_start_end($self->chr_name(),
						   $self->chr_start(),
						   $self->chr_end());
  return @bands; 
}


sub get_Chromosome {
  my $self = shift @_;

  my $ca =  $self->adaptor->db->get_ChromosomeAdaptor();


  return $ca->fetch_by_chr_name($self->chr_name());
}

sub get_repeatmasked_seq {
    my ($self) = @_;

    $self->warn("Slice: get_repeatmasked_seq\n");

    my @repeats = $self->get_all_RepeatFeatures();
    my $dna = $self->seq();
    my $masked_dna = $self->mask_features($dna, @repeats);
    my $masked_seq = Bio::PrimarySeq->new(   '-seq'        => $masked_dna,
                                             '-display_id' => $self->id,
                                             '-primary_id' => $self->id,
                                             '-moltype' => 'dna',
					     );
    return $masked_seq;
}


sub mask_features {
    my ($self, $dnastr,@repeats) = @_;

   $self->warn("Slice: mask_features\n");


    my $dnalen = length($dnastr);
    #print "there are ".@repeats."\n";
  REP:foreach my $f (@repeats) {
      #print $f->start." ".$f->end." ".$f->repeat_id."\n";
      my $start    = $f->start;
      my $end	   = $f->end;
      my $length = ($end - $start) + 1;
      
      if ($start < 0 || $start > $dnalen || $end < 0 || $end > $dnalen) {
	  print STDERR "Eeek! Coordinate mismatch - start $start or  end $end not within $dnalen\n";
	  next REP;
      }
      
      $start--;
      
      my $padstr = 'N' x $length;
      
      substr ($dnastr,$start,$length) = $padstr;
  }
    return $dnastr;
} 


sub length {
  my ($self) = @_;

  return $self->chr_end() - $self->chr_start() + 1;
}


sub get_all_MapFrags {
    my $self = shift;
    my $mapset = shift;
    return $self->adaptor->db->get_MapFragAdaptor->fetch_mapset_chr_start_end( 
        $mapset, $self->chr_name, $self->chr_start, $self->chr_end
    );
}    

sub has_MapSet {
    my( $self, $mapset_name ) = @_;
    return $self->adaptor()->db()->get_MapFragAdaptor->has_mapset( $mapset_name );
}


sub get_tiling_path {
  my ($self) = @_;

  my $mapper = $self->adaptor()->db->get_AssemblyMapperAdaptor()->
    fetch_by_type($self->assembly_type());


  # Get the ids of the raw_contigs in this region specified in chrmsml coords 
  $mapper->register_region( $self->chr_name, $self->chr_start(),
			     $self->chr_end() );
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
  my $raw_contigs = 
    $self->adaptor->db->get_RawContigAdaptor()->
      fetch_filled_by_dbIDs(@raw_contig_ids);

  my @tiling_path;
  my $current_start = 1;

  foreach my $coord ( @mapped ) {
    my $length = $coord->end() - $coord->start() + 1; 

    if ( $coord->isa("Bio::EnsEMBL::Mapper::Coordinate" ) ) {
      # this is a contig, create a tiling path piece from it
      my $tile = {};
      $tile->{'start'} = $current_start;
      $tile->{'end'} = ($current_start + $length-1);
      $tile->{'contig'} = $raw_contigs->{ $coord->id() };
      $tile->{'strand'} = $coord->strand();
      
      $current_start += $length;

      push(@tiling_path, $tile);
    } else {
      # this is a gap, just add the length and discard it
      $current_start += $length;
    }
  }
  return @tiling_path;
}
  

sub get_landmark_MarkerFeatures {
  my $self = shift;

  my $lma = $self->adaptor()->db()->get_LandmarkMarkerAdaptor();
  if( ! defined $lma ) {
    return ();
  } else {
    return $lma->fetch_by_Slice( $self );
  }

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

   my $mapper = $self->adaptor()->db->get_AssemblyMapperAdaptor()->
     fetch_by_type($self->assembly_type());

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

   foreach my $extf ( $self->adaptor()->db()->_each_DASFeatureFactory ) {
       
       if( $extf->can('get_Ensembl_SeqFeatures_DAS') ) {
	       foreach my $sf (
                $extf->get_Ensembl_SeqFeatures_DAS(
                    $self->chr_name,$self->chr_start,$self->chr_end,
                    \@fpccontigs, \@clones,\@raw_contig_names, $chr_length)
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
#                    if(my $contig_from_clone = $self->contig_from_clone($sf->seqname()) ) {
#                        print STDERR "CONTIG NAME FROM CLONE >$contig_from_clone<\n";
#                        $sf->seqname($contig_from_clone);
# 	                    push(@contig_features, $sf);
#                    }
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
	        warn "Slice: Got a DAS feature factory that can't do get_Ensembl_SeqFeatures_DAS\n";
	        #$self->throw("Got a DAS feature factory that can't do get_Ensembl_SeqFeatures_DAS");
       }
   }
   
   my $chr_start = $self->chr_start();
   my $chr_end   = $self->chr_end();
   
   foreach my $sf ( @contig_features ) {

#            print STDERR "SEG ID: ",         $sf->seqname(), "\t";
#            print STDERR "ID: ",             $sf->das_id(), "\t";
#            print STDERR "DSN: ",            $sf->das_dsn(), "\t";
#            print STDERR "FEATURE START: ",  $sf->das_start(), "\t";
#            print STDERR "FEATURE END: ",    $sf->das_end(), "\t";
#            print STDERR "FEATURE STRAND: ", $sf->das_strand(), "\t";
#            print STDERR "FEATURE TYPE: ",   $sf->das_type_id(), "\n";
     
     # map to a chromosomal coordinate
     my @coord_list = $mapper->map_coordinates_to_assembly
       ( $contig_name_hash{ $sf->seqname() }->dbID(), 
	 $sf->das_start, $sf->das_end, $sf->das_strand );
     # this should work with one coordinate
     my $coord = shift( @coord_list );
     
     # if its not mappable than ignore the feature
     if( $coord->isa( "Bio::EnsEMBL::Mapper::Gap" )) {
       next;
     }
     
     $sf->das_start( $coord->start() - $self->chr_start + 1 );
     $sf->das_end( $coord->end() - $self->chr_start + 1 );
     $sf->das_strand( $coord->strand() * $self->strand() );

     if($sf->das_start <= $self->length && $sf->das_end >= 1) {
       push(@genomic_features, $sf);
     }

   }
   
   

   foreach my $sf ( @chr_features ) {
     
     # chromosome to slice coords mapping
     $sf->das_start( $sf->das_start() - $self->chr_start() + 1 );
     $sf->das_end( $sf->das_end() - $self->chr_start + 1 );
     $sf->das_strand( $sf->das_strand() * $self->strand() );

     if($sf->das_start <= $self->length && $sf->das_end >= 1) {
       push(@genomic_features, $sf);
     }
   }
   $self->{'_das_cached_features'} = \@genomic_features;

   return @genomic_features;
}









=head2 Backward Compatibility functions

=cut

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

   my ($p,$f,$l) = caller;
   $self->warn("$f:$l get_all_Genes_exononly has been deprecated. get_all_Genes called");

   return $self->get_all_Genes();
}


sub get_all_SangerGenes_startend_lite {
  my $self = shift;

  $self->warn("Slice->get_all_SangerGenes_startend_lite deprecated" . 
	      " use get_allGenes() instead\n");
  
  return $self->get_Genes_by_source('sanger');
}
  
sub get_all_VirtualGenes_startend_lite {
  my $self = shift;

  $self->warn("Slice->get_all_VirtualGenes_startend_lite deprecated" .
	      " use get_all_Genes() instead\n");

  return $self->get_all_Genes();
}


sub get_all_EMBLGenes_startend_lite {
  my $self = shift;

  $self->warn("Slice->get_all_EMBLGenes_startend_lite deprecated" .
	      " use get_Genes_by_source() instead\n");

  return $self->get_Genes_by_source('embl');
}



sub fetch_chromosome_length {
  my ($self) = @_;

  $self->warn( "Call to deprecated method fetch_chromosome_length\n" .
	       "use \$slice->get_Chromosome()->length(); instead.\n" .
	       $self->stack_trace_dump());

  return $self->get_Chromosome()->length();
}


        
sub fetch_karyotype_band_start_end {
   my ($self,@args) = @_;

   $self->warn( "Call to deprecated method fetch_karyotype_band_start_end\n" .
		"use \$slice->get_KaryotypeBands(); instead.\n" .
	       $self->stack_trace_dump());

   return $self->get_KaryotypeBands();
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
   my ($self) = @_;

  $self->warn("Slice->get_all_ExternalFeatures is not currently implemented, and will either be implemented or removed completely in the future\n");
  
   return ();

   #return $self->_get_all_SeqFeatures_type('external');
}

=head2 _get_all_SeqFeatures_type

 Title   : _get_all_SeqFeatures_type
 Usage   : Internal function which encapsulates getting
           features of a particular type and returning
           them in the slice coordinates.
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _get_all_SeqFeatures_type {
   my ($self,$type) = @_;
   $self->throw('interface fault') if @_ != 2;

  $self->warn("Slice->get_all_SeqFeatures_type is not currently correctly implemented.  This may arrive in the future, or may be removed completely");
 
   return ();
 
#   my $mapper = $self->adaptor->db->get_AssemblyMapperAdaptor->fetch_by_type
#                                    ( $self->assembly_type );

#   # register the VC
#   $mapper->register_region( $self->chr_name,
#			     $self->chr_start,
#			     $self->chr_end );
  
#   # get contig IDs for the VC
#   my @cids = $mapper->list_contig_ids( $self->chr_name,
#				        $self->chr_start,
#				        $self->chr_end );
   
#   my $rca = $self->adaptor->db->get_RawContigAdaptor;
#   my @vcsf = ();
#   foreach my $id (@cids) {
#     my $c = $rca->fetch_by_dbID($id);

#     # get start and end of the golden path fragment of the raw contig
#     my ($chr, $start, $end) = $self->adaptor->get_chr_start_end_of_contig(
#                                                                  $c->name);

#     if ($type eq 'external') {
#       foreach my $f ($c->get_all_ExternalFeatures) {
#         my @feature_mapped_to_assembly = $mapper->map_coordinates_to_assembly
#                         ($id, $f->start, $f->end, $f->strand);
#         if($feature_mapped_to_assembly[0]->isa("Bio::EnsEMBL::Mapper::Gap")) {
#           next;
#	 }
#	 my $newstrand = $feature_mapped_to_assembly[0]->strand
#	                                       * $self->strand;
#	 my $newstart = $feature_mapped_to_assembly[0]->start
#	                                  - $self->chr_start + 1;
#	 my $newend = $newstart + $f->end - $f->start;
#	 my $newf = Bio::EnsEMBL::SeqFeature->new();
#	 %$newf = %$f;
#	 $newf->start($newstart);
#         $newf->end($newend);
#         $newf->strand($newstrand);
#	 push @vcsf, $newf;
#       }
#     } else {
#       $self->throw("Type $type not recognised");
#     }
#   }

#   return @vcsf;
}



=head2 get_all_SimilarityFeatures_above_score

  Args      : $logic_name, $score
  Function  : Retrieves all of the DnaDnaAlignFeatures and all of the
              ProteinDnaAlignFeatures on this slice.
  Returntype: @Bio::EnsEMBL::BaseAlignFeature
  Exceptions: none
  Caller    : general

=cut

sub get_all_SimilarityFeatures_above_score {
  my ($self, $logic_name, $score) = @_;

  #
  # To deprecate, or not deprecate...
  # It seems to me that this function isn't very useful in that it's task
  # could be performed much faster through a call to either
  # get_DnaAlignFeatures_above_score or
  # get_ProteinAlignFeatures_above_score.  
  # Reminder Warning follows:
  #
  $self->warn("Call to Slice->get_all_SimilarityFeatures_above_score." .
	     "logic name = $logic_name.  Should this be deprecated??");

  my @prot_feats = 
    $self->get_all_ProteinAlignFeatures_above_score($logic_name, $score);
  my @dna_feats = 
    $self->get_all_DnaAlignFeatures_above_score($logic_name, $score);

  return (@prot_feats, @dna_feats);
}


1;
