
#
# BioPerl module for Contig
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::AceDB::Contig - Handle onto a database stored contig

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Handle onto a database stored contig for use with camace.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::AceDB::Contig;
use vars qw(@ISA);
use strict;

use Bio::Seq;
use Bio::Root::Object;
use Bio::EnsEMBL::AceDB::Obj;
use Bio::EnsEMBL::ContigOverlapHelper;
use Bio::EnsEMBL::FeatureFactory;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Repeat;
use Bio::EnsEMBL::DB::RawContigI;
use Bio::SeqFeature::Generic;

# Object preamble - inheriets from Bio::Root::Object



@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::RawContigI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);
  my ($dbobj, $id, $clone) = $self->_rearrange([qw(DBOBJ
					  ID CLONE
					  )],@args);

  $id || $self->throw("Cannot make AceDB Contig object without id");
  $dbobj || $self->throw("Cannot make AceDB Contig object without db object");
  $clone || $self->throw("Cannot make AceDB Contig object without clone");
  $dbobj->isa('Bio::EnsEMBL::AceDB::Obj') || $self->throw("Cannot make contig db object with a $dbobj object");  
    # set stuff in self from @args
  $self->id($id); 
  $self->dbobj($dbobj);
  $self->clone($clone);
  return $make; # success - we hope!
}


=head2 chromosome

 Title   : chromosome
 Usage   : $obj->chromosome()
 Function: 
 Example : 
 Returns : The chromosome object
 Args    : 


=cut

sub chromosome {
   my ($self) = @_;
   
    ## Get the sequence object
    #my $seq = $self->ace_seq();
    #my $chromosome_name = $seq->at('Origin.Chromosome[1]');  
     
    # Create a new human Species object
    my $species = new Bio::EnsEMBL::Species("human");
    return $species->chromosome_by_name("human", "unknown");
}


=head2 dbobj

 Title   : dbobj
 Usage   : $obj->dbobj($newval)
 Function: 
 Example : 
 Returns : value of _dbobj
 Args    : newvalue (optional)


=cut

sub dbobj {
   my ($self,$value) = @_;
   if (defined $value) {
      $self->{'_dbobj'} = $value;
    }
    return $self->{'_dbobj'};
}


=head2 embl_offset

 Title   : embl_offset
 Usage   : $offset = $contig->embl_offset()
 Function: Returns the embl offset which is always 1 for AceDB cotigs
 Example :
 Returns : integer
 Args    : none


=cut

sub embl_offset {
   my ($self) = @_;

   return 1;
}


=head2 embl_order

 Title   : embl_order
 Usage   : $order = $contig->embl_order()
 Function: Provides the order of the contig, starting at
         : zero.
 Example :
 Returns : integer
 Args    : none


=cut

sub embl_order {
   my ($self) = @_;

   # All AceDB contigs are single contigs in a single clone.
   return 0; 

}   


=head2 embl_accession

 Title   : embl_accession
 Usage   : $accession = $contig->embl_accession()
 Function: Returns the embl accession
 Example :
 Returns : 
 Args    : 


=cut

sub embl_accession {
   my ($self) = @_;

    unless (defined $self->{'embl_accession'}) {
        my $seq = $self->ace_seq();
        $self->{'embl_accession'} = $seq->at('DB_info.Database.EMBL[2]')->name;
    }
    return $self->{'embl_accession'};

} 


=head2 embl_id

 Title   : embl_id
 Usage   : $id = $contig->embl_id()
 Function: Returns the embl accession
 Example :
 Returns : 
 Args    : 


=cut

sub embl_id {
   my ($self) = @_;

    unless (defined $self->{'embl_id'}) {
        my $seq = $self->ace_seq();
        $self->{'embl_id'} = $seq->at('DB_info.Database.EMBL[1]')->name;
    }
    return $self->{'embl_id'};

} 


=head2 get_all_ExternalFeatures

 Title   : get_all_ExternalFeatures
 Usage   : foreach my $sf ( $contig->get_all_ExternalFeatures )
 Function: Gets all the external features on a contig which is just an 
            empty list for AceDB 
 Example :
 Returns : Array of Bio::EnsEMBL::FeaturePair objects
 Args    : none


=cut


sub get_all_ExternalFeatures {
    my ($self) = @_;
    return;     
}


=head2 get_all_RepeatFeatures

 Title   : get_all_RepeatFeatures
 Usage   : foreach my $sf ( $contig->get_all_RepeatFeatures )
 Function: Gets all the motif homols found by the RepeatMaster and RepeatMaster_SINE methods
            and the tandem repeat features for the contig.
 Example :
 Returns : Array of Bio::EnsEMBL::FeaturePair objects
 Args    : none


=cut

sub get_all_RepeatFeatures {
    my ($self) = @_;

    # Get tandem features
    print STDERR "About get tandems ...\n";
    my @repeat_features = $self->_get_tandems();
    print STDERR "About to get motifs...\n";

    # Append motif_homols with RepeatMaster, RepeatMaster_SINE, hmmfs.3 or scan
    my @methods = qw/RepeatMaster RepeatMaster_SINE hmmfs.3 scan/;
    push(@repeat_features, $self->_get_homols('Homol.Motif_homol[1]', @methods));
     
    print STDERR "otta here!\n";

    # Return all the sequence features 
    return @repeat_features;
}

 
=head2 get_all_PredictionFeatures

 Title   : get_all_PredictionFeatures
 Usage   : foreach my $sf ( $contig->get_all_PredictionFeatures ) 
 Function: 
 Example :
 Returns : Array of Bio::EnsEMBL::Feature objects
 Args    :


=cut

sub get_all_PredictionFeatures {
    my ($self) = @_;
               
    # Create a hash of methods we're interested in
    my %methods = map{$_, 1} ('FGENES', 'FGENESH', 'GENESCAN');
    
    my $id = $self->id();        
    my $seq = $self->ace_seq();
    my $phase = $self->_get_phase();
    my @prediction_features;
    
    # Loop through the subsequences
    foreach my $sub ($seq->at('Structure.Subsequence')) {
        
        # Fetch the method and check we're interested in this subsequence 
        my $method = map($_->name, $sub->fetch->at("Method[1]") );
                  
        if ($methods{ $method }) {
             
            my ($start, $end) = map($_->name(), $sub->row(1));
            my $strand = ($start < $end) ? 1 : -1;
            my $seqname = $sub->name;

            my $analysis = new Bio::EnsEMBL::Analysis(
                                            -db_version      => "NULL",					    
					    -program         => $method,
					    -program_version => 1,
					    -gff_source      => $method,
					    -gff_feature     => "exon"
                                            );    
             
            my $feature = Bio::EnsEMBL::FeatureFactory->new_feature();
            $feature->start($start);
            $feature->end($end);
            $feature->strand($strand);
            $feature->seqname($id);
            $feature->score(0);
            $feature->source_tag($method);
            $feature->primary_tag("Tag");
            $feature->analysis($analysis);
             
            
            # Fetch all the exons on the sequence, create subfeatures and add them to the feature 
            my $subseq = $sub->fetch();           
            foreach my $exon ($subseq->at('Structure.From.Source_Exons[1]')) {  
                                         
	        my ($starte, $ende) = map($_->name(), $exon->row());
 
                my $subFeature = Bio::EnsEMBL::FeatureFactory->new_feature();
                
                # We have to map acedb coordinates which are relative to the
                # start/end in the subsequence to the exon coordinates, which
                # are absolute.
                
                if( $strand == 1 ) {
                        $subFeature->start($start + $starte - 1);
                        $subFeature->end($start + $ende - 1);
                }
                else {
                        $subFeature->start($start - $ende + 1);
                        $subFeature->end($start - $starte + 1);
                }                
                
                $subFeature->strand($strand);
                $subFeature->seqname($id);
                $subFeature->score(0);
                $subFeature->source_tag($method);
                $subFeature->primary_tag("similarity");
                $subFeature->analysis($analysis);
                
                $feature->add_sub_SeqFeature($subFeature);    
            }            
            
            $feature->validate;
            # Add the feature to the prediction_features array
            push(@prediction_features, $feature);
        }
    }

    return @prediction_features; 
} 
      
             
=head2 get_all_SeqFeatures

 Title   : get_all_SeqFeatures
 Usage   : foreach my $sf ( $contig->get_all_SeqFeatures ) 
 Function: Finds all the sequence features for the contig by calling get_all_RepeatFeatures(), 
            get_all_ExternalFeatures() and get_all_SimilarityFeatures().
 Example :
 Returns : Array of Bio::EnsEMBL::FeaturePair objects
 Args    :


=cut

sub get_all_SeqFeatures {
    my ($self,@args) = @_;
           
    # Get repeat, external similarity and prediction features
    my @seq_features = $self->get_all_RepeatFeatures();
    push(@seq_features, $self->get_all_ExternalFeatures());
    push(@seq_features, $self->get_all_SimilarityFeatures());
    push(@seq_features, $self->get_all_PredictionFeatures());
#    print(STDERR "Fetched all features\n");
    
    # Return all the sequence features
    return @seq_features; 
}


=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   : @genes = $contig->get_all_Genes
 Function:
 Example : 
 Returns : 
 Args    : none
 Note    : WARNING. At the moment this just
           gets genes from the Sequence object, not
           any potential genes in link objects above
           this sequence object. To be fixed!

=cut  

sub get_all_Genes {
    my ($self) = @_;
    
    # Create a hash of methods we're interested in
    my %methods = map{$_, 1} ('supported_CDS', 'curated', 'GD_in_progress');
    
    my $id = $self->id();        
    my $seq = $self->ace_seq();
    my $phase = $self->_get_phase();
    my $annotated = map($_->name, $seq->at('Properties.Status.Annotated'));
    my @genes;
    
    # Loop through the subsequences
    foreach my $sub ($seq->at('Structure.Subsequence')) {
    
        my $genename = "$sub";        
        if ($genename =~ /^.TRIMMED/) {
            next;
        }
        
        # Fetch the method and check we're interested in it            
        if ($methods{$sub->fetch->at("Method[1]")}) {
             
                                              
            my ($start, $end) = map($_->name(), $sub->row(1));
            my $strand = ($start < $end) ? 1 : -1;
            my $subseq = $sub->fetch();             
            my $index = 1;
            my @exons;
            
            # Fetch all the exons            
            foreach my $hit ($subseq->at('Structure.From.Source_Exons[1]')) {  
                                         
	        my ($starte, $ende) = map($_->name(), $hit->row());                              	                        
                my $exon = $self->_create_exon($id, $strand, $start, 
                    $starte, $end, $ende, $genename, $index, $phase, $annotated); 
                
                # Set index and phase for the next exon
                $index++;                
	        $phase = $exon->end_phase(); 
                               	                        
                push(@exons, $exon);
            }            
            
            # Create a new Translation object with the start and end exon IDs
            my $translation = new Bio::EnsEMBL::Translation();
           
            $translation->id($genename); 
            $translation->version(1);
            
            if ($strand) {                   
	        $translation->start($start);
	        $translation->end($end);
                $translation->start_exon_id($exons[0]->id());
	        $translation->end_exon_id($exons[$#exons]->id());
            }
            else {            
                $translation->start($end);
	        $translation->end($start);
                $translation->start_exon_id($exons[$#exons]->id());
	        $translation->end_exon_id($exons[0]->id());
            } 
           
	    # Create a new Transcript object from the exons and the Translation
            my $transcript = new Bio::EnsEMBL::Transcript(@exons);
            $transcript->id($genename); 
            $transcript->version(1);   
	    $transcript->translation($translation);

            # Create a new Gene object and add the Transcript 
            my $gene = new Bio::EnsEMBL::Gene();
            $gene->id($genename);
            $gene->version(1); 
            $gene->add_Transcript($transcript);
            $gene->add_cloneid_neighbourhood($self->id);
            
            # Add the gene to the genes array
            push(@genes, $gene);
        }
    }

    # Return all the genes
    return @genes;
}


=head2 clone

 Title   : clone
 Usage   : $obj->clone($newval)
 Function: 
 Example : 
 Returns : value of clone
 Args    : newvalue (optional)


=cut

sub clone {
    my ($self,$value) = @_;
    if (defined $value) {
	$self->{'clone'} = $value;
    }
    return $self->{'clone'};
}


=head2 get_all_SimilarityFeatures

 Title   : get_all_SimilarityFeatures
 Usage   : foreach my $sf ( $contig->get_all_SimilarityFeatures )
 Function: Gets all the DNA and Pep homols for the contig.
 Example :
 Returns : Array of Bio::EnsEMBL::FeaturePair objects
 Args    : none


=cut

sub get_all_SimilarityFeatures {
     my ($self) = @_;
 
    # Get DNA, Pep, EST and STS homols
    my @similarity_features = $self->_get_homols('Homol.DNA_homol[1]');
    push (@similarity_features, $self->_get_homols('Homol.Pep_homol[1]'));
    push (@similarity_features, $self->_get_homols('Homol.EST_homol[1]'));
    push (@similarity_features, $self->_get_homols('Homol.STS_homol[1]'));
    
    # Append motif_homols with Queryprosite
    my @methods = qw/Queryprosite/;
    push(@similarity_features, $self->_get_homols('Homol.Motif_homol[1]', @methods));
     
    # Return all the sequence features 
    return @similarity_features;    
}


=head2 get_left_overlap

 Title   : get_left_overlap
 Usage   : $overlap = $contig->get_left_overlap()
 Function: Checks whether left_overlap has already been defined.
            If it hasn't it calls _get_overlap with the parameter "left"
 Example :
 Returns : Bio::EnsEMBL::ContigOverlapHelper
 Args    : none


=cut

sub get_left_overlap {
    my ($self) = @_;
    
    unless (defined $self->{'left_overlap'}) {
	$self->{'left_overlap'} = $self->_create_overlap("left");
    }
    return $self->{'left_overlap'};
}


=head2 get_right_overlap

 Title   : get_right_overlap
 Usage   : $overlap = $contig->get_right_overlap()
 Function: Checks whether right_overlap has already been defined.
            If it hasn't it calls _get_overlap with the parameter "right"
 Example :
 Returns : Bio::EnsEMBL::ContigOverlapHelper
 Args    : none


=cut

sub get_right_overlap {   
    my ($self) = @_;
    
    unless (defined $self->{'right_overlap'}) {
	$self->{'right_overlap'} = $self->_create_overlap("right");
    }
    return $self->{'right_overlap'};
}


=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Example : 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id {
    my ($self,$value) = @_;
    if (defined $value) {
	$self->{'id'} = $value;
    }
    return $self->{'id'};
}

 
=head2 internal_id

 Title   : internal_id
 Usage   : $obj->internal_id($newval)
 Function: 
 Example : 
 Returns : value of internal_id
 Args    : newvalue (optional)


=cut

sub internal_id {
    my ($self,$value) = @_;

    return $self->id($value);
}
 
  
=head2 length

 Title   : length
 Usage   : $len = $contig->length
 Function: Provides the length of the contig
 Example :
 Returns : integer
 Args    : none


=cut

sub length {
    my ($self) = @_;
    return $self->primary_seq->length;
}


=head2 orientation

 Title   : orientation
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub orientation {
   my ($self,@args) = @_;   
   return 1;
}


=head2 primary_seq

 Title   : primary_seq
 Usage   : $seq = $contig->primary_seq();
 Function: Gets a Bio::PrimarySeq object out from the contig
 Example :
 Returns : Bio::PrimarySeq object
 Args    : None


=cut

sub primary_seq {
   my ($self) = @_;
   my $id = $self->id();

    unless ($self->{'_primary_seq'}) {
       my $seq = $self->ace_seq;

       my $dna = $seq->asDNA || $self->throw("Could not retrieve DNA from $id");

       $dna =~ s/^>.*\n//g;
       $dna =~ s/\s//g;
       $dna =~ tr/[a-z]/[A-Z]/;
       $self->{'_primary_seq'} = Bio::PrimarySeq->new ( -seq => $dna , '-id' => $id, -type => 'DNA' ) ;
    }
   
   return $self->{'_primary_seq'};
}


=head2 seq_date

 Title   : seq_date
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub seq_date {
    my ($self) = @_;
    if (my $date = $self->ace_seq->at('Properties.Status.Finished[1]')) {
        return $date->name;
    }
    return 0; 
}

=head2 seq_version

 Title   : seq_version
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub seq_version {
    my ($self) = @_;
    if (my $version = $self->ace_seq->at('DB_info.Sequence_version[1]')) {
        return $version->name;
    }
    return 0; 
}


=head2 ace_seq

 Title   : ace_seq
 Usage   : $seq = $contig->acce_seq();
 Function: 
 Example :
 Returns : The Ace::Object for this sequence
 Args    : None

=cut

sub ace_seq {
   my ($self) = @_;
   
   unless($self->{'_ace_seq'}) {
        my $id = $self->id();
        $self->{'_ace_seq'} = $self->dbobj->fetch('Genome_sequence', $id) || 
           $self->throw("Could not retrieve $id from acedb" . Ace->error());
   }
   return $self->{'_ace_seq'};
}


=head2 _create_exon

 Title   : _create_exon
 Usage   : 
 Function: 
 Example :
 Returns : 
 Args    : 

=cut

sub _create_exon {
    my ($self, $id, $strand, $start, $starte, $end, $ende, $genename, $index, $phase, $annotated) = @_;
     
    my $bioseq = $self->primary_seq();

    my $exon = new Bio::EnsEMBL::Exon();        
    $exon->clone_id($id);
    $exon->contig_id($id);

    $exon->strand($strand);
    $exon->phase($phase);
    $exon->seqname($genename);
    $exon->version(1);
    $exon->modified($annotated);
    $exon->created($annotated);
    
    # We have to map acedb coordinates which are relative to the
    # start/end in the subsequence to the exon coordinates, which
    # are absolute.
    if( $strand == 1 ) {
            $exon->start($start+$starte-1);
            $exon->end($start+$ende-1);
    }
    else {
	    $exon->start($start-$ende+1);
	    $exon->end($start-$starte+1);
    } 
    $exon->attach_seq($bioseq);
#    print STDERR "Exon " . $exon->start . "\t" . $exon-> end . "\t" . $exon->seq->seq . "\n";
    my $exonid;
    if( $self->dbobj->_exon_id_start() ) {
	     $exonid = $self->dbobj->_exon_id_start();
	     my $nexte = $exonid++;
	     $self->dbobj->_exon_id_start($exonid);
    } 
    else {
	     $exonid = "dummy_exon_id.$genename.$index";
    }
    $exon->id($exonid);
    
   
    return $exon;
}


=head2 _create_feature_pair

 Title   : _create_feature_pair
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _create_feature_pair {
    my ($self, $meth, $score, $start, $end, $hstart, $hend, $hid) = @_;        
    my $strand = 1;
    
    if ($start > $end) {
        ($start, $end) = ($end, $start);
        $strand = -1;        
    }         
    if ($hstart > $hend) {                
        ($hstart, $hend) = ($hend, $hstart); 
         $strand = -$strand;
    }
 
    my ($program, $source);   
    if ($meth =~ /^TBLASTX/) {
        $program = "TBLASTX";
        $source = substr($meth, 8);     
    }
    else {
        $program = $source = $meth;
    }   

    my $feature = Bio::EnsEMBL::FeatureFactory->new_feature_pair();   
    # Set the fields with the set_all_fields method of FeaturePair 
    $feature->set_all_fields($start,      # start
                            $end,         # end,
                            $strand,      # strand,
                            $score,       # score,
                            $source,      # source,
                            "similarity", # primary,
                            $self->id,    # seqname,
                            $hstart,      # hstart,
                            $hend,        # hend,
                            1,            # hstrand is always 1
                            $score,       # hscore is the same as score
                            $source,      # hsource is the same as source
                            "similarity", # hprimary,
                            $hid);        # hseqname
                            
    my $analysis = new Bio::EnsEMBL::Analysis(					    
					    -program         => $program,
					    -program_version => 1,
					    -gff_source      => $source,
					    -gff_feature     => "similarity"
                                            );
    $feature->analysis($analysis);                            
    $feature->validate();
    return $feature;
}


=head2 _create_repeat

 Title   : _create_repeat
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _create_repeat {
    my ($self, $meth, $score, $start, $end, $hstart, $hend, $hid) = @_;        
    my $strand = 1;
     
    if ($start > $end) {
        ($start, $end) = ($end, $start);
        $strand = -1;        
    }         
    if ($hstart > $hend) {                
        ($hstart, $hend) = ($hend, $hstart); 
         $strand = -$strand;
    }
    
    # Create a new Bio::EnsEMBL::Repeat 
    my $repeat = Bio::EnsEMBL::FeatureFactory->new_repeat();
   
    # Set the Repeat features with the set_all_fields method of FeaturePair 
    $repeat->set_all_fields($start,       # start
                            $end,         # end,
                            $strand,      # strand,
                            $score,       # score,
                            $meth,        # source,
                            "similarity", # primary,
                            $self->id,    # seqname,
                            $hstart,      # hstart,
                            $hend,        # hend,
                            1,            # hstrand is always 1
                            $score,       # hscore is the same as score
                            $meth,        # hsource is the same as source
                            "similarity", # hprimary is the same as primary
                            $hid);        # hseqname
                            
    my $analysis = new Bio::EnsEMBL::Analysis(					    
					    -program         => $meth,
					    -program_version => 1,
					    -gff_source      => $meth,
					    -gff_feature     => "repeat"
                                            );
    $repeat->analysis($analysis);                            
    $repeat->validate();
    return $repeat;
}


=head2 _create_overlap

 Title   : _create_overlap
 Usage   :  $overlap = $contig->_create_overlap("left" or "right")
 Function: Called by get_left_overlap or get_right_overlap to search
            through the link the contig is part of and find the appropriate
            overlap which is returned as a ContigOverlapHelper
 Example :
 Returns : Bio::EnsEMBL::ContigOverlapHelper
 Args    : "left" or "right"


=cut

sub _create_overlap {   
    my ($self, $side) = @_;

    # Get the sequence object
    my $seq = $self->ace_seq();
    # Get the link object (chromosome) it's part of
    my $link = $seq->at('Structure.From.Source[1]');
    # If Structure.From.Source hasn't been defined return nothing
    return unless ($link);
    $link = $link->fetch;
    
    # Get an array of sequnces to be tested for overlaps
    my @testSequences = $link->at('Structure.Subsequence');    
    # Get a tag to this contig from the link
    my $seq_tag = $link->fetch->at("Structure.Subsequence.\Q$seq")
        or $self->throw("$seq is not in the link object $link");           
    # Find the start and end of this contig from this tag    
    my ($seq_start, $seq_end) = map $_->name, $seq_tag->row(1);
    unless (defined $seq_end) {
        $self->throw("End is not defined for contig $seq");
    }
    
    # Find the sequence polarity from the start and end
    my $seqPolarity = ($seq_start < $seq_end) ? 1 : -1;
    
    # Define the overlap point and the selfPosition depending on the side choosen 
    my $overlapPoint;
    my $selfPosition = 1;
    if ($side eq "left") {
        $overlapPoint = $seq_start;
    }
    elsif ($side eq "right") {
        $overlapPoint = $seq_end;
        $selfPosition += ($seq_end - $seq_start);
    }
    
    # Find the largest overlap of the overlapPoint from the testSequences    
    my ($overlapID, $overlapSize, $overlapPosition, $overlapPolarity) = 
        $self->_get_largest_overlap($seqPolarity, $overlapPoint, @testSequences);
    
    # If no overlap was found return undef to indicate there is no overlap
    unless (defined $overlapID) {
        return undef;
    }
    # Create a new Contig and return the ContigOverlapHelper
    my $overlapContig = new Bio::EnsEMBL::AceDB::Contig(-'dbobj' => $self->dbobj(),
                                                        -'id' => $overlapID,
                                                        -'clone' => $self->clone);     
    return new Bio::EnsEMBL::ContigOverlapHelper(-sister => $overlapContig, 
                                                -sisterposition => $overlapPosition, 
                                                -sisterpolarity => $overlapPolarity, 
                                                -source => $self, 
                                                -distance => 1,        # There is no distance between them
                                                -selfposition => $selfPosition);   
}


=head2 _get_phase

 Title   : _get_phase
 Usage   : $phase = get_phase()
 Function: Gets the phase of the contig from AceDB.
 Example :
 Returns : integer
 Args    : none


=cut

sub _get_phase {
    my ($self) = @_;
    
    # Get the sequence object
    my $seq = $self->ace_seq();
    
    my $phase = $seq->at('Properties.Coding.CDS[1]');
    
    # If the phase is defined set it to phase - 1 
    if ($phase) {
        return $phase->name() - 1;
    } 
    
    # Otherwise it's 0
    return 0;
}


=head2 _from_ace_seqfeature

 Title   : _from_ace_seqfeature
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _from_ace_seqfeature {
    my ($score,$start,$end,$tstart,$tend) = @_;
    my $strand;

    if ( $start !~ /^\d+$/ || $end !~ /^\d+$/ || $score !~ /^\d+\.?\d*$/ ) {
	&Bio::Root::Object::throw("start $start, end $end and score $score look dodgy");
    }

    if( $start > $end ) {
	my $temp =$end;
	$end = $start;
	$start = $temp;
	$strand = -1;
    } else {
	$strand = +1;
    }
    
    my $out = new Bio::SeqFeature::Generic;
    $out->start($start);
    $out->end($end);
    $out->strand($strand);
    $out->score($score);
    
    return $out;
}


=head2 _get_largest_overlap

 Title   : _get_largest_overlap
 Usage   :  $overlap = $contig->_get_largest_overlap()
 Function: Called by _create_overlap
 Example :
 Returns : ($overlapID, $overlapSize, $overlapPolarity)
 Args    : $seqPolarity, $overlapPoint, @testSequences


=cut

sub _get_largest_overlap {   
    my ($self, $seqPolarity, $overlapPoint, @testSequences) = @_;
       
    my $overlapSize = 0;   
    my $overlapPolarity;
    my $overlapPosition;
    my $overlapID;
    # Loop through the test sequences and find any that overlap the overlap point. 
    # If they do check that they are the largest one so far and if so keep.    
    foreach my $testSeq (@testSequences) {
    
        my ($start, $end) = map($_->name(), $testSeq->row(1));
        my $testPolarity = ($start < $end) ? 1 : -1;     
        # Test if the overlapPoint is between start and end for either polarity 
        if ((($testPolarity) && ($start < $overlapPoint) && ($overlapPoint < $end)) 
            || ((! $testPolarity) && ($end < $overlapPoint) && ($overlapPoint < $start))) {
        
            if ($end - $start > $overlapSize) {
                $overlapSize = $end - $start;
                $overlapPolarity = ($testPolarity == $seqPolarity) ? 1 : -1;
                $overlapID = $testSeq->name();
                $overlapPosition = ($testPolarity) ? $overlapPoint - $start : $end - $overlapPoint;
            }
        }
    }
    
    return ($overlapID, $overlapSize, $overlapPosition, $overlapPolarity); 
}


=head2 _get_homols

 Title   : _get_homols
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _get_homols {
    my ($self, $type, @test_methods) = @_;        
    my %reqired_methods = map {$_, 1} @test_methods;
    my @seq_features;
   
    # Loop through the homols of the appropriate type
    print STDERR "calling with $type\n";
    
    my @homols = $self->ace_seq->at($type);    
    foreach my $hid (@homols) {
        print STDERR "got a homol with ",$hid->name,"\n";
        
        # Loop through the different methods
        my @methods = $hid->col(1);        
        foreach my $meth (@methods) {
        
            # Test if the no test methods were given as parameters 
            # or if they were whether this is one of them.           
            if ((! @test_methods) || ($reqired_methods{$meth})) {
            
                # Loop through the scores
                foreach my $score ($meth->col(1)) {
                
                    # Loop through the column to the right of $score
                    foreach my $pos ($score->col(1)) {
                    
                        # Initalise $start etc from the row of data
                        my ($start, $end, $hstart, $hend) = $pos->row();
                         
                        # Create a repeat or a feature pair and add it to the array 
                        if (($type eq "Homol.Motif_homol[1]") && ($meth ne "Queryprosite")) {
                            my $repeat = $self->_create_repeat
                                ($meth, $score, $start, $end, $hstart, $hend, $hid);
                            push(@seq_features, $repeat);
                        } 
                        else {
                            my $feature = $self->_create_feature_pair
                                ($meth, $score, $start, $end, $hstart, $hend, $hid);
                            push(@seq_features, $feature);
                        }                                                                           
                    }
                }
            }
        }
    }
    
    # Return all the sequence features
    return @seq_features;
}


=head2 _get_tandems

 Title   : _get_tandems
 Usage   :
 Function:
 Example :
 Returns : array of Bio::EnsEMBL::Features
 Args    : 


=cut

sub _get_tandems {
    my ($self) = @_;
    my @seq_features;

    print STDERR "got into get_tandems\n";

    my @tandems = $self->ace_seq->at('Feature.Tandem[1]');

    my $analysis = new Bio::EnsEMBL::Analysis(					    
					    -program         => "tandem",
					    -program_version => 1,
					    -gff_source      => "tandem",
					    -gff_feature     => "repeat"
                                            );
    
    # Create a new feature object for each of the tandems found
    foreach my $tandem (@tandems) {
        print STDERR "Looking at a tandem...\n";

        my($start, $end, $score, $remark) = $tandem->row();               
        my $feature = Bio::EnsEMBL::FeatureFactory->new_repeat();
        
        my $strand = 1;
        if ($start > $end) {
            $strand = -1;
        }
        # Set the Repeat features with the set_all_fields method of FeaturePair 
        $feature->set_all_fields($start,       # start
                            $end,         # end,
                            $strand,      # strand,
                            $score,       # score,
                            "tandem",      # source,
                            "similarity",      # primary,
                            $self->id,    # seqname,
                            $start,      # hstart,
                            $end,        # hend,
                            -$strand,     # hstrand 
                            $score,       # hscore is the same as score
                            "tandem",      # hsource is the same as source
                            "similarity",      # hprimary is the same as primary
                            $self->id);        # hseqname        
        
        $feature->analysis($analysis);
        $feature->validate();
        push(@seq_features, $feature);
    }
    
    # Return all the sequence features
    return @seq_features;
} 

1;
