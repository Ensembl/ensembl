
#
# BioPerl module for Bio::EnsEMBL::Virtual::Contig.pm
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Virtual::Contig - View a region of a genome as continuous DNA with features

=head1 SYNOPSIS

    # get a raw contig in the genome somehow.
    $rc = $db->get_Contig('AP000012.00001');

    # make a virtual contig from this position
    $vc = Bio::EnsEMBL::Virtual::Contig->new (
					      -focuscontig => $rc,
					      -focusposition => 100,
					      -orientation => 1,
					      -left => 1000000,
					      -right => 100000
					      );

    # this is now a ContigI compliant object, so the following calls work

    # This will be a long line!
    print "Sequence is ",$vc->primary_seq->seq,"\n";

    # if one part of a gene is on this then the gene is returned
    $vc->get_all_Genes();
   
    # this gets genes, but with start/end points in this virtual contig
    $vc->get_all_VirtualGenes();

    # can get out sequence features
    $vc->get_all_RepeatFeatures();
    $vc->get_all_SimilarityFeatures();
    $vc->get_all_PredictionFeatures();

    # from a virtual contig you can build either translations or 
    # sub virtual contigs

    # moves the VC 100,000 bp to the left
    $new_vc = $vc->extend(-100000,-100000);

    # builds a smaller vc, focused at 200,000 in this VC coordinates with
    # 1000 on each side
    $windowed_vc = $vc->(200000,1000,1000);


=head1 DESCRIPTION

This object provides a key abstraction for Ensembl: it presents a
series of underlying overlapping and ordered fragments which make up a
genome as a single piece of DNA with features. All the sequence,
sequence features and genes which are stored in the underlying
coordinates on the fragments (called RawContigs in Ensembl
terminology) can be "viewed" on this object, with coordinates all making sense
internally (though of course, comparing things between virtualcontigs is
simply crazy)

NB: To display a possible misunderstanding here, the fragments or
RawContigs are the local assemblies made by, for example, Phrap of
clones. In the future, as things move to larger assemblies, the
RawContigs are likely to the unit used in the computational
procedures.

This object is where the Ensembl software provides programmers with a
crucial abstraction which saves most people an very frustrating
programming task

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Virtual::Contig;
use vars qw(@ISA);
use strict;
# ContigI uses RootI.
use Bio::EnsEMBL::DB::ContigI;
use Bio::EnsEMBL::Virtual::Map;
use Bio::EnsEMBL::Virtual::PrimarySeq;
use Bio::EnsEMBL::Utils::Eprof qw(eprof_start eprof_end);


@ISA = qw(Bio::EnsEMBL::DB::ContigI);

sub _make_datastructures {
    my $self = shift;

    # data structures for caching and coordinating
    # top_SeqFeature call with EMBL dumping (ie, repressing BLAST hits in EMBL)
    $self->{'_sf_cache'} = {};
    $self->{'_skip_feature'} = {};
    my $vmap=Bio::EnsEMBL::Virtual::Map->new();
    $self->_vmap($vmap);

    return $self;
}

=head2 new_from_one

 Title   : new_from_one
 Usage   :
 Function: makes a virtual contig out of a single contig, including the
           bits that are not on the GoldenPath
 Example :
 Returns : 
 Args    : a Contig (virtual or otherwise)

=cut

sub new_from_one {
    my ($class,$contig) = @_;

    my $self = {};
    bless $self,$class;
    $self->_make_datastructures();

    if (! $contig->isa('Bio::EnsEMBL::DB::ContigI') ) {
	$self->throw("$contig is not a Bio::EnsEMBL::DB::ContigI object, cannot make Virtual Contig!");
    }
    
    $self->_vmap->create_MapContig($contig,1,$contig->length,1,1);
    $self->id($contig->id);
    return $self;
}



=head2 new

 Title   : new
 Args    : Now does nothing excepts throws an exception

=cut

sub new {
    my ($class,@args) = @_;
    
    my $self = {};
    bless $self,$class;
    $self->_make_datastructures();
    $self->{'_unmapped_exons'}=[];
    $self->_old_exon_call(0);
    $self->throw("Bare new now useless. Should be static or something similar");

    return $self;
}

=head2 new_from_inverted

 Title   : new_from_inverted
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub new_from_inverted{
   my ($class,$vc) = @_;

   my $self = {};
   bless $self,$class;
   $self->_make_datastructures();

   if( !ref $vc || !$vc->isa("Bio::EnsEMBL::Virtual::Contig") ) {
       $self->throw("no virtual contig provided to new from inverted");
   }

   my $length = $vc->length;
   foreach my $mc ( $vc->_vmap->each_MapContig ) {
       # this looks too easy ;)
       $self->_vmap->create_MapContig($mc->contig,$length-$mc->end+1,$length-$mc->start+1,$mc->rawcontig_start,$mc->orientation*-1);
   }


   # these need to set separately to deal with overhangs
   $self->_vmap->length($length);
   $self->length($length);

   $self->id($vc->id.".inverted");

   return $self;
}

=head2 invert

 Title   : invert
 Usage   : $newvc = $vc->invert;
 Function: Makes a new vc which is the reverse complement of this
           vc. New vc is completely ok for coordinate transformation
           etc
 Example :
 Returns : 
 Args    :


=cut

sub invert{
   my ($self) = @_;

   return Bio::EnsEMBL::Virtual::Contig->new_from_inverted($self);
}

		 
		    
=head1 Virtual::Contig functions.

=head2 primary_seq

 Title   : primary_seq
 Usage   : $seq = $contig->primary_seq();
 Function: Gets a Bio::EnsEMBL::Virtual::PrimarySeq object out from the contig
 Returns : Bio::EnsEMBL::Virtual::PrimarySeq object
 Args    : none

=cut

sub primary_seq {
    my ($self) = @_;

    if( $self->{'_virtual_primary_seq'} ){
	return $self->{'_virtual_primary_seq'};
    }

    my $vseq = Bio::EnsEMBL::Virtual::PrimarySeq->new( -vmap =>$self->_vmap, 
						       -id =>$self->id						        );
    $self->{'_virtual_primary_seq'} = $vseq;
    return $vseq;
}

=head1 Functions implementing the Bio::SeqI interface inherieted by ContigI

=head2 id
    
 Title   : id
 Usage   : 
 Function: 
 Example :
 Returns : 
 Args    :


=cut

sub id {
    my $obj = shift;
    
    if( @_ ) {
	my $value = shift;
	$obj->{'id'} = $value;
    }
    return $obj->{'id'};
}

sub moltype { 
    return 'dna';
}

=head2 top_SeqFeatures
    
 Title   : top_SeqFeatures
 Usage   :
 function:
 Example :
 Returns : 
 Args    :


=cut
    
sub top_SeqFeatures {
    my ($self,@args) = @_;
    my (@f);

    
    if( !$self->skip_SeqFeature('similarity')  ) { 
	push(@f,$self->get_all_SimilarityFeatures());
    } 
    
    if( !$self->skip_SeqFeature('repeat')  ) { 
	push(@f,$self->get_all_RepeatFeatures());
    } 
    
    if( !$self->skip_SeqFeature('external')  ) { 
	push(@f,$self->get_all_ExternalFeatures());
    } 

    if ( !$self->skip_SeqFeature('prediction') ) {
	push(@f,$self->get_all_PredictionFeatures());
    } 

    
    if( !$self->skip_SeqFeature('gene') ) {
	foreach my $gene ( $self->get_all_Genes()) {
	    my $vg = Bio::EnsEMBL::VirtualGene->new(-gene => $gene,-contig => $self);
	    push(@f,$vg);
	}
    }

    if( !$self->skip_SeqFeature('contig') ) {
	push(@f,$self->_vmap->each_MapContig);
    }

    if( !$self->skip_SeqFeature('meta') ) {
	my $sf = Bio::SeqFeature::Generic->new();
	$sf->start(1);
	$sf->end($self->length());
	$sf->strand(1);
	$sf->primary_tag('source');
	$sf->add_tag_value('organism',$self->species->binomial);
	push(@f,$sf);
    }

    return @f;
}

=head2 length
    
 Title   : length
 Usage   : 
 Function: Provides the length of the contig
 Example :
 Returns : 
 Args    :


=cut

sub length {
    my $obj = shift;
    
    if( @_ ) {
	my $value = shift;
	$obj->{'length'} = $value;
    }
    return $obj->{'length'};
}


=head2 embl_accession

 Title   : embl_accession
 Usage   : $obj->embl_accession($newval)
 Function: 
 Returns : value of embl_accession
 Args    : newvalue (optional)


=cut

sub embl_accession {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'embl_accession'} = $value;
    }
    return $obj->{'embl_accession'};

}

=head1 Functions to get sequence features and genes

=head2 get_all_VirtualGenes
    
 Title   : get_all_VirtualGenes
 Usage   : foreach my $virtualgene ( $contig->get_all_VirtualGenes ) 
 Function: Gets all the genes on this VC as VirtualGene objects
 Example : 
 Returns : array of Bio::EnsEMBL::VirtualGene objects
 Args    : none


=cut

sub get_all_VirtualGenes {
    my ($self,$supporting) = @_;
        
    my @out;
    foreach my $gene ( $self->get_all_Genes($supporting)) {
	my $vg = Bio::EnsEMBL::VirtualGene->new(-gene => $gene,-contig => $self);
	push(@out,$vg);
    }
    return @out;
}



=head2 get_all_VirtualGenes_startend

 Title   : get_all_VirtualGenes_startend
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_VirtualGenes_startend{
   my ($self,@args) = @_;

   my (%trans,%exon,%exonconverted);
   
   my %gene;
   my @ret;
   
   foreach my $contig ($self->_vmap->get_all_RawContigs) {
       foreach my $gene ( $contig->get_all_Genes() ) {      
	   $gene{$gene->id()} = $gene;
       }
   }
   
   
 GENE:
   foreach my $gene ( values %gene ) {
       
       my $internalExon=0;
       my $genestart = $self->length;
       my $geneend   = 1;
       my $genestr;
       foreach my $trans ( $gene->each_Transcript ) {
	   foreach my $exon ( $trans->get_all_Exons ) {
	       
	       my ($st,$en,$str) = $self->_convert_start_end_strand_vc($exon->contig_id,$exon->start,$exon->end,$exon->strand);
               if( !defined $st ) { 
		   next;
               }
	       if( $st < $genestart ) {
		   $genestart = $st;
	       }
	       if( $en > $geneend ) {
		   $geneend = $en;
	       }
	       $genestr = $str;
	       

	       # hack to get things to behave
	       $exon->seqname($exon->contig_id);
	       $exon{$exon->id} = $exon;
	       
	       ### got to treat sticky exons separately.
	       if( $exon->isa('Bio::EnsEMBL::StickyExon') ) {
		   my @stickies = $exon->each_component_Exon();
		   # sort them by start-end
		   @stickies = sort { $a->start <=> $b->start } @stickies;
		   my $st_start;
		   my $st_end;
		   my $st_strand;
		   my $current_end;
		   my $mapped_sticky = 1;
		   
		   foreach my $sticky ( @stickies ) {
		       if( $self->_convert_seqfeature_to_vc_coords($sticky) == 0 ) {
			   $mapped_sticky = 0;
			   last;
		       } else {
			   if( defined $current_end ) {
			       if( $sticky->start-1 != $current_end ) {
				   $mapped_sticky = 0;
				   last; 
			       }
			}
			   if( !defined $st_start ) {
			       $st_start = $sticky->start;
			   }
			   # at the end of this loop, will be the last one
			   $st_end = $sticky->end;
			   $st_strand = $sticky->strand;
		       }
		   }
		   
		   if( $mapped_sticky == 1 ) {	
		    $exonconverted{$exon->dbID} = 1;
		    $internalExon = 1;
		} else {
		    # do nothing
		}
		   
	       } else {
		   if ($self->_convert_seqfeature_to_vc_coords($exon)) {
		       $internalExon = 1;
		       $exonconverted{$exon->dbID} = 1;
		   } else {$internalExon=0;}               
		   
	       }   
	   }
       }
       
       if ($internalExon) {    
	   my $vg = Bio::EnsEMBL::VirtualGene->new(-gene => $gene,
                                                   -contig => $self, 
						   -start => $genestart, 
                                                   -end => $geneend, 
						   -strand => $genestr);
	   push(@ret,$vg);
 
       }
   }
   
   return @ret;
}


=head2 get_all_SeqFeatures

 Title   : get_all_SeqFeatures
 Usage   : foreach my $sf ( $contig->get_all_SeqFeatures ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SeqFeatures {
   my ($self) = @_;
   my @out;
   push(@out,$self->get_all_SimilarityFeatures());
   push(@out,$self->get_all_RepeatFeatures());

   return @out;

}




=head2 get_all_SimilarityFeatures_above_score

 Title   : get_all_SimilarityFeatures_above_score
 Usage   : foreach my $sf ( $contig->get_all_SimilarityFeatures_above_score(analysis_type, score) ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SimilarityFeatures_above_score{

    my ($self, $analysis_type, $score) = @_;
    
    my $sf = [];
    
    foreach my $c ($self->_vmap->get_all_RawContigs) {
	
	
	push(@$sf,$c->get_all_SimilarityFeatures_above_score($analysis_type,$score));
	
	# Need to clip seq features to fit the boundaries of
	# our v/c so displays don't break
	my @vcsf = ();
	my $count = 0;
	foreach $sf ( @$sf ) {
	    $sf = $self->_convert_seqfeature_to_vc_coords($sf);
	    
	    if( !defined $sf ) {next;}        
	    
	    if (($sf->start < 0 ) || ($sf->end >$self->length)) {$count++;}
	    
	    else{push (@vcsf, $sf);}
	}
	
	return @vcsf;	
    }
    
}








=head2 get_all_SimilarityFeatures

 Title   : get_all_SimilarityFeatures
 Usage   : foreach my $sf ( $contig->get_all_SimilarityFeatures ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SimilarityFeatures {
   my ($self) = @_;
   
   return $self->_get_all_SeqFeatures_type('similarity');

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
   my ($self) = @_;
   
   return $self->_get_all_SeqFeatures_type('repeat');

}




=head2 get_all_ExternalFeatures

 Title   : get_all_ExternalFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut




sub get_all_ExternalFeatures {
   my ($self) = @_;

   return $self->_get_all_SeqFeatures_type('external');
}

=head2 get_all_PredictionFeatures

 Title   : get_all_PredictionFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_PredictionFeatures {
   my ($self) = @_;

   return $self->_get_all_SeqFeatures_type('prediction');
}


sub get_MarkerFeatures {

    my ($self)=@_;
     return $self->_get_all_SeqFeatures_type('marker');
}



=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Genes {
    my ($self) = @_;
    my (%gene,%trans,%exon,%exonconverted);
    
    if (defined $self->{'_all_genes'}) {return @{$self->{'_all_genes'}}}

    my @internal_id;
    foreach my $contig ( $self->_vmap->get_all_RawContigs) {
	push(@internal_id,$contig->internal_id);
    }

    if( scalar(@internal_id) == 0 ) { return (); }

    my $idlist = join(',',@internal_id);

    my $sth = $self->dbobj->prepare("
        SELECT distinct(t.gene_id)
        FROM transcript t
          , exon_transcript et
          , exon e
        WHERE e.contig_id IN ($idlist)
          AND et.exon_id = e.exon_id
          AND et.transcript_id = t.transcript_id
        ");
    my $res = $sth->execute;

    my @geneid;
    while( my ($gene) = $sth->fetchrow_array ) {
	push(@geneid,$gene);
    }

    if( scalar(@geneid) == 0 ) {
	return ();
    }

    my $gadp = $self->dbobj->get_GeneAdaptor();
    my @gene;

    foreach my $geneid ( @geneid ) {
	push(@gene,$gadp->fetch_by_dbID($geneid));
    }

    foreach my $gene ( @gene ) {
	$gene{$gene->dbID()}= $gene;
    }

    my @genes=$self->_gene_query(%gene);

    $self->{'_all_genes'}=\@genes;
    return @{$self->{'_all_genes'}};



}

=head2 get_all_ExternalGenes

 Title   : get_all_ExternalGenes
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_ExternalGenes {
   my ($self) = @_;
   my (%gene);

   my @rawids;


   if( exists $self->{'_external_gene_cache'} ) {
     return @{$self->{'_external_gene_cache'}};
   }

   foreach my $rc ( $self->_vmap->get_all_RawContigs ) {
       push(@rawids,$rc->id);
   }

   &eprof_start("external_gene_retrieve");
   foreach my $extf ( $self->dbobj->_each_ExternalFeatureFactory ) {
     if( $extf->can('get_Ensembl_Genes_contig_list')) {
	     foreach my $gene ( $extf->get_Ensembl_Genes_contig_list(@rawids) ) {
	       #print STDERR "Retrieved gene with ",$gene->id,"\n";
	       $gene{$gene->dbID()} = $gene;
	       #$gene{$gene->stable_id()} = $gene;
	     }
       }
   }
   &eprof_end("external_gene_retrieve");
   
   # foreach my $contig ($self->_vmap->get_all_RawContigs) {
   # foreach my $gene ( $contig->get_all_ExternalGenes() ) {
   #    $gene{$gene->id()} = $gene;
   #}
   #}
   &eprof_start("external_gene_lift");
   my @array=$self->_gene_query(%gene);
   &eprof_end("external_gene_lift");

   $self->{'_external_gene_cache'} = \@array;

   return @array;
}



=head2 get_old_Exons

 Title   : get_old_exons
 Usage   : foreach my $exon ( $contig->get_old_Exons ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_old_Exons{
    my ($self,$logfile,$mapref) = @_;
    
    $self->_old_exon_call(1);
    
    my @exons;
    foreach my $c ($self->_vmap->get_all_RawContigs) {
	push(@exons,$c->get_old_Exons($logfile,$mapref));
	push(@{$self->{'_unmapped_exons'}},$c->unmapped_exons);
    }
    print "fetched ".scalar(@exons)."\n";
    
    my @vcexons = ();
    foreach my $exon ( @exons ) {
	my ($st,$en,$str) = $self->_convert_start_end_strand_vc($exon->contig_id,$exon->start,$exon->end,$exon->strand);
	
	if( !defined $st ) {next;}        
	if (($st < 0 ) || ($en >$self->length)) {next;}
	
	else 
	{
	    $exon->start($st);
	    $exon->end($en);
	    $exon->strand($str);   
	    
	    push @vcexons,$exon;
	}
    }
    print "converted ".scalar(@vcexons)."\n";

    #Make sure it's not redundant
    my @unique=();
    my %seen=();
    foreach my $exon (@vcexons) {
	push (@unique,$exon) unless $seen{$exon->id}++;
    }
    print "unique ".scalar(@unique)."\n";
    
    return @unique;	  
}

=head2 get_old_Genes

 Title   : get_old_Genes
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_old_Genes {
    my ($self,$mapref) = @_;
    my (%gene,%trans,%exon,%exonconverted);
    my @genes;
    foreach my $contig ($self->_vmap->get_all_RawContigs) {
	push (@genes,$contig->get_old_Genes($mapref));
    }
	#foreach my $gene ( $contig->get_old_Genes($mapref) ) {      
	#    $gene{$gene->id()} = $gene;
	#}
    #}
    #return $self->_gene_query(%gene);
    return @genes;
}

=head2 unmapped_exons

 Title   : unmapped_exons
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub unmapped_exons{
   my ($self,@args) = @_;

   if (! $self->old_exon_call) {
       $self->throw("Cannot ask for unmapped exons if get_old_Exons is not called first");
   }
   else {
       return @{$self->{'_unmapped_exons'}};
   }
}


=head2 get_all_Exons

 Title   : get_all_Exons
 Usage   : foreach my $exon ( $contig->get_all_Exons ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Exons{

    my ($self) = @_;
    
    my @exons;


    foreach my $c ($self->_vmap->get_all_RawContigs) {push(@exons,$c->get_all_Exons);}
    
    my @vcexons = ();
    foreach my $exon ( @exons ) {
	my ($st,$en,$str) = $self->_convert_start_end_strand_vc($exon->contig_id,$exon->start,$exon->end,$exon->strand);
	
	if( !defined $st ) {next;}        
	if (($st < 0 ) || ($en >$self->length)) {next;}
	
	else 
	{
	    $exon->start($st);
	    $exon->end($en);
	    $exon->strand($str);   
	    
	    push @vcexons,$exon;
	}
    }
    #Make sure it's not redundant
    my @unique=();
    my %seen=();
    foreach my $exon (@vcexons) {
	push (@unique,$exon) unless $seen{$exon->id}++;
    }
    
    return @unique;	
}



=head2 get_Genes_by_Type

 Title   : get_Genes_by_Type
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Genes_by_Type {
    my ($self,$type,$supporting) = @_;
    my (%gene,%trans,%exon,%exonconverted);
    
    foreach my $contig ($self->_vmap->get_all_RawContigs) {
	foreach my $gene ( $contig->get_Genes_by_Type($type,$supporting) ) {      
	    $gene{$gene->dbID()} = $gene;
	}
    }

    return $self->_gene_query(%gene);

}

# (PL: the naming of this function is, shall we say, a bit unobvious :-)
# Furthermore, it's not entirely clear which coordinate system is assumed,
# and if nothing is lost during the transfer)
sub _gene_query{
    
    my ($self,%gene) = @_;
    my (%trans,%exon,%exonconverted);
        
    my $internalExon = 0;

    foreach my $gene ( values %gene ) {

        $internalExon =0;

	foreach my $exon ( $gene->get_all_Exons() ) {
	    # hack to get things to behave
	    $exon->seqname($exon->contig->id);
	    $exon{$exon->dbID} = $exon;
	    
     #      print STDERR "Exon for gene ",$gene->dbID," is on ",$exon->seqname," ",$exon->start,":",$exon->end,"\n";


	    ### got to treat sticky exons separately.
	    if( $exon->isa('Bio::EnsEMBL::StickyExon') ) {
		my @stickies = $exon->each_component_Exon();

                # walk along stickies in order. If any sticky does
                # not convert, then abort. Otherwise remember start-end
                # strand. If strand disagrees, abort as well.
                # Have not implemented component exons being in sync with assembly

                my $mapped_sticky = 1;
                my $vc_start;
                my $vc_end;
                my $vc_strand = undef;

                foreach my $sticky ( @stickies ) {
		    print STDERR "Mapping sticky $sticky ",$sticky->start," ",$sticky->end,"\n";
                    unless ( $self->_convert_seqfeature_to_vc_coords($sticky) ) {
			print STDERR "Unmappable!\n";
                        # unmappable component exon, abort.
                        $mapped_sticky = 0;
                        last;
                    } else {
                        # handle start end points.
			print STDERR "In sticky, seen ",$sticky->start," ",$sticky->end,"\n";
			if( !defined $vc_strand ) {
                            $vc_start = $sticky->start;
                            $vc_end   = $sticky->end;
                            $vc_strand = $sticky->strand;
                        } else {
                            if( $vc_strand != $sticky->strand ) {
                                $self->warn("sticky exon mappable but strand switching");
                                $mapped_sticky = 0;
                                last; # end of foreach my $sticky
                            }
                            # could test for contigous with assembly
                            if( $vc_start > $sticky->start ) {
                                $vc_start = $sticky->start;
                            } 
                            if( $vc_end  < $sticky->end ) {
                                $vc_end   = $sticky->end;
                            }
                        }

                    }
                }

                if( $mapped_sticky == 0 ) {
                    # sticky was not mapped...
                    # no need to do anything - current exon object
                    # is valid
                } else {
                    # sticky exon mapped. Reset sequnece and
                    # coordinates of the StickyExon
		    $exon->attach_seq($self->primary_seq);
		    $exon->seqname($self->id);
		    $exon->start($vc_start);
		    $exon->end($vc_end);
		    $exon->strand($vc_strand);
		    $exonconverted{$exon->dbID} = 1;
                    $internalExon = 1;
                }

	    } else {                    # not a Sticky
		# soooooo much simpler
	        
		if ($self->_convert_seqfeature_to_vc_coords($exon)) {
		    $internalExon = 1;
		    $exonconverted{$exon->dbID} = 1;
		}               
	    }

	}                               # foreach exon
        
        if ($internalExon == 0) { 
# Utterly weird perl behaviour on acari: using unless breaks this 
#        unless  ($internalExon) {
        #    my $geneid = $gene->id;
            delete $gene{$gene->dbID};
        } 
    }                                   # foreach gene

# PL: there used (rev. 1.17) to be code dealing with converting raw to
# virtual coords for Translation objs; that is not needed anymore. 

    
    return values %gene;
}


=head2 _get_all_SeqFeatures_type

 Title   : _get_all_SeqFeatures_type
 Usage   : Internal function which encapsulates getting
           features of a particular type and returning
           them in the VC coordinates, optionally cacheing
           them.
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _get_all_SeqFeatures_type {
   my ($self,$type) = @_;

   if( $self->_cache_seqfeatures() && $self->_has_cached_type($type) ) {
       return $self->_get_cache($type);
   }

   # ok - build the sequence feature list...

   my $sf;
   
   if( $self->_cache_seqfeatures() ) {
       $sf = $self->_make_cache($type);
   } else {
       $sf = []; 
   }
   
   foreach my $c ($self->_vmap->get_all_RawContigs) {
       if( $type eq 'repeat' ) {
	   push(@$sf,$c->get_all_RepeatFeatures());
       } elsif ( $type eq 'similarity' ) {
	   push(@$sf,$c->get_all_SimilarityFeatures());
       } elsif ( $type eq 'prediction' ) {
	   push(@$sf,$c->get_all_PredictionFeatures());
       } elsif ( $type eq 'external' ) {
	   push(@$sf,$c->get_all_ExternalFeatures());
       } elsif ( $type eq 'marker' ) {
	   push(@$sf,$c->get_MarkerFeatures());
       } else {
	   $self->throw("Type $type not recognised");
       }
   }

   #print STDERR "before clipping ",scalar(@$sf)," for $type\n";

   my @vcsf = ();


   # need to clip seq features to fit the boundaries of
   # our v/c so displays don't break

   my $count = 0;
   foreach $sf ( @$sf ) {
       $sf = $self->_convert_seqfeature_to_vc_coords($sf);
       
       if( !defined $sf ) {      
	   next;
       }
       
       
       if($sf->start < 0 ){
	   $count++;
       }
       elsif ($sf->end > $self->length){
	   $count++;
       }
       else{
	   
	   if( $type eq 'prediction' ) {
	       foreach my $sub ( $sf->sub_SeqFeature ) {
	       }
	   }
	   push (@vcsf, $sf);
       }
   }
   
   return @vcsf;
}


=head2 _convert_seqfeature_to_vc_coords

 Title   : _convert_seqfeature_to_vc_coords
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _convert_seqfeature_to_vc_coords {
    my ($self,$sf) = @_;
    
    my $cid = $sf->seqname();
    my $mc;
    if( !defined $cid ) {
    	$self->throw("sequence feature [$sf] has no seqname!");
    }
    # potentially we could be asked to convert something
    # that wasn't on this VC at all, eg, an exon from a distant contig
    eval {
	    $mc = $self->_vmap->get_MapContig_by_id($cid);
    };
    if ($@ || !ref $mc) { 
	    return undef;
    }

    # if this is something with subfeatures, then this is much more complex
    my @sub = $sf->sub_SeqFeature();
    
    if( $#sub >=  0 ) {
	    # chain to constructor of the object. Not pretty this.
	    $sf->flush_sub_SeqFeature();
	    my $seen = 0;
	    my $strand;
	    foreach my $sub ( @sub ) {
	        #print STDOUT "Converting sub ",$sub->id,":",$sub->seqname,":",$sub->start,":",$sub->end,":",$sub->strand,"\n";
	        $sub = $self->_convert_seqfeature_to_vc_coords($sub);
	        if( !defined $sub ) {        
		        next;
	        }
	        if( $seen == 0 ){
		        $sf->start($sub->start);
		        $sf->end($sub->end);
	        }
	        $seen =1;
	        $strand = $sub->strand;
	        $sf->add_sub_SeqFeature($sub,'EXPAND');
	    }
	    if( $seen == 1 ) {       
	        # we assumme that the mapping was unambiguous wrt to the strand
	        $sf->strand($strand);
	        return $sf;
	    } else {
	        return undef;
	    }
    }

    # might be clipped left/right
    #print ("Leftmost " . $mc->leftmost . " " . $mc->orientation . " " . $mc->start_in . " " . $mc->end_in  . " " . $sf->start . " " . $sf->end . "\n");
    # Could be clipped on ANY contig  
    if ($sf->start < $mc->rawcontig_start) {
        # print STDERR "Binning $cid\n";
	    return undef;              
    }
    if ($sf->end >  $mc->rawcontig_end) {  
        # print STDERR "Binning $cid\n";
	    return undef;              
    }
    my ($rstart,$rend,$rstrand) = $self->_convert_start_end_strand_vc($cid,$sf->start,$sf->end,$sf->strand);

    if( $sf->can('attach_seq') ) {
	    if (!$self->noseq) {
	        $sf->attach_seq($self->primary_seq);
	    }
    }
    $sf->start ($rstart);
    $sf->end   ($rend);
    $sf->strand($rstrand);
    $sf->seqname($self->id);
    return $sf;

}

=head2 _convert_start_end_strand_vc

 Title   : _convert_start_end_strand_vc
 Usage   : Essentially an internal for _convert_seqfeature,
           but sometimes we have coordinates not on seq features

 Function: convert RawContig coordinates of a given RawContig to
           the coords of the VC.
 Example : ($start,$end,$strand) = $self->_convert_start_end_strand_vc($contigid,$start,$end,$strand)
 Returns : A list of start,end,strand in  VC coords.
 Args    : RawContig display_id, start, end, strand in RC coords

=cut

sub _convert_start_end_strand_vc {
    my ($self,$contig,$start,$end,$strand) = @_;
    my ($rstart,$rend,$rstrand);

    my $mc;
    eval {
	$mc=$self->_vmap->get_MapContig_by_id($contig);
    };
 
    if ($@ || !ref $mc) { 
	return undef;
    }

    if( $mc->orientation == 1 ) {
       
        # ok - forward with respect to vc. Only need to add offset
	my $offset = $mc->start - $mc->rawcontig_start;
	$rstart = $start + $offset;
	$rend   = $end + $offset;
	$rstrand = $strand;
    } else {

	# flip strand
	$rstrand = $strand * -1;
	
	# yup. A number of different off-by-one errors possible here

	$rstart = $mc->end   - ($end   - $mc->rawcontig_start);
	$rend   = $mc->end   - ($start - $mc->rawcontig_start);

    }
    return ($rstart,$rend,$rstrand);
}


=head1 Helper functions for Virtual::Contig

=head2 raw_contig_position 
 Title   : raw_contig_position
 Usage   : my ($map_contig,$rc_position,$rc_strand) = $vmap->raw_contig_position($vc_pos,$vc_strand)
 Function: Maps a VirtualContig position to a RawContig + RawContig position
 Returns : The underlying RawContig and a position on it (in RC coords),
           and optionally the RC strandedness
 Args   : position on VirtualContig (in VC coords), and optionally
          VirtualContig strand.
=cut

sub raw_contig_position {
    my ($self, $vcpos, $vcstrand)=@_;

    return $self->_vmap->raw_contig_position($vcpos,$vcstrand);
}

=head2 skip_SeqFeature

 Title   : skip_SeqFeature
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub skip_SeqFeature {
   my ($self,$tag,$value) = @_;

   if( defined $value ) {
       $self->{'feature_skip'}->{$tag} = $value;
   }

   return $self->{'feature_skip'}->{$tag};
}

=head2 rawcontig_ids

 Title   : rawcontig_ids
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub list_rawcontig_ids {
   my ($self) = @_;
   
   my @contig_ids;
   for my $contig ( $self->_vmap->get_all_RawContigs() ) {
     push ( @contig_ids, $contig->id() );
   }
   return @contig_ids;
}

=head2 found_left_end

 Title   : found_left_end
 Usage   : 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub found_left_end {
    my ($self) = @_;
    return $self->_vmap->found_left_end;
}


=head2 found_right_end

 Title   : found_right_end
 Usage   : 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub found_right_end {
    my ($self) = @_;
    
    return $self->_vmap->found_right_end;
}

=head2 is_truncated

 Title   : is_truncated
 Usage   : 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub is_truncated {
    my ($self) = @_;
    
    my $flag = 0;
    
    if ($self->found_right_end == 1) {
	$flag = 1;
    } 
    if ($self->found_left_end == 1) {
	$flag = 1;
    }
    
    return $flag;
}

=head2 _dump_map

 Title   : _dump_map
 Usage   : Produces a dumped map for debugging purposes
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _dump_map {
   my ($self,$fh) = @_;

   ! defined $fh && do { $fh = \*STDERR};
   print $fh "Contig Map Dump: \n";
   
   foreach my $mc ($self->_vmap->each_MapContig) {
       print $fh "Contig ".$mc->contig->id." starts:",$mc->start," ends:".$mc->end." start in contig ",$mc->rawcontig_start," end in contig ".$mc->rawcontig_end." orientation ",$mc->orientation,"\n";
   }
}


=head2 get_all_RawContigs

 Title   : get_all_RawContigs
 Usage   : $obj->get_all_RawContigs()
 Function: 
 Example : 
 Returns : array of raw contigs sorted by starting position in vc
 Args    : 


=cut

sub get_all_RawContigs {
    my ($self) = @_;

    my @contigs = ();

#    if( !ref $self || ! $self->isa('Bio::EnsEMBL::DB::VirtualContigI') ) {
#        $self->throw("Must supply a VirtualContig to get_all_RawContigs: Bailing out...");
#    }
    
    return $self->_vmap->get_all_RawContigs;
}

=head2 get_rawcontig_by_position

 Title   : get_rawcontig_by_position
 Usage   : $obj->get_rawcontig_by_position($position)
 Function: 
 Example : 
 Returns : returns a raw contig object or undef on error
 Args    : 


=cut

sub get_rawcontig_by_position {

    my ($self, $pos) = @_;
    
    if( !ref $self || ! $self->isa('Bio::EnsEMBL::DB::VirtualContigI') ) {
        $self->throw("Must supply a VirtualContig to get_all_RawContigs: Bailing out...");
    }
    
    if ($pos < 1 || $pos > $self->length){
        $self->throw("get_rawcontig_by_position error: Position must be > 0 and < vc->length");
    }

    foreach my $mc ($self->_vmap->get_all_MapContigs) {
        if ($pos > $mc->end) {
            next;
        } else {
            return $mc->contig;
        }
    }

    return (undef);
}

=head1 Convertible Virtual::Contig stuff

=head2 convert_Gene_to_raw_contig

 Title   : convert_Gene_to_raw_contig
 Usage   : $newgene = $cvc->convert_Gene_to_raw_contig($gene)
 Function: Converts a gene built on this VirtualContig back to being
           a gene built on raw contig positions, read to be written back 
           

           Internally this builds a copy of the genes,transcripts and translations
           in RC coordinate space, making heavy use of VirtualMap raw_contig_position
           function. Exons could be split into sub exons across boundaries.

 Example :
 Returns : 
 Args    :


=cut

sub convert_Gene_to_raw_contig {
   my ($self,$gene) = @_;

   if( !ref $gene || ! $gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("Got to write a gene, not a [$gene]");
   }

   $self->_sanity_check($gene);

   # we need to map the exons back into RC coordinates.
   
   # in some cases the mapping will cause us to get a large number
   # of sticky exons.

   # the basically means we have to clone all the objects. :(

   my $clonedgene = Bio::EnsEMBL::Gene->new();
   my %translation;
   $clonedgene->type($gene->type);
   $clonedgene->analysis($gene->analysis);

   foreach my $dbl ( $gene->each_DBLink() ) {
       $clonedgene->add_DBLink($dbl);
   }

   #
   #
   #


   # 
   # convert exons first, as unique exons. In particular this is to 
   # handle alternative splicing. Exons crossing boundaries become
   # sticky exons
   #
   my %convertedexon;

   foreach my $exon ( $gene->get_all_Exons ) {
       $convertedexon{$exon} = $self->_reverse_map_Exon($exon);
   }



   foreach my $trans ( $gene->each_Transcript ) {
       my $clonedtrans = Bio::EnsEMBL::Transcript->new();
       #print STDERR "Reverse mapping ",$trans->id,"\n";
       foreach my $dbl ( $trans->each_DBLink() ) {
	   $clonedtrans->add_DBLink($dbl);
       }

       $clonedgene->add_Transcript($clonedtrans);

       foreach my $exon ( $trans->get_all_Exons ) {
	   $clonedtrans->add_Exon($convertedexon{$exon});

	   # translations.
           # (PL: looks like a 'deep copy' is being made;
           # can't we reuse the existing object, and forget about deep
           # copying ?)
	   if( exists $translation{$trans->translation} ) {
#	     print STDERR "Translation already exists " . $trans->id . " " . $trans->translation->id . "\n";
               # (PL: looks like this is cached; should it? )
	       $clonedtrans->translation($translation{$trans->translation});
	   } else {
	#     print STDERR "Making new translation " . $trans->id . " " . $trans->translation->id . "\n";
	       my $trl = $trans->translation(); 

	       my $clonedtrl = Bio::EnsEMBL::Translation->new();
#	       $clonedtrl->id($trl->id);
	       $clonedtrl->start_exon($convertedexon{$trl->start_exon});
               $clonedtrl->start( $trl->start );
	       $clonedtrl->end_exon($convertedexon{$trl->end_exon});
               $clonedtrl->end( $trl->end );
	       $clonedtrl->version($trl->version);

               ## code converting Translation obj. coordinates to VC
               ## coords is now gone (since rev. 1.17)
	       
	       $translation{$trans->translation} = $clonedtrl;
	       $clonedtrans->translation($clonedtrl);
	   }
       }
   }

#    $self->_sanity_check($clonedgene);

   return $clonedgene;

}                                       # convert_Gene_to_raw_contig

=head2 _reverse_map_Exon

 Title   : _reverse_map_Exon
 Usage   : $exon = $self->_reverse_map_Exon($exon)
 Function: Makes exons in RawContig coordinates from exon in VC coordinates.
           When an exon crosses a contig boundary, it makes a sticky exon
           
 Example :
 Returns : 
 Args    :


=cut

sub _reverse_map_Exon {
   my ($self,$exon) = @_;

   if( !ref $exon || !$exon->isa('Bio::EnsEMBL::Exon') ) {
       $self->throw("Must supply reverse map an exon not an [$exon]");
   }

   # 
   # this maps a virtual contig exon to a raw contig exon
   # for exons that do not cross contig boundaries in the golden
   # path this is easy - we just map coordinates and transfer ids.
   
   # for exons which do cross boundaries, this is where the magic
   # happens. A StickyExon object, which holds a number of component
   # exons with ascending sticky_rank numbers holds the exon across
   # the join coordinates.

   # for stickyexons on the reverse strand, the sticky_rank ordering
   # has to be right to left along the golden path, not left to right
   # as the strand informaiton of the component exons do not contribute
   # to the final strandness of the assembled sticky exon. 

   # yes - this confuses us as well regularly.



#   print STDERR "Reverse mapping $exon ",$exon->start,":",$exon->end,"\n";

   my ($scontig,$start,$sstrand) = 
     $self->_vmap->raw_contig_position($exon->start,$exon->strand);

   my ($econtig,$end,$estrand)   = 
     $self->_vmap->raw_contig_position($exon->end  ,$exon->strand);

   # print STDERR "Got $scontig ",$start," to $econtig ",$end,"\n";

   if ($scontig eq 'gapcontig' || $econtig eq 'gapcontig') {
       $self->throw("gap contigs not allowed: exon " . $exon->id );
   }

   my $exon_to_return;
  
   if( $scontig->id eq $econtig->id ) {
       if( $sstrand != $estrand ) {
	   $self->throw("Bad internal error. Exon mapped to same contig but different strands!");
       }

#       print STDERR "Straight forward mapping\n";

       my $rmexon = Bio::EnsEMBL::Exon->new(); # the re-mapped exon
#       $rmexon->id($exon->id);
       $rmexon->created($exon->created);
       $rmexon->modified($exon->modified);
       $rmexon->version($exon->version);
       $rmexon->phase($exon->phase);
       $rmexon->sticky_rank(1);


       # we could test on strand changes. This just assummes everything works
       # as it says on the tin ;)
       if( $start < $end ) {
	   $rmexon->start($start);
	   $rmexon->end($end);
       } else {
	   $rmexon->start($end);
	   $rmexon->end($start);
       }
       $rmexon->strand($sstrand);
       $rmexon->contig($scontig);
       $rmexon->seqname($scontig->id);
       $rmexon->attach_seq($scontig->primary_seq);
       $exon_to_return = $rmexon;
   } else {
       # we are in the world of sticky-ness....
#       print STDERR "Into sticky exon\n";

       my @mapcontigs = $self->_vmap->each_MapContig();

       my $found=0;
       my $mc;
       while ( $mc = shift @mapcontigs ) { 
           if( $mc->contig->id eq $scontig->id ) {
# print STDERR "Unshifting ",$mc->contig->id,"\n";
               unshift(@mapcontigs,$mc);
               $found = 1;
               last;
           }
       }

       if ( $found == 0 ) {
	   $self->throw("Internal error - unable to find map contig with this id");
       }

       my $vcstart = $exon->start;
       #print STDERR "Looking from exon-wise",$exon->start,":",$exon->end,"\n";

       # ok. Move from start towards end, after we hit end.
       my @exported_exons;
       my $sticky = 1;

       foreach my $c ( @mapcontigs ) {	   
	   my $vcend;
	   #print STDERR "***Looking at",$c->contig->id," - $vcstart...\n";

	   if( $c->contig->id eq $econtig->id ) {
	       # go to end position
	       #print STDERR "Going for end...",$econtig->id,"\n";
	       $vcend = $exon->end();
	   } else {
	       #print STDERR "Going for end of contig\n";
	       $vcend = $c->end();
	   }

	   #print STDERR "....Going to call with $vcstart:$vcend\n";
	   #$self->_dump_map(\*STDERR);

	   my ($iscontig,$istart,$isstrand) = $self->_vmap->raw_contig_position($vcstart,$exon->strand);
	   my ($iecontig,$iend,$iestrand)   = $self->_vmap->raw_contig_position($vcend  ,$exon->strand);
  
	   if( $iscontig->id ne $iecontig->id || $isstrand != $iestrand) {
	       $self->throw("Bad internal error. Sticky Exon mapped to different contig/strand for a correct contig placement ".$iscontig->id.":".$iecontig->id);
	   }

	   my $rmexon = Bio::EnsEMBL::Exon->new();
	   $rmexon->dbID($exon->dbID);
	   $rmexon->modified($exon->modified);
	   $rmexon->version($exon->version);
	   $rmexon->phase($exon->phase);
	   $rmexon->sticky_rank($sticky++);
	   $rmexon->attach_seq($c->contig->primary_seq);

	   # could have flipped around, in which case, flip again.
	   # Again we assumme everything works as advertised
	   if( $istart > $iend ) {
	       $rmexon->start($iend);
	       $rmexon->end($istart);
	   } else {
	       $rmexon->start($istart);
	       $rmexon->end($iend);
	   }

	   $rmexon->strand($isstrand);
	   $rmexon->contig($c->contig);
	   $rmexon->seqname($c->contig->id);
	   push(@exported_exons,$rmexon);

	   if( $c->contig->id eq $econtig->id ) {
	       last;
	   }
	   $vcstart = $vcend+1;
       }
       my $sticky_exon = Bio::EnsEMBL::StickyExon->new();
       $sticky_exon->dbID($exon->dbID);
       # for reverse strand exons, we need to reverse the
       # order of the components and renumber
       if( $exon->strand == -1) {
           @exported_exons = reverse (@exported_exons);
           $sticky = 1;
           foreach my $e ( @exported_exons ) {
               $e->sticky_rank($sticky++);
           }
       }
       
       # add exons to final StickyExon object
       foreach my $e ( @exported_exons) {
           $sticky_exon->add_component_Exon($e);
       }
       
       # give it back
       $exon_to_return = $sticky_exon;
   }

   #
   # here we handle supporting features
   #

   foreach my $se ( $exon->each_Supporting_Feature ) {
     # we only map featurepairs
     if( ! $se->isa('Bio::EnsEMBL::FeaturePairI') ) {
       $self->warn("In reverse map exon, cannot map supporting feature $se");
       next;
     }

     my @res = $self->_vmap->raw_contig_interval($se->start,$se->end,$se->strand);

     my $hstart = $se->hstart;
     foreach my $res ( @res ) {
       if( $res->{'gap_start'} ) {
	 $hstart = $hstart + $res->{'gap_end'} - $res->{'gap_start'} +1;
	 next; # no features in gaps
       }

       $se->validate();
       print STDERR "SE validation done!";

       my $new_feature = Bio::EnsEMBL::FeatureFactory->new_feature_pair();
       $new_feature->start($res->{'raw_start'});
       $new_feature->end($res->{'raw_end'});
       print STDERR "Setting strand to ",$res->{'raw_strand'},"\n";
       $new_feature->strand($res->{'raw_strand'});
       $new_feature->seqname($res->{'raw_contig_id'});

       $new_feature->hstart($hstart);
       $new_feature->hend($hstart + $res->{'raw_end'} - $res->{'raw_start'});
       $hstart =$hstart + $res->{'raw_end'} - $res->{'raw_start'}+1;
       $new_feature->hstrand($se->hstrand);
       $new_feature->score($se->score);
       $new_feature->hscore($se->score);
       $new_feature->hseqname($se->hseqname);
       
       $new_feature->analysis($se->analysis);
       $new_feature->source_tag($se->source_tag);
       $new_feature->primary_tag($se->primary_tag);

       $new_feature->hsource_tag($se->hsource_tag);
       $new_feature->hprimary_tag($se->hprimary_tag);

       print STDERR "Adding feature to exon ",$exon_to_return->dbID,"\n";
       $new_feature->validate();

       $exon_to_return->add_Supporting_Feature($new_feature);
     }
   }

   return $exon_to_return;

}

# internal function used by _sanity_check; returns undef if all OK, error
# string otherwise
sub _check_exon_start_end {
    my ($transl, $exon, $which_one) = @_;
    my $message = undef;
    my $to_test;

    # depending on the 'which_one' arg, check beginning or end of it. 
    if ($which_one eq 'start' ) { 
        $to_test =  $transl->start;
    } elsif ($which_one eq 'end' ) {
        $to_test =  $transl->end;
    } else { 
        return "Internal error: call with 'start' or 'end'";
    }

    if ( $to_test < 1 ) { 
        $message .= "Translation's $which_one < 1: " 
          . (($which_one eq 'start')?$transl->start :$transl->end)
          . " (translation:". $transl->stable_id()  . ")";
    }
    
    if ( $to_test > $exon->length) { 
        $message .= "Translation's $which_one ("
          . (($which_one eq 'start')?$transl->start :$transl->end)
          .") > exon length (".$exon->length.")"
            . " (translation:". $transl->stable_id() . ",exon:".$exon->stable_id().")";
    }

    return $message;
}                                       # _check_exon_start_end


=head2 _sanity_check

 Title   : _sanity_check
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _sanity_check{
   my ($self,$gene) = @_;
   my $error =0;
   my $message;

   if( !defined $gene->analysis ) {
     $error = 1;
     $message .= "Gene has no analysis object";
   }

   foreach my $transc ( $gene->each_Transcript ) {

       if( !defined $transc->translation || !ref $transc->translation) {
	   $error = 1;
	   $message .= "Transcript has no translation;";
       } else {
	   if( !defined $transc->translation->start ) {
	       $error = 1;
	       $message .= "Translation has no start";
	   } 
	   if( !defined $transc->translation->start_exon ) {
	       $error = 1;
	       $message .= "Translation has no start exon id";
	   } 
	   if( !defined $transc->translation->end ) {
	       $error = 1;
	       $message .= "Translation has no end";
	   } 
	   if( !defined $transc->translation->end_exon ) {
	       $error = 1;
	       $message .= "Translation has no end exon id";
	   } 
       }
       foreach my $exon ( $transc->get_all_Exons ) {

	   if( !defined $exon->contig_id  ) {
	       $error = 1;
	       $message .= "Exon has no contig id";
	   } else {
	 
	       if( $exon->contig_id ne $self->id ) {
	#	   $error = 1;
	#	   $message .= "Exon [".$exon->id."] does not seem to be on this VirtualContig";
	       }
	   }
	   if( !defined $exon->start || !defined $exon->end) {
	       $error = 1;
	       $message .= "Exon has error in start/end";
	   } 
	   if (!defined $exon->strand) {
	       $error = 1;
	       $message .= "Exon does not have a strand";
	   } 
	   if (!defined $exon->phase) {
	       $error = 1;
	       $message .= "Exon does not have a phase";
	   }           
       }

       # now see if the changes (using exon coords rather than VC)
       # have worked:
       # start exon:
       my $t = $transc->translation;
       if( $t->start > $t->start_exon->length ) {
	   $message .= "Translation start ".$t->start." is greater than start exon length ".$t->start_exon->length;
	   $error = 1;
       }
       if( $t->end > $t->end_exon->length ) {
	   $message .= "Translation end ".$t->end." is greater than end exon length ".$t->end_exon->length;
	   $error = 1;
       }

   }                                    # each_Transcript

   if( $error ) {
       $self->throw("Cannot write gene due to: $message");
   }
}                                       # _sanity_check


=head1 Get/Set for Virtual::Contig orientated datastructures.

=head2 _cache_seqfeatures

 Title   : _cache_seqfeatures
 Usage   : $obj->_cache_seqfeatures($newval)
 Function: 
 Returns : value of _cache_seqfeatures
 Args    : newvalue (optional)


=cut

sub _cache_seqfeatures {
    my $obj = shift;
    if( @_ ) {
	my $value = shift;
	$obj->{'_cache_seqfeatures'} = $value;
    }
    return $obj->{'_cache_seqfeatures'};
}

=head2 _has_cached_type

 Title   : _has_cached_type
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _has_cached_type {
    my ($self,$type) = @_;
    
    if ( exists $self->{'_sf_cache'}->{$type} ) {
	return 1;
    } else {
	return 0;
    }
}

=head2 _make_cache

 Title   : _make_cache
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _make_cache {
    my ($self,$type) = @_;
    
    if( $self->_has_cached_type($type) == 1) {
	$self->throw("Already got a cache for $type! Error in logic here");
    }
    
    $self->{'_sf_cache'}->{$type} = [];
    
    return $self->{'_sf_cache'}->{$type};
}

=head2 _get_cache

 Title   : _get_cache
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _get_cache {
    my ($self,$type) = @_;
    
    return $self->{'_sf_cache'}->{$type};
    
}

=head2 _clear_vc_cache

 Title   : _clear_vc_cache
 Usage   : $virtual_contig->_clear_vc_cache
 Function: clears a v/c internal cache (use when extending a virtual contig)
 Example :
 Returns : 
 Args    :


=cut

sub _clear_vc_cache {
    my $self = shift;
    $self->{'_seq_cache'} = undef;
    foreach my $c (keys %{$self->{'_sf_cache'}}){
        $self->{'_sf_cache'}->{$c} = undef;
    }   
}

=head2 _seq_cache

 Title   : _seq_cache
 Usage   : $obj->_seq_cache($newval)
 Function: 
 Returns : value of _seq_cache
 Args    : newvalue (optional)


=cut

sub _seq_cache {
    my $obj = shift;
    if( @_ ) {
	my $value = shift;
	$obj->{'_seq_cache'} = $value;
    }
    return $obj->{'_seq_cache'};
}

=head2 _unique_number

 Title   : _unique_number
 Usage   : $obj->_unique_number($newval)
 Function: 
 Example : 
 Returns : value of _unique_number
 Args    : newvalue (optional)


=cut

sub _unique_number{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_unique_number'} = $value;
    }
    return $obj->{'_unique_number'};

}

=head2 _vmap

 Title   : _vmap
 Usage   : $obj->_vmap($newval)
 Function: 
 Example : 
 Returns : value of _vmap
 Args    : newvalue (optional)


=cut

sub _vmap{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_vmap'} = $value;
    }
    return $obj->{'_vmap'};

}


=head2 noseq

 Title   : noseq
 Usage   : $vc->noseq
 Function: If set to 1, vc does not attach a primary seq to its features (exons) 
 Example : $vc->noseq(1)
 Returns : nothing
 Args    : 1/0


=cut

sub noseq{
    my ($obj,$value) = @_;
    
    if( defined $value) {
	$obj->{'_noseq'} = $value;
    }
    return $obj->{'_noseq'};    
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

=head2 sv

 Title   : sv
 Usage   : $obj->sv($newval)
 Function: 
 Returns : value of sv
 Args    : newvalue (optional)


=cut

sub sv{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'sv'} = $value;
    }
    return $obj->{'sv'};

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
   my ($self,$value) = @_;

   if( defined $value ) {
       $self->throw("Can't set species any more - comes from database");
   }

   if( !defined $self->{'_species_cache'} ) {
       $self->{'_species_cache'} = $self->dbobj->get_MetaContainer->get_Species();
   }
   return $self->{'_species_cache'};
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
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'annotation'} = $value;
    } else {
	if( !defined $self->{'annotation'} ) {
	    my $annotation = Bio::Annotation->new();
	    $self->{'annotation'} = $annotation;
	}
    }

    return $self->{'annotation'};
}

=head2 _old_exon_call

 Title   : _old_exon_call
 Usage   : $obj->_old_exon_call($newval)
 Function: Getset for _old_exon_call value
 Returns : value of _old_exon_call
 Args    : newvalue (optional)


=cut

sub _old_exon_call{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_old_exon_call'} = $value;
    }
    return $obj->{'_old_exon_call'};

}


=head2 dbobj

 Title   : dbobj
 Usage   : $obj->dbobj($newval)
 Function: 
 Example : 
 Returns : value of dbobj
 Args    : newvalue (optional)


=cut

sub dbobj{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_dbobj'} = $value;
    }
   return $obj->{'_dbobj'};

}

1;





