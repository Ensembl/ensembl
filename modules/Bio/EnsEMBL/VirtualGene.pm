
#
# BioPerl module for Bio::EnsEMBL::VirtualGene
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::VirtualGene - A Gene viewed from one Contig's perspective

=head1 SYNOPSIS

   $vg = Bio::EnsEMBL::VirtualGene->new( -gene => $gene, -contig => $contig );
 
   print "The gene ",$vg->gene->id," has it is first exon at ", $vg->start," and last at ",$vg->end,
         "from the perspective of ",$contig->id,"\n";

   # valid methods - this is a Bio::SeqFeatureI compliant object
   $vg->start();
   $vg->end();
   $vg->strand(); 
   $vg->primary_tag(); # returns 'genefragment';
   $vg->source_tag();  # returns 'ensembl'

   # you can get out GFF if you really want to
   print $vg->gff_string,"\n";

   # you can add it to Bio::Seq objects  
   $seq->add_SeqFeature($vg);

   # you can get the original gene object
   $vg->gene

   # test whether there are exons elsewhere or not
   if( $vg->is_complete ) {
      print "The entire gene",$vg->gene->id," is on ",$contig->id,"\n";
   }



=head1 DESCRIPTION

VirtualGene provides a view of a Gene from the perspective of a
contig. In this contig's perspective, the gene has a start and end, being
the first and last exon on the contig respectively. The strand is taken to
be arbitarily the first exon it encounters on the call to each_unique_Exon
on this contig. If the gene is jumping strand (a possibility due to software
issues, not biologically sane of course) then this is not indicated.

VirtualGene, by having a start,end,strand is-a seqfeature

=head1 CONTACT

Ewan Birney

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::VirtualGene;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::SeqFeatureI;
use Bio::Root::Object

@ISA = qw(Bio::Root::Object Bio::SeqFeatureI);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);

  my ($gene,$contig) = $self->_rearrange(['GENE','CONTIG'],@args);
  if( !defined $gene ) {
      $self->throw("No gene in virtualgene object");
  }
  if( !defined $contig || ! ref $contig || ! $contig->isa('Bio::EnsEMBL::DB::ContigI') ) {
      $self->throw("you have to have a virtual gene on a particular contig");
  }

  $self->gene($gene);
  $self->dbobj($contig->dbobj);
  $self->_calculate_coordinates($gene,$contig);

  return $make; # success - we hope!
}

=head2 start

 Title   : start
 Usage   : $obj->start($newval)
 Function: 
 Returns : value of start
 Args    : newvalue (optional)


=cut

sub start{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'start'} = $value;
    }
    return $obj->{'start'};

}

=head2 end

 Title   : end
 Usage   : $obj->end($newval)
 Function: 
 Returns : value of end
 Args    : newvalue (optional)


=cut

sub end{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'end'} = $value;
    }
    return $obj->{'end'};

}

=head2 strand

 Title   : strand
 Usage   : $obj->strand($newval)
 Function: 
 Returns : value of strand
 Args    : newvalue (optional)


=cut

sub strand{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'strand'} = $value;
    }
    return $obj->{'strand'};

}

=head2 seqname

 Title   : seqname
 Usage   : $obj->seqname($newval)
 Function: 
 Returns : value of seqname
 Args    : newvalue (optional)


=cut

sub seqname{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'seqname'} = $value;
    }
    return $obj->{'seqname'};

}

=head2 score

 Title   : score
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub score{
   my ($self,@args) = @_;
   return undef;
}

=head2 frame

 Title   : frame
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub frame{
   my ($self,@args) = @_;

   return undef;
}

=head2 dbobj

 Title   : dbobj
 Usage   : $obj->dbobj($newval)
 Function: 
 Returns : value of dbobj
 Args    : newvalue (optional)


=cut

sub dbobj{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'dbobj'} = $value;
    }
    return $obj->{'dbobj'};

}

=head2 _calculate_coordinates

 Title   : _calculate_coordinates
 Usage   : internal function to fill in start,end,strand,seqname
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _calculate_coordinates{
   my ($self,$gene,$contig) = @_;

   if( !defined $contig || ! ref $contig || ! $contig->isa('Bio::EnsEMBL::DB::ContigI') ) {
       $self->throw("_calculate_coordinates(gene,contig)");
   }

   if( ! ref $gene || ! $gene->isa('Bio::EnsEMBL::Gene') ) {
       $self->throw("_calculate_coordinates(gene,contig)");
   }

   my @exons = $gene->each_unique_Exon();
   my $cid = $contig->id;
   my $outside_exon = 0;
   my $inside_exon = 0;
   my ($start,$end,$strand);
   foreach my $exon ( @exons ) {
       if( $cid = $exon->contig_id ) {
	   if( $inside_exon == 0 ) {
	       $start = $exon->start();
	       $end   = $exon->end();
	       $strand = $exon->strand();
	       $inside_exon = 1;
	   } else {
	       if( $start > $exon->start ) {
		   $start = $exon->start;
	       } 
	       if( $end <  $exon->end ) {
		   $end = $exon->end;
	       } 
	   }

	   $self->add_contained_Exon($exon);
       } else {
	   $outside_exon = 1;
       }
   }

   if( $inside_exon == 0 ) {
       $self->throw("trying to make a virtualgene on a contig which does not contain the gene. Not possible");
   }

   if( $outside_exon == 0 ) {
       $self->is_complete(1);
   } else {
       $self->is_complete(0);
   }

   $self->start($start);
   $self->end($end);
   $self->seqname($cid);
   $self->strand($strand);

}

=head2 primary_tag

 Title   : primary_tag
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub primary_tag{
   my ($self,@args) = @_;

   return "genefragment";
}

=head2 source_tag

 Title   : source_tag
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub source_tag{
   my ($self,@args) = @_;

   return "ensembl";
}

=head2 has_tag

 Title   : has_tag
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub has_tag{
   my ($self,@args) = @_;
   
   return 0;
}

=head2 all_tags

 Title   : all_tags
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub all_tags{
   my ($self,@args) = @_;

   return;
}


=head2 each_tag_value

 Title   : each_tag_value
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_tag_value{
   my ($self,@args) = @_;

   $self->throw("Calling each tag value on a VirtualGene. Not possible");

}

=head2 gene

 Title   : gene
 Usage   : $obj->gene($newval)
 Function: 
 Returns : value of gene
 Args    : newvalue (optional)


=cut

sub gene{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      if( ! ref $value || ! $value->isa("Bio::EnsEMBL::Gene") ) {
	  $obj->throw("Gene object must inheriet from Gene...");
      }

      $obj->{'gene'} = $value;
    }
    return $obj->{'gene'};

}

=head2 is_complete

 Title   : is_complete
 Usage   : $obj->is_complete($newval)
 Function: 
 Returns : value of is_complete
 Args    : newvalue (optional)


=cut

sub is_complete{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'is_complete'} = $value;
    }
    return $obj->{'is_complete'};

}

=head2 add_contained_Exon

 Title   : add_contained_Exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_contained_Exon{
   my ($self,$exon) = @_;

   if( !ref $exon || !$exon->isa("Bio::EnsEMBL::Exon") ) {
       $self->throw("add_contained_Exon $exon");
   }

   push(@{$self->{'_contained_exon'}},$exon);
}

=head2 all_contained_Exons

 Title   : all_contained_Exons
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub all_contained_Exons {
   my ($self) = @_;

   return @{$self->{'_contained_exon'}};
}

=head2 sub_SeqFeature

 Title   : sub_SeqFeature
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub sub_SeqFeature{
   my ($self) = @_;

   return $self->all_contained_Exons();
}


#
# This is the magic of the FT EMBL file production.
#

sub to_FTHelper {
    my ($self) = @_;
    my (@out);

    my %contig;
    my $id = $self->id();

    foreach my $trans ( $self->each_Transcript ) {
	foreach my $ptrans ( $trans->split_Transcript_to_Partial ) {
	    my $trans = $ptrans->translate();
	    my $join = "";
	    foreach my $exon ( $ptrans->each_Exon ) {
		if( $exon->seqname eq $id ) {
		    # in this contigs coordinate systems - fine
		    if( $exon->strand == 1 ) {
			$join .= $exon->start."..".$exon->end.",";
		    } else {
			$join .= "complement(".$exon->start."..".$exon->end."),";
		    }
		} else {
		    # in someone else's coordinate system. Yuk.
		    if( !defined $contig{$exon->contig_id} ) {
			$contig{$exon->contig_id} = $self->dbobj->get_Contig($exon->contig_id);
		    }
		    my $tstart = $exon->start + $contig{$exon->contig_id}->embl_offset;
		    my $tend   = $exon->end   + $contig{$exon->contig_id}->embl_offset;
		    my $acc = $contig{$exon->contig_id}->embl_accession;

		    if( $exon->strand == 1 ) {
			$join .= "$acc:".$exon->start."..".$exon->end.",";
		    } else {
			$join .= "complement($acc:".$exon->start."..".$exon->end."),";
		    }
		}
	    }

	    # strip off trailing comma

	    $join =~ s/\,$//g;
	    # build FTHelper object

	    my $ft = Bio::SeqIO::FTHelper->new();
	    $ft->loc($join);
	    $ft->key('CDS');
	    $ft->add_field('translate',$trans->seq);
	    $ft->add_field('cds',$trans->translation->id);
	    push(@out,$ft);
	}
    }

    foreach my $exon ( $self->each_contained_Exon() ) {
	my $ft = Bio::SeqIO::FTHelper->new();

	if( $exon->strand == 1 ) {
	    $ft->loc($exon->start."..".$exon->end);
	} else {
	    $ft->loc("complement(".$exon->start."..".$exon->end.")");
	}

	$ft->key("exon");
	# add other stuff to Exon?
	if ($self->strict_EMBL_dumping) {
	    $ft->add_field('db_xref', 'ENSEMBL:HUMAN-Exon-'. $exon->id);
	} else {
	    $ft->add_field('created',     scalar(gmtime($exon->created())));
	    $ft->add_field('modified',    scalar(gmtime($exon->modified())));
	    $ft->add_field('exon_id',     $exon->id());
	    $ft->add_field('start_phase', $exon->phase());
	    $ft->add_field('end_phase',   $exon->end_phase());
	}

	push(@out,$ft);
    }

    return @out;
}

=head2 strict_EMBL_dumping

 Title   : strict_EMBL_dumping
 Usage   : $obj->strict_EMBL_dumping($newval)
 Function: 
 Returns : value of strict_EMBL_dumping
 Args    : newvalue (optional)


=cut

sub strict_EMBL_dumping{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'strict_EMBL_dumping'} = $value;
    }
    return $obj->{'strict_EMBL_dumping'};

}


		    
1;





