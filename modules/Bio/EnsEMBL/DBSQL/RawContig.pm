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

Bio::EnsEMBL::DB::RawContig - Handle onto a database stored raw contiguous DNA 

=head1 SYNOPSIS

    # get a contig object somehow, eg from an DB::Obj

    @genes = $contig->get_all_Genes();
    @sf    = $contig->get_all_RepeatFeatures();
    @sf    = $contig->get_all_SimilarityFeatures();

    $contig->id();
    $contig->length();
    $primary_seq = $contig->primary_seq();

=head1 DESCRIPTION

A RawContig is physical piece of DNA coming out a sequencing project,
ie a single product of an assembly process. A RawContig defines an atomic
coordinate system on which features and genes are placed (remember that
genes can cross between atomic coordinate systems).

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::RawContig;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Object;

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DB::RawContigI;

use Bio::EnsEMBL::Repeat;
use Bio::EnsEMBL::ContigOverlap;
use Bio::EnsEMBL::FeatureFactory;
use Bio::PrimarySeq;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::RawContigI);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  my ($dbobj,$id) = $self->_rearrange([qw(DBOBJ
					  ID
					  )],@args);

  $id    || $self->throw("Cannot make contig db object without id");
  $dbobj || $self->throw("Cannot make contig db object without db object");
  $dbobj->isa('Bio::EnsEMBL::DBSQL::Obj') || $self->throw("Cannot make contig db object with a $dbobj object");

  $self->id($id);
  $self->_dbobj($dbobj);
  $self->_got_overlaps(0);

# set stuff in self from @args
  return $make; # success - we hope!
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
   my ($self,$supporting) = @_;
   my @out;
   my $contig_id = $self->internal_id();
   # prepare the SQL statement
   my %got;
   my $gene;
   
   print STDERR "Getting into get_all_Genes";

   my $sth = $self->_dbobj->prepare("select p3.gene from transcript as p3, exon_transcript as p1, exon as p2 where p2.contig = '$contig_id' and p1.exon = p2.id and p3.id = p1.transcript");

   my $res = $sth->execute();
   while( my $rowhash = $sth->fetchrow_hashref) {
       if( ! exists $got{$rowhash->{'gene'}}  ) {
	   if ($supporting && $supporting eq 'evidence') {
	       $gene = $self->_dbobj->get_Gene($rowhash->{'gene'},'evidence');
	   }
	   else {
	       $gene = $self->_dbobj->get_Gene($rowhash->{'gene'});
	   }
	   print STDERR "Got a $gene in get all Genes\n";

	   push(@out,$gene);
	   $got{$rowhash->{'gene'}} = 1;
       }
       
   }
   

   return @out;

}


=head2 primary_seq

 Title   : seq
 Usage   : $seq = $contig->primary_seq();
 Function: Gets a Bio::PrimarySeqI object out from the contig
 Example :
 Returns : Bio::PrimarySeqI object
 Args    :


=cut

sub primary_seq {
   my ($self) = @_;
   my $id = $self->internal_id();

   if( $self->_seq_cache() ) {
       return $self->_seq_cache();
   }

   my $sth = $self->_dbobj->prepare("select d.sequence from dna as d,contig as c where c.internal_id = $id and c.dna = d.id");
   my $res = $sth->execute();

   # should be a better way of doing this
   while(my $rowhash = $sth->fetchrow_hashref) {
     my $str = $rowhash->{sequence};

     if( ! $str) {
       $self->throw("No DNA sequence in contig " . $self->id . " " . $id);
     } 

     $str =~ /[^ATGCNRY]/ && $self->warn("Got some non standard DNA characters here! Yuk!");
     $str =~ s/\s//g;
     $str =~ s/[^ATGCNRY]/N/g;

     my $ret =Bio::PrimarySeq->new ( -seq => $str, -id => $id, -moltype => 'dna' );
     $self->_seq_cache($ret);
     
     return $ret;
   }


   $self->throw("No dna sequence associated with $id!");
   
}

=head2 _seq_cache

 Title   : _seq_cache
 Usage   : $obj->_seq_cache($newval)
 Function: 
 Returns : value of _seq_cache
 Args    : newvalue (optional)


=cut

sub _seq_cache{
   my $obj = shift;
   if( @_ ) {
       my $value = shift;
       $obj->{'_seq_cache'} = $value;
   }
   return $obj->{'_seq_cache'};

}

=head2 get_all_SeqFeatures

 Title   : get_all_SeqFeatures
 Usage   : foreach my $sf ( $contig->get_all_SeqFeatures
 Function: Gets all the sequence features on the whole contig
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SeqFeatures {
    my ($self) = @_;

    my @out;

    push(@out,$self->get_all_SimilarityFeatures);
    push(@out,$self->get_all_RepeatFeatures);
#   push(@out,$self->get_all_PredictionFeatures);

    print(STDERR "Fetched all features\n");
    return @out;
}

=head2 get_all_SimilarityFeatures

 Title   : get_all_SimilarityFeatures
 Usage   : foreach my $sf ( $contig->get_all_SimilarityFeatures($start,$end) ) 
 Function: Gets all the sequence similarity features on the whole contig
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SimilarityFeatures{
   my ($self) = @_;

   my @array;

   my $id     = $self->internal_id();
   my $length = $self->length();

   my %analhash;

   #First of all, get all features that are part of a feature set

   my $sth = $self->_dbobj->prepare("select  p1.id, p1.seq_start, p1.seq_end, " . 
				    "p1.strand,p1.score,p1.analysis,p1.name,  " .
				    "p1.hstart,p1.hend,p1.hid,p2.fset,p2.rank " . 
				    "from feature as p1, fset_feature as p2 where " .
				    "p1.contig ='$id' and p2.feature = p1.id order by p2.fset");
   $sth->execute();

   my ($fid,$start,$end,$strand,$score,$analysisid,$name,$hstart,$hend,$hid,$fset,$rank);
   my $seen = 0;
   
   # bind the columns
   $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$score,\$analysisid,\$name,\$hstart,\$hend,\$hid,\$fset,\$rank);

   my $out;
   
   my $fset_id_str = "";

   while($sth->fetch) {

       my $analysis;

       if (!$analhash{$analysisid}) {
	   $analysis = $self->_dbobj->get_Analysis($analysisid);
	   $analhash{$analysisid} = $analysis;
	   
       } else {
	   $analysis = $analhash{$analysisid};
       }
       
       if( !defined $name ) {
	   $name = 'no_source';
       }
       
       #Build fset feature object if new fset found
       if ($fset != $seen) {
	   print(STDERR "Making new fset feature $fset\n");
	   $out =  new Bio::EnsEMBL::SeqFeature;
	   $out->id($fset);
	   $out->analysis($analysis);
	   $seen = $fset;
	   push(@array,$out);
       }
       $fset_id_str = $fset_id_str . $fid . ",";       
       #Build Feature Object
       my $feature = new Bio::EnsEMBL::SeqFeature;
       $feature->seqname   ($id);
       $feature->start     ($start);
       $feature->end       ($end);
       $feature->strand    ($strand);
       $feature->source_tag($name);
       $feature->primary_tag('similarity');
       $feature->id         ($fid);
       
       if( defined $score ) {
	   $feature->score($score);
       }
       
       $feature->analysis($analysis);
       
       # Final check that everything is ok.
       $feature->validate();

       #Add this feature to the fset
       $out->add_sub_SeqFeature($feature,'EXPAND');

   }
   
   #Then get the rest of the features, i.e. featurepairs and single features that are not part of a fset
   $fset_id_str =~ s/\,$//;

   if ($fset_id_str) {
       $sth = $self->_dbobj->prepare("select id,seq_start,seq_end,strand,score,analysis,name,hstart,hend,hid " .
				     "from feature where id not in (" . $fset_id_str . ") and contig = \"$id\"");
   } else {
       $sth = $self->_dbobj->prepare("select id,seq_start,seq_end,strand,score,analysis,name,hstart,hend,hid ".
				     "from feature where contig = \"$id\"");
   }

   $sth->execute();

   # bind the columns
   $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$score,\$analysisid,\$name,\$hstart,\$hend,\$hid);
   
   while($sth->fetch) {
       my $out;
       my $analysis;
              
       if (!$analhash{$analysisid}) {
	   $analysis = $self->_dbobj->get_Analysis($analysisid);
	   $analhash{$analysisid} = $analysis;
	   
       } else {
	   $analysis = $analhash{$analysisid};
       }
       
       if( !defined $name ) {
	   $name = 'no_source';
       }
       
       if( $hid ne '__NONE__' ) {
	   # is a paired feature
	   # build EnsEMBL features and make the FeaturePair
	 
	   $out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();


	   $out->set_all_fields($start,$end,$strand,$score,$self->id(),'similarity',$id,
				$hstart,$hend,1,$score,$name,'similarity',$hid);

	   $out->analysis    ($analysis);
       } else {
	   $out = new Bio::EnsEMBL::SeqFeature;
	   $out->seqname   ($id);
	   $out->start     ($start);
	   $out->end       ($end);
	   $out->strand    ($strand);
	   $out->source_tag($name);
	   $out->primary_tag('similarity');
	   $out->id         ($fid);

	   if( defined $score ) {
	       $out->score($score);
	   }
	   $out->analysis($analysis);
       }
       # Final check that everything is ok.
       $out->validate();
       
      push(@array,$out);
      
   }
   
   return @array;
}

=head2 get_all_RepeatFeatures

 Title   : get_all_RepeatFeatures
 Usage   : foreach my $sf ( $contig->get_all_RepeatFeatures )
 Function: Gets all the repeat features on a contig.
 Example :
 Returns : 
 Args    :


=cut

sub get_all_RepeatFeatures {
   my ($self) = @_;

   my @array;

   my $id     = $self->internal_id();
   my $length = $self->length();

   my %analhash;

   # make the SQL query

   my $sth = $self->_dbobj->prepare("select id,seq_start,seq_end,strand,score,analysis,hstart,hend,hid " . 
				    "from repeat_feature where contig = '$id'");

   $sth->execute();

   my ($fid,$start,$end,$strand,$score,$analysisid,$hstart,$hend,$hid);

   # bind the columns
   $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$score,\$analysisid,\$hstart,\$hend,\$hid);

   while( $sth->fetch ) {
       my $out;
       my $analysis;

       if (!$analhash{$analysisid}) {
	   $analysis = $self->_dbobj->get_Analysis($analysisid);
	   $analhash{$analysisid} = $analysis;

       } else {
	   $analysis = $analhash{$analysisid};
       }


       if( $hid ne '__NONE__' ) {
	   # is a paired feature
	   # build EnsEMBL features and make the FeaturePair

	   $out = Bio::EnsEMBL::FeatureFactory->new_repeat();

	   $out->set_all_fields($start,$end,$strand,$score,'repeatmasker','repeat',$id,
				$hstart,$hend,1,$score,'repeatmasker','repeat',$hid);

	   $out->analysis($analysis);

       } else {
	   $self->warn("Repeat feature does not have a hid. bad news....");
       }
       
       $out->validate();

      push(@array,$out);
  }
 
   return @array;
}
=head2 get_all_RepeatFeatures

 Title   : get_all_RepeatFeatures
 Usage   : foreach my $sf ( $contig->get_all_RepeatFeatures )
 Function: Gets all the repeat features on a contig.
 Example :
 Returns : 
 Args    :


=cut

sub get_all_PredictionFeatures {
   my ($self) = @_;

   my @array;

   my $id     = $self->internal_id();
   my $length = $self->length();

   my %analhash;

   # make the SQL query

   my $sth = $self->_dbobj->prepare("select id,seq_start,seq_end,strand,score,analysis" . 
				    "from feature where contig = '$id' and name = 'genscan'");
   
   $sth->execute();
   
   my ($fid,$start,$end,$strand,$score,$analysisid);
   
   # bind the columns
   $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$score,\$analysisid);
   
   while( $sth->fetch ) {
       my $out;
       my $analysis;
       
       if (!$analhash{$analysisid}) {
	   $analysis = $self->_dbobj->get_Analysis($analysisid);
	   $analhash{$analysisid} = $analysis;
	   
       } else {
	   $analysis = $analhash{$analysisid};
       }


       $out = new Bio::EnsEMBL::SeqFeature;
       
       $out->seqname   ($id);
       $out->start     ($start);
       $out->end       ($end);
       $out->strand    ($strand);

       $out->source_tag('genscan');
       $out->primary_tag('prediction');
       
       if( defined $score ) {
	   $out->score($score);
       }

       $out->analysis($analysis);

       # Final check that everything is ok.
       
       $out->validate();

      push(@array,$out);
  }
 
   return @array;
}

=head2 length

 Title   : length
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub length{
   my ($self,@args) = @_;
   my $id= $self->internal_id();

   if (!defined($self->{_length})) {
       my $sth = $self->_dbobj->prepare("select length from contig where id = \"$id\" ");
       $sth->execute();
       
       my $rowhash = $sth->fetchrow_hashref();
       
       $self->{_length} = $rowhash->{'length'};
   }

   return $self->{_length};
       
}


=head2 seq_version

 Title   : seq_version
 Usage   : $obj->seq_version($newval)
 Function: 
 Example : 
 Returns : value of seq_version
 Args    : newvalue (optional)


=cut

sub seq_version{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'seq_version'} = $value;
    }
    return $obj->{'seq_version'};

}

=head2 embl_order

 Title   : order
 Usage   : $obj->embl_order
 Function: 
 Returns : 
 Args    : 


=cut

sub embl_order{
   my $self = shift;
   my $id = $self->id();
   my $sth = $self->_dbobj->prepare("select corder from contig where id = \"$id\" ");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref();
   return $rowhash->{'corder'};
   
}




=head2 embl_offset

 Title   : embl_offset
 Usage   : 
 Returns : 
 Args    :


=cut

sub embl_offset{
   my $self = shift;
   my $id = $self->id();


   my $sth = $self->_dbobj->prepare("select offset from contig where id = \"$id\" ");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref();
   return $rowhash->{'offset'};

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
   if( defined $value) {
      $self->{'id'} = $value;
    }
    return $self->{'id'};

}

=head2 internal_id

 Title   : internal_id
 Usage   : $obj->internal_id($newval)
 Function: 
 Example : 
 Returns : value of database internal id
 Args    : newvalue (optional)


=cut

sub internal_id {
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'internal_id'} = $value;
    }
    return $self->{'internal_id'};

}

=head2 seq_date

 Title   : seq_date
 Usage   : $contig->seq_date()
 Function: Gives the unix time value of the dna table created datetime field, which indicates
           the original time of the dna sequence data
 Example : $contig->seq_date()
 Returns : unix time
 Args    : none


=cut

sub seq_date{
   my ($self) = @_; 

   my $id = $self->internal_id();

   my $sth = $self->_dbobj->prepare("select UNIX_TIMESTAMP(created) from dna where sequence = \"$id\" ");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref(); 
   return $rowhash->{'UNIX_TIMESTAMP(created)'};
}


=head2 get_left_overlap

 Title   : get_left_overlap
 Usage   : $overlap_object = $contig->get_left_overlap();
 Function: Returns the overlap object of contig to the left.
           This could be undef, indicating no overlap
 Returns : A Bio::EnsEMBL::ContigOverlap object
 Args    : None

=cut

sub get_left_overlap{
   my ($self,@args) = @_;
   if( $self->_got_overlaps == 0 ) {
       $self->_load_overlaps() ;
   }

   return $self->_left_overlap();
}


=head2 get_right_overlap

 Title   : get_right_overlap
 Usage   : $overlap_object = $contig->get_right_overlap();
 Function: Returns the overlap object of contig to the left.
           This could be undef, indicating no overlap
 Returns : A Bio::EnsEMBL::ContigOverlap object
 Args    : None

=cut

sub get_right_overlap{
   my ($self,@args) = @_;

   if( $self->_got_overlaps == 0 ) {
       $self->_load_overlaps() ;
   }

   return $self->_right_overlap();
}

=head2 _got_overlaps

 Title   : _got_overlaps
 Usage   : $obj->_got_overlaps($newval)
 Function: 
 Returns : value of _got_overlaps
 Args    : newvalue (optional)


=cut

sub _got_overlaps{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_got_overlaps'} = $value;
    }
    return $obj->{'_got_overlaps'};

}

=head2 _load_overlaps

 Title   : _load_overlaps
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _load_overlaps{
   my ($self,@args) = @_;

   my $id = $self->internal_id();
   my $version = $self->seq_version();

   print STDERR "select contig_a,contig_b,contig_a_position,contig_b_position,overlap_type from contigoverlap where (contig_a = '$id' and contig_a_version = $version ) or (contig_b = '$id' and contig_b_version = $version )\n";

   my $sth = $self->_dbobj->prepare("select contig_a,contig_b,contig_a_position,contig_b_position,overlap_type from contigoverlap where (contig_a = '$id' and contig_a_version = $version ) or (contig_b = '$id' and contig_b_version = $version )");
   
   $sth->execute();

   #
   # Certainly worth explaining here.
   #
   # The overlap type indicates which end on the two contigs this overlap is.
   # left means 5', right means 3'. There are four options. From these four
   # options we can figure out
   #    a) whether this overlap is on the 5' or the 3' of our contig
   #    b) which polarity the overlap on the next contig is 
   #
   # Polarity indicates whether the sequence is being read in the same 
   # direction as this contig. 
   # 
   # The sequence has to be appropiately versioned otherwise this gets complicated
   # in the update scheme.

   # 
   # start by switching on whether things are in the a or b contig
   # positions, then builds left/right and polarity variables.
   #

   while( my $rowhash = $sth->fetchrow_hashref ) {

       my $contigid = $rowhash->{'id'};

       if( $rowhash->{'contig_a'} eq $id) {
	   my $t = $rowhash->{'overlap_type'};
	   my ($selflr,$sisterpol);
	   if( $t eq 'right2left' ) {
	       $selflr = 'right';
	       $sisterpol = 1;
	   } elsif( $t eq 'right2right' ) {
	       $selflr = 'right';
	       $sisterpol = -1;
	   } elsif( $t eq 'left2right' ) {
	       $selflr = 'left';
	       $sisterpol = 1;
	   } elsif ( $t eq 'left2left' ) {
	       $selflr = 'left';
	       $sisterpol = -1;
	   } else {
	       $self->throw("Impossible type position $t\n");
	   }

	   my $sis = $self->_dbobj->get_Contig($contigid);
	   my $co = Bio::EnsEMBL::ContigOverlap->new(
						     -sister => $sis,
						     -sisterposition => $rowhash->{'contig_b_position'}, 
						     -selfposition => $rowhash->{'contig_a_position'},
						     -sisterpolarity => $sisterpol );
	   if( $selflr eq 'left' ) {
	       $self->_left_overlap($co);
	   } else {
	       $self->_right_overlap($co);
	   }

       } else {
	   my $t = $rowhash->{'overlap_type'};
	   my ($selflr,$sisterpol);
	   if( $t eq 'right2left' ) {
	       $selflr = 'left';
	       $sisterpol = 1;
	   } elsif( $t eq 'right2right' ) {
	       $selflr = 'right';
	       $sisterpol = -1;
	   } elsif( $t eq 'left2right' ) {
	       $selflr = 'right';
	       $sisterpol = 1;
	   } elsif ( $t eq 'left2left' ) {
	       $selflr = 'left';
	       $sisterpol = -1;
	   } else {
	       $self->throw("Impossible type position $t\n");
	   }
	   my $sis = $self->_dbobj->get_Contig($contigid);

	   my $co = Bio::EnsEMBL::ContigOverlap->new(
						     -sister => $sis,
						     -sisterposition => $rowhash->{'contig_a_position'}, 
						     -selfposition => $rowhash->{'contig_b_position'},
						     -sisterpolarity => $sisterpol );
	   if( $selflr eq 'left' ) {
	       $self->_left_overlap($co);
	   } else {
	       $self->_right_overlap($co);
	   }
       }
   }

   $self->_got_overlaps(1);

}

=head2 _right_overlap

 Title   : _right_overlap
 Usage   : $obj->_right_overlap($newval)
 Function: 
 Example : 
 Returns : value of _right_overlap
 Args    : newvalue (optional)


=cut

sub _right_overlap{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_right_overlap'} = $value;
    }
    return $obj->{'_right_overlap'};

}

=head2 _left_overlap

 Title   : _left_overlap
 Usage   : $obj->_left_overlap($newval)
 Function: 
 Example : 
 Returns : value of _left_overlap
 Args    : newvalue (optional)


=cut

sub _left_overlap{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'_left_overlap'} = $value;
    }
    return $obj->{'_left_overlap'};

}


=head2 _dbobj

 Title   : _dbobj
 Usage   : $obj->_dbobj($newval)
 Function: 
 Example : 
 Returns : value of _dbobj
 Args    : newvalue (optional)


=cut

sub _dbobj{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_dbobj'} = $value;
    }
    return $self->{'_dbobj'};

}



1;
