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

Bio::EnsEMBL::DB::Contig - Handle onto a database stored contig

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::Contig;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::Object;

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DB::ContigI;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Repeat;

@ISA = qw(Bio::Root::Object Bio::EnsEMBL::DB::ContigI);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  my ($dbobj,$id) = $self->_rearrange([qw(DBOBJ
					  ID
					  )],@args);

  $id || $self->throw("Cannot make contig db object without id");
  $dbobj || $self->throw("Cannot make contig db object without db object");
  $dbobj->isa('Bio::EnsEMBL::DBSQL::Obj') || $self->throw("Cannot make contig db object with a $dbobj object");

  $self->id($id);
  $self->_dbobj($dbobj);

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
   my $contig_id = $self->id();
   # prepare the SQL statement
   my %got;
   my $gene;

   my $sth = $self->_dbobj->prepare("select p3.gene from transcript as p3, exon_transcript as p1, exon as p2 where p2.contig = '$contig_id' and p1.exon = p2.id and p3.id = p1.transcript");

   my $res = $sth->execute();
   while( my $rowhash = $sth->fetchrow_hashref) {
       if( $got{$rowhash->{'gene'}} != 1 ) {
	   if ($supporting && $supporting eq 'evidence') {
	       $gene = $self->_dbobj->get_Gene($rowhash->{'gene'},'evidence');
	   }
	   else {
	       $gene = $self->_dbobj->get_Gene($rowhash->{'gene'});
	   }
	   push(@out,$gene);
	   $got{$rowhash->{'gene'}} = 1;
       }
       
   }
   

   return @out;

}


=head2 seq

 Title   : seq
 Usage   : $seq = $contig->seq();
 Function: Gets a Bio::Seq object out from the contig
 Example :
 Returns : Bio::Seq object
 Args    :


=cut

sub seq{
   my ($self) = @_;
   my $id = $self->id();

   if( $self->_seq_cache() ) {
       return $self->_seq_cache();
   }

   my $sth = $self->_dbobj->prepare("select sequence from dna where contig = \"$id\"");
   my $res = $sth->execute();

   # should be a better way of doing this
   while(my $rowhash = $sth->fetchrow_hashref) {
     my $str = $rowhash->{sequence};

     if( ! $str) {
       $self->throw("No DNA sequence in contig $id");
     } 
     $str =~ /[^ATGCNRY]/ && $self->warn("Got some non standard DNA characters here! Yuk!");
     $str =~ s/\s//g;
     $str =~ s/[^ATGCNRY]/N/g;

     my $ret =Bio::Seq->new ( -seq => $str, -id => $id, -type => 'Dna' );
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
    push(@out,$self->get_all_PredictionFeatures);

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

   my $id     = $self->id();
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
       # skip if this is actually a prediction feature

       # FIXME: this is not good.

       if( $analysis->program eq 'Genscan' ) {
	   next;
       }

       print STDERR "Found analysis program type of ",$analysis->program,"\n";

       
       if( !defined $name ) {
	   $name = 'no_source';
       }
       
       #Build fset feature object if new fset found
       if ($fset != $seen) {
	   #print("Making new fset feature $fset\n");
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

       if( $analysis->program eq 'Genscan' ) {
	   next;
       }
       
       if( !defined $name ) {
	   $name = 'no_source';
       }
       
       if( $hid ne '__NONE__' ) {
	   # is a paired feature
	   # build EnsEMBL features and make the FeaturePair
	   my $feature1 = new Bio::EnsEMBL::SeqFeature;
	   my $feature2 = new Bio::EnsEMBL::SeqFeature;
	   
	   $out = Bio::EnsEMBL::FeaturePair->new( -feature1 => $feature1, 
						  -feature2 => $feature2);

	   $out->hstart     ($hstart);
	   $out->hend       ($hend);
	   $out->hseqname   ($hid);
	   $out->hsource_tag($name);
	   $out->hprimary_tag('similarity');
	   $out->hstrand     ($strand);
	   $out->analysis    ($analysis);

	   if( defined $score ) {
	       $out->hscore($score);
	   }

       } else {
	   $out = new Bio::EnsEMBL::SeqFeature;
       }

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

   my $id     = $self->id();
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
	   my $feature1 = new Bio::EnsEMBL::SeqFeature;
	   my $feature2 = new Bio::EnsEMBL::SeqFeature;

	   $out = Bio::EnsEMBL::Repeat->new( -feature1 => $feature1, 
					     -feature2 => $feature2);

	   $out->hstart      ($hstart);
	   $out->hend       ($hend);
	   $out->hseqname   ($hid);
	   $out->hsource_tag('repeat');
	   $out->hprimary_tag('similarity');
	   $out->hstrand     ($strand);
	   $out->analysis    ($analysis);

	   if( defined $score ) {
	       $out->hscore($score);
	   }

       } else {
#	   $out = new Bio::EnsEMBL::SeqFeature;
	   $self->warn("Repeat feature doesn't have hid. Skipping\n");
       }

      
       $out->seqname   ($id);
       $out->start     ($start);
       $out->end       ($end);
       $out->strand    ($strand);
       $out->source_tag('repeat');
       $out->primary_tag('similarity');

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

   my $id     = $self->id();
   my $length = $self->length();

   my %analhash;

   # make the SQL query

   my $sth = $self->_dbobj->prepare("select  p1.id, p1.seq_start, p1.seq_end, " . 
				    "p1.strand,p1.score,p1.analysis,p1.name,  " .
				    "p1.hstart,p1.hend,p1.hid,p2.fset,p2.rank " . 
				    "from feature as p1, fset_feature as p2 analysis as p3 where " .
				    "p1.contig ='$id' and p2.feature = p1.id and p3.id = p1.analysis ".
                                    "and p3.program = 'Genscan' order by p2.fset");
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
       # skip if this is actually a prediction feature

       # FIXME: this is not good.

       if( $analysis->program ne 'Genscan' ) {
	   next;
       }

       
       if( !defined $name ) {
	   $name = 'no_source';
       }
       
       #Build fset feature object if new fset found
       if ($fset != $seen) {
	   #print("Making new fset feature $fset\n");
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
   my $id= $self->id();

   if (!defined($self->{_length})) {
       my $sth = $self->_dbobj->prepare("select length from contig where id = \"$id\" ");
       $sth->execute();
       
       my $rowhash = $sth->fetchrow_hashref();
       
       $self->{_length} = $rowhash->{'length'};
   }

   return $self->{_length};
       
}


=head2 order

 Title   : order
 Usage   : $obj->order($newval)
 Function: 
 Returns : value of order
 Args    : newvalue (optional)


=cut

sub order{
   my $self = shift;
   my $id = $self->id();
   my $sth = $self->_dbobj->prepare("select corder from contig where id = \"$id\" ");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref();
   return $rowhash->{'corder'};
   
}

=head2 offset

 Title   : offset
 Usage   : 
 Returns : 
 Args    :


=cut

sub offset{
   my $self = shift;
   my $id = $self->id();

   my $sth = $self->_dbobj->prepare("select offset from contig where id = \"$id\" ");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref();
   return $rowhash->{'offset'};

}


=head2 orientation

 Title   : orientation
 Usage   : 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub orientation{
   my ($self) = @_;
   my $id = $self->id();

   my $sth = $self->_dbobj->prepare("select orientation from contig where id = \"$id\" ");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref();
   return $rowhash->{'orientation'};
}


=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Example : 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'id'} = $value;
    }
    return $self->{'id'};

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

   my $id = $self->id();

   my $sth = $self->_dbobj->prepare("select UNIX_TIMESTAMP(created) from dna where contig = \"$id\" ");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref(); 
   return $rowhash->{'UNIX_TIMESTAMP(created)'};
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



