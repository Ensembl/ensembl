
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

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;
use Bio::SeqFeature::Generic;
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
   my ($self,@args) = @_;
   my @out;
   my $contig_id = $self->id();
   # prepare the SQL statement
   my %got;

   my $sth = $self->_dbobj->prepare("select p3.gene from transcript as p3, exon_transcript as p1, exon as p2 where p2.contig = '$contig_id' and p1.exon = p2.id and p3.id = p1.transcript");

   my $res = $sth->execute();
   while( my $rowhash = $sth->fetchrow_hashref) {
       if( $got{$rowhash->{'gene'}} != 1 ) {
          my $gene = $self->_dbobj->get_Gene($rowhash->{'gene'});
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
 Usage   : foreach my $sf ( $contig->get_all_SeqFeatures($start,$end) ) 
 Function: Gets all the sequence features that lie between the coordinates
           specified or on the whole contig.
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SeqFeatures{
   my ($self) = @_;

   my @array;

   my $id     = $self->id();
   my $length = $self->length();

   my %analhash;

   # make the SQL query

   my $sth = $self->_dbobj->prepare("select id,seq_start,seq_end,strand,score,analysis,name,hstart,hend,hid " . 
				    "from feature where contig = '$id'");
   $sth->execute();

   my ($fid,$start,$end,$strand,$score,$analysisid,$name,$hstart,$hend,$hid);

   # bind the columns
   $sth->bind_columns(undef,\$fid,\$start,\$end,\$strand,\$score,\$analysisid,\$name,\$hstart,\$hend,\$hid);

   while( $sth->fetch ) {
       my $out;

       if( $hid ne '__NONE__' ) {
	   # is a paired feature
	   # build EnsEMBL features and make the FeaturePair
	   my $feature1 = new Bio::EnsEMBL::SeqFeature;
	   my $feature2 = new Bio::EnsEMBL::SeqFeature;

	   $out = Bio::EnsEMBL::FeaturePair->new( -feature1 => $feature1, -feature2 => $feature2);

	   $out->hstart     ($hstart);
	   $out->hend       ($hend);
	   $out->hseqname   ($hid);
	   $out->hsource_tag($name);
	   $out->hprimary_tag('similarity');
	   #$out->strand    ($strand);

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

       if( defined $score ) {
	   $out->score($score);
       }

       my $analysis;

       if (!$analhash{$analysisid}) {
	   $analysis = $self->_dbobj->get_Analysis($analysisid);
	   
	   $analhash{$analysisid} = $analysis;

       } else {
	   $analysis = $analhash{$analysisid};
       }

       $out->analysis($analysis);

       # downcast to repeat for repeats. Not pretty.
       #
       if( $out->source_tag() =~ /Repeat/ ) {
	   bless $out, "Bio::EnsEMBL::Repeat";
       }


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

   my $sth = $self->_dbobj->prepare("select created from dna where contig = \"$id\" ");
   $sth->execute();
   my $rowhash = $sth->fetchrow_hashref(); 
   my $datetime = $rowhash->{'created'};
   $sth = $self->_dbobj->prepare("select UNIX_TIMESTAMP('".$datetime."')");
   $sth->execute();
   $rowhash = $sth->fetchrow_arrayref();
   return $rowhash->[0];
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



