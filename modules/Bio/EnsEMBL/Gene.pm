
#
# BioPerl module for Gene
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Gene - Object for confirmed Genes

=head1 SYNOPSIS

Confirmed genes. Basically has a set of transcripts

=head1 DESCRIPTION

Needs more description.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Gene;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::SeqFeature::Generic

use Bio::Root::RootI;
use Bio::EnsEMBL::Transcript;
use Bio::DBLinkContainerI;
use Bio::Annotation::DBLink;


@ISA = qw(Bio::Root::RootI Bio::DBLinkContainerI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub new {
  my($class,@args) = @_;

  my $self = bless {}, $class;

  $self->{'_transcript_array'} = [];
  $self->{'_clone_neighbourhood'} = [];
  $self->{'_db_link'} = [];
# set stuff in self from @args
  return $self; # success - we hope!
}

=head2 is_known

 Title   : is_known
 Usage   : if( $gene->is_known ) 
 Function: returns true if there are any dblinks on the gene or transcript objects
 Example :
 Returns : 
 Args    :


=cut

sub is_known{
   my ($self) = @_;
   my @array;
   @array = $self->each_DBLink();
   if( scalar(@array) > 0 ) {
       return 1;
   }
   foreach my $trans ( $self->each_Transcript ) {
       @array = $trans->each_DBLink();
       if( scalar(@array) > 0 ) {
	   return 1;
       }
   }
       
   return 0;
}

=head2 description

 Title   : description
 Usage   : $gene->description
 Function: get/set description
 Example :
 Returns : a string
 Args    : none

=cut

sub description {
    my ($self, $value) = @_;

    if (defined $value) {
        $self->{description} = $value;
    }
    $self->{description};    
}

=head2 each_DBLink

 Title   : each_DBLink
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_DBLink {
   my ($self,@args) = @_;

   return @{$self->{'_db_link'}}
}


=head2 add_DBLink

 Title   : add_DBLink
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_DBLink{
   my ($self,$value) = @_;

   if( !defined $value || !ref $value || ! $value->isa('Bio::Annotation::DBLink') ) {
       $self->throw("This [$value] is not a DBLink");
   }

   push(@{$self->{'_db_link'}},$value);
}



=head2 each_unique_Exon

 Title   : each_unique_Exon
 Usage   : foreach my $exon ( $gene->each_unique_Exon )
 Function: retrieves an array of exons associated with this
           gene, made nonredudant on the basis of $exon->id
 Example :
 Returns : 
 Args    :


=cut

sub each_unique_Exon{
   my ($self) = @_;

   $self->{_unique_exons}=undef;

   foreach my $trans ( $self->each_Transcript ) {
#       print STDERR "Transcript " . $trans->id . "\n";
       foreach my $exon ( $trans->each_Exon ) {
#	   print STDERR "Found exon $exon " . $exon->id . "\t" . $exon->start . "\t" . $exon->end . "\n";
	   $self->{_unique_exons}{$exon->id()} = $exon;
       }
   }

   return values %{$self->{_unique_exons}};
}

=head2 get_Exon_by_id

 Title   : get_Exon_by_id
 Usage   : $gene->get_Exon($exon_id);
 Function: resolve the id into an exon (or return undef)
 Example :
 Returns : an Exon or undef
 Args    :

=cut

sub get_Exon_by_id {
    my ($self, $id) = @_;

    if (! defined($self->{_unique_exons}) ) {
        $self->each_unique_Exon;
    }
    return $self->{_unique_exons}{$id};
}

=head2 type

 Title   : type
 Usage   : $obj->type($newval)
 Function: 
 Returns : value of type
 Args    : newvalue (optional)


=cut

sub type {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'type'} = $value;
    }
    return $obj->{'type'};
}

=head2 all_Exon_objects

 Title   : all_Exon_objects
 Usage   : foreach $e ( $gene->all_Exon_objects() ) 
 Function: Gives an array of all the exon objects, with each
           object represented only once. This is non redundant
           on the basis of object location, not on id
           (see each_unique_Exon for a difference)
 Example :
 Returns : 
 Args    :


=cut

sub all_Exon_objects{
   my ($self,@args) = @_;
   my %h;

   foreach my $trans ( $self->each_Transcript ) {
       foreach my $exon ( $trans->each_Exon ) {
	   $h{"$exon"} = $exon;
       }
   }

   return values %h;


}

=head2 unique_contig_ids

 Title   : unique_contig_ids
 Usage   : foreach $id ( $gene->unique_contig_ids ) 
 Function: returns an array of contig ids made unique linked
           to this gene
 Example :
 Returns : 
 Args    :


=cut

sub unique_contig_ids{
   my ($self) = @_;
   my %h;

   foreach my $exon ( $self->all_Exon_objects ) {
       if( $exon->isa('Bio::EnsEMBL::StickyExon') ) {
	   foreach my $se ( $exon->each_component_Exon ) {
	       $h{$se->contig_id()} = 1;
	   }
       } else {
	   $h{$exon->contig_id()} = 1;
       }
   }

   return keys %h;
}


=head2 add_Transcript

 Title   : add_Transcript
 Usage   : $gene->add_Transcript($tr)
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Transcript{
   my ($self,$trans) = @_;

   if( !ref $trans || ! $trans->isa("Bio::EnsEMBL::Transcript") ) {
       $self->throw("$trans is not a Bio::EnsEMBL::Transcript!");
   }

   # at the moment, use the SeqFeature sub hash. But in the future,
   # possibly do something better?

   push(@{$self->{'_transcript_array'}},$trans);
}

=head2 each_Transcript

 Title   : each_Transcript
 Usage   : foreach $trans ( $gene->each_Transcript)
 Function:
 Example :
 Returns : An array of Transcript objects
 Args    :


=cut

sub each_Transcript {
   my ($self) = @_;
   my @sub;
   my @ret;

   return @{$self->{'_transcript_array'}};   

}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'id'} = $value;
    }
    return $obj->{'id'};

}


=head2 each_cloneid_neighbourhood

 Title   : each_cloneid_neighbourhood
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_cloneid_neighbourhood{
   my ($self) = @_;

   return @{$self->{'_clone_neighbourhood'}};
}

=head2 add_cloneid_neighbourhood

 Title   : add_cloneid_neighbourhood
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_cloneid_neighbourhood{
   my ($self,$value) = @_;

   if( !defined $value || ref $value ) {
       $self->throw("Value [$value] does not look good for clone neighbourhood id!");
   }

   push(@{$self->{'_clone_neighbourhood'}},$value);
}

=head2 created

 Title   : created
 Usage   : $obj->created($newval)
 Function: 
 Returns : value of created
 Args    : newvalue (optional)


=cut

sub created{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'created'} = $value;
    }
    return $obj->{'created'};

}

=head2 modified

 Title   : modified
 Usage   : $obj->modified($newval)
 Function: 
 Returns : value of modified
 Args    : newvalue (optional)


=cut

sub modified{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'modified'} = $value;
    }
    return $obj->{'modified'};

}

=head2 _stored

 Title   : _stored
 Usage   : $obj->_stored($newval)
 Function: Internal method should not really be needed
           stores the time of storage of the deleted object
 Returns : value of _stored
 Args    : newvalue (optional)


=cut

sub _stored{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_stored'} = $value;
    }
    return $obj->{'_stored'};
}

=head2 version

 Title   : version
 Usage   : $obj->version($newval)
 Function: 
 Returns : value of version
 Args    : newvalue (optional)


=cut

sub version{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'version'} = $value;
    }
    return $obj->{'version'};

}

=head2 analysis

 Title   : analysis
 Usage   : $gene->analysis($analysisObject)
 Function: get/set this genes analysis object
 Returns : on get the analysis object
 Args    : newvalue (optional)


=cut

sub analysis {
  my ($self,$value) = @_;
  if( defined $value ) {
    $self->{'analysis'} = $value;
  }
  return $self->{'analysis'};
}


=head2 _dump

 Title   : _dump
 Usage   : dump data structure for debugging
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _dump{
   my ($self,$fh) = @_;

   if( ! $fh ) {
       $fh = \*STDOUT;
   }

   print $fh "Gene ", $self->id(), "\n";
   foreach my $t ( $self->each_Transcript ) {
       print $fh "  Trans ", $t->id(), " :";
       foreach my $e ( $t->each_Exon ) {
	   print $fh " ",$e->id(),",";
       }
       print "\n";
   }


}


1;
