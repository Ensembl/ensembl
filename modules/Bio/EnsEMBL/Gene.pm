
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

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Transcript;
use Bio::DBLinkContainerI;
use Bio::Annotation::DBLink;
use Bio::EnsEMBL::DBEntry;


@ISA = qw(Bio::EnsEMBL::Root Bio::DBLinkContainerI);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub new {
  my($class,@args) = @_;

  my $self = bless {}, $class;

  $self->{'_transcript_array'} = [];
#  $self->{'_db_link'} = [];
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

=head2 adaptor

 Title   : adaptor
 Usage   :
 Function: give this genes GeneAdaptor if known
 Example :
 Returns : 
 Args    :


=cut

sub adaptor {
   my ($self, $arg) = @_;

   if ( defined $arg ) {
      $self->{'_adaptor'} = $arg ;
   }
   return $self->{'_adaptor'};
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

=head2 dbID

 Title   : dbID
 Usage   :
 Function: internal db id if available
 Example :
 Returns : 
 Args    :


=cut

sub dbID {
   my ($self, $arg) = @_;

   if ( defined $arg ) {
      $self->{'_dbID'} = $arg ;
   }
   return $self->{'_dbID'};
}


=head2 description

 Title   : description
 Usage   : $gene->description
 Function: gets the gene description line. Setting is not allowed
 Example :
 Returns : a string
 Args    : none

=cut

sub description {
    my ($self) = @_;
    
    if( exists $self->{'_description'} ) {
      return $self->{'_description'};
    }
    $self->{'_description'} = $self->adaptor->get_description($self->dbID);
    return $self->{'_description'};
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

   if( !defined $self->{'_db_link'} ) {
       $self->{'_db_link'} = [];
       if( defined $self->adaptor ) {
	 $self->adaptor->db->get_DBEntryAdaptor->fetch_by_gene($self);
       }
   } 


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

   if( !defined $self->{'_db_link'} ) {
       $self->{'_db_link'} = [];
   }

   push(@{$self->{'_db_link'}},$value);
}




sub each_unique_Exon{
   my ($self) = @_;

   my ($p,$f,$l) = caller;
   $self->warn("$f:$l each_unique_Exon deprecated. use get_all_Exons instead. Exon objects should be unique memory locations");

   return $self->get_all_Exons;
}


sub all_Exon_objects{

   my ($self) = @_;

   my ($p,$f,$l) = caller;
   $self->warn("$f:$l all_Exon_objects deprecated. use get_all_Exons instead. Exon objects should be unique memory locations");

   return $self->get_all_Exons;
}



=head2 get_all_Exons

 Title   : get_all_Exons
 Usage   : foreach my $exon ( $gene->each_unique_Exon )
 Function: retrieves an array of exons associated with this
           gene, guarenteed to be nonredundant
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Exons {
   my ($self,@args) = @_;
   my %h;

   foreach my $trans ( $self->each_Transcript ) {
       foreach my $exon ( $trans->get_all_Exons ) {
	   $h{"$exon"} = $exon;
       }
   }

   return values %h;
}

=head2 refresh

 Title   : refresh
 Usage   :
 Function: This function is for cacheing of external genes. It
           refreshs the coordinate system of the underlying exons
           to be what they were retrieved with

           This function may be obselete with new exons. EB and Elia
 Example :
 Returns : 
 Args    :


=cut

sub refresh {
   my ($self) = @_;

   foreach my $e ($self->get_all_Exons) {
       $e->start($e->ori_start);
       $e->end($e->ori_end);
       $e->strand($e->ori_strand);
   }
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

    # perhaps not ideal
    foreach my $exon ( $self->get_all_Exons ) {
      # should this be stable_id
      if( $exon->dbID eq $id ) {
	return $exon;
      }
    }

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
  my $self = shift;
  my $value = shift;

   my ($p,$f,$l) = caller;
   $self->warn("$f:$l id deprecated. Please choose from stable_id or dbID");

  if( defined $value ) {
    $self->warn("$f:$l stable ids are loaded separately and dbIDs are generated on writing. Ignoring set value $value");
    return;
  }


   if( defined $self->stable_id ) {
     return $self->stable_id();
   } else {
     return $self->dbID;
   }

}

=head2 Stable id 

Stable id information is fetched on demand from stable tables

=head2 created

 Title   : created
 Usage   : $obj->created()
 Function: 
 Returns : value of created
 Args    :


=cut

sub created{
    my ($self,$value) = @_;

    if(defined $value ) {
      my ($p,$f,$l) = caller;
      $self->warn("$f $l  created dates are loaded. Ignoring set value $value");
      return;
    }


    if( exists $self->{'_created'} ) {
      return $self->{'_created'};
    }

    $self->_get_stable_entry_info();

    return $self->{'_created'};

}

=head2 modified

 Title   : modified
 Usage   : $obj->modified()
 Function: 
 Returns : value of modified
 Args    : 


=cut

sub modified{
    my ($self,$value) = @_;
    

    if( defined $value ) {
      my ($p,$f,$l) = caller;
      $self->warn("$f $l  modified dates are loaded. Ignoring set value $value");
      return;
    }

    if( exists $self->{'_modified'} ) {
      return $self->{'_modified'};
    }

    $self->_get_stable_entry_info();

    return $self->{'_modified'};
}


=head2 version

 Title   : version
 Usage   : $obj->version()
 Function: 
 Returns : value of version
 Args    : 

=cut

sub version{

    my ($self,$value) = @_;
    

    if( defined $value ) {
      my ($p,$f,$l) = caller;
      $self->warn("$f $l  modified dates are loaded. Ignoring set value $value");
      return;
    }

    if( exists $self->{'_version'} ) {
      return $self->{'_version'};
    }

    $self->_get_stable_entry_info();

    return $self->{'_version'};

}


=head2 stable_id

 Title   : stable_id
 Usage   : $obj->stable_id
 Function: 
 Returns : value of stable_id
 Args    : 


=cut

sub stable_id{

    my ($self,$value) = @_;
    

    if( defined $value ) {
      $self->{'_stable_id'} = $value;
      return;
    }

    if( exists $self->{'_stable_id'} ) {
      return $self->{'_stable_id'};
    }

    $self->_get_stable_entry_info();

    return $self->{'_stable_id'};

}

sub _get_stable_entry_info {
   my $self = shift;

   if( !defined $self->adaptor ) {
     return undef;
   }

   $self->adaptor->get_stable_entry_info($self);

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

=head2 temporary_id

 Title   : temporary_id
 Usage   : $obj->temporary_id($newval)
 Function: Temporary ids are used for Genscan predictions - which should probably
           be moved over to being stored inside the gene tables anyway. Bio::EnsEMBL::TranscriptFactory use this.
           MC Over my dead body they will.  Unless you can speed up the database by a couple of orders of magnitude.
 Example : 
 Returns : value of temporary_id
 Args    : newvalue (optional)


=cut

sub temporary_id {
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'temporary_id'} = $value;
    }
    return $obj->{'temporary_id'};

}


1;
