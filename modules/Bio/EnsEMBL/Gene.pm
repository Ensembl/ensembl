
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
use Bio::EnsEMBL::TranscriptI;
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


=head2 start
  Title    : start
  Usage    : $start = $gene->start()
  Function : Gets/Sets the lowest start coordinate in of this genes exons 
             in slice coordinates. No consistancy check is performed if this 
             is used as a setter and potentially the start could be set to 
             a value which does not correspond to the lowest  exon start.
             Will throw an exception if called on a gene containing exons
             in raw contig coordinates.
  Returns  : The start of this gene in slice coordinates
  Args     : none, or the start coordinate of this gene
=cut

sub start {
  my($self, $start) = @_;

  if($start) {
    $self->{start} = $start;   
  } elsif(!defined $self->{start}) {

    # calculate the start of this gene
    $self->{start} = -1;
    foreach my $exon ($self->get_all_Exons()) {
      
      unless(defined $exon->contig() && 
             $exon->contig()->isa('Bio::EnsEMBL::Slice')) {
        $self->throw('Gene start cannot be calculated in slice coordinates' .
                     'because exons are in raw contig coordinates');
      }

      if($self->{start} == -1 || $exon->start() < $self->{start}) {
        $self->{start} = $exon->start();
      }

      if($self->{start} == -1) {
        $self->throw('Gene start cannot be calculated because gene does not' .
                     'contain exons with start coordinates');
      }
    }
  }
  return $self->{start};
}


=head2 end
  Title    : end
  Usage    : $end = $gene->end()
  Function : Gets/Sets the highest end coordinate in of this genes exons 
             in slice coordinates. No consistancy check is performed if this 
             is used as a setter and potentially the end could be set to 
             a value which does not correspond to the highest exon end.
             Will throw an exception if called on a gene containing exons
             in raw contig coordinates.
  Returns  : The end of this gene in slice coordinates
  Args     : none, or the end coordinate of this gene
=cut

sub end {
  my($self, $end) = @_;

  if($end) {
    $self->{end} = $end;   
  } elsif(!defined $self->{end}) {

    # calculate the end of this gene
    $self->{end} = -1;
    foreach my $exon ($self->get_all_Exons()) {
      
      unless(defined $exon->contig() && 
             $exon->contig()->isa('Bio::EnsEMBL::Slice')) {
        $self->throw('Gene end cannot be calculated in slice coordinates' .
                     'because exons are in raw contig coordinates');
      }

      if($self->{end} == -1 || $exon->end() > $self->{end}) {
        $self->{end} = $exon->end();
      }

      if($self->{end} == -1) {
        $self->throw('Gene end cannot be calculated because gene does not' .
                     'contain exons with end coordinates');
      }
    }
  }

  return $self->{end};
}


sub strand {
  my $self = shift;
  my $arg = shift;

  if( defined $arg ) {
    $self->{'strand'} = $arg;
  } elsif( ! defined $self->{strand} ) {
    $self->warn( "Gene strand not set, difficult to calculate..." );
  }
  return $self->{'strand'};

}

sub source {
  my ($self, $source) = @_;

  if(defined $source) {
    $self->{'_source'} = $source;
  }

  return $self->{'_source'};
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
  foreach my $trans ( $self->get_all_Transcripts ) {
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
 Usage   : $id = $obj->dbID();
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


=head2 external_name

 Title   : external_name
 Usage   : $ext_name = $obj->external_name();
 Function: external_name if available
 Example : 
 Returns : the external name of this gene
 Args    : new external name (optional)

=cut

sub external_name {
  my ($self, $arg ) = @_;

  if( defined $arg ) {
    $self->{'_external_name'} = $arg;
  }

  return $self->{'_external_name'};
}

sub external_db {
  my ($self, $arg ) = @_;

  if( defined $arg ) {
    $self->{'_external_db'} = $arg;
  }

  return $self->{'_external_db'};
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





=head2 get_all_Exons

 Title   : get_all_Exons
 Usage   : foreach my $exon ( $gene->each_unique_Exon )
 Function: retrieves an array of exons associated with this
           gene, guaranteed to be nonredundant
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Exons {
   my ($self,@args) = @_;
   my %h;

   foreach my $trans ( $self->get_all_Transcripts ) {
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

   if( !ref $trans || ! $trans->isa("Bio::EnsEMBL::TranscriptI") ) {
       $self->throw("$trans is not a Bio::EnsEMBL::TranscriptI!");
   }

   #invalidate the start and end since they may need to be recalculated
   $self->{start} = undef;
   $self->{end} = undef;
   $self->{strand} = undef;

   push(@{$self->{'_transcript_array'}},$trans);
}


=head2 get_all_Transcripts

 Title   : get_all_Transcripts
 Usage   : foreach $trans ( $gene->get_all_Transcripts)
 Function:
 Example :
 Returns : An array of Transcript objects
 Args    :

=cut

sub get_all_Transcripts {
  my ($self) = @_;

  return @{$self->{'_transcript_array'}};
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

   print $fh "Gene ", $self->dbID(), "\n";
   foreach my $t ( $self->get_all_Transcripts() ) {
       print $fh "  Trans ", $t->dbID(), " :";
       foreach my $e ( $t->get_all_Exons ) {
	   print $fh " ",$e->dbID(),",";
       }
       print "\n";
   }

}


=head2 transform

  Arg  1    : Bio::EnsEMBL::Slice $slice
              make this Gene slice coords
              Without a $slice, transform to rawContig Coords
  Function  : make slice coords from raw contig coords or vice versa
              changes gene in place and returns itself.
  Returntype: Bio::EnsEMBL::Gene
  Exceptions: none
  Caller    : object::methodname or just methodname

=cut


sub transform {
  my $self = shift;
  my $slice = shift;

  # hash arrray to store the refs of transformed exons
  my %exon_transforms;

  # transform Exons
  for my $exon ( $self->get_all_Exons() ) {
    my $newExon = $exon->transform( $slice );
    $exon_transforms{ $exon } = $newExon;
  }

  # now need to re-jiggle the transcripts and their
  # translations to account for the re-mapping process

  for my $transcript ( $self->get_all_Transcripts() ) {

    # need to grab the translation before starting to 
    # re-jiggle the exons

    $transcript->transform( \%exon_transforms );
    
  }

  #unset the start, end, and strand - they need to be recalculated
  $self->{start} = undef;
  $self->{end} = undef;
  $self->{strand} = undef;

  return $self;
}


=head2 _set_adaptors

  Arg  1    : Bio::EnsEMBL::DBSQL::DBAdaptor $dba
  Function  : Set adaptors for Genes and Exons, so they can transform

  Returntype: Bio::EnsEMBL::Gene
  Exceptions: none
  Caller    : object::methodname or just methodname

=cut

sub _set_adaptors {
  my $self = shift;
  my $dba = shift;

  if( defined $self->adaptor() ) {
    return;
  }

  # set geneadaptor and exonadaptors,
  # this is to get mappers otherwise where should they come from
  # ? maybe we could get them from the slice

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





#############################
#
# DEPRECATED METHODS FOLLOW
#
#############################

=head2 each_Transcript

 Title   : each_Transcript
 Usage   : DEPRECATED foreach $trans ( $gene->each_Transcript)
 Function: DEPRECATED
 Example : DEPRECATED
 Returns : DEPRECATED An array of Transcript objects
 Args    : DEPRECATED


=cut

sub each_Transcript {
   my ($self) = @_;

   $self->warn("Gene->each_Transcript is deprecated.  " .
	     "Use get_all_Transcripts().\n" . $self->stack_trace_dump() ."\n");
   
   

   return $self->get_all_Transcripts();   
}

=head2 id

 Title   : id
 Usage   : DEPRECATED $obj->id($newval)
 Function: DEPRECATED
 Returns : DEPRECATED value of id
 Args    : DEPRECATED newvalue (optional)


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


1;
