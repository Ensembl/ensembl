#
# BioPerl module for Gene
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

A representation of a Gene within the ensembl system. A gene is basically a 
set of one or more alternative transcripts.

=head1 CONTACT

Contact the EnsEMBL development mailing list for info <ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::Gene;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);


@ISA = qw(Bio::EnsEMBL::Feature);


=head2 new

  Arg [start]  : int $start
  Arg [end]    : int $end
  Arg [strand] : 1,-1 $strand
  Arg [slice]  : Bio::EnsEMBL::Slice $slice
  Example    : $gene = Bio::EnsEMBL::Gene->new();
  Description: Creates a new gene object
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

  return $self;
}



=head2 chr_name

  Arg [1]    : (optional) string $chr_name
  Example    : $chr_name = $gene->chr_name
  Description: Getter/Setter for the name of the chromosome that this
               Gene is on.  This is really just a shortcut to the slice
               attached this genes exons, but the value can also be set, which 
               is useful for use with the lite database and web code.
               This function will return undef if this gene is not attached
               to a slice and the chr_name attribute has not already been set. 
  Returntype : string
  Exceptions : warning if chr_name is not defined and Gene is in RawContig 
               coordinates
  Caller     : Lite GeneAdaptor, domainview

=cut

sub chr_name {
  my $self = shift;

  deprecate( "Use project() to obtain other coordinate systems" );

  my $gene_slice = $self->slice();
  if( $gene_slice->coord_system()->name eq "chromosome" ) {
    return $gene_slice->seq_region_name();
  }

  my $sa = $self->slice->adaptor();
  throw( "need db connection for chr_name call" ) unless $sa;

  my $ca = $sa->db()->get_CoordSystemAdaptor();
  my $coord_system = $ca->fetch_by_name( "chromosome" );
  if( ! $coord_system ) {
    throw( "Chromosome coordinate system not available" );
  }
  my $coords = $self->project( $coord_system );

  if( @$coords ) {
    return $coords->[0]->[2]->seq_region_name();
  }
}



=head2 is_known

  Args       : none
  Example    : none
  Description: returns true if the Gene or one of its Transcripts have
               DBLinks
  Returntype : 0,1
  Exceptions : none
  Caller     : general

=cut


sub is_known{
  my ($self) = @_;

  for my $entry ( @{$self->get_all_DBLinks()} ) {
    return 1 if $entry->status =~ /KNOWN/;
  }

  foreach my $trans ( @{$self->get_all_Transcripts} ) {
    for my $entry ( @{$trans->get_all_DBLinks()} ) {
      return 1 if $entry->status =~ /KNOWN/;
    }
  }
  
  return 0;
}


=head2 adaptor

  Arg [1]    : Bio::EnsEMBL::DBSQL::GeneAdaptor $adaptor
  Example    : none
  Description: get/set for attribute adaptor
  Returntype : Bio::EnsEMBL::DBSQL::GeneAdaptor
  Exceptions : none
  Caller     : set only used by adaptor on store or retrieve

=cut


sub adaptor {
   my $self = shift;

   $self->{'adaptor'} = shift if( @_ );

   return $self->{'adaptor'};
}




=head2 analysis

  Arg [1]    : Bio::EnsEMBL::Analysis $analysis
  Example    : none
  Description: get/set for attribute analysis
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : none
  Caller     : general

=cut


sub analysis {
  my $self = shift;

  $self->{'analysis'} = shift if(@_);

  return $self->{'analysis'};
}




=head2 dbID

  Arg [1]    : int $dbID
  Example    : none
  Description: get/set for attribute dbID
  Returntype : int
  Exceptions : none
  Caller     : set only by adaptor on store or retrieve

=cut


sub dbID {
   my $self = shift;

   $self->{'dbID'} = shift if( @_ );
   return $self->{'dbID'};
}



=head2 external_name

  Arg [1]    : string $external_name
  Example    : none
  Description: get/set for attribute external_name.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub external_name {
  my ($self, $ext_name) = @_;

  $self->{'external_name'} = shift if( @_ );

  if( exists $self->{'external_name'} ) {
    return $self->{'external_name'};
  }

  my $display_xref = $self->display_xref();

  if( defined $display_xref ) {
    return $display_xref->display_id()
  } else {
    return undef;
  }
}


=head2 external_db	

  Arg [1]    : string $external_db
  Example    : none
  Description: get/set for attribute external_db. The db is the one that 
               belongs to the external_name.  
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub external_db {
  my $self = shift;

  $self->{'external_db'} = shift if( @_ );

  if( exists $self->{'external_db'} ) {
    return $self->{'external_db'};
  }

  my $display_xref = $self->display_xref();

  if( defined $display_xref ) {
    return $display_xref->dbname()
  } else {
    return undef;
  }
}

=head2 external_status

  Arg [1]    : string $external_status
  Example    : none
  Description: get/set for attribute external_status. The status of
               the external db of the one that belongs to the external_name.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub external_status {
  my ( $self, $ext_status ) = @_;

  return $self->{'_ext_status'} = $ext_status if defined $ext_status;
  return $self->{'_ext_status'} if exists $self->{'_ext_status'};

  my $display_xref = $self->display_xref();

  if( defined $display_xref ) {
    return $display_xref->status()
  } else {
    return undef;
  }
}



=head2 description

  Arg [1]    : (optional) string $description
  Example    : none
  Description: you can set get this argument. If not set its lazy loaded
               from attached adaptor.
  Returntype : string
  Exceptions : if no GeneAdaptor is set and no description is there
  Caller     : general

=cut


sub description {
    my ($self, $arg) = @_;
    
    if( defined $arg ) {
      $self->{'_description'} = $arg;
      return $arg;
    }

    if( exists $self->{'_description'} ) {
      return $self->{'_description'};
    }
    $self->{'_description'} = $self->adaptor->get_description($self->dbID);
    return $self->{'_description'};
}



=head2 add_DBEntry

  Arg [1]    : Bio::EnsEMBL::DBEntry $dbe
               The dbEntry to be added
  Example    : @dbentries = @{$gene->get_all_DBEntries()};
  Description: Associates a DBEntry with this gene. Note that adding DBEntries
               will prevent future lazy-loading of DBEntries for this gene
               (see get_all_DBEntries).
  Returntype : none
  Exceptions : thrown on incorrect argument type
  Caller     : general

=cut

sub add_DBEntry {
  my $self = shift;
  my $dbe = shift;

  unless($dbe && ref($dbe) && $dbe->isa('Bio::EnsEMBL::DBEntry')) {
    $self->throw('Expected DBEntry argument');
  }

  $self->{'dbentries'} ||= [];
  push @{$self->{'dbentries'}}, $dbe;
}


=head2 get_all_DBEntries

  Arg [1]    : none
  Example    : @dbentries = @{$gene->get_all_DBEntries()};
  Description: Retrieves DBEntries (xrefs) for this gene.  This does _not_ 
               include DBEntries that are associated with the transcripts and
               corresponding translations of this gene (see get_all_DBLinks).

               This method will attempt to lazy-load DBEntries from a
               database if an adaptor is available and no DBEntries are present
               on the gene (i.e. they have not already been added or loaded).
  Returntype : list reference to Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : get_all_DBLinks, GeneAdaptor::store

=cut

sub get_all_DBEntries {
  my $self = shift;

  #if not cached, retrieve all of the xrefs for this gene
  if(!defined $self->{'dbentries'} && $self->adaptor()) {
    $self->{'dbentries'} = 
      $self->adaptor->db->get_DBEntryAdaptor->fetch_all_by_Gene($self);
  }

  $self->{'dbentries'} ||= [];

  return $self->{'dbentries'};
}


=head2 get_all_DBLinks

  Arg [1]    : none
  Example    : @dblinks = @{$gene->get_all_DBLinks()};
  Description: Retrieves _all_ related DBEntries for this gene.  This includes
               all DBEntries that are associated with the transcripts and
               corresponding translations of this gene.

               If you only want to retrieve the DBEntries associated with the
               gene (and not the transcript and translations) then you should
               use the get_all_DBEntries call instead.
  Returntype : list reference to Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : general

=cut

sub get_all_DBLinks {
   my $self = shift;

   my @links = @{$self->get_all_DBEntries()};

   # add all of the transcript and translation xrefs to the return list
   foreach my $transc (@{$self->get_all_Transcripts()}) {
     push @links, @{$transc->get_all_DBEntries};

     my $transl = $transc->translation();
     push @links, @{$transl->get_all_DBEntries} if($transl);
   }

   return \@links;
}


=head2 get_all_Exons

  Args       : none
  Example    : none
  Description: a set off all the exons associated with this gene.
  Returntype : listref Bio::EnsEMBL::Exon
  Exceptions : none
  Caller     : general

=cut


sub get_all_Exons {
   my ($self,@args) = @_;
   my %h;

   my @out = ();

   foreach my $trans ( @{$self->get_all_Transcripts} ) {
       foreach my $e ( @{$trans->get_all_Exons} ) {
	   $h{$e->start()."-".$e->end()."-".$e->strand()."-".$e->phase()} = $e;
       }
   }

   push @out, values %h;

   return \@out;
}



=head2 type

  Arg [1]    : string $type
  Example    : none
  Description: get/set for attribute type
  Returntype : string
  Exceptions : none
  Caller     : general

=cut


sub type {
   my $self = shift;

   $self->{'type'} = shift if( @_ );

   return $self->{'type'};
}



=head2 add_Transcript

  Arg  1     : Bio::EnsEMBL::Transcript $transcript
  Example    : none
  Description: adds another Transcript to the set of alternativly
               spliced Transcripts off this gene. If it shares exons 
               with another Transcript, these should be object-identical
  Returntype : none
  Exceptions : none
  Caller     : general

=cut


sub add_Transcript{
   my ($self,$trans) = @_;

   if( !ref $trans || ! $trans->isa("Bio::EnsEMBL::Transcript") ) {
       $self->throw("$trans is not a Bio::EnsEMBL::Transcript!");
   }

   if( defined $self->{'start'} ) {
     if( $self->{'start'} > $trans->start() ) {
       $self->start( $trans->start() );
     }
   } else {
     $self->start( $trans->start() );
   }

   if( defined $self->{'end'} ) {
     if( $self->{'end'} < $trans->end() ) {
       $self->end( $trans->end() );
     }
   } else {
     $self->end( $trans->end() );
   }

   $self->{strand} = $trans->strand();;

   $self->{'_transcript_array'} ||= [];
   push(@{$self->{'_transcript_array'}},$trans);
}



=head2 get_all_Transcripts

  Args       : none
  Example    : none
  Description: return the Transcripts in this gene
  Returntype : listref Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general

=cut


sub get_all_Transcripts {
  my ($self) = @_;

  if( ! exists $self->{'_transcript_array'} ) {
    if( defined $self->adaptor() ) {
      my $ta = $self->adaptor()->db()->get_TranscriptAdaptor();
      my $transcripts = $ta->fetch_by_Gene( $self );
      $self->{'_transcript_array'} = $transcripts;
    }
  }
  return $self->{'_transcript_array'};
}



=head2 created

  Arg [1]    : string $created
               The time the stable id for this gene was created. Not very well
               maintained data (at release 9)
  Example    : none
  Description: get/set/lazy_load for the created timestamp
  Returntype : string
  Exceptions : none
  Caller     : general

=cut


sub created{
    my ($self,$value) = @_;

    deprecated( "The created attribute isnt available any more" );

    if(defined $value ) {
      $self->{'_created'} = $value;
    }


    if( exists $self->{'_created'} ) {
      return $self->{'_created'};
    }

    $self->_get_stable_entry_info();

    return $self->{'_created'};

}


=head2 modified

  Arg [1]    : string $modified
               The time the gene with this stable_id was last modified.
               Not well maintained data (release 9)
  Example    : none
  Description: get/set/lazy_load of modified timestamp
  Returntype : string
  Exceptions : none
  Caller     : general

=cut


sub modified{
    my ($self,$value) = @_;

    deprecate( "The modified item isnt available any more" );

    if( defined $value ) {
      $self->{'_modified'} = $value;
    }

    if( exists $self->{'_modified'} ) {
      return $self->{'_modified'};
    }

    $self->_get_stable_entry_info();

    return $self->{'_modified'};
}



=head2 version

  Arg [1]    : int $version
               A version number for the stable_id
  Example    : nonen
  Description: get/set/lazy_load for the version number
  Returntype : int
  Exceptions : none
  Caller     : general

=cut


sub version{

    my ($self,$value) = @_;
    

    if( defined $value ) {
      $self->{'_version'} = $value;
    }

    if( exists $self->{'_version'} ) {
      return $self->{'_version'};
    }

    $self->_get_stable_entry_info();

    return $self->{'_version'};

}



=head2 stable_id

  Arg [1]    : string $stable_id
  Example    : ("ENSG0000000001")
  Description: get/set/lazy_loaded stable id for this gene
  Returntype : string
  Exceptions : none
  Caller     : general

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



=head2 _get_stable_entry_info

  Args       : none
  Example    : none
  Description: does the lazy loading for all stable id related information
  Returntype : none
  Exceptions : none
  Caller     : internal

=cut


sub _get_stable_entry_info {
   my $self = shift;

   if( !defined $self->adaptor ) {
     return undef;
   }

   $self->adaptor->get_stable_entry_info($self);

}



sub _dump{
   my ($self,$fh) = @_;

   if( ! $fh ) {
       $fh = \*STDOUT;
   }

   print $fh "Gene ", $self->dbID(), "\n";
   foreach my $t ( @{$self->get_all_Transcripts()} ) {
       print $fh "  Trans ", $t->dbID(), " :";
       foreach my $e ( @{$t->get_all_Exons} ) {
	   print $fh " ",$e->dbID(),",";
       }
       print "\n";
   }

}


=head2 transform

  Arg  1     : String $coordinate_system_name
  Arg [2]    : String $coordinate_system_version
  Description: moves this gene to the given coordinate system. If this gene has Transcripts
               attached, they move as well.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : wrong parameters
  Caller     : general

=cut


sub transform {
  my $self = shift;

  # catch for old style transform calls
  if( ref $_[0] && $_[0]->isa( "Bio::EnsEMBL::Slice" )) {
    throw( "transform needs coordinate systems details now, please use transfer" );
  }

  my $new_gene = $self->SUPER::transform( @_ );

  if( exists $self->{'_transcript_array'} ) {
    my @new_transcript_array;
    for my $old_transcript ( @{$self->{'_transcript_array'}} ) {
      my $new_transcript = $old_transcript->transform( @_ );
      push( @{$new_gene->{'_transcript_array'}}, $new_transcript );
    }
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
   deprecated( "I cant see what a temporary_id is good for, please use dbID or stableID or\n"
	       ."try without an id." );

   if( defined $value) {
      $obj->{'temporary_id'} = $value;
    }
    return $obj->{'temporary_id'};

}


=head2 display_xref

  Arg [1]    : Bio::EnsEMBL::DBEntry $display_xref
  Example    : $gene->display_xref($db_entry);
  Description: get/set/lazy_loaded display_xref for this gene
  Returntype : Bio::EnsEMBL::DBEntry
  Exceptions : none
  Caller     : general

=cut

sub display_xref {

    my $self = shift;
    if( @_ ) {
      $self->{'display_xref'} = shift;
    } elsif( exists $self->{'display_xref'} ) {
      return $self->{'display_xref'};
    } elsif ( defined $self->adaptor() ) {
      $self->{'display_xref'} = $self->adaptor->get_display_xref( $self );
    }

    return $self->{'display_xref'};
}



=head2 DEPRECATED add_DBLink

  Arg [1]    : DEPRECATED Bio::Annotation::DBLink $link
               a link is a database entry somewhere else.
               Usually this is a Bio::EnsEMBL::DBEntry.
  Example    : DEPRECATED 
  Description: This method has been deprecated in favor of the add_DBEntry
               method.  Objects are responible for holding only xrefs directly
               associated with themselves now.
  Returntype : none
  Exceptions : none
  Caller     : general

=cut


sub add_DBLink{
  my ($self,$value) = @_;

  $self->throw("add_DBLink is deprecated.  You probably want add_DBEntry.");

#  unless(defined $value && ref $value 
#	 && $value->isa('Bio::Annotation::DBLink') ) {
#    $self->throw("This [$value] is not a DBLink");
#  }
  
#  if( !defined $self->{'_db_link'} ) {
#    $self->{'_db_link'} = [];
#  }

#  push(@{$self->{'_db_link'}},$value);
}

1;
