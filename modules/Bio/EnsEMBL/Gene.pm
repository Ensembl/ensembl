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

  my ( $stable_id, $version, $external_name, $type, $external_db, 
       $external_status, $display_xref ) = 
    rearrange( [ 'STABLE_ID', 'VERSION', 'EXTERNAL_NAME', 'TYPE',
		 'EXTERNAL_DB', 'EXTERNAL_STATUS', 'DISPLAY_XREF' ], @_ );
  
  $self->stable_id( $stable_id );
  $self->version( $version );
  $self->external_name( $external_name ) if( defined $external_name );
  $self->external_db( $external_db ) if( defined $external_db );
  $self->external_status( $external_status ) if( defined $external_status );
  $self->display_xref( $display_xref ) if( defined $display_xref );
  $self->type( $type ) if( defined $type );
  return $self;
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


=head2 external_name

  Arg [1]    : string $external_name
  Example    : none
  Description: get/set for attribute external_name.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub external_name {
  my  $self  = shift;

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
  Description: getter setter for gene description
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub description {
    my $self = shift;
    $self->{'description'} = shift if( @_ );
    return $self->{'description'};
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
    throw('Expected DBEntry argument');
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
       throw("$trans is not a Bio::EnsEMBL::Transcript!");
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
   $self->recalculate_coordinates();
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
      my $transcripts = $ta->fetch_all_by_Gene( $self );
      $self->{'_transcript_array'} = $transcripts;
    }
  }
  return $self->{'_transcript_array'};
}




=head2 version

  Arg [1]    : int $version
               A version number for the stable_id
  Example    : nonen
  Description: getter/setter for version number
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub version{
  my $self = shift;
  $self->{'version'} = shift if(@_);
  return $self->{'version'};
}


=head2 stable_id

  Arg [1]    : string $stable_id
  Example    : ("ENSG0000000001")
  Description: getter/setter for stable id for this gene
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub stable_id{
  my $self = shift;
  $self->{'stable_id'} = shift if(@_);
  return $self->{'stable_id'};
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
  if( !@_ || ( ref $_[0] && $_[0]->isa( "Bio::EnsEMBL::Slice" ))) {
    throw( "transform needs coordinate systems details now, please use transfer" );
  }

  my $new_gene = $self->SUPER::transform( @_ );
  return undef unless $new_gene;

  if( exists $self->{'_transcript_array'} ) {
    my @new_transcripts;
    for my $old_transcript ( @{$self->{'_transcript_array'}} ) {
      my $new_transcript = $old_transcript->transform( @_ );
      push( @new_transcripts, $new_transcript );
    }
    $new_gene->{'_transcript_array'} = \@new_transcripts;
  }
  return $new_gene;
}



=head2 transfer

  Arg [1]    : Bio::EnsEMBL::Slice $destination_slice
  Example    : none
  Description: MOves this Gene to given target slice coordinates. If Transcripts
               are attached they are moved as well. Returns a new gene.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general

=cut

sub transfer {
  my $self  = shift;
  
  my $new_gene = $self->SUPER::transfer( @_ );
  return undef unless $new_gene;

  if( exists $self->{'_transcript_array'} ) {
    my @new_transcripts;
    for my $old_transcript ( @{$self->{'_transcript_array'}} ) {
      my $new_transcript = $old_transcript->transfer( @_ );
      push( @new_transcripts, $new_transcript );
    }
    $new_gene->{'_transcript_array'} = \@new_transcripts;
  }
  return $new_gene;
}



=head2 display_xref

  Arg [1]    : Bio::EnsEMBL::DBEntry $display_xref
  Example    : $gene->display_xref($db_entry);
  Description: get/set display_xref for this gene
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
    } 

    return $self->{'display_xref'};
}


=head2 recalculate_coordinates

  Args       : none
  Example    : none
  Description: called when transcript added to the gene
               tries to set coords for the gene.
  Returntype : none
  Exceptions : none
  Caller     : internal

=cut

sub recalculate_coordinates {
  my $self = shift;

  if( ! defined $self->{'_transcript_array'} ) {
    warning( "Cant recalculate position without transcripts" );
    return;
  }
  
  my $transcripts = $self->{'_transcript_array'};

  my ( $slice, $start, $end, $strand );
  $slice = $transcripts->[0]->slice();
  $strand = $transcripts->[0]->strand();
  $start = $transcripts->[0]->start();
  $end = $transcripts->[0]->end();

  my $transsplicing = 0;

  for my $t ( @$transcripts ) {
    if( $t->start() < $start ) {
      $start = $t->start();
    }
  
    if( $t->end() < $end ) {
      $end = $t->end();
    }
  
    if( $t->slice()->name() ne $slice->name() ) {
      throw( "Transcripts with different slices not allowed on one Gene" );
    }
    
    if( $t->strand() != $strand ) {
      $transsplicing = 1;
    }
  }
  if( $transsplicing ) {
    warning( "Gene contained trans splicing event" );
  }

  $self->start( $start );
  $self->end( $end );
  $self->strand( $strand );
  $self->slice( $slice );
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


###########################
# DEPRECATED METHODS FOLLOW
###########################

=head2 DEPRECATED add_DBLink

  Description: DEPRECATED This method has been deprecated in favor of the
               add_DBEntry method.  Objects are responible for holding only
               xrefs directly associated with themselves now.

=cut


sub add_DBLink{
  my ($self,$value) = @_;

  throw("add_DBLink is deprecated.  You probably want add_DBEntry.");

#  unless(defined $value && ref $value 
#	 && $value->isa('Bio::Annotation::DBLink') ) {
#    throw("This [$value] is not a DBLink");
#  }
  
#  if( !defined $self->{'_db_link'} ) {
#    $self->{'_db_link'} = [];
#  }

#  push(@{$self->{'_db_link'}},$value);
}




=head2 temporary_id

 Function: DEPRECATED:  Use dbID or stable_id or something else instead

=cut

sub temporary_id {
   my ($obj,$value) = @_;
   deprecate( "I cant see what a temporary_id is good for, please use " .
               "dbID or stableID or\n try without an id." );
   if( defined $value) {
      $obj->{'temporary_id'} = $value;
    }
    return $obj->{'temporary_id'};
}

=head2 created

  Description: DEPRECATED - Transcripts no longer have a created attribute

=cut

sub created{
    my ($self,$value) = @_;
    deprecated( "The created attribute isnt available any more" );
    if(defined $value ) {
      $self->{'created'} = $value;
    }
    return $self->{'created'};
}

=head2 modified

  Description: DEPRECATED - Transcripts no longer have a modified attribute

=cut

sub modified {
    my ($self,$value) = @_;
    deprecate( "The modified item isnt available any more" );
    if( defined $value ) {
      $self->{'modified'} = $value;
    }
    return $self->{'modified'};
}

=head2 chr_name

  Description: DEPRECATED.  Use project, tranform, or transfer to obtain this
               gene in another coord system.  Use $gene->slice->seq_region_name
               to get the name of the underlying coord system. Or
               $gene->slice->name().

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

1;
