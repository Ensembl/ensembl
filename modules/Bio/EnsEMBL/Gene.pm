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

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::TranscriptI;
use Bio::Annotation::DBLink;
use Bio::EnsEMBL::DBEntry;

@ISA = qw(Bio::EnsEMBL::Root);


=head2 new

  Arg [1]    : none
  Example    : $gene = Bio::EnsEMBL::Gene->new();
  Description: Creates a new gene object
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general

=cut

sub new {
  my($class,@args) = @_;

  my $self = bless {}, $class;
  $self->{'_transcript_array'} = [];

  return $self;
}



=head2 start

  Arg [1]    : (optional) int $start
  Example    : $start = $gene->start;
  Description: This is a convenience method.  It may be better to calculate
               your own start by looking at the exons yourself.
               Gets/sets the lowest start coordinate of this genes exons.
               No consistancy check is performed and if this is used as a
               setter potentially the start could be set to a value which
               does not correspond to the lowest exon start.  If this
               gene is in RawContig coordinates and its exons span multiple
               contigs the lowest value is still returned and a warning is
               issued.
  Returntype : int
  Exceptions : warning if gene spans multiple contigs
  Caller     : general, contigview

=cut

sub start {
  my($self, $start) = @_;

  my $multi_flag = 0;

  if($start) {
    $self->{start} = $start;   
  } elsif(!defined $self->{start}) {
    my $last_contig;
    foreach my $exon (@{$self->get_all_Exons}) {
      if(!defined($self->{start}) || $exon->start() < $self->{start}) {
        $self->{start} = $exon->start();
      }
      $multi_flag = 1 if($last_contig && $last_contig ne $exon->contig->name);
      $last_contig = $exon->contig->name;
    }
  }

  if($multi_flag) {
    $self->warn("Bio::EnsEMBL::Gene::start - Gene spans multiple contigs." .
		"The return value from start may not be what you want");
  }    
  
  return $self->{start};
}



=head2 end

  Arg [1]    : (optional) int $end
  Example    : $end = $gene->end;
  Description: This is a convenience method.  It may be better to calculate
               your own end by looking at the exons yourself.
               Gets/sets the highest end coordinate of this genes exons.
               No consistancy check is performed and if this is used as a
               setter potentially the end could be set to a value which
               does not correspond to the highest exon end.  If this
               gene is in RawContig coordinates and its exons span multiple
               contigs the highest value is still returned and a warning is
               issued.
  Returntype : int
  Exceptions : warning if gene spans multiple contigs
  Caller     : general, contigview

=cut

sub end {
  my($self, $end) = @_;

  my $multi_flag = 0;

  if($end) {
    $self->{end} = $end;   
  } elsif(!defined $self->{end}) {
    my $last_contig;
    foreach my $exon (@{$self->get_all_Exons()}) {
      if(!defined($self->{end}) || $exon->end() > $self->{end}) {
        $self->{end} = $exon->end();
      }
      $multi_flag = 1 if($last_contig && $last_contig ne $exon->contig->name);
      $last_contig = $exon->contig->name;
    }
  }

  if($multi_flag) {
    $self->warn("Bio::EnsEMBL::Gene::end - Gene spans multiple contigs." .
		"The return value from end may not be what you want");
  }

  return $self->{end};
}



=head2 strand

  Arg [1]    : (optional) int strand 
  Example    : $strand = $gene->strand;
  Description: This is a convenience method. It may be better just to
               get the strand from this genes exons yourself.  
               Gets/Sets the strand of this gene. No consistancy check is
               performed and if used as a setter the strand can be set 
               incorrectly.  If this gene is in RawContig coords and spans 
               multiple contigs it is not possible to calculate the strand
               correctly, and a warning is returned.  
  Returntype : int
  Exceptions : Warning if strand is not defined and Gene is in RawContig coords
               so strand cannot be calculated
  Caller     : general

=cut

sub strand {
  my $self = shift;
  my $arg = shift;

  if( defined $arg ) {
    $self->{'strand'} = $arg;
  } elsif( ! defined $self->{strand} ) {
    my $exons = $self->get_all_Exons();
    if(@$exons) {
      if($exons->[0]->contig && 
	 $exons->[0]->contig->isa("Bio::EnsEMBL::RawContig")) {
	my $last_contig;
	foreach my $exon (@$exons) {
	  if($last_contig && $last_contig ne $exon->contig->name) {
	    $self->warn("Bio::EnsEMBL::Gene::strand - strand can not be " .
			"calculated for a Gene in RawContig coordinates that "
			. "spans multiple contigs");
	    return 0;
	  }
	  $last_contig = $exon->contig->name;
	}
      }
      $self->{'strand'} = $exons->[0]->strand();
    }
  }
  return $self->{'strand'};

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
  my ($self, $chr_name) = @_;

  if(defined $chr_name) { 
    $self->{'_chr_name'} = $chr_name;
  } elsif(!defined $self->{'_chr_name'}) {
    #attempt to get the chr_name from the contig attached to the exons
    my ($exon, $contig);
    ($exon) = @{$self->get_all_Exons()};
    if($exon && ($contig = $exon->contig())) {
      if(ref $contig && $contig->isa('Bio::EnsEMBL::Slice')) {
        $self->{'_chr_name'} = $contig->chr_name();
      } else {
	$self->warn('Gene::chr_name - Gene is in RawContig coords, and must '
                  . 'be in Slice coords to have a valid chr_name');
      }
    }
  } 

  return $self->{'_chr_name'};
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
    return 1 if uc($entry->status) eq "KNOWN";
  }

  foreach my $trans ( @{$self->get_all_Transcripts} ) {
    for my $entry ( @{$trans->get_all_DBLinks()} ) {
      return 1 if uc($entry->status) eq "KNOWN";
    }
  }
  
  return 0;
}


=head2 adaptor

  Arg [1]    : Bio::EnsEMBL::DBSQL::GeneAdaptor $adaptor
  Example    : $gene->adaptor($gene_adaptor);
               $gene->adaptor(undef);   # to drop adaptor
  Description: get/set for attribute adaptor
  Returntype : Bio::EnsEMBL::DBSQL::GeneAdaptor
  Exceptions : none
  Caller     : set only used by adaptor on store or retrieve

=cut


sub adaptor {
    my $self = shift;
    
    if (@_) {
        # Testing for any arguments allows undef to be
        # passed as an argument to unset the adaptor
        $self->{'_adaptor'} = shift;
    }
    return $self->{'_adaptor'};
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
  my ($self,$value) = @_;
  if( defined $value ) {
    $self->{'analysis'} = $value;
  }
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
   my ($self, $arg) = @_;

   if ( defined $arg ) {
      $self->{'_dbID'} = $arg ;
   }
   return $self->{'_dbID'};
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

  if(defined $ext_name) { 
    return ( $self->{'_ext_name'} = $ext_name );
  } 

  if( exists $self->{'_ext_name'} ) {
    return $self->{'_ext_name'};
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
  my ( $self, $ext_dbname ) = @_;

  if(defined $ext_dbname) { 
    return ( $self->{'_ext_dbname'} = $ext_dbname );
  } 

  if( exists $self->{'_ext_dbname'} ) {
    return $self->{'_ext_dbname'};
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
    elsif (my $aptr = $self->adaptor) {
        $self->{'_description'} = $aptr->get_description($self->dbID);
    }
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
       foreach my $exon ( @{$trans->get_all_Exons} ) {
	   $h{"$exon"} = $exon;
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
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'type'} = $value;
    }
    return $obj->{'type'};
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

  Args       : none
  Example    : none
  Description: return the Transcripts in this gene
  Returntype : listref Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general

=cut


sub get_all_Transcripts {
  my ($self) = @_;

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

  Arg  1     : (optional) Bio::EnsEMBL::Slice $slice
  Description: when passed a Slice as argument,
               it will transform this Gene to the Slice coordinate system.
               Without an argument it  transforms the Gene (which should be
               in a slice) to a RawContig
               coordinate system.
               The method changes the Gene in place and returns itself.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : object::methodname or just methodname

=cut


sub transform {
  my $self = shift;
  my $slice = shift;

  # hash arrray to store the refs of transformed exons
  my %exon_transforms;

  # transform Exons
  for my $exon ( @{$self->get_all_Exons()} ) {
    my $newExon = $exon->transform( $slice );
    $exon_transforms{ $exon } = $newExon;
  }

  # now need to re-jiggle the transcripts and their
  # translations to account for the re-mapping process

  foreach my $transcript ( @{$self->get_all_Transcripts()} ) {

    # need to grab the translation before starting to 
    # re-jiggle the exons

    $transcript->transform( \%exon_transforms );
  }

  #unset the start, end, and strand - they need to be recalculated
  $self->{start} = undef;
  $self->{end} = undef;
  $self->{strand} = undef;
  $self->{_chr_name} = undef;

  return $self;
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

=head2 get_all_DASFactories

  Arg [1]   : none
  Function  : Retrieves a listref of registered DAS objects
              TODO: Abstract to a DBLinkContainer obj
  Returntype: [ DAS_objects ]
  Exceptions:
  Caller    :
  Example   : $dasref = $prot->get_all_DASFactories

=cut

sub get_all_DASFactories {
   my $self = shift;
   return [ $self->adaptor()->db()->_each_DASFeatureFactory ];
}

=head2 get_all_DASFeatures

  Arg [1]    : none
  Example    : $features = $prot->get_all_DASFeatures;
  Description: Retreives a hash reference to a hash of DAS feature
               sets, keyed by the DNS, NOTE the values of this hash
               are an anonymous array containing:
                (1) a pointer to an array of features;
                (2) a pointer to the DAS stylesheet
              TODO: Abstract to a DBLinkContainer obj
  Returntype : hashref of Bio::SeqFeatures
  Exceptions : ?
  Caller     : webcode

=cut

sub get_all_DASFeatures{
  my ($self,@args) = @_;
  $self->{_das_features} ||= {}; # Cache
  my %das_features;
  foreach my $dasfact( @{$self->get_all_DASFactories} ){
    my $dsn  = $dasfact->adaptor->dsn;
    my $name = $dasfact->adaptor->name;
    $name ||= $dasfact->adaptor->url .'/'. $dsn;
    if( $self->{_das_features}->{$name} ){ # Use cached
      $das_features{$name} = $self->{_das_features}->{$name};
      next;
    }
    else{ # Get fresh data
      my @featref = $dasfact->fetch_all_by_DBLink_Container( $self );
      $self->{_das_features}->{$name} = [@featref];
      $das_features{$name} = [@featref];
    }
  }
  return \%das_features;
}



1;
