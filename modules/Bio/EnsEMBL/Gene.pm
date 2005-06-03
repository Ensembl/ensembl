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
       $external_status, $display_xref, $description, $transcripts,
       $created_date, $modified_date, $confidence, $biotype, $source ) = 
    rearrange( [ 'STABLE_ID', 'VERSION', 'EXTERNAL_NAME', 'TYPE',
		 'EXTERNAL_DB', 'EXTERNAL_STATUS', 'DISPLAY_XREF', 'DESCRIPTION',
                 'TRANSCRIPTS', 'CREATED_DATE', 'MODIFIED_DATE', 
	         'CONFIDENCE', 'BIOTYPE', 'SOURCE'], @_ );

  if ($transcripts) {
    $self->{'_transcript_array'} = $transcripts;
    $self->recalculate_coordinates();
  }

  $self->stable_id( $stable_id );
  $self->version( $version );
  $self->{'created_date'} = $created_date;
  $self->{'modified_date'} = $modified_date;

  $self->external_name( $external_name ) if( defined $external_name );
  $self->external_db( $external_db ) if( defined $external_db );
  $self->external_status( $external_status ) if( defined $external_status );
  $self->display_xref( $display_xref ) if( defined $display_xref );
  $self->biotype( $type ) if( defined $type );
  $self->biotype( $biotype ) if( defined $biotype );
  $self->description($description);
  $self->confidence( $confidence );
  $self->source( $source );
  return $self;
}



=head2 is_known

  Args       : none
  Example    : none
  Description: returns true if this gene has a display_xref
  Returntype : 0,1
  Exceptions : none
  Caller     : general

=cut


sub is_known{
  my $self = shift;
  return ( $self->{'confidence'} eq "KNOWN" );
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


=head2 confidence

  Arg [1]    : string $confidence
  Example    : none
  Description: get/set for attribute confidence
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub confidence {
   my $self = shift;
  $self->{'confidence'} = shift if( @_ );
  return $self->{'confidence'};
}



=head2 source

  Arg [1]    : string $source
  Example    : none
  Description: get/set for attribute source
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub source {
  my $self = shift;
  $self->{'source'} = shift if( @_ );
  return ( $self->{'source'} || "ensembl" );
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
  my $self = shift;

  $self->{'_ext_status'} = shift if ( @_ );
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

=head2 get_all_homologous_Genes

  Args       : none
  Example    : 
  Description: Queries the Ensembl Compara database and retrieves all
               Genes from other species that are orthologous.
               REQUIRES properly setup Registry conf file.
  Returntype : listref [
                        Bio::EnsEMBL::Gene,
                        Bio::EnsEMBL::Compara::Homology,
                        string $species, # Need as cannot get spp from Gene 
                       ]
  Exceptions : none
  Caller     : general

=cut

sub get_all_homologous_Genes{
  my $self = shift;

  if( exists( $self->{'homologues'} ) ){
    return $self->{'homologues'};
  }
  $self->{'homologues'} = [];

  # TODO: Find a robust way of retrieving compara dba directly.
  # For now look through all DBAs
  my $compara_dba;
  foreach my $dba( Bio::EnsEMBL::Registry->get_all_DBAdaptors ){
    if( $dba->isa('Bio::EnsEMBL::Compara::DBSQL::DBAdaptor') ){
      $compara_dba = $dba;
      last;
    }
  }
  unless( $compara_dba ){
    warning("No compara in Bio::EnsEMBL::Registry");
    return $self->{'homologues'};
  }

  # Get the compara 'member' corresponding to self
  my $member_adaptor   = $compara_dba->get_adaptor('Member');
  my $query_member = $member_adaptor->fetch_by_source_stable_id
      ("ENSEMBLGENE",$self->stable_id);
  unless( $query_member ){ return $self->{'homologues'} };

  # Get the compara 'homologies' corresponding to 'member'
  my $homology_adaptor = $compara_dba->get_adaptor('Homology');
  my @homolos = @{$homology_adaptor->fetch_by_Member($query_member)};
  unless( scalar(@homolos) ){ return $self->{'homologues'} };

  # Get the ensembl 'genes' corresponding to 'homologies'
  foreach my $homolo( @homolos ){
    foreach my $member_attrib( @{$homolo->get_all_Member_Attribute} ){
      my ($member, $attrib) = @{$member_attrib};
      my $hstable_id = $member->stable_id;
      next if ($hstable_id eq $query_member->stable_id); # Ignore self     
      my $hgene = undef;
      eval { $hgene = $member->get_Gene;} ;
      unless( $hgene ){
        # Something up with DB. Create a new gene is best we can do
        $hgene = Bio::EnsEMBL::Gene->new
            ( -stable_id=>$hstable_id,
              -description=>$member->description, );
      }
      my $hspecies = $member->genome_db->name;
      push @{$self->{'homologues'}}, [$hgene,$homolo,$hspecies];
    }
  }
  return $self->{'homologues'};
}

=head2 type

  Arg [1]    : string $type
  Example    : none
  Description: get/set for attribute type
               This function is going to be deprecated soon, use biotype instead
  Returntype : string
  Exceptions : none
  Caller     : general

=cut


sub type {
   my $self = shift;
   
   $self->{'biotype'} = shift if( @_ );

   return ( $self->{'biotype'} || "protein_coding" );
}



=head2 biotype

  Arg [1]    : string $biotype
  Example    : $gene->biotype( "protein_coding" );
  Description: get/set for the biotype attribute
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub biotype {
  my $self = shift;

  $self->{'biotype'} = shift if( @_ );
  return ( $self->{'biotype'} || "protein_coding" );
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


sub add_Transcript {
   my ($self,$trans) = @_;

   if( !ref $trans || ! $trans->isa("Bio::EnsEMBL::Transcript") ) {
       throw("$trans is not a Bio::EnsEMBL::Transcript!");
   }

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



=head2 get_all_alt_alleles

  Arg [1]    : none
  Example    :  ( optional )
  Description: Return a listref of Gene objects that represent this Gene on
               an alternative haplotype. Empty list if there is no such
               Gene. (eg there is no overlapping haplotype)
  Returntype : listref of Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general

=cut

sub get_all_alt_alleles {
  my $self = shift;
  my $result = $self->adaptor()->fetch_all_alt_alleles( $self );
  return $result;
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


sub created_date {
  my $self = shift;
  $self->{'created_date'} = shift if ( @_ );
  return $self->{'created_date'};
}


sub modified_date {
  my $self = shift;
  $self->{'modified_date'} = shift if ( @_ );
  return $self->{'modified_date'};
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
  if( !@_  || ( ref $_[0] && $_[0]->isa( "Bio::EnsEMBL::Slice" ))) {
    deprecate('Calling transform without a coord system name is deprecated.');
    return $self->_deprecated_transform(@_);
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
  $self->{'display_xref'} = shift if(@_);
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

  my $transcripts = $self->get_all_Transcripts();

  return if(!$transcripts || !@$transcripts);

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

    if( $t->end() > $end ) {
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



=head2 display_id

  Arg [1]    : none
  Example    : print $gene->display_id();
  Description: This method returns a string that is considered to be
               the 'display' identifier.  For genes this is the stable id if
               it is available otherwise it is an empty string.
  Returntype : string
  Exceptions : none
  Caller     : web drawing code

=cut

sub display_id {
  my $self = shift;
  return $self->{'stable_id'} || '';
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

sub get_all_DASFeatures_by_slice{
  my ($self, $slice, @args) = @_;
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
      my @featref = ("${name}_ENSMUSG", ($dasfact->fetch_all_by_Slice( $slice ))[0]);
      $self->{_das_features}->{$name} = [@featref];
      $das_features{$name} = [@featref];
    }
  }
  return \%das_features;
}


=head2 get_all_DAS_Features

  Arg [1]    : none
  Example    : $features = $prot->get_all_DAS_Features;
  Description: Retreives a hash reference to a hash of DAS feature
               sets, keyed by the DNS, NOTE the values of this hash
               are an anonymous array containing:
                (1) a pointer to an array of features;
                (2) a pointer to the DAS stylesheet
  Returntype : hashref of Bio::SeqFeatures
  Exceptions : ?
  Caller     : webcode


=cut

sub get_all_DAS_Features{
  my ($self,@args) = @_;
  $self->{_das_features} ||= {}; # Cache
  my %das_features;

  my $slice = $self->feature_Slice;

  foreach my $dasfact( @{$self->get_all_DASFactories} ){
    my $dsn = $dasfact->adaptor->dsn;
    my $name = $dasfact->adaptor->name;
    my $type = $dasfact->adaptor->type;
    my $key = defined($dasfact->adaptor->url) ? $dasfact->adaptor->url .'/'. $dsn : $dasfact->adaptor->protocol .'://'.$dasfact->adaptor->domain.'/'. $dsn;

    $name ||= $key;

    if( $self->{_das_features}->{$key} ){ # Use cached
	$das_features{$key} = $self->{_das_features}->{$key};
	next;
    } else{ # Get fresh data
	my @featref = ($type eq 'ensembl_location') ?  ($name, ($dasfact->fetch_all_by_Slice( $slice ))[0]) : $dasfact->fetch_all_by_DBLink_Container( $self );
	$self->{_das_features}->{$key} = [@featref];
	$das_features{$key} = [@featref];
    }
  }
  return \%das_features;
}


=head2 get_all_regulatory_features

  Arg [1]    : If set, regulatory features on transcripts belonging to this gene
               are returned as well.
  Example    : @features = $gene->gene_all_regulatory_features(1);
  Description: Gets all the regulatory features associated with a
               particular gene, and (optionally) its transcripts.
               Each feature only appears once.
  Returntype : Listref of Bio::EnsEMBL::RegulatoryFeature
  Exceptions : If arg is not of correct type.
  Caller     : ?

=cut

sub get_all_regulatory_features {

   my ($self, $recursive) = @_;

   my $rfa = $self->adaptor->db->get_RegulatoryFeatureAdaptor();

   return $rfa->fetch_all_by_gene($self, $recursive);

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

  my $coords = $self->project( "toplevel" );

  if( @$coords ) {
    return $coords->[0]->[2]->seq_region_name();
  }
}


1;
