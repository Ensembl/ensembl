#
# Ensembl module for Transcript
#
# Copyright (c) 1999-2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#

=head1 NAME

Transcript - gene transcript object

=head1 SYNOPSIS

Creation:

     my $tran = new Bio::EnsEMBL::Transcript();
     my $tran = new Bio::EnsEMBL::Transcript(-EXONS => \@exons);

Manipulation:

     # Returns an array of Exon objects
     my @exons = @{$tran->get_all_Exons()};

     # Returns the peptide translation of the exons as a Bio::Seq
     if($tran->translation() {
       my $pep   = $tran->translate();
     } else {
       print "Transcript ", $tran->stable_id(), " is non-coding\n";
     }

=head1 DESCRIPTION

A representation of a transcript within the ensembl system.  A transcript
consists of a set of Exons and (possibly) a Translation which defines the
coding and non-coding regions of the exons.

=head1 CONTACT

Email questions to the ensembl developer mailing list <ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut

package Bio::EnsEMBL::Transcript;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Intron;
use Bio::EnsEMBL::TranscriptMapper;
use Bio::EnsEMBL::Utils::TranscriptSNPs;
use Bio::EnsEMBL::SeqEdit;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( deprecate warning throw );


@ISA = qw(Bio::EnsEMBL::Feature);



=head2 new

  Arg [-EXONS] :
        reference to list of Bio::EnsEMBL::Exon objects - exons which make up 
        this transcript
  Arg [-STABLE_ID] :
        string - the stable identifier of this transcript
  Arg [-VERSION] :
        int - the version of the stable identifier of this transcript
  Arg [-EXTERNAL_NAME] :
        string - the external database name associated with this transcript
  Arg [-EXTERNAL_DB] :
        string - the name of the database the external name is from
  Arg [-EXTERNAL_STATUS]:
        string - the status of the external identifier
  Arg [-DISPLAY_XREF]:
        Bio::EnsEMBL::DBEntry - The external database entry that is used
        to label this transcript when it is displayed.
  Example    : $tran = new Bio::EnsEMBL::Transcript(-EXONS => \@exons);
  Description: Constructor. Instantiates a Transcript object.
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : throw on bad arguments
  Caller     : general

=cut



sub new {
  my($class) = shift;

  if( ref $class ) { 
      $class = ref $class;
  }

  my $self = $class->SUPER::new(@_);

  my ( $exons, $stable_id, $version, $external_name, $external_db,
       $external_status, $display_xref, $created_date, $modified_date );

  #catch for old style constructor calling:
  if((@_ > 0) && ref($_[0])) {
    $exons = [@_];
    deprecate("Transcript constructor should use named arguments.\n" .
              'Use Bio::EnsEMBL::Transcript->new(-EXONS => \@exons);' .
              "\ninstead of Bio::EnsEMBL::Transcript->new(\@exons);");
  }
  else {
    ( $exons, $stable_id, $version, $external_name, $external_db,
      $external_status, $display_xref, $created_date, $modified_date ) = 
        rearrange( [ "EXONS", 'STABLE_ID', 'VERSION', 'EXTERNAL_NAME', 
                     'EXTERNAL_DB', 'EXTERNAL_STATUS', 'DISPLAY_XREF',
		     'CREATED_DATE', 'MODIFIED_DATE' ], @_ );
  }

  if( $exons ) {
    $self->{'_trans_exon_array'} = $exons;
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
  $self->edits_enabled(1);


  return $self;
}


=head2 get_all_DBLinks

  Arg [1]    : none
  Example    : @dblinks = @{$transcript->get_all_DBLinks()};
  Description: Retrieves _all_ related DBEntries for this transcript.  
               This includes all DBEntries that are associated with the
               corresponding translation.

               If you only want to retrieve the DBEntries associated with the
               transcript then you should use the get_all_DBEntries call 
               instead.
  Returntype : list reference to Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : general

=cut

sub get_all_DBLinks {
  my $self = shift;

  my @links;

  push @links, @{$self->get_all_DBEntries};

  my $transl = $self->translation();
  push @links, @{$transl->get_all_DBEntries} if($transl);

  return \@links;
}


=head2 get_all_DBEntries

  Arg [1]    : none
  Example    : @dbentries = @{$gene->get_all_DBEntries()};
  Description: Retrieves DBEntries (xrefs) for this transcript.  
               This does _not_ include the corresponding translations 
               DBEntries (see get_all_DBLinks).

               This method will attempt to lazy-load DBEntries from a
               database if an adaptor is available and no DBEntries are present
               on the transcript (i.e. they have not already been added or 
               loaded).
  Returntype : list reference to Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : get_all_DBLinks, TranscriptAdaptor::store

=cut

sub get_all_DBEntries {
  my $self = shift;

  #if not cached, retrieve all of the xrefs for this gene
  if(!defined $self->{'dbentries'} && $self->adaptor()) {
    $self->{'dbentries'} = 
      $self->adaptor->db->get_DBEntryAdaptor->fetch_all_by_Transcript($self);
  }

  $self->{'dbentries'} ||= [];

  return $self->{'dbentries'};
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



=head2 get_all_supporting_features

  Arg [1]    : none
  Example    : @evidence = @{$transcript->get_all_supporting_features()};
  Description: Retreives any supporting features added manually by 
               calls to add_supporting_features.
  Returntype : listreference of Bio::EnsEMBL::FeaturePair
  Exceptions : none
  Caller     : general

=cut

sub get_all_supporting_features {
  my $self = shift;
  
  return $self->{_supporting_evidence} || [];
}



=head2 add_supporting_features

  Arg [1]    : Bio::EnsEMBL::FeaturePair $feature
  Example    : $exon->add_supporting_features(@features);
  Description: Adds a list of supporting features to this Transcript.
               The added features can be retieved by get_all_supporting_features
  Returntype : none
  Exceptions : throw if any of the features are not FeaturePairs
               throw if any of the features are not in the same coordinate
               system as the Transcript
  Caller     : general
 
=cut
 
sub add_supporting_features {
  my ($self,@features) = @_;

  return unless @features;
 
  $self->{_supporting_evidence} ||= [];
  
  # check whether this feature object has been added already
  FEATURE: foreach my $feature (@features) {

    unless($feature && $feature->isa("Bio::EnsEMBL::FeaturePair")) {
      throw("Supporting feat [$feature] not a " .
            "Bio::EnsEMBL::FeaturePair");
    } 
    
    if ((defined $self->slice() && defined $feature->slice())&&
	    ( $self->slice()->name() ne $feature->slice()->name())){
      throw("Supporting feat not in same coord system as exon\n" .
            "exon is attached to [".$self->slice()->name()."]\n" .
            "feat is attached to [".$feature->slice()->name()."]");
    }

    foreach my $added_feature ( @{ $self->{_supporting_evidence} } ){
      # compare objects
      if ( $feature == $added_feature ){
	#this feature has already been added
	next FEATURE;
      }
    }
    
    #no duplicate was found, add the feature
    push(@{$self->{_supporting_evidence}},$feature);
  }
}


=head2 external_db

 Title   : external_db
 Usage   : $ext_db = $obj->external_db();
 Function: external_name if available
 Returns : the external db link for this transcript
 Args    : new external db (optional)

=cut

sub external_db {
  my ( $self, $ext_dbname ) = @_;

  if(defined $ext_dbname) { 
    return ( $self->{'external_db'} = $ext_dbname );
  }

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

 Title   : external_status
 Usage   : $ext_db = $obj->external_status();
 Function: external_name if available
 Returns : the external db link for this transcript
 Args    : new external db (optional)

=cut

sub external_status { 
  my ( $self, $ext_status ) = @_;

  if(defined $ext_status) {
    return ( $self->{'external_status'} = $ext_status );
  }

  if( exists $self->{'external_status'} ) {
    return $self->{'external_status'};
  }

  my $display_xref = $self->display_xref();

  if( defined $display_xref ) {
    return $display_xref->status()
  } else {
    return undef;
  }
}



=head2 external_name

 Title   : external_name
 Usage   : $ext_name = $obj->external_name();
 Function: external_name if available
 Example : 
 Returns : the external name of this transcript
 Args    : new external name (optional)

=cut


sub external_name {
  my ($self, $ext_name) = @_;

  if(defined $ext_name) { 
    return ( $self->{'external_name'} = $ext_name );
  }

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


=head2 is_known

  Args       : none
  Example    : none
  Description: returns true if this transcript has a display_xref
  Returntype : 0,1
  Exceptions : none
  Caller     : general

=cut

sub is_known {
  my $self = shift;
  return ($self->{'display_xref'}) ? 1 : 0;
}


sub type {
  my $self = shift;

  $self->{'type'} = shift if( @_ );
  return $self->{'type'};
}


=head2 display_xref

  Arg [1]    : Bio::EnsEMBL::DBEntry $display_xref
  Example    : $transcript->display_xref(Bio::EnsEMBL::DBEntry->new(...));
  Description: getter setter for display_xref for this transcript
  Returntype : Bio::EnsEMBL::DBEntry
  Exceptions : none
  Caller     : general

=cut

sub display_xref {
  my $self = shift;
  $self->{'display_xref'} = shift if(@_);
  return $self->{'display_xref'};
}


=head2 translation

  Arg [1]    : Bio::EnsEMBL::Translation
  Example    : if($transcript->translation()) {
                 print $translation->stable_id(), "\n";
               } else {
                 print "Pseudogene\n";
               }
  Description: Getter/setter for the Translation object which defines the
               CDS (and as a result the peptide encoded by) this transcript.
               This function will return undef if this Transcript is a
               pseudogene - i.e. a non-translating transcript such as an
               ncRNA.  This is the accepted method of determining whether
               a transcript is a pseudogene or not.
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general

=cut

sub translation {
  my $self = shift;
  if( @_ ) {
    my $value = shift;
    if( defined($value) &&
        (!ref($value) || !$value->isa('Bio::EnsEMBL::Translation'))) {
      throw("Bio::EnsEMBL::Translation argument expected.");
    }
    $self->{'translation'} = $value;
  } elsif( !exists($self->{'translation'}) and defined($self->adaptor())) {
    $self->{'translation'} =
      $self->adaptor()->db()->get_TranslationAdaptor()->
        fetch_by_Transcript( $self );
  }
  return $self->{'translation'};
}



=head2 spliced_seq

  Args       : none
  Example    : none
  Description: Retrieves all Exon sequences and concats them together.
               No phase padding magic is done, even if phases do not align.
  Returntype : txt
  Exceptions : none
  Caller     : general

=cut

sub spliced_seq {
  my ( $self ) = @_;

  my $seq_string = "";
  for my $ex ( @{$self->get_all_Exons()} ) {
    my $seq = $ex->seq();

    if(!$seq) {
      warning("Could not obtain seq for exon.  Transcript sequence may not " .
              "be correct.");
      $seq_string .= 'N' x $ex->length();
    } else {
      $seq_string .= $seq->seq();
    }
  }

  # apply post transcriptional edits
  if($self->edits_enabled()) {
    my @seqeds = @{$self->get_all_SeqEdits()};

    # sort edits in reverse order to remove complication of
    # adjusting downstream edits
    @seqeds = sort {$b->start() <=> $a->start()} @seqeds;

    foreach my $se (@seqeds) {
      $se->apply_edit(\$seq_string);
    }
  }

  return $seq_string;
}


=head2 translateable_seq

  Args       : none
  Example    : print $transcript->translateable_seq(), "\n";
  Description: Returns a sequence string which is the the translateable part
               of the transcripts sequence.  This is formed by splicing all
               Exon sequences together and apply all defined RNA edits.
               Then the coding part of the sequence is extracted and returned.
               The code will not support monkey exons any more. If you want to
               have non phase matching exons, defined appropriate _rna_edit
               attributes!

               An empty string is returned if this transcript is a pseudogene
               (i.e. is non-translateable).
  Returntype : txt
  Exceptions : none
  Caller     : general

=cut

sub translateable_seq {
  my ( $self ) = @_;

  if(!$self->translation()) {
    return '';
  }

  my $mrna = $self->spliced_seq();
  my $start = $self->cdna_coding_start();
  my $end = $self->cdna_coding_end();

  $mrna = substr( $mrna, $start-1, $end-$start+1 );

  my $start_phase = $self->translation->start_Exon->phase();
  if( $start_phase > 0 ) {
    $mrna = "N"x$start_phase . $mrna;
  }
  if( ! $start || ! $end ) {
    return "";
  }

  return $mrna;
}



=head2 cdna_coding_start

  Arg [1]    : (optional) $value
  Example    : $relative_coding_start = $transcript->cdna_coding_start;
  Description: Retrieves the position of the coding start of this transcript
               in cdna coordinates (relative to the start of the 5prime end of
               the transcript, excluding introns, including utrs).

               This will return undef if this is a pseudogene (i.e. a
               transcript with no translation).
  Returntype : int
  Exceptions : none
  Caller     : five_prime_utr, get_all_snps, general

=cut

sub cdna_coding_start {
  my $self = shift;

  if( @_ ) {
    $self->{'cdna_coding_start'} = shift;
  }

  if(!defined $self->{'cdna_coding_start'} && defined $self->translation){
    # calc coding start relative from the start of translation (in cdna coords)
    my $start = 0;

    my @exons = @{$self->get_all_Exons};
    my $exon;

    while($exon = shift @exons) {
      if($exon == $self->translation->start_Exon) {
        #add the utr portion of the start exon
        $start += $self->translation->start;
        last;
      } else {
        #add the entire length of this non-coding exon
        $start += $exon->length;
      }
    }

    # adjust cdna coords if sequence edits are enabled
    if($self->edits_enabled()) {
      my @seqeds = @{$self->get_all_SeqEdits()};
      # sort in reverse order to avoid adjustment of downstream edits
      @seqeds = sort {$b->start() <=> $a->start()} @seqeds;

      foreach my $se (@seqeds) {
        # use less than start so that start of CDS can be extended
        if($se->start() < $start) {
          $start += $se->length_diff();
        }
      }
    }

    $self->{'cdna_coding_start'} = $start;
  }

  return $self->{'cdna_coding_start'};
}



=head2 cdna_coding_end

  Arg [1]    : (optional) $value
  Example    : $cdna_coding_end = $transcript->cdna_coding_end;
  Description: Retrieves the end of the coding region of this transcript in
               cdna coordinates (relative to the five prime end of the
               transcript, excluding introns, including utrs).

               This will return undef if this transcript is a pseudogene
               (i.e. a transcript with no translation and therefor no CDS).
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub cdna_coding_end {
  my $self = shift;

  if( @_ ) {
    $self->{'cdna_coding_end'} = shift;
  }

  if(!defined $self->{'cdna_coding_end'} && defined $self->translation) {
    my @exons = @{$self->get_all_Exons};

    my $end = 0;
    while(my $exon = shift @exons) {
      if($exon == $self->translation->end_Exon) {
        # add coding portion of the final coding exon
        $end += $self->translation->end;
        last;
      } else {
        # add entire exon
        $end += $exon->length;
      }
    }

    # adjust cdna coords if sequence edits are enabled
    if($self->edits_enabled()) {
      my @seqeds = @{$self->get_all_SeqEdits()};
      # sort in reverse order to avoid adjustment of downstream edits
      @seqeds = sort {$b->start() <=> $a->start()} @seqeds;

      foreach my $se (@seqeds) {
        # use less than or equal to end+1 so end of the CDS can be extended
        if($se->start() <= $end + 1) {
          $end += $se->length_diff();
        }
      }
    }

    $self->{'cdna_coding_end'} = $end;
  }

  return $self->{'cdna_coding_end'};
}



=head2 coding_region_start

  Arg [1]    : (optional) $value
  Example    : $coding_region_start = $transcript->coding_region_start
  Description: Retrieves the start of the coding region of this transcript
               in genomic coordinates (i.e. in either slice or contig coords).
               By convention, the coding_region_start is always lower than
               the value returned by the coding_end method.
               The value returned by this function is NOT the biological
               coding start since on the reverse strand the biological coding
               start would be the higher genomic value.

               This function will return undef if this is a pseudogene
               (a non-translated transcript).
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub coding_region_start {
  my ($self, $value) = @_;

  if( defined $value ) {
    $self->{'coding_region_start'} = $value;
  } elsif(!defined $self->{'coding_region_start'} &&
	  defined $self->translation) {
    #calculate the coding start from the translation
    my $start;
    my $strand = $self->translation()->start_Exon->strand();
    if( $strand == 1 ) {
      $start = $self->translation()->start_Exon->start();
      $start += ( $self->translation()->start() - 1 );
    } else {
      $start = $self->translation()->end_Exon->end();
      $start -= ( $self->translation()->end() - 1 );
    }
    $self->{'coding_region_start'} = $start;
  }

  return $self->{'coding_region_start'};
}



=head2 coding_region_end

  Arg [1]    : (optional) $value
  Example    : $coding_region_end = $transcript->coding_region_end
  Description: Retrieves the end of the coding region of this transcript
               in genomic coordinates (i.e. in either slice or contig coords).
               By convention, the coding_region_end is always higher than the
               value returned by the coding_region_start method.
               The value returned by this function is NOT the biological
               coding end since on the reverse strand the biological coding
               end would be the lower genomic value.

               This function will return undef if this is a pseudogene
               (a non-translated transcript).
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub coding_region_end {
  my ($self, $value ) = @_;

  my $strand;
  my $end;

  if( defined $value ) {
    $self->{'coding_region_end'} = $value;
  } elsif( ! defined $self->{'coding_region_end'}
	   && defined $self->translation() ) {
    $strand = $self->translation()->start_Exon->strand();
    if( $strand == 1 ) {
      $end = $self->translation()->end_Exon->start();
      $end += ( $self->translation()->end() - 1 );
    } else {
      $end = $self->translation()->start_Exon->end();
      $end -= ( $self->translation()->start() - 1 );
    }
    $self->{'coding_region_end'} = $end;
  }

  return $self->{'coding_region_end'};
}



=head2 edits_enabled

  Arg [1]    : (optional) boolean $newval
  Example    : $transcript->edits_enabled(1);
  Description: Enables/Disables the application of SeqEdits to this transcript.
               Edits are enabled by default, and affect the cdna/mrna
               sequences coordinates and the resultant translation.
  Returntype : boolean - the current value of the edits
  Exceptions : none
  Caller     : general, cdna_coding_start, cdna_coding_end, length

=cut

sub edits_enabled {
  my $self = shift;

  if(@_) {
    $self->{'edits_enabled'} = shift;
    # flush cached values that will be different with/without edits
    $self->{'cdna_coding_start'} = undef;
    $self->{'cdna_coding_end'}   = undef;
    $self->{'transcript_mapper'} = undef;
  }

  return $self->{'edits_enabled'};
}


=head2 get_all_SeqEdits

  Arg [1]    : none
  Example    : my @seqeds = @{$transcript->get_all_SeqEdits()};
  Description: Retrieves all post transcriptional sequence modifications for
               this transcript.
  Returntype : Bio::EnsEMBL::SeqEdit
  Exceptions : none
  Caller     : spliced_seq()

=cut

sub get_all_SeqEdits {
  my $self = shift;

  my @seqeds;

  my $attribs = $self->get_all_Attributes('_rna_edit');

  # convert attributes to SeqEdit objects
  foreach my $a (@$attribs) {
    push @seqeds, Bio::EnsEMBL::SeqEdit->new(-ATTRIB => $a);
  }

  return \@seqeds;
}


=head2 get_all_Attributes

  Arg [1]    : optional string $attrib_code
               The code of the attribute type to retrieve values for.
  Example    : ($rna_edits) = @{$transcript->get_all_Attributes('_rna_edit')};
               @transc_attributes    = @{$transcript->get_all_Attributes()};
  Description: Gets a list of Attributes of this transcript.
               Optionally just get Attrubutes for given code.
  Returntype : listref Bio::EnsEMBL::Attribute
  Exceptions : warning if transcript does not have attached adaptor and 
               attempts lazy load.
  Caller     : general

=cut

sub get_all_Attributes {
  my $self = shift;
  my $attrib_code = shift;

  if( ! exists $self->{'attributes' } ) {
    if(!$self->adaptor() ) {
#      warning('Cannot get attributes without an adaptor.');
      return [];
    }

    my $attribute_adaptor = $self->adaptor->db->get_AttributeAdaptor();
    $self->{'attributes'} = $attribute_adaptor->fetch_all_by_Transcript($self);
  }

  if( defined $attrib_code ) {
    my @results = grep { uc($_->code()) eq uc($attrib_code) }
    @{$self->{'attributes'}};
    return \@results;
  } else {
    return $self->{'attributes'};
  }
}


=head2 add_Attributes

  Arg [1...] : Bio::EnsEMBL::Attribute $attribute
               You can have more Attributes as arguments, all will be added.
  Example    : $transcript->add_Attributes($rna_edit_attribute);
  Description: Adds an Attribute to the Transcript. Usefull to do _rna_edits.
               If you add an attribute before you retrieve any from database, 
               lazy load will be disabled.
  Returntype : none
  Exceptions : throw on incorrect arguments
  Caller     : general

=cut

sub add_Attributes {
  my $self = shift;
  my @attribs = @_;

  if( ! exists $self->{'attributes'} ) {
    $self->{'attributes'} = [];
  }

  for my $attrib ( @attribs ) {
    if( ! $attrib->isa( "Bio::EnsEMBL::Attribute" )) {
     throw( "Argument to add_Attribute has to be an Bio::EnsEMBL::Attribute" );
    }
    push( @{$self->{'attributes'}}, $attrib );
  }

  # flush cdna coord cache b/c we may have added a SeqEdit
  $self->{'cdna_coding_start'} = undef;
  $self->{'cdna_coding_end'} = undef;
  $self->{'transcript_mapper'} = undef;

  return;
}


=head2 add_Exon

 Title   : add_Exon
 Usage   : $trans->add_Exon($exon)
 Returns : Nothing
 Args    :

=cut

sub add_Exon{
  my ($self,$exon, $rank) = @_;

  #yup - we are going to be picky here...
  unless(defined $exon && ref $exon && $exon->isa("Bio::EnsEMBL::Exon") ) {
    throw("[$exon] is not a Bio::EnsEMBL::Exon!");
  }

  $self->{'_trans_exon_array'} ||= [];

  if(defined($rank)) {
    $self->{'_trans_exon_array'}->[$rank-1] = $exon;
    return;
  }

  my $was_added = 0;

  my $ea = $self->{'_trans_exon_array'};
  if( @$ea ) {
    if( $exon->strand() == 1 ) {
      if( $exon->start() > $ea->[$#$ea]->end() ) {
        push(@{$self->{'_trans_exon_array'}},$exon);
        $was_added = 1;
      } else {
        # insert it at correct place
        for( my $i=0; $i <= $#$ea; $i++ ) {
          if( $exon->end() < $ea->[$i]->start() ) {
            splice( @$ea, $i, 0, $exon );
            $was_added = 1;
            last;
          }
        }
      }
    } else {
      if( $exon->end() < $ea->[$#$ea]->start() ) {
        push(@{$self->{'_trans_exon_array'}},$exon);
        $was_added = 1;
      } else {
        # insert it at correct place
        for( my $i=0; $i <= $#$ea; $i++ ) {
          if( $exon->start() > $ea->[$i]->end() ) {
            splice( @$ea, $i, 0, $exon );
            $was_added = 1;
            last;
          }
        }
      }
    }
  } else {
    push( @$ea, $exon );
    $was_added = 1;
  }

  # sanity check:
  if(!$was_added) {
    # exon was not added because it has same end coord as start
    # of another exon
    my $all_str = '';
    foreach my $e (@$ea) {
      $all_str .= '  '.$e->start .'-'.$e->end.' ('.$e->strand.') ' .
        ($e->stable_id || '') . "\n";
    }
    my $cur_str = '  '.$exon->start.'-'.$exon->end. ' ('.$exon->strand.') '.
      ($exon->stable_id || '')."\n";
    throw("Exon overlaps with other exon in same transcript.\n" .
          "Transcript Exons:\n$all_str\n" .
          "This Exon:\n$cur_str");
  }

  # recalculate start, end, slice, strand
  $self->recalculate_coordinates();
}



=head2 get_all_Exons

  Arg [1]    : none
  Example    : my @exons = @{$transcript->get_all_Exons()};
  Description: Returns an listref of the exons in this transcipr in order.
               i.e. the first exon in the listref is the 5prime most exon in 
               the transcript.
  Returntype : a list reference to Bio::EnsEMBL::Exon objects
  Exceptions : none
  Caller     : general

=cut

sub get_all_Exons {
   my ($self) = @_;
   if( ! defined $self->{'_trans_exon_array'} && defined $self->adaptor() ) {
     $self->{'_trans_exon_array'} = $self->adaptor()->db()->
       get_ExonAdaptor()->fetch_all_by_Transcript( $self );
   }
   return $self->{'_trans_exon_array'};
}

=head2 get_all_Introns

  Arg [1]    : none
  Example    : my @introns = @{$transcript->get_all_Introns()};
  Description: Returns an listref of the introns in this transcipr in order.
               i.e. the first intron in the listref is the 5prime most exon in 
               the transcript.
  Returntype : a list reference to Bio::EnsEMBL::Intron objects
  Exceptions : none
  Caller     : general

=cut

sub get_all_Introns {
   my ($self) = @_;
   if( ! defined $self->{'_trans_exon_array'} && defined $self->adaptor() ) {
     $self->{'_trans_exon_array'} = $self->adaptor()->db()->
       get_ExonAdaptor()->fetch_all_by_Transcript( $self );
   }

   my @introns=();
   my @exons = @{$self->{'_trans_exon_array'}};
   for(my $i=0; $i < scalar(@exons)-1; $i++){
     my $intron = new Bio::EnsEMBL::Intron($exons[$i],$exons[$i+1]);
     push(@introns, $intron)
   }
   return \@introns;
}



=head2 length


    my $t_length = $transcript->length

Returns the sum of the length of all the exons in
the transcript.

=cut

sub length {
  my( $self ) = @_;

  my $length = 0;
  foreach my $ex (@{$self->get_all_Exons}) {
    $length += $ex->length;
  }

  # adjust the length if post transcriptional edits are enabled
  if($self->edits_enabled()) {
    foreach my $se (@{$self->get_all_SeqEdits()}) {
      $length += $se->length_diff();
    }
  }

  return $length;
}





=head2 flush_Exons

  Arg [1]    : none
  Example    : $transcript->flush_Exons();
  Description: Removes all Exons from this transcript and flushes related
               internal caches.
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub flush_Exons{
   my ($self,@args) = @_;
   $self->{'transcript_mapper'} = undef;
   $self->{'coding_region_start'} = undef;
   $self->{'coding_region_end'} = undef;
   $self->{'cdna_coding_start'} = undef;
   $self->{'cdna_coding_end'} = undef;
   $self->{'start'} = undef;
   $self->{'end'} = undef;
   $self->{'strand'} = undef;

   $self->{'_trans_exon_array'} = [];
}



=head2 five_prime_utr

  Arg [1]    : none
  Example    : my $five_prime  = $transcrpt->five_prime_utr
                 or warn "No five prime UTR";
  Description: Obtains a Bio::Seq object of the five prime UTR of this
               transcript.  If this transcript is a pseudogene
               (i.e. non-translating) or has no five prime UTR undef is
               returned instead.
  Returntype : Bio::Seq or undef
  Exceptions : none
  Caller     : general

=cut

sub five_prime_utr {
  my $self = shift;

  my $cdna_coding_start  = $self->cdna_coding_start();

  return undef if(!$cdna_coding_start);

  my $seq = substr($self->spliced_seq, 0, $cdna_coding_start - 1);

  return undef if(!$seq);

  return Bio::Seq->new(
	       -DISPLAY_ID => $self->stable_id,
	       -MOLTYPE    => 'dna',
	       -SEQ        => $seq);
}



=head2 three_prime_utr

  Arg [1]    : none
  Example    : my $three_prime  = $transcrpt->three_prime_utr
                 or warn "No five prime UTR";
  Description: Obtains a Bio::Seq object of the three prime UTR of this
               transcript.  If this transcript is a pseudogene
               (i.e. non-translating) or has no three prime UTR,
               undef is returned instead.
  Returntype : Bio::Seq or undef
  Exceptions : none
  Caller     : general

=cut

sub three_prime_utr {
  my $self = shift;

  my $cdna_coding_end = $self->cdna_coding_end();

  return undef if(!$cdna_coding_end);

  my $seq = substr($self->spliced_seq, $cdna_coding_end);

  return undef if(!$seq);

  return Bio::Seq->new(
	       -DISPLAY_ID => $self->stable_id,
	       -MOLTYPE    => 'dna',
	       -SEQ        => $seq);
}


=head2 get_all_translateable_Exons

  Args       : none
  Example    : none
  Description: Returns a list of exons that translate with the
               start and end exons truncated to the CDS regions.
               This function does not take into account any SeqEdits
               (post transcriptional RNA modifictions) when constructing the
               the 'translateable' exons, and it does not update the phase
               information of the created 'translateable' exons.

               If this transcript is a pseudogene (i.e. non-translateable)
               a reference to an empty list is returned.

  Returntype : listref Bio::EnsEMBL::Exon
  Exceptions : throw if translation has invalid information
  Caller     : Genebuild

=cut


sub get_all_translateable_Exons {
  my ( $self ) = @_;

  #return an empty list if there is no translation (i.e. pseudogene)
  my $translation = $self->translation or return [];
  my $start_exon      = $translation->start_Exon;
  my $end_exon        = $translation->end_Exon;
  my $t_start         = $translation->start;
  my $t_end           = $translation->end;

  my( @translateable );

  foreach my $ex (@{$self->get_all_Exons}) {

    if ($ex ne $start_exon and ! @translateable) {
      next;   # Not yet in translated region
    }

    my $length  = $ex->length;

    my $adjust_start = 0;
    my $adjust_end = 0;
    # Adjust to translation start if this is the start exon
    if ($ex == $start_exon ) {
      if ($t_start < 1 or $t_start > $length) {
        throw("Translation start '$t_start' is outside exon $ex length=$length");
      }
      $adjust_start = $t_start - 1;
    }

    # Adjust to translation end if this is the end exon
    if ($ex == $end_exon) {
#      if ($t_end < 1 or $t_end > $length) {
#        throw("Translation end '$t_end' is outside exon $ex length=$length");
#      }
      $adjust_end = $t_end - $length;
    }

    # Make a truncated exon if the translation start or
    # end causes the coordinates to be altered.
    if ($adjust_end || $adjust_start) {
      my $newex = $ex->adjust_start_end( $adjust_start, $adjust_end );

      push( @translateable, $newex );
    } else {
      push(@translateable, $ex);
    }

    # Exit the loop when we've found the last exon
    last if $ex eq $end_exon;
  }
  return \@translateable;
}


=head2 translate

  Args       : none
  Example    : none
  Description: return the peptide (plus eventuel stop codon) for this 
               transcript. Does N padding of non phase matching exons. 
               It uses translateable_seq internally.
               Returns undef if this Transcript does not have a translation
               (i.e. pseudogene).
  Returntype : Bio::Seq or undef
  Exceptions : If no Translation is set in this Transcript
  Caller     : general

=cut

sub translate {
  my ($self) = @_;

  if(!$self->translation()) {
    return undef;
  }

  my $mrna = $self->translateable_seq();
  my $display_id;
  if( defined $self->translation->stable_id ) {
    $display_id = $self->translation->stable_id;
  } elsif ( defined $self->translation->dbID ) {
    $display_id = $self->translation->dbID();
  } else {
    #use memory location as temp id
    $display_id = scalar($self->translation());
  }

  # remove final stop codon from the mrna if it is present
  # produced peptides will not have '*' at end
  # if terminal stop codon is desired call translatable_seq directly
  # and produce a translation from it

  if( CORE::length( $mrna ) % 3 == 0 ) {
    $mrna =~ s/TAG$|TGA$|TAA$//i;
  }

  my $peptide = Bio::Seq->new( -seq      => $mrna,
                               -moltype  => "dna",
                               -alphabet => 'dna',
                               -id       => $display_id );

  # Alternative codon tables (such as the mitochondrial codon table) can
  # be sepcified for a sequence region via the seq_region_attrib table.
  # A list of codon tables and their codes is at:
  # http://www.ncbi.nlm.nih.gov/htbin-post/Taxonomy/wprintgc?mode=c

  my $codon_table;
  if($self->slice()) {
    my ($attrib) = @{$self->slice()->get_all_Attributes('codon_table')};
    $codon_table = $attrib->value() if($attrib);
  }

  $codon_table ||= 1; # default vertebrate codon table

  my $translation = $peptide->translate(undef,undef,undef,$codon_table);

  if($self->edits_enabled()) {
    $self->translation()->modify_translation( $translation );
  }

  return $translation;
}


=head2 seq

Returns a Bio::Seq object which consists of just
the sequence of the exons concatenated together,
without messing about with padding with N\'s from
Exon phases like B<dna_seq> does.

=cut

sub seq {
  my( $self ) = @_;
  return Bio::Seq->new
    (-DISPLAY_ID => $self->stable_id,
     -MOLTYPE    => 'dna',
     -SEQ        => $self->spliced_seq);
}


=head2 pep2genomic

  Description: See Bio::EnsEMBL::TranscriptMapper::pep2genomic

=cut

sub pep2genomic {
  my $self = shift;
  return $self->get_TranscriptMapper()->pep2genomic(@_);
}


=head2 genomic2pep

  Description: See Bio::EnsEMBL::TranscriptMapper::genomic2pep

=cut

sub genomic2pep {
  my $self = shift;
  return $self->get_TranscriptMapper()->genomic2pep(@_);
}


=head2 cdna2genomic

  Description: See Bio::EnsEMBL::TranscriptMapper::cdna2genomic

=cut

sub cdna2genomic {
  my $self = shift;
  return $self->get_TranscriptMapper()->cdna2genomic(@_);
}

=head2 genomic2cdna

  Description: See Bio::EnsEMBL::TranscriptMapper::genomic2cdna

=cut

sub genomic2cdna {
  my $self = shift;
  return $self->get_TranscriptMapper->genomic2cdna(@_);
}

=head2 get_TranscriptMapper

  Args       : none
  Example    : my $trans_mapper = $transcript->get_TranscriptMapper();
  Description: Gets a TranscriptMapper object which can be used to perform
               a variety of coordinate conversions relating this transcript,
               genomic sequence and peptide resulting from this transcripts
               translation.
  Returntype : Bio::EnsEMBL::TranscriptMapper
  Exceptions : none
  Caller     : cdna2genomic, pep2genomic, genomic2cdna, cdna2genomic

=cut

sub get_TranscriptMapper {
  my ( $self ) = @_;
  return $self->{'transcript_mapper'} ||=
    Bio::EnsEMBL::TranscriptMapper->new($self);
}



=head2 start_Exon

 Title   : start_Exon
 Usage   : $start_exon = $transcript->start_Exon;
 Returns : The first exon in the transcript.
 Args    : NONE

=cut

sub start_Exon{
  my $self = shift;
  return $self->get_all_Exons()->[0];
}

=head2 end_Exon

 Title   : end_exon
 Usage   : $end_exon = $transcript->end_Exon;
 Returns : The last exon in the transcript.
 Args    : NONE

=cut

sub end_Exon{
   my $self = shift;
   return $self->get_all_Exons()->[-1];
}



=head2 description

 Title   : description
 Usage   : $obj->description($newval)
 Function: 
 Returns : value of description
 Args    : newvalue (optional)


=cut

sub description{
  my $self = shift;
  $self->{'description'} = shift if( @_ );
  return $self->{'description'};
}


=head2 version

 Title   : version
 Usage   : $obj->version()
 Function: 
 Returns : value of version
 Args    : 

=cut

sub version{
  my $self = shift;
  $self->{'version'} = shift if( @_ );
  return $self->{'version'};
}


=head2 stable_id

 Title   : stable_id
 Usage   : $obj->stable_id
 Function: 
 Returns : value of stable_id
 Args    : 


=cut

sub stable_id{
  my $self = shift;
  $self->{'stable_id'} = shift if( @_ );
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



=head2 swap_exons

  Arg [1]    : Bio::EnsEMBL::Exon $old_Exon
               An exon that should be replaced
  Arg [2]    : Bio::EnsEMBL::Exon $new_Exon
               The replacement Exon
  Example    : none
  Description: exchange an exon in the current Exon list with a given one.
               Usually done before storing of Gene, so the Exons can
               be shared between Transcripts.
  Returntype : none
  Exceptions : none
  Caller     : GeneAdaptor->store()

=cut

sub swap_exons {
  my ( $self, $old_exon, $new_exon ) = @_;
  
  my $arref = $self->{'_trans_exon_array'};
  for(my $i = 0; $i < @$arref; $i++) {
    if($arref->[$i] == $old_exon) {
      $arref->[$i] = $new_exon;
      last;
    }
  }

  if( defined $self->{'translation'} ) {
     if( $self->translation()->start_Exon() == $old_exon ) {
      $self->translation()->start_Exon( $new_exon );
    }
    if( $self->translation()->end_Exon() == $old_exon ) {
      $self->translation()->end_Exon( $new_exon );
    }
  }
}

=head2 transform

  Arg  1     : String $coordinate_system_name
  Arg [2]    : String $coordinate_system_version
  Example    : $transcript = $transcript->transform('contig');
               $transcript = $transcript->transform('chromosome', 'NCBI33');
  Description: Moves this Transcript to the given coordinate system.
               If this Transcript has Exons attached, they move as well.
               A new Transcript is returned. If the transcript cannot be
               transformed to the destination coordinate system undef is
               returned instead.
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : wrong parameters
  Caller     : general

=cut


sub transform {
  my $self = shift;

  # catch for old style transform calls
  if( ref $_[0] eq 'HASH') {
    deprecate("Calling transform with a hashref is deprecate.\n" .
              'Use $trans->transfer($slice) or ' .
              '$trans->transform("coordsysname") instead.');
    my (undef, $new_ex) = each(%{$_[0]});
    return $self->transfer($new_ex->slice);
  }

  my $new_transcript = $self->SUPER::transform( @_ );
  return undef unless $new_transcript;

  if( defined $self->{'translation'} ) {
    my $new_translation;
    %$new_translation = %{$self->{'translation'}};;
    bless $new_translation, ref( $self->{'translation'} );
    $new_transcript->{'translation'} = $new_translation;
  }

  if( exists $self->{'_trans_exon_array'} ) {
    my @new_exons;
    for my $old_exon ( @{$self->{'_trans_exon_array'}} ) {
      my $new_exon = $old_exon->transform( @_ );
      if( defined $new_transcript->{'translation'} ) {
        if( $new_transcript->translation()->start_Exon() == $old_exon ) {
          $new_transcript->translation()->start_Exon( $new_exon );
        }
        if( $new_transcript->translation()->end_Exon() == $old_exon ) {
          $new_transcript->translation()->end_Exon( $new_exon );
        }
      }
      push( @new_exons, $new_exon );
    }
    $new_transcript->{'_trans_exon_array'} = \@new_exons;
  }

  # flush cached internal values that depend on the exon coords
  $new_transcript->{'transcript_mapper'} = undef;
  $new_transcript->{'coding_region_start'} = undef;
  $new_transcript->{'coding_region_end'} = undef;
  $new_transcript->{'cdna_coding_start'} = undef;
  $new_transcript->{'cdna_coding_end'} = undef;

  return $new_transcript;
}


=head2 transfer

  Arg  1     : Bio::EnsEMBL::Slice $destination_slice
  Example    : $transcript = $transcript->transfer($slice);
  Description: Moves this transcript to the given slice.
               If this Transcripts has Exons attached, they move as well.
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general

=cut


sub transfer {
  my $self = shift;

  my $new_transcript = $self->SUPER::transfer( @_ );
  return undef unless $new_transcript;

  if( defined $self->{'translation'} ) {
    my $new_translation;
    %$new_translation = %{$self->{'translation'}};;
    bless $new_translation, ref( $self->{'translation'} );
    $new_transcript->{'translation'} = $new_translation;
  }

  if( exists $self->{'_trans_exon_array'} ) {
    my @new_exons;
    for my $old_exon ( @{$self->{'_trans_exon_array'}} ) {
      my $new_exon = $old_exon->transfer( @_ );
      if( defined $new_transcript->{'translation'} ) {
        if( $new_transcript->translation()->start_Exon() == $old_exon ) {
          $new_transcript->translation()->start_Exon( $new_exon );
        }
        if( $new_transcript->translation()->end_Exon() == $old_exon ) {
          $new_transcript->translation()->end_Exon( $new_exon );
        }
      }
      push( @new_exons, $new_exon );
    }

    $new_transcript->{'_trans_exon_array'} = \@new_exons;
  }

  # flush cached internal values that depend on the exon coords
  $new_transcript->{'transcript_mapper'} = undef;
  $new_transcript->{'coding_region_start'} = undef;
  $new_transcript->{'coding_region_end'} = undef;
  $new_transcript->{'cdna_coding_start'} = undef;
  $new_transcript->{'cdna_coding_end'} = undef;

  return $new_transcript;
}




=head recalculate_coordinates

  Args       : none
  Example    : none
  Description: called when exon coordinate change happened to recalculate the
               coords of the transcript.  This method should be called if one
               of the exons has been changed.
  Returntype : none
  Exceptions : none
  Caller     : internal

=cut

sub recalculate_coordinates {
  my $self = shift;

  my $exons = $self->get_all_Exons();

  return if(!$exons || !@$exons);

  my ( $slice, $start, $end, $strand );
  my $e_index;
  for ($e_index=0; $e_index<@$exons; $e_index++) {
    my $e = $exons->[$e_index];
    next if (!defined($e) or !defined($e->start)); # Skip missing or unmapped exons!
    $slice = $e->slice();
    $strand = $e->strand();
    $start = $e->start();
    $end = $e->end();
    last;
  }

  my $transsplicing = 0;

  # Start loop after first exon with coordinates
  for (; $e_index<@$exons; $e_index++) {
    my $e = $exons->[$e_index];
    next if (!defined($e) or !defined($e->start)); # Skip missing or unmapped exons!
    if( $e->start() < $start ) {
      $start = $e->start();
    }

    if( $e->end() > $end ) {
      $end = $e->end();
    }

    if( $slice && $e->slice() && $e->slice()->name() ne $slice->name() ) {
      throw( "Exons with different slices not allowed on one Transcript" );
    }

    if( $e->strand() != $strand ) {
      $transsplicing = 1;
    }
  }
  if( $transsplicing ) {
    warning( "Transcript contained trans splicing event" );
  }

  $self->start( $start );
  $self->end( $end );
  $self->strand( $strand );
  $self->slice( $slice );

  # flush cached internal values that depend on the exon coords
  $self->{'transcript_mapper'} = undef;
  $self->{'coding_region_start'} = undef;
  $self->{'coding_region_end'} = undef;
  $self->{'cdna_coding_start'} = undef;
  $self->{'cdna_coding_end'} = undef;
}


=head2 display_id

  Arg [1]    : none
  Example    : print $transcript->display_id();
  Description: This method returns a string that is considered to be
               the 'display' identifier.  For transcripts this is the 
               stable id if it is available otherwise it is an empty string.
  Returntype : string
  Exceptions : none
  Caller     : web drawing code

=cut

sub display_id {
  my $self = shift;
  return $self->{'stable_id'} || '';
}


=head2 get_all_peptide_variations

Description: See Bio::EnsEMBL::Utils::TranscriptSNPs::get_all_peptide_variations

=cut

sub get_all_peptide_variations {
  my ($self, $source, $snps) = @_;

  if(!$snps) {
    my $shash = Bio::EnsEMBL::Utils::TranscriptSNPs::get_all_cdna_SNPs($self, $source);
    $snps = $shash->{'coding'};
  }

  return Bio::EnsEMBL::Utils::TranscriptSNPs::get_all_peptide_variations($self,
                                                                        $snps);
}

=head2 get_all_SNPs

Description: See Bio::EnsEMBL::Utils::TranscriptSNPs::get_all_SNPs

=cut

sub get_all_SNPs {
  return Bio::EnsEMBL::Utils::TranscriptSNPs::get_all_SNPs(@_);
}

=head2 get_all_peptide_variations

Description: See Bio::EnsEMBL::Utils::TranscriptSNPs::get_all_cdna_SNPs

=cut

sub get_all_cdna_SNPs {
  return Bio::EnsEMBL::Utils::TranscriptSNPs::get_all_cdna_SNPs(@_);
}

###########################
# DEPRECATED METHODS FOLLOW
###########################

=head2 sort

  Description: DEPRECATED.  This method is no longer needed.  Exons are sorted
               automatically when added to the transcript.

=cut

sub sort {
  my $self = shift;

  deprecate( "Exons are kept sorted, you dont have to call sort any more" );
  # Fetch all the features
  my @exons = @{$self->get_all_Exons()};
  
  # Empty the feature table
  $self->flush_Exons();

  # Now sort the exons and put back in the feature table
  my $strand = $exons[0]->strand;

  if ($strand == 1) {
    @exons = sort { $a->start <=> $b->start } @exons;
  } elsif ($strand == -1) {
    @exons = sort { $b->start <=> $a->start } @exons;
  }

  foreach my $e (@exons) {
    $self->add_Exon($e);
  }
}


# _translation_id
# Usage   : DEPRECATED - not needed anymore

sub _translation_id {
   my $self = shift;
   deprecate( "This method shouldnt be necessary any more" );
   if( @_ ) {
      my $value = shift;
      $self->{'_translation_id'} = $value;
    }
    return $self->{'_translation_id'};

}

=head2 created

 Description: DEPRECATED - this attribute is not part of transcript anymore

=cut

sub created{
   my $obj = shift;
   deprecate( "This attribute is no longer supported" );
   if( @_ ) {
      my $value = shift;
      $obj->{'created'} = $value;
    }
    return $obj->{'created'};
}


=head2 modified

  Description: DEPRECATED - this attribute is not part of transcript anymore

=cut

sub modified{
   my $obj = shift;
   deprecate( "This attribute is no longer supported" );
   if( @_ ) {
      my $value = shift;
      $obj->{'modified'} = $value;
    }
    return $obj->{'modified'};
}

=head2 temporary_id

 Function: DEPRECATED: Use dbID or stable_id or something else instead

=cut

sub temporary_id{
   my ($obj,$value) = @_;
   deprecate( "I cant see what a temporary_id is good for, please use dbID" .
               "or stableID or\ntry without an id." );
   if( defined $value) {
      $obj->{'temporary_id'} = $value;
    }
    return $obj->{'temporary_id'};
}


=head2 get_all_DASFactories

  Arg [1]   : none
  Function  : Retrieves a listref of registered DAS objects
  Returntype: [ DAS_objects ]
  Exceptions:
  Caller    :
  Example   : $dasref = $prot->get_all_DASFactories

=cut

sub get_all_DASFactories {
   my $self = shift;
   return [ $self->adaptor()->db()->_each_DASFeatureFactory ];
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

  my $db = $self->adaptor->db;
  my $GeneAdaptor = $db->get_GeneAdaptor;
  my $Gene = $GeneAdaptor->fetch_by_transcript_stable_id($self->stable_id);	
  my $slice = $Gene->feature_Slice;

  foreach my $dasfact( @{$self->get_all_DASFactories} ){
    my $dsn = $dasfact->adaptor->dsn;
    my $type = $dasfact->adaptor->type;
    my $key = defined($dasfact->adaptor->url) ? $dasfact->adaptor->url .'/'. $dsn : $dasfact->adaptor->protocol .'://'.$dasfact->adaptor->domain.'/'. $dsn;
    if( $self->{_das_features}->{$key} ){ # Use cached
		  $das_features{$key} = $self->{_das_features}->{$key};
		  next;
    } else{ # Get fresh data
		  my @featref = ($type eq 'ensembl_location') ?  ($key, ($dasfact->fetch_all_by_Slice( $slice ))[0]) : $dasfact->fetch_all_by_DBLink_Container( $self );
		  $self->{_das_features}->{$key} = [@featref];
		  $das_features{$key} = [@featref];
	 }
  }
  return \%das_features;
}



1;
