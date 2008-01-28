package Bio::EnsEMBL::Translation;

=head1 NAME

Bio::EnsEMBL::Translation - A class representing the translation of a
transcript

=head1 SYNOPSIS

  my $translation = Bio::EnsEMBL::Translation->new(
    -START_EXON => $exon1,
    -END_EXON   => $exon2,
    -SEQ_START  => 98,
    -SEQ_END    => 39
  );

  # stable ID setter
  $translation->stable_id('ENSP00053458');

  # get start and end position in start/end exons
  my $start = $translation->start;
  my $end = $translation->end;

=head1 DESCRIPTION

A Translation object defines the CDS and UTR regions of a Transcript
through the use of start_Exon/end_Exon, and start/end attributes.

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut


use vars qw($AUTOLOAD @ISA);
use strict;

use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use Bio::EnsEMBL::Storable;

@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

  Arg [-START_EXON] : The Exon object in which the translation (CDS) starts
  Arg [-END_EXON]   : The Exon object in which the translation (CDS) ends
  Arg [-SEQ_START]  : The offset in the start_Exon indicating the start
                      position of the CDS.
  Arg [-SEQ_END]    : The offset in the end_Exon indicating the end
                      position of the CDS.
  Arg [-STABLE_ID]  : The stable identifier for this Translation
  Arg [-VERSION]    : The version of the stable identifier
  Arg [-DBID]       : The internal identifier of this Translation
  Arg [-ADAPTOR]    : The TranslationAdaptor for this Translation
  Arg [-SEQ]        : Manually sets the peptide sequence of this translation.
                      May be useful if this translation is not stored in
                      a database.
  Arg [-CREATED_DATE]: the date the translation was created
  Arg [-MODIFIED_DATE]: the date the translation was modified
  Example    : my $tl = Bio::EnsEMBL::Translation->new
                   (-START_EXON => $ex1,
                    -END_EXON   => $ex2,
                    -SEQ_START  => 98,
                    -SEQ_END    => 39);
  Description: Constructor.  Creates a new Translation object
  Returntype : Bio::EnsEMBL::Translation
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my ( $start_exon, $end_exon, $seq_start, $seq_end,
       $stable_id, $version, $dbID, $adaptor, $seq,
       $created_date, $modified_date ) = 
    rearrange( [ "START_EXON", "END_EXON", "SEQ_START", "SEQ_END",
                 "STABLE_ID", "VERSION", "DBID", "ADAPTOR",
                 "SEQ", "CREATED_DATE", "MODIFIED_DATE" ], @_ );

  my $self = bless {
		    'start_exon' => $start_exon,
		    'end_exon'   => $end_exon,
		    'adaptor'    => $adaptor,
		    'dbID'       => $dbID,
		    'start'      => $seq_start,
		    'end'        => $seq_end,
		    'stable_id'  => $stable_id,
		    'version'    => $version,
		    'created_date' => $created_date,
		    'modified_date' => $modified_date,
        'seq'        => $seq
		   }, $class;

  return $self;
}


=head2 start

  Arg [1]    : (optional) int $start - start position to set
  Example    : $translation->start(17);
  Description: Getter/setter for the value of start, which is a position within
               the exon given by start_Exon.

               If you need genomic coordinates, use the genomic_start()
               method.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub start{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      
      $obj->{'start'} = $value;
    }
    return $obj->{'start'};

}


=head2 end

  Arg [1]    : (optional) int $end - end position to set
  Example    : $translation->end(8);
  Description: Getter/setter for the value of end, which is a position within
               the exon given by end_Exon.

               If you need genomic coordinates, use the genomic_end()
               method.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub end {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      
      $self->{'end'} = $value;
    }
    return $self->{'end'};

}


=head2 start_Exon

  Arg [1]    : (optional) Bio::EnsEMBL::Exon - start exon to assign
  Example    : $translation->start_Exon($exon1);
  Description: Getter/setter for the value of start_Exon, which denotes the
               exon at which translation starts (and within this exon, at the
               position indicated by start, see above).
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : thrown on wrong argument type
  Caller     : general
  Status     : Stable

=cut

sub start_Exon {
   my $self = shift;

   if( @_ ) {
      my $value = shift;
      if( !ref $value || !$value->isa('Bio::EnsEMBL::Exon') ) {
         throw("Got to have an Exon object, not a $value");
      }
      $self->{'start_exon'} = $value;
    }
   return $self->{'start_exon'};
}


=head2 end_Exon

  Arg [1]    : (optional) Bio::EnsEMBL::Exon - start exon to assign
  Example    : $translation->start_Exon($exon1);
  Description: Getter/setter for the value of end_Exon, which denotes the
               exon at which translation ends (and within this exon, at the
               position indicated by end, see above).
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : thrown on wrong argument type
  Caller     : general
  Status     : Stable

=cut

sub end_Exon {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      if( !ref $value || !$value->isa('Bio::EnsEMBL::Exon') ) {
         throw("Got to have an Exon object, not a $value");
      }
      $self->{'end_exon'} = $value;
    } 

    return $self->{'end_exon'};
}

=head2 cdna_start

    Arg  [1]    : (optional) Bio::EnsEMBL::Transcript $transcript
                  The transcript which this is a translation of.
    Example     : $translation_cdna_start = $translation->cdna_start();
    Description : Returns the start position of the translation in cDNA
                  coordinates.
                  If no transcript is given, the method will use
                  TranscriptAdaptor->fetch_by_translation_id() to locate
                  the correct transcript.
    Return type : Integer
    Exceptions  : Throws if the given (optional) argument is not a
                  transcript.
    Caller      : General
    Status      : At Risk (Under Development)

=cut

sub cdna_start {
    my $self = shift;
    my ($transcript) = @_;

    if ( defined $transcript
         && (    !ref $transcript
              || !$transcript->isa('Bio::EnsEMBL::Transcript') ) )
    {
        throw("Argument is not a transcript");
    }

    if ( !exists $self->{'cdna_start'} ) {
        if ( !defined $transcript ) {
            # We were not given a transcript, get the transcript out of
            # the database.

            my $transcript_adaptor =
              $self->adaptor()->db()->get_TranscriptAdaptor();

            $transcript =
              $transcript_adaptor->fetch_by_translation_id(
                                                        $self->dbID() );
        }

        $self->{'cdna_start'} =
          $self->start_Exon()->cdna_coding_start($transcript);
    }

    return $self->{'cdna_start'};
} ## end sub cdna_start

=head2 cdna_end

    Arg  [1]    : (optional) Bio::EnsEMBL::Transcript $transcript
                  The transcript which this is a translation of.
    Example     : $translation_cdna_end = $translation->cdna_end();
    Description : Returns the end position of the translation in cDNA
                  coordinates.
                  If no transcript is given, the method will use
                  TranscriptAdaptor->fetch_by_translation_id() to locate
                  the correct transcript.
    Return type : Integer
    Exceptions  : Throws if the given (optional) argument is not a
                  transcript.
    Caller      : General
    Status      : At Risk (Under Development)

=cut

sub cdna_end {
    my $self = shift;
    my ($transcript) = @_;

    if ( defined $transcript
         && (    !ref $transcript
              || !$transcript->isa('Bio::EnsEMBL::Transcript') ) )
    {
        throw("Argument is not a transcript");
    }

    if ( !exists $self->{'cdna_end'} ) {
        if ( !defined $transcript ) {
            # We were not given a transcript, get the transcript out of
            # the database.

            my $transcript_adaptor =
              $self->adaptor()->db()->get_TranscriptAdaptor();

            $transcript =
              $transcript_adaptor->fetch_by_translation_id(
                                                        $self->dbID() );
        }

        $self->{'cdna_end'} =
          $self->end_Exon()->cdna_coding_end($transcript);
    }

    return $self->{'cdna_end'};
} ## end sub cdna_end

=head2 genomic_start

    Args        : None
    Example     : $translation_genomic_start =
                      $translation->genomic_start();
    Description : Returns the start position of the translation in
                  genomic coordinates on the forward strand.
    Return type : Integer
    Exceptions  : None
    Caller      : General
    Status      : At Risk (Under Development)

=cut

sub genomic_start {
    my $self = shift;

    if ( !exists $self->{'genomic_start'} ) {
        if ( $self->start_Exon()->strand() >= 0 ) {
            $self->{'genomic_start'} =
              $self->start_Exon()->start() + ( $self->start() - 1 );
        } else {
            $self->{'genomic_start'} =
              $self->end_Exon()->end() - ( $self->end() - 1 );
        }
    }

    return $self->{'genomic_start'};
}

=head2 genomic_end

    Args        : None
    Example     : $translation_genomic_end = $translation->genomic_end();
    Description : Returns the end position of the translation in genomic
                  coordinates on the forward strand.
    Return type : Integer
    Exceptions  : None
    Caller      : General
    Status      : At Risk (Under Development)

=cut

sub genomic_end {
    my $self = shift;

    if ( !exists $self->{'genomic_end'} ) {
        if ( $self->end_Exon()->strand() >= 0 ) {
            $self->{'genomic_end'} =
              $self->end_Exon()->start() + ( $self->end() - 1 );
        } else {
            $self->{'genomic_end'} =
              $self->start_Exon()->end() - ( $self->start() - 1 );
        }
    }

    return $self->{'genomic_end'};
}

=head2 version

  Arg [1]    : (optional) string $version - version to set
  Example    : $translation->version(2);
  Description: Getter/setter for attribute version
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub version {
   my $self = shift;
  $self->{'version'} = shift if( @_ );
  return $self->{'version'};
}


=head2 stable_id

  Arg [1]    : (optional) string $stable_id - stable ID to set
  Example    : $translation->stable_id('ENSP0059890');
  Description: Getter/setter for attribute stable_id
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub stable_id {
   my $self = shift;
  $self->{'stable_id'} = shift if( @_ );
  return $self->{'stable_id'};
}

=head2 created_date

  Arg [1]    : (optional) string $created_date - created date to set
  Example    : $translation->created_date('2007-01-10 20:52:00');
  Description: Getter/setter for attribute created date
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut 

sub created_date {
  my $self = shift;
  $self->{'created_date'} = shift if ( @_ );
  return $self->{'created_date'};
}


=head2 modified_date

  Arg [1]    : (optional) string $modified_date - modification date to set
  Example    : $translation->modified_date('2007-01-10 20:52:00');
  Description: Getter/setter for attribute modified date
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut 

sub modified_date {
  my $self = shift;
  $self->{'modified_date'} = shift if ( @_ );
  return $self->{'modified_date'};
}



=head2 transform

  Arg [1]    : hashref $old_new_exon_map
               a hash that maps old to new exons for a whole gene
  Description: maps start end end exon according to mapping table.
              If an exon is not mapped, just keep the old one.
  Returntype: none
  Exceptions : none
  Caller     : Transcript->transform() 
  Status     : Stable

=cut

sub transform {
  my $self = shift;
  my $href_exons = shift;

  my $start_exon = $self->start_Exon();
  my $end_exon = $self->end_Exon();

  if ( exists $href_exons->{$start_exon} ) {
    $self->start_Exon($href_exons->{$start_exon});
  } else {
    # do nothing, the start exon wasnt mapped
  }

  if ( exists $href_exons->{$end_exon} ) {
    $self->end_Exon($href_exons->{$end_exon});
  } else { 
    # do nothing, the end exon wasnt mapped
  }
}


=head2 get_all_DBEntries

  Arg [1]    : (optional) $ex_db_exp - external db name
  Example    : @dbentries = @{$translation->get_all_DBEntries()};
  Description: Retrieves DBEntries (xrefs) for this translation.  

               This method will attempt to lazy-load DBEntries from a
               database if an adaptor is available and no DBEntries are present
               on the translation (i.e. they have not already been added or 
               loaded).
  Returntype : list reference to Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : get_all_DBLinks, TranslationAdaptor::store
  Status     : Stable

=cut

sub get_all_DBEntries {
  my $self = shift;
  my $ex_db_exp = shift;
  my $ex_db_type = shift;

  my $cache_name = "dbentries";

  if(defined($ex_db_exp)){
    $cache_name .= $ex_db_exp;
  }
  if(defined($ex_db_type)){
    $cache_name .= $ex_db_type;
  }

  # if not cached, retrieve all of the xrefs for this gene
  if(!defined $self->{$cache_name}) {
    my $adaptor = $self->adaptor();
    my $dbID    = $self->dbID();

    return [] if(!$adaptor || !$dbID);
    $self->{$cache_name} =
      $self->adaptor->db->get_DBEntryAdaptor->fetch_all_by_Translation($self, $ex_db_exp, $ex_db_type);
  }

  $self->{$cache_name} ||= [];

  return $self->{$cache_name};
}


=head2 add_DBEntry

  Arg [1]    : Bio::EnsEMBL::DBEntry $dbe
               The dbEntry to be added
  Example    : $translation->add_DBEntry($xref);
  Description: Associates a DBEntry with this translation. Note that adding
               DBEntries will prevent future lazy-loading of DBEntries for this
               translation (see get_all_DBEntries).
  Returntype : none
  Exceptions : thrown on incorrect argument type
  Caller     : general
  Status     : Stable

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


=head2 get_all_DBLinks

  Arg [1]    : see get_all_DBEntries
  Example    : see get_all_DBEntries
  Description: This is here for consistancy with the Transcript and Gene 
               classes. It is a synonym for the get_all_DBEntries method.
  Returntype : see get_all_DBEntries
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_DBLinks {
  my $self = shift;
  return $self->get_all_DBEntries(@_);
}


=head2 get_all_ProteinFeatures

  Arg [1]    : (optional) string $logic_name
               The analysis logic_name of the features to retrieve.  If not
               specified, all features are retrieved instead.
  Example    : $features = $self->get_all_ProteinFeatures('PFam');
  Description: Retrieves all ProteinFeatures associated with this 
               Translation. If a logic_name is specified, only features with 
               that logic_name are returned.  If no logic_name is provided all
               associated protein_features are returned.

               ProteinFeatures are lazy-loaded from the database unless they
               added manually to the Translation or had already been loaded.
  Returntype : Bio::EnsEMBL::ProteinFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_ProteinFeatures {
  my $self = shift;
  my $logic_name = shift;

  if(!$self->{'protein_features'}) {
    my $adaptor = $self->adaptor();
    my $dbID    = $self->dbID();

    return [] if (!$adaptor || !$dbID);

    my %hash;
    $self->{'protein_features'} = \%hash;

    my $pfa = $adaptor->db()->get_ProteinFeatureAdaptor();
    my $name;
    foreach my $f (@{$pfa->fetch_all_by_translation_id($dbID)}) {
      my $analysis = $f->analysis();
      if($analysis) {
        $name = lc($f->analysis->logic_name());
      } else {
        warning("ProteinFeature has no attached analysis\n");
        $name = '';
      }
      $hash{$name} ||= [];
      push @{$hash{$name}}, $f;
    }
  }

  # a specific type of protein feature was requested
  if(defined($logic_name)) {
    $logic_name = lc($logic_name);
    return $self->{'protein_features'}->{$logic_name} || [];
  }

  my @features = ();

  # all protein features were requested
  foreach my $type (keys %{$self->{'protein_features'}}) {
    push @features, @{$self->{'protein_features'}->{$type}};
  }

  return \@features;    
}


=head2 get_all_DomainFeatures

  Example    : @domain_feats = @{$translation->get_all_DomainFeatures};
  Description: A convenience method which retrieves all protein features
               that are considered to be 'Domain' features.  Features which
               are 'domain' features are those with analysis logic names:
               'pfscan', 'scanprosite', 'superfamily', 'pfam', 'prints',
               'smart', 'pirsf', 'tigrfam'.
  Returntype : listref of Bio::EnsEMBL::ProteinFeatures
  Exceptions : none
  Caller     : webcode (protview)
  Status     : Stable

=cut

sub get_all_DomainFeatures{
 my ($self) = @_;

 my @features;

 my @types = ('pfscan',      #profile (prosite or pfam motifs)
              'scanprosite', #prosite
              'superfamily',
              'pfam',
              'smart',
              'tigrfam',
              'pirsf',
              'prints');

 foreach my $type (@types) {
   push @features, @{$self->get_all_ProteinFeatures($type)};
 }

 return \@features;
}


=head2 add_ProteinFeature

  Arg [1]    : Bio::EnsEMBL::ProteinFeature $pf
               The ProteinFeature to be added
  Example    : $translation->add_ProteinFeature($pf);
  Description: Associates a ProteinFeature with this translation. Note that
               adding ProteinFeatures will prevent future lazy-loading of
               ProteinFeatures for this translation (see
               get_all_ProteinFeatures).
  Returntype : none
  Exceptions : thrown on incorrect argument type
  Caller     : general
  Status     : Stable

=cut

sub add_ProteinFeature {
  my $self = shift;
  my $pf = shift;

  unless ($pf && ref($pf) && $pf->isa('Bio::EnsEMBL::ProteinFeature')) {
    throw('Expected ProteinFeature argument');
  }

  my $analysis = $pf->analysis;
  throw("ProteinFeature has no attached Analysis.") unless $analysis;

  push @{ $self->{'protein_features'}->{$analysis->logic_name} }, $pf;
}


=head2 display_id

  Example    : print $translation->display_id();
  Description: This method returns a string that is considered to be
               the 'display' identifier. For translations this is (depending on
               availability and in this order) the stable Id, the dbID or an
               empty string.
  Returntype : string
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub display_id {
  my $self = shift;
  return $self->{'stable_id'} || $self->dbID || '';
}


=head2 length

  Example    : print "Peptide length =", $translation->length();
  Description: Retrieves the length of the peptide sequence (i.e. number of
               amino acids) represented by this Translation object.
  Returntype : int
  Exceptions : none
  Caller     : webcode (protview etc.)
  Status     : Stable

=cut

sub length {
  my $self = shift;
  my $seq = $self->seq();

  return ($seq) ? CORE::length($seq) : 0;
}


=head2 seq

  Example    : print $translation->seq();
  Description: Retrieves a string representation of the peptide sequence
               of this Translation.  This retrieves the transcript from the
               database and gets its sequence, or retrieves the sequence which
               was set via the constructor.
  Returntype : string
  Exceptions : warning if the sequence is not set and cannot be retrieved from
               the database.
  Caller     : webcode (protview etc.)
  Status     : Stable

=cut

sub seq {
  my $self = shift;

  if(@_) {
    $self->{'seq'} = shift;
    return $self->{'seq'};
  }

  return $self->{'seq'} if($self->{'seq'});

  my $adaptor = $self->{'adaptor'};
  if(!$adaptor) {
    warning("Cannot retrieve sequence from Translation - adaptor is not set.");
  }

  my $dbID = $self->{'dbID'};
  if(!$dbID) {
    warning("Cannot retrieve sequence from Translation - dbID is not set.");
  }
  
  my $tr_adaptor = $self->{'adaptor'}->db()->get_TranscriptAdaptor;

  my $seq = $tr_adaptor->fetch_by_translation_id($dbID)->translate();
  if($seq){
    $self->{'seq'} = $seq->seq();
    return $self->{'seq'};
  }
  else{
    return ''; #empty string
  }
}


=head2 get_all_Attributes

  Arg [1]    : optional string $attrib_code
               The code of the attribute type to retrieve values for.
  Example    : ($sc_attr) = @{$tl->get_all_Attributes('_selenocysteine')};
               @tl_attributes = @{$translation->get_all_Attributes()};
  Description: Gets a list of Attributes of this translation.
               Optionally just get Attrubutes for given code.
               Recognized attribute "_selenocysteine"
  Returntype : listref Bio::EnsEMBL::Attribute
  Exceptions : warning if translation does not have attached adaptor and 
               attempts lazy load.
  Caller     : general, modify_translation
  Status     : Stable

=cut

sub get_all_Attributes {
  my $self = shift;
  my $attrib_code = shift;

  if( ! exists $self->{'attributes' } ) {
    if(!$self->adaptor() ) {
#      warning('Cannot get attributes without an adaptor.');
      return [];
    }

    my $aa = $self->adaptor->db->get_AttributeAdaptor();
    $self->{'attributes'} = $aa->fetch_all_by_Translation( $self );
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

  Arg [1..N] : Bio::EnsEMBL::Attribute $attribute
               Attributes to add.
  Example    : $translation->add_Attributes($selenocysteine_attribute);
  Description: Adds an Attribute to the Translation. Usefull to 
               do _selenocysteine.
               If you add an attribute before you retrieve any from database, 
               lazy load will be disabled.
  Returntype : none
  Exceptions : throw on incorrect arguments
  Caller     : general
  Status     : Stable

=cut

sub add_Attributes {
  my $self = shift;
  my @attribs = @_;

  if( ! exists $self->{'attributes'} ) {
    $self->{'attributes'} = [];
  }

  for my $attrib ( @attribs ) {
    if( ! $attrib->isa( "Bio::EnsEMBL::Attribute" )) {
      throw( "Argument to add_Attribute must be a Bio::EnsEMBL::Attribute" );
    }
    push( @{$self->{'attributes'}}, $attrib );
    $self->{seq}=undef;
  }
}


=head2 get_all_SeqEdits

  Example    : my @seqeds = @{$transcript->get_all_SeqEdits()};
  Description: Retrieves all post transcriptional sequence modifications for
               this transcript.
  Returntype : Bio::EnsEMBL::SeqEdit
  Exceptions : none
  Caller     : spliced_seq()
  Status     : Stable

=cut

sub get_all_SeqEdits {
  my $self = shift;

  my @seqeds;

  my $attribs;
  
  my @edits = ('initial_met', '_selenocysteine', 'amino_acid_sub');
  

  foreach my $edit(@edits){
    $attribs = $self->get_all_Attributes($edit);

    # convert attributes to SeqEdit objects
    foreach my $a (@$attribs) {
      push @seqeds, Bio::EnsEMBL::SeqEdit->new(-ATTRIB => $a);
    }
  }
  return \@seqeds;
}
  

=head2 modify_translation

  Arg [1]    : Bio::Seq $peptide 
  Example    : my $seq = Bio::Seq->new(-SEQ => $dna)->translate();
               $translation->modify_translation($seq);
  Description: Applies sequence edits such as selenocysteines to the Bio::Seq 
               peptide thats passed in
  Returntype : Bio::Seq
  Exceptions : none
  Caller     : Bio::EnsEMBL::Transcript->translate
  Status     : Stable

=cut

sub modify_translation {
  my ($self, $seq) = @_;

  my @seqeds = @{$self->get_all_SeqEdits()};

  # sort in reverse order to avoid complication of adjusting downstream edits
  @seqeds = sort {$b <=> $a} @seqeds;

  # apply all edits
  my $peptide = $seq->seq();
  foreach my $se (@seqeds) {
    $se->apply_edit(\$peptide);
  }
  $seq->seq($peptide);

  return $seq;
}


=head2 temporary_id

  Description: DEPRECATED This method should not be needed. Use dbID,
               stable_id or something else.

=cut

sub temporary_id {
   my $self = shift;
   deprecate( "I cant see what a temporary_id is good for, please use " .
               "dbID or stableID or\n try without an id." );
  $self->{'temporary_id'} = shift if( @_ );
  return $self->{'temporary_id'};
}


=head2 get_all_DASFactories

  Function  : Retrieves a listref of registered DAS objects
  Returntype: DAS objects
  Exceptions: none
  Caller    : webcode
  Example   : $dasref = $prot->get_all_DASFactories;
  Status    : Stable

=cut

sub get_all_DASFactories {
   my $self = shift;
   return [ $self->adaptor()->db()->_each_DASFeatureFactory ];
}


=head2 get_all_DAS_Features

  Example    : $features = $prot->get_all_DAS_Features;
  Description: Retreives a hash reference to a hash of DAS feature
               sets, keyed by the DNS, NOTE the values of this hash
               are an anonymous array containing:
                (1) a pointer to an array of features;
                (2) a pointer to the DAS stylesheet
  Returntype : hashref of Bio::SeqFeatures
  Exceptions : none
  Caller     : webcode
  Status     : Stable

=cut

sub get_all_DAS_Features{
  my $self = shift;

  my $db = $self->adaptor->db;
  my $GeneAdaptor = $db->get_GeneAdaptor;
  my $Gene = $GeneAdaptor->fetch_by_translation_stable_id($self->stable_id) || return;
  my $slice = $Gene->feature_Slice;
 
  return $self->SUPER::get_all_DAS_Features($slice);
}


1;
