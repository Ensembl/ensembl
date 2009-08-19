=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Transcript - object representing an Ensembl transcript

=head1 SYNOPSIS

Creation:

  my $tran = new Bio::EnsEMBL::Transcript();
  my $tran = new Bio::EnsEMBL::Transcript( -EXONS => \@exons );

Manipulation:

  # Returns an array of Exon objects
  my @exons = @{ $tran->get_all_Exons() };

  # Returns the peptide translation of the exons as a Bio::Seq
  if ( $tran->translation() ) {
    my $pep = $tran->translate();
  } else {
    print "Transcript ", $tran->stable_id(), " is non-coding\n";
  }

=head1 DESCRIPTION

A representation of a transcript within the Ensembl system.  A transcript
consists of a set of Exons and (possibly) a Translation which defines the
coding and non-coding regions of the exons.

=head1 METHODS

=cut

package Bio::EnsEMBL::Transcript;

use strict;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Intron;
use Bio::EnsEMBL::TranscriptMapper;
use Bio::EnsEMBL::Utils::TranscriptSNPs;
use Bio::EnsEMBL::SeqEdit;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( deprecate warning throw );

use vars qw(@ISA);
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
  Arg [-CREATED_DATE]:
        string - the date the transcript was created
  Arg [-MODIFIED_DATE]:
        string - the date the transcript was last modified
  Arg [-DESCRIPTION]:
        string - the transcipts description
  Arg [-BIOTYPE]: 
        string - the biotype e.g. "protein_coding"
  Arg [-STATUS]:
        string - the transcripts status i.e. "KNOWN","NOVEL"
  Arg [-IS_CURRENT]:
        Boolean - specifies if this is the current version of the transcript
  Example    : $tran = new Bio::EnsEMBL::Transcript(-EXONS => \@exons);
  Description: Constructor. Instantiates a Transcript object.
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : throw on bad arguments
  Caller     : general
  Status     : Stable

=cut

sub new {
  my ($class) = shift;

  if (ref $class) { 
      $class = ref $class;
  }

  my $self = $class->SUPER::new(@_);

  my ( $exons, $stable_id, $version, $external_name, $external_db,
       $external_status, $display_xref, $created_date, $modified_date,
       $description, $biotype, $confidence, $external_db_name, $status,
       $is_current );

  #catch for old style constructor calling:
  if((@_ > 0) && ref($_[0])) {
    $exons = [@_];
    deprecate("Transcript constructor should use named arguments.\n" .
              'Use Bio::EnsEMBL::Transcript->new(-EXONS => \@exons);' .
              "\ninstead of Bio::EnsEMBL::Transcript->new(\@exons);");
}
  else {
      ( $exons, $stable_id, $version, $external_name, $external_db,
	$external_status, $display_xref, $created_date, $modified_date,
	$description, $biotype, $confidence, $external_db_name, $status,
	$is_current ) = 
	    rearrange( [ "EXONS", 'STABLE_ID', 'VERSION', 'EXTERNAL_NAME', 
			 'EXTERNAL_DB', 'EXTERNAL_STATUS', 'DISPLAY_XREF',
			 'CREATED_DATE', 'MODIFIED_DATE', 'DESCRIPTION',
			 'BIOTYPE', 'CONFIDENCE', 'EXTERNAL_DB_NAME', 'STATUS',
			 'IS_CURRENT' ], @_ );
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

  $self->description( $description );
  $self->status( $confidence );  # old style name
  $self->status( $status );      # new style name
  $self->biotype( $biotype );

  # default is_current
  $is_current = 1 unless (defined($is_current));
  $self->{'is_current'} = $is_current;

  return $self;
}


=head2 get_all_DBLinks

  Example    : my @dblinks = @{ $transcript->get_all_DBLinks };
  Description: Retrieves _all_ related DBEntries for this transcript.  
               This includes all DBEntries that are associated with the
               corresponding translation.

               If you only want to retrieve the DBEntries associated with the
               transcript then you should use the get_all_DBEntries call 
               instead.
  Returntype : Listref of Bio::EnsEMBL::DBEntry objects, sorted by
               priority (desc), external db name (asc), display_id (asc)
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_DBLinks {
  my $self = shift;
  my $ex_db_exp = shift;
  my $ex_db_type = shift;

  my @links;

  push @links, @{$self->get_all_DBEntries($ex_db_exp, $ex_db_type)};

  my $transl = $self->translation();
  push @links, @{$transl->get_all_DBEntries($ex_db_exp, $ex_db_type)} if($transl);

  @links = sort {_compare_xrefs()} @links;

  return \@links;
}


=head2 get_all_DBEntries

  Example    : my @dbentries = @{ $gene->get_all_DBEntries };
  Description: Retrieves DBEntries (xrefs) for this transcript.  
               This does _not_ include the corresponding translations 
               DBEntries (see get_all_DBLinks).

               This method will attempt to lazy-load DBEntries from a
               database if an adaptor is available and no DBEntries are present
               on the transcript (i.e. they have not already been added or 
               loaded).
  Returntype : Listref of Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : get_all_DBLinks, TranscriptAdaptor::store
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
  if(!defined $self->{$cache_name} && $self->adaptor()) {
    $self->{$cache_name} = 
      $self->adaptor->db->get_DBEntryAdaptor->fetch_all_by_Transcript($self, $ex_db_exp, $ex_db_type);
  }

  $self->{$cache_name} ||= [];

  return $self->{$cache_name};
}


=head2 add_DBEntry

  Arg [1]    : Bio::EnsEMBL::DBEntry $dbe
               The dbEntry to be added
  Example    : my $dbe = Bio::EnsEMBL::DBEntery->new(...);
               $transcript->add_DBEntry($dbe);
  Description: Associates a DBEntry with this transcript. Note that adding
               DBEntries will prevent future lazy-loading of DBEntries for this
               gene (see get_all_DBEntries).
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


=head2 get_all_supporting_features

  Example    : my @evidence = @{ $transcript->get_all_supporting_features };
  Description: Retreives any supporting features added manually by 
               calls to add_supporting_features.
  Returntype : Listref of Bio::EnsEMBL::FeaturePair objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_supporting_features {
  my $self = shift;

  if( !exists  $self->{_supporting_evidence} )  {
    if($self->adaptor) {
      my $tsfa = $self->adaptor->db->get_TranscriptSupportingFeatureAdaptor();
      $self->{_supporting_evidence} = $tsfa->fetch_all_by_Transcript($self);
    }
  }
  
  return $self->{_supporting_evidence} || [];
}


=head2 add_supporting_features

  Arg [1-N]  : Bio::EnsEMBL::FeaturePair $feature
               The supporting features to add
  Example    : $transcript->add_supporting_features(@features);
  Description: Adds a list of supporting features to this Transcript.
               The added features can be retieved by
               get_all_supporting_features().
  Returntype : none
  Exceptions : throw if any of the features are not FeaturePairs
               throw if any of the features are not in the same coordinate
               system as the Transcript
  Caller     : general
  Status     : Stable
 
=cut
 
sub add_supporting_features {
  my ($self, @features) = @_;

  return unless @features;
 
  $self->{_supporting_evidence} ||= [];
  
  # check whether this feature object has been added already
  FEATURE: foreach my $feature (@features) {

    if (!defined($feature) || ref($feature) eq "ARRAY") {
      throw("Element in transcript supporting features array is undefined or is an ARRAY for " . $self->dbID);
    }
    if (!$feature || !$feature->isa("Bio::EnsEMBL::FeaturePair")) {
      print "feature = " . $feature . "\n";
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
    push(@{$self->{_supporting_evidence}}, $feature);
  }
}


=head2 flush_supporting_features

  Example     : $transcript->flush_supporting_features;
  Description : Removes all supporting evidence from the transcript.
  Return type : (Empty) listref
  Exceptions  : none
  Caller      : general
  Status      : Stable

=cut

sub flush_supporting_features {
  my $self = shift;
  $self->{'_supporting_evidence'} = [];
}


=head2 external_db

  Arg [1]    : (optional) String - name of external db to set
  Example    : $transcript->external_db('HGNC');
  Description: Getter/setter for attribute external_db. The db is the one that 
               belongs to the external_name.  
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

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

  Arg [1]    : (optional) String - status of the external db
  Example    : $transcript->external_status('KNOWNXREF');
  Description: Getter/setter for attribute external_status. The status of
               the external db of the one that belongs to the external_name.
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

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

  Arg [1]    : (optional) String - the external name to set
  Example    : $transcript->external_name('BRCA2-001');
  Description: Getter/setter for attribute external_name.
  Returntype : String or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

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

  Example    : print "Transcript ".$transcript->stable_id." is KNOWN\n" if
                  $transcript->is_known;
  Description: Returns TRUE if this gene has a status of 'KNOWN'
  Returntype : TRUE if known, FALSE otherwise
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub is_known {
  my $self = shift;
  return ( $self->{'status'} eq "KNOWN" );
}


=head2 status

  Arg [1]    : string $status
  Example    : none
  Description: get/set for attribute status
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Medium Risk

=cut

sub status {
   my $self = shift;
  $self->{'status'} = shift if( @_ );
  return $self->{'status'};
}

=head2 biotype

  Arg [1]    : string $biotype
  Example    : none
  Description: get/set for attribute biotype
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub biotype {
   my $self = shift;
  $self->{'biotype'} = shift if( @_ );
  return ( $self->{'biotype'} || "protein_coding" );
}


=head2 display_xref

  Arg [1]    : (optional) Bio::EnsEMBL::DBEntry - the display xref to set
  Example    : $transcript->display_xref($db_entry);
  Description: Getter/setter for display_xref for this transcript.
  Returntype : Bio::EnsEMBL::DBEntry
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub display_xref {
  my $self = shift;
  $self->{'display_xref'} = shift if(@_);
  return $self->{'display_xref'};
}


=head2 translation

  Args       : None
  Example    : if ( $transcript->translation() ) {
                 print( $transcript->translation()->stable_id(), "\n" );
               } else {
                 print("Pseudogene\n");
               }
  Description: Getter/setter for the Translation object which
               defines the CDS (and as a result the peptide encoded
               by) this transcript.  This function will return
               undef if this transcript is a pseudogene, i.e. a
               non-translating transcript such as an ncRNA.  This
               is the accepted method of determining whether a
               transcript is a pseudogene or not.
  Returntype : Bio::EnsEMBL::Translation
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub translation {
  my ( $self, $value ) = @_;

  if ( defined($value) ) {
    if ( !ref($value) || !$value->isa('Bio::EnsEMBL::Translation') ) {
      throw("Bio::EnsEMBL::Translation argument expected.");
    }

    $self->{'translation'} = $value;

  } elsif ( !exists( $self->{'translation'} )
    && defined( $self->adaptor() ) )
  {
    $self->{'translation'} =
      $self->adaptor()->db()->get_TranslationAdaptor()
      ->fetch_by_Transcript($self);
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
  Status     : Stable

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
  Status     : Stable

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
  Status     : Stable

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
  Status     : Stable

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
  Status     : Stable

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
  Status     : Stable

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
  Status     : Stable

=cut

sub edits_enabled {
  my ( $self, $boolean ) = @_;

  if ( defined($boolean) ) {
    $self->{'edits_enabled'} = $boolean;

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
  Status     : Stable

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
  Status     : Stable

=cut

sub get_all_Attributes {
  my $self = shift;
  my $attrib_code = shift;

  if( ! exists $self->{'attributes' } ) {
    if(!$self->adaptor() ) {
      return [];
    }

    my $attribute_adaptor = $self->adaptor->db->get_AttributeAdaptor();
    $self->{'attributes'} = $attribute_adaptor->fetch_all_by_Transcript($self);
  }

  if( defined $attrib_code) {
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
  Status     : Stable

=cut

sub add_Attributes {
  my $self    = shift;
  my @attribs = @_;

  if ( !exists $self->{'attributes'} ) {
    $self->{'attributes'} = [];
  }

  my $seq_change = 0;
  for my $attrib (@attribs) {
    if ( !$attrib->isa("Bio::EnsEMBL::Attribute") ) {
      throw("Argument to add_Attribute "
          . "has to be an Bio::EnsEMBL::Attribute" );
    }
    push( @{ $self->{'attributes'} }, $attrib );
    if ( $attrib->code eq "_rna_edit" ) {
      $seq_change = 1;
    }
  }
  if ($seq_change) {
    foreach my $ex ( @{ $self->get_all_Exons() } ) {
      $ex->{'_trans_exon_array'} = undef;
      $ex->{'_seq_cache'}        = undef;
    }
    my $translation = $self->translation;
    if ( defined($translation) ) {
      $translation->{seq} = undef;
    }
  }

  # flush cdna coord cache b/c we may have added a SeqEdit
  $self->{'cdna_coding_start'} = undef;
  $self->{'cdna_coding_end'}   = undef;
  $self->{'transcript_mapper'} = undef;
} ## end sub add_Attributes


=head2 add_Exon

 Title   : add_Exon
 Usage   : $trans->add_Exon($exon)
 Returns : Nothing
 Args [1]: Bio::EnsEMBL::Exon object to add
 Args [2]: rank
 Exceptions: throws if not a valid Bio::EnsEMBL::Exon
           : or exon clasjes with another one
 Status  : Stable

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

  Arg [CONSTITUTIVE]    : Boolean
                          Only return constitutive exons if true (non-zero)

  Examples  :   my @exons = @{ $transcript->get_all_Exons() };

                my @exons =
                  @{ $transcript->get_all_Exons( -constitutive => 1 ) };

  Description: Returns an listref of the exons in this transcript
               in order, i.e. the first exon in the listref is the
               5prime most exon in the transcript.  Only returns
               constitutive exons if the CONSTITUTIVE argument is
               true.

  Returntype : a list reference to Bio::EnsEMBL::Exon objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Exons {
  my ( $self, @args ) = @_;

  my $constitutive;
  if (@args) {
    $constitutive = rearrange( 'CONSTITUTIVE', @args );
  }

  if (!defined( $self->{'_trans_exon_array'} )
    && defined( $self->adaptor() ) )
  {
    $self->{'_trans_exon_array'} =
      $self->adaptor()->db()->get_ExonAdaptor()
      ->fetch_all_by_Transcript($self);
  }

  my @result;
  if ( defined($constitutive) && $constitutive != 0 ) {
    foreach my $exon ( @{ $self->{'_trans_exon_array'} } ) {
      if ( $exon->is_constitutive() ) {
        push( @result, $exon );
      }
    }
  } else {
    @result = @{ $self->{'_trans_exon_array'} };
  }

  return \@result;
} ## end sub get_all_Exons

=head2 get_all_constitutive_Exons

  Arg        :  None

  Examples   :  my @exons = @{ $transcript->get_all_constitutive_Exons() };

  Description:  Returns an listref of the constitutive exons in this
                transcript in order, i.e. the first exon in the
                listref is the 5prime most exon in the transcript.

  Returntype : a list reference to Bio::EnsEMBL::Exon objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_constitutive_Exons {
  my ($self) = @_;
  return $self->get_all_Exons( '-constitutive' => 1 );
}

=head2 get_all_Introns

  Arg [1]    : none
  Example    : my @introns = @{$transcript->get_all_Introns()};
  Description: Returns an listref of the introns in this transcript in order.
               i.e. the first intron in the listref is the 5prime most exon in 
               the transcript.
  Returntype : a list reference to Bio::EnsEMBL::Intron objects
  Exceptions : none
  Caller     : general
  Status     : Stable

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

  Args       : none
  Example    : my $t_length = $transcript->length
  Description: Returns the sum of the length of all the exons in the transcript.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

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
  Status     : Stable

=cut

sub flush_Exons {
  my ($self) = @_;

  $self->{'transcript_mapper'}   = undef;
  $self->{'coding_region_start'} = undef;
  $self->{'coding_region_end'}   = undef;
  $self->{'cdna_coding_start'}   = undef;
  $self->{'cdna_coding_end'}     = undef;
  $self->{'start'}               = undef;
  $self->{'end'}                 = undef;
  $self->{'strand'}              = undef;

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
  Status     : Stable

=cut

sub five_prime_utr {
  my $self = shift;

  my $cdna_coding_start  = $self->cdna_coding_start();

  return undef if(!$cdna_coding_start);

  my $seq = substr($self->spliced_seq, 0, $cdna_coding_start - 1);

  return undef if(!$seq);

  return Bio::Seq->new(
	       -DISPLAY_ID => $self->display_id,
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
  Status     : Stable

=cut

sub three_prime_utr {
  my $self = shift;

  my $cdna_coding_end = $self->cdna_coding_end();

  return undef if(!$cdna_coding_end);

  my $seq = substr($self->spliced_seq, $cdna_coding_end);

  return undef if(!$seq);

  return Bio::Seq->new(
	       -DISPLAY_ID => $self->display_id,
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
  Status     : Stable

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
        warning("WARN: Translation start '$t_start' is outside exon $ex length=$length");
	return [];
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
  Description: Return the peptide (plus eventual stop codon) for
               this transcript.  Does N-padding of non-phase
               matching exons.  It uses translateable_seq
               internally.  Returns undef if this Transcript does
               not have a translation (i.e. pseudogene).
  Returntype : Bio::Seq or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub translate {
  my ($self) = @_;

  if ( !defined( $self->translation() ) ) { return undef }

  my $mrna = $self->translateable_seq();

  # Alternative codon tables (such as the mitochondrial codon table)
  # can be specified for a sequence region via the seq_region_attrib
  # table.  A list of codon tables and their codes is at:
  # http://www.ncbi.nlm.nih.gov/htbin-post/Taxonomy/wprintgc?mode=c

  my $codon_table_id;
  my ( $complete5, $complete3 );
  if ( defined( $self->slice() ) ) {
    my $attrib;

    ($attrib) = @{ $self->slice()->get_all_Attributes('codon_table') };
    if ( defined($attrib) ) {
      $codon_table_id = $attrib->value();
    }

    ($attrib) = @{ $self->slice()->get_all_Attributes('complete5') };
    if ( defined($attrib) ) {
      $complete5 = $attrib->value();
    }

    ($attrib) = @{ $self->slice()->get_all_Attributes('complete3') };
    if ( defined($attrib) ) {
      $complete3 = $attrib->value();
    }
  }
  $codon_table_id ||= 1;    # default vertebrate codon table

  # Remove final stop codon from the mrna if it is present.  Produced
  # peptides will not have '*' at end.  If terminal stop codon is
  # desired call translatable_seq directly and produce a translation
  # from it.

  if ( CORE::length($mrna) % 3 == 0 ) {
    my $codon_table =
      Bio::Tools::CodonTable->new( -id => $codon_table_id );

    if ( $codon_table->is_ter_codon( substr( $mrna, -3, 3 ) ) ) {
      substr( $mrna, -3, 3, '' );
    }
  }

  if ( CORE::length($mrna) < 1 ) { return undef }

  my $display_id = $self->translation->display_id()
    || scalar( $self->translation() );

  my $peptide = Bio::Seq->new( -seq      => $mrna,
                               -moltype  => 'dna',
                               -alphabet => 'dna',
                               -id       => $display_id );

  my $translation =
    $peptide->translate( undef, undef, undef, $codon_table_id, undef,
                         undef, $complete5, $complete3 );

  if ( $self->edits_enabled() ) {
    $self->translation()->modify_translation($translation);
  }

  return $translation;
} ## end sub translate


=head2 seq

  Description: Returns a Bio::Seq object which consists of just
             : the sequence of the exons concatenated together,
             : without messing about with padding with N\'s from
             : Exon phases like B<dna_seq> does.
  Args       : none
  Example    : none
  Returntype : Bio::Seq
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seq {
  my( $self ) = @_;
  return Bio::Seq->new
    (-DISPLAY_ID => $self->display_id,
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
  Status     : Stable

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
 Status  : Stable

=cut

sub start_Exon {
  my $self = shift;
  return $self->get_all_Exons()->[0];
}


=head2 end_Exon

 Title   : end_exon
 Usage   : $end_exon = $transcript->end_Exon;
 Returns : The last exon in the transcript.
 Args    : NONE
 Status  : Stable

=cut

sub end_Exon {
   my $self = shift;
   return $self->get_all_Exons()->[-1];
}


=head2 description

 Title   : description
 Usage   : $obj->description($newval)
 Function: 
 Returns : value of description
 Args    : newvalue (optional)
 Status  : Stable

=cut

sub description {
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
 Status  : Stable

=cut

sub version {
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
 Status  : Stable

=cut

sub stable_id {
  my $self = shift;
  $self->{'stable_id'} = shift if( @_ );
  return $self->{'stable_id'};
}


=head2 is_current

  Arg [1]    : Boolean $is_current
  Example    : $transcript->is_current(1)
  Description: Getter/setter for is_current state of this transcript.
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub is_current {
  my $self = shift;
  $self->{'is_current'} = shift if (@_);
  return $self->{'is_current'};
}


=head2 created_date

  Arg [1]    : (optional) string to be used for the created date
  Example    : none
  Description: get/set for attribute created date
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

  Arg [1]    : (optional) string to be used for the modified date
  Example    : none
  Description: get/set for attribute modified date
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
  Status     : Stable

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
  Status     : Medium Risk
             : deprecation needs to be removed at some time

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
  if( ! defined $new_transcript ) {
    my @segments = $self->project( @_ );
    # if it projects, maybe the exons transform well?
    # lazy load them here
    return undef if( ! @segments );  
    $self->get_all_Exons();
  }


  if( exists $self->{'_trans_exon_array'} ) {
    my @new_exons;
    my ( $low_start, $hi_end, $slice );
    # we want to check whether the transform preserved 5prime 3prime
    # ordering. This assumes 5->3 order. No complaints on transsplicing.

    my ( $last_new_start, $last_old_strand, 
	 $last_new_strand, $start_exon, $end_exon,
	$last_seq_region_name );
    my $first = 1;
    my $ignore_order = 0;
    my $order_broken = 0;

    for my $old_exon ( @{$self->{'_trans_exon_array'}} ) {      
      my $new_exon = $old_exon->transform( @_ );
      return undef if( !defined $new_exon );
      if( ! defined $new_transcript ) {
	if( !$first ) {
	  if( $old_exon->strand() != $last_old_strand ) {
	    # transsplicing, ignore ordering
	    $ignore_order = 1;  
	  }

	  if( $new_exon->slice()->seq_region_name() ne 
	      $last_seq_region_name ) {
	    return undef;
	  }

	  if( $last_new_strand == 1 and 
	      $new_exon->start() < $last_new_start ) {
	    $order_broken = 1;
	  }

	  if( $last_new_strand == -1 and
	      $new_exon->start() > $last_new_start ) {
	    $order_broken = 1;
	  }

	  if( $new_exon->start() < $low_start ) {
	    $low_start = $new_exon->start();
	  }
	  if( $new_exon->end() > $hi_end ) {
	    $hi_end = $new_exon->end();
	  }
	} else {
	  $first = 0;
	  $low_start = $new_exon->start();
	  $hi_end = $new_exon->end();
	}

	$last_seq_region_name = $new_exon->slice()->seq_region_name();
	$last_old_strand = $old_exon->strand();
	$last_new_start = $new_exon->start();
	$last_new_strand = $new_exon->strand();
      }

      if( defined $self->{'translation'} ) {
        if( $self->translation()->start_Exon() == $old_exon ) {
          $start_exon = $new_exon;
        }
        if( $self->translation()->end_Exon() == $old_exon ) {
          $end_exon = $new_exon;
        }
      }
      push( @new_exons, $new_exon );
    }

    if( $order_broken && !$ignore_order ) {
      warning( "Order of exons broken in transform of ".$self->dbID() ); 
      return undef;
    }

    if( !defined $new_transcript ) {
      %$new_transcript = %$self;
      bless $new_transcript, ref( $self );
      $new_transcript->start( $low_start );
      $new_transcript->end( $hi_end );
      $new_transcript->slice( $new_exons[0]->slice() );
      $new_transcript->strand( $new_exons[0]->strand() );
    }

    $new_transcript->{'_trans_exon_array'} = \@new_exons;

    # should be ok to do inside exon array loop
    # translations only exist together with the exons ...

    if( defined $self->{'translation'} ) {
      my $new_translation;
      %$new_translation = %{$self->{'translation'}};;
      bless $new_translation, ref( $self->{'translation'} );
      $new_transcript->{'translation'} = $new_translation;
      $new_translation->start_Exon( $start_exon );
      $new_translation->end_Exon( $end_exon );
    }
  }

  if( exists $self->{'_supporting_evidence'} ) {
    my @new_features;
    for my $old_feature ( @{$self->{'_supporting_evidence'}} ) {
      my $new_feature = $old_feature->transform( @_ );
      if (defined $new_feature) { 
        push @new_features, $new_feature;
      }
    }
    $new_transcript->{'_supporting_evidence'} = \@new_features;
  }


  # flush cached internal values that depend on the exon coords
  $new_transcript->{'transcript_mapper'}   = undef;
  $new_transcript->{'coding_region_start'} = undef;
  $new_transcript->{'coding_region_end'}   = undef;
  $new_transcript->{'cdna_coding_start'}   = undef;
  $new_transcript->{'cdna_coding_end'}     = undef;

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
  Status     : Stable

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

  if( exists $self->{'_supporting_evidence'} ) {
    my @new_features;
    for my $old_feature ( @{$self->{'_supporting_evidence'}} ) {
      my $new_feature = $old_feature->transfer( @_ );
      push( @new_features, $new_feature );
    }
    $new_transcript->{'_supporting_evidence'} = \@new_features;
  }


  # flush cached internal values that depend on the exon coords
  $new_transcript->{'transcript_mapper'}   = undef;
  $new_transcript->{'coding_region_start'} = undef;
  $new_transcript->{'coding_region_end'}   = undef;
  $new_transcript->{'cdna_coding_start'}   = undef;
  $new_transcript->{'cdna_coding_end'}     = undef;

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
  Status     : Stable

=cut

sub recalculate_coordinates {
  my ($self) = @_;

  my $exons = $self->get_all_Exons();

  if ( !$exons || !@$exons ) { return }

  my ( $slice, $start, $end, $strand );

  my $e_index;
  for ( $e_index = 0; $e_index < @$exons; $e_index++ ) {
    my $e = $exons->[$e_index];
    # Skip missing or unmapped exons!
    next if ( !defined($e) or !defined( $e->start ) );
    $slice  = $e->slice();
    $strand = $e->strand();
    $start  = $e->start();
    $end    = $e->end();
    last;
  }

  my $transsplicing = 0;

  # Start loop after first exon with coordinates
  for ( ; $e_index < @$exons; $e_index++ ) {
    my $e = $exons->[$e_index];

    # Skip missing or unmapped exons!
    if ( !defined($e) or !defined( $e->start ) ) { next }

    if ( $e->start() < $start ) {
      $start = $e->start();
    }

    if ( $e->end() > $end ) {
      $end = $e->end();
    }

    if ( $slice
      && $e->slice()
      && $e->slice()->name() ne $slice->name() )
    {
      throw("Exons with different slices "
          . "are not allowed on one Transcript" );
    }

    if ( $e->strand() != $strand ) {
      $transsplicing = 1;
    }
  }
  if ($transsplicing) {
    warning("Transcript contained trans splicing event");
  }

  $self->start($start);
  $self->end($end);
  $self->strand($strand);
  $self->slice($slice);

  # flush cached internal values that depend on the exon coords
  $self->{'transcript_mapper'}   = undef;
  $self->{'coding_region_start'} = undef;
  $self->{'coding_region_end'}   = undef;
  $self->{'cdna_coding_start'}   = undef;
  $self->{'cdna_coding_end'}     = undef;
} ## end sub recalculate_coordinates


=head2 display_id

  Arg [1]    : none
  Example    : print $transcript->display_id();
  Description: This method returns a string that is considered to be
               the 'display' identifier. For transcripts this is (depending on
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


=head2 get_all_peptide_variations

  Description: See Bio::EnsEMBL::Utils::TranscriptSNPs::get_all_peptide_variations
  Status  : At Risk
          : Will be replaced with modules from the ensembl-variation package


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

  Status  : At Risk
          : Will be replaced with modules from the ensembl-variation package

=cut

sub get_all_SNPs {
  return Bio::EnsEMBL::Utils::TranscriptSNPs::get_all_SNPs(@_);
}


=head2 get_all_cdna_SNPs

  Description: See Bio::EnsEMBL::Utils::TranscriptSNPs::get_all_cdna_SNPs
 
  Status  : At Risk
          : Will be replaced with modules from the ensembl-variation package

=cut

sub get_all_cdna_SNPs {
  return Bio::EnsEMBL::Utils::TranscriptSNPs::get_all_cdna_SNPs(@_);
}


=head2 get_all_DASFactories

  Arg [1]   : none
  Function  : Retrieves a listref of registered DAS objects
  Returntype: [ DAS_objects ]
  Exceptions:
  Caller    :
  Example   : $dasref = $prot->get_all_DASFactories
  Status    : Stable

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
  Status    : Stable


=cut

sub get_all_DAS_Features {
  my ($self,@args) = @_;

  my $db = $self->adaptor->db;
  my $GeneAdaptor = $db->get_GeneAdaptor;
  my $Gene = $GeneAdaptor->fetch_by_transcript_stable_id($self->stable_id);	
  my $slice = $Gene->feature_Slice;
  return $self->SUPER::get_all_DAS_Features($slice);
}



=head2 _compare_xrefs

  Description: compare xrefs based on priority (descending), then name (ascending),
               then display_label (ascending)

=cut

sub _compare_xrefs {
  # compare on priority first (descending)
  if ($a->priority() != $b->priority()) {
    return $b->priority() <=> $a->priority();
  } else { # equal priorities, compare on external_db name
    if ($a->dbname() ne $b->dbname()) {
      return $a->dbname() cmp $b->dbname();
    } else { # equal priorities and names, compare on display_label
      return $a->display_id() cmp $b->display_id();
    }
  }
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


=head2 type

  Description: DEPRECATED. Use biotype() instead.

=cut

sub type {
  deprecate("Use biotype() instead");
  biotype(@_);
}


=head2 confidence

  Description: DEPRECATED. Use status() instead.

=cut

sub confidence {
  deprecate("Use status() instead");
  status(@_);
}


1;

