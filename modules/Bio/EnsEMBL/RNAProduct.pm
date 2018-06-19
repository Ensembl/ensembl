=head1 LICENSE

Copyright [2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::RNAProduct - A class representing the mature RNA product
of a transcript

=head1 DESCRIPTION

TODO

=head1 SYNOPSIS

  my $rnaproduct = Bio::EnsEMBL::RNAProduct->new(
    -SEQ_START => 36,
    -SEQ_END   => 58
  );

  # Stable-ID setter
  $rnaproduct->stable_id('ENSM00090210');

  # Get start and end position in the precursor transcript
  my $start = $rnaproduct->start();
  my $end = $rnaproduct->end();

=cut


package Bio::EnsEMBL::RNAProduct;

use vars qw($AUTOLOAD);
use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref wrap_array );
use Scalar::Util qw(weaken);

use Bio::EnsEMBL::Storable;

use parent qw(Bio::EnsEMBL::Storable);


=head2 new

  Arg [-SEQ_START]    : The offset in the Translation indicating the start
                        position of the product sequence.
  Arg [-SEQ_END]      : The offset in the Translation indicating the end
                        position of the product sequence.
  Arg [-STABLE_ID]    : The stable identifier for this RNAPRoduct
  Arg [-VERSION]      : The version of the stable identifier
  Arg [-DBID]         : The internal identifier of this RNAProduct
  Arg [-ADAPTOR]      : The TranslationAdaptor for this RNAProduct
  Arg [-SEQ]          : Manually sets the nucleotide sequence of this
                        rnaproduct. May be useful if this rnaproduct is not
                        stored in a database.
  Arg [-CREATED_DATE] : the date the rnaproduct was created
  Arg [-MODIFIED_DATE]: the date the rnaproduct was modified
  Example    : my $rp = Bio::EnsEMBL::RNAProduct->new(
                 -SEQ_START => 36,
                 -SEQ_END   => 58
               );
  Description: Constructor.  Creates a new RNAProduct object
  Returntype : Bio::EnsEMBL::RNAProduct
  Exceptions : none
  Caller     : general
  Status     : In Development

=cut

sub new { ## no critic (Subroutines::RequireArgUnpacking)
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my ($seq_start, $seq_end, $stable_id, $version, $dbID, $adaptor, $seq,
      $created_date, $modified_date ) =
	rearrange(["SEQ_START", "SEQ_END", "STABLE_ID", "VERSION", "DBID",
		   "ADAPTOR", "SEQ", "CREATED_DATE", "MODIFIED_DATE"], @_);

  # For consistency between stable_id() and stable_id_version()
  $stable_id //= '';

  # Default version
  $version //= 1;

  my $self = bless {
    'start'      => $seq_start,
    'end'        => $seq_end,
    'stable_id'  => $stable_id,
    'version'    => $version,
    'dbID'       => $dbID,
    'seq'        => $seq,
    'created_date' => $created_date,
    'modified_date' => $modified_date
  }, $class;

  $self->adaptor($adaptor);

  return $self;
}


=head2 add_Attributes

  Arg [1..N] : Bio::EnsEMBL::Attribute $attribute
               Attributes to add.
  Example    : $rnaproduct->add_Attributes($selenocysteine_attribute);
  Description: Adds an Attribute to the RNAProduct.
               If you add an attribute before you retrieve any from database,
               lazy load will be disabled.
  Returntype : none
  Exceptions : throw on incorrect arguments
  Caller     : general
  Status     : Stable

=cut

sub add_Attributes {
  my ($self, @attribs) = @_;

  if (! exists $self->{'attributes'}) {
    $self->{'attributes'} = [];
  }

  for my $attrib (@attribs) {
    if (! $attrib->isa("Bio::EnsEMBL::Attribute")) {
      throw("Argument to add_Attribute must be a Bio::EnsEMBL::Attribute");
    }
    push (@{$self->{'attributes'}}, $attrib);
    # FIXME: why is this here???
    $self->{seq}=undef;
  }

  return;
}


=head2 add_DBEntry

  Arg [1]    : Bio::EnsEMBL::DBEntry $dbe
               The dbEntry to be added
  Example    : $rnaproduct->add_DBEntry($xref);
  Description: Associates a DBEntry with this rnaproduct. Note that adding
               DBEntries will prevent future lazy-loading of DBEntries for this
               rnaproduct (see get_all_DBEntries).
  Returntype : none
  Exceptions : thrown on incorrect argument type
  Caller     : general
  Status     : Stable

=cut

sub add_DBEntry {
  my ($self, $dbe) = @_;

  if (!$dbe || !ref($dbe) || !$dbe->isa('Bio::EnsEMBL::DBEntry')) {
    throw('Expected DBEntry argument');
  }

  $self->{'dbentries'} ||= [];
  push @{$self->{'dbentries'}}, $dbe;
}


=head2 cdna_end

    Example     : $rnaproduct_cdna_end = $rnaproduct->cdna_end();
    Description : Returns the end position of the rnaproduct in cDNA
                  coordinates.
                  Since rnaproducts do not span multiple exons, this is
                  simply an alias for end().
    Return type : Integer
    Caller      : General
    Status      : Stable

=cut

sub cdna_end {
  my $self = shift;

  return $self->end();
}


=head2 cdna_start

    Example     : $rnaproduct_cdna_start = $rnaproduct->cdna_start();
    Description : Returns the start position of the rnaproduct in cDNA
                  coordinates.
                  Since rnaproducts do not span multiple exons, this is
                  simply an alias for start().
    Return type : Integer
    Caller      : General
    Status      : Stable

=cut

sub cdna_start {
  my $self = shift;

  return $self->start();
}


=head2 created_date

  Arg [1]    : (optional) string $created_date - created date to set
  Example    : $rnaproduct->created_date('2007-01-10 20:52:00');
  Description: Getter/setter for attribute created_date
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub created_date {
  my $self = shift;
  $self->{'created_date'} = shift if (@_);
  return $self->{'created_date'};
}


=head2 display_id

  Example    : print $rnaproduct->display_id();
  Description: This method returns a string that is considered to be
               the 'display' identifier. For rnaproducts this is (depending on
               availability and in this order) the stable ID, the dbID or an
               empty string.
  Returntype : string
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub display_id {
  my $self = shift;
  return $self->stable_id() || $self->dbID() || '';
}


=head2 end

  Arg [1]    : (optional) int $end - end position to set
  Example    : $rnaproduct->end(39);
  Description: Getter/setter for the value of end, which is a position within
               the precursor Transcript.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub end {
  my $self = shift;
  $self->{'end'} = shift if (@_);
  return $self->{'end'};
}


=head2 genomic_end

    Args        : None
    Example     : $rnaproduct_genomic_end = $rnaproduct->genomic_end();
    Description : Returns the end position of the rnaproduct in genomic
                  coordinates on the forward strand.
    Return type : Integer
    Exceptions  : None
    Caller      : General
    Status      : At Risk (Under Development)

=cut

sub genomic_end {
  my $self = shift;

  if (!exists $self->{'genomic_end'}) {
    my $transcript = $self->transcript();

    if ($transcript->strand() >= 0) {
      $self->{'genomic_end'} =
	$transcript->start() + ($self->end() - 1);
    } else {
      $self->{'genomic_end'} =
	$transcript->end() - ($self->start() - 1);
    }
  }

  return $self->{'genomic_end'};
}


=head2 genomic_start

    Args        : None
    Example     : $rnaproduct_genomic_start = $rnaproduct->genomic_start();
    Description : Returns the start position of the rnaproduct in
                  genomic coordinates on the forward strand.
    Return type : Integer
    Exceptions  : None
    Caller      : General
    Status      : At Risk (Under Development)

=cut

sub genomic_start {
  my $self = shift;

  if (!exists $self->{'genomic_start'}) {
    my $transcript = $self->transcript();

    if ($transcript->strand() >= 0) {
      $self->{'genomic_start'} =
	$transcript->start() + ($self->start() - 1);
    } else {
      $self->{'genomic_start'} =
	$transcript->end() - ($self->end() - 1);
    }
  }

  return $self->{'genomic_start'};
}


=head2 get_all_Attributes

  Arg [1]    : optional string $attrib_code
               The code of the attribute type to retrieve values for.
  Example    : ($n_attr) = @{$tl->get_all_Attributes('note')};
               @rp_attributes = @{$rnaproduct->get_all_Attributes()};
  Description: Gets a list of Attributes of this rnaproduct.
               Optionally just get Attributes for given code.
  Returntype : listref Bio::EnsEMBL::Attribute
  Exceptions : Throws if there is no adaptor or defined for the rnaproduct
               object.
  Caller     : general
  Status     : Stable

=cut

sub get_all_Attributes {
  my ($self, $attrib_code) = @_;

  if (! exists $self->{'attributes'}) {
    if (!$self->adaptor()) {
      throw("Adaptor not set for rnaproduct, cannot fetch its attributes");
    }

    my $aa = $self->adaptor->db->get_AttributeAdaptor();
    $self->{'attributes'} = $aa->fetch_all_by_RNAProduct($self);
  }

  if (defined $attrib_code) {
    my @results = grep { uc($_->code()) eq uc($attrib_code) }
      @{$self->{'attributes'}};
    return \@results;
  } else {
    return $self->{'attributes'};
  }
}


=head2 get_all_DBEntries

  Arg [1]    : (optional) String, external database name,
               SQL wildcard characters (_ and %) can be used to
               specify patterns.

  Arg [2]    : (optional) String, external_db type,
               ('ARRAY','ALT_TRANS','ALT_GENE','MISC','LIT','PRIMARY_DB_SYNONYM','ENSEMBL'),
               SQL wildcard characters (_ and %) can be used to
               specify patterns.

  Example    : my @dbentries = @{ $rnaproduct->get_all_DBEntries() };
               @dbentries = @{ $rnaproduct->get_all_DBEntries('Uniprot%') };
               @dbentries = @{ $rnaproduct->get_all_DBEntries('%', 'ENSEMBL') };

  Description: Retrieves DBEntries (xrefs) for this rnaproduct.

               This method will attempt to lazy-load DBEntries
               from a database if an adaptor is available and no
               DBEntries are present on the rnaproduct (i.e. they
               have not already been added or loaded).

  Returntype : Listref to Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub get_all_DBEntries {
  my ($self, $ex_db_exp, $ex_db_type) = @_;

  my $cache_name = 'dbentries';

  if (defined($ex_db_exp)) {
    $cache_name .= $ex_db_exp;
  }

  if (defined($ex_db_type)) {
    $cache_name .= $ex_db_type;
  }

  # if not cached, retrieve all of the xrefs for this rnaproduct
  if (!defined($self->{$cache_name}) && defined($self->adaptor())) {
    $self->{$cache_name} = $self->adaptor()->db()->get_DBEntryAdaptor()->
      fetch_all_by_RNAProduct( $self, $ex_db_exp, $ex_db_type );
  }

  $self->{$cache_name} ||= [];

  return $self->{$cache_name};
}





=head2 get_all_DBLinks

  Arg [1]    : (optional) String, database name
               SQL wildcard characters (_ and %) can be used to
               specify patterns.

  Arg [2]    : (optional) String, external database type, can be one of
               ('ARRAY','ALT_TRANS','ALT_GENE','MISC','LIT','PRIMARY_DB_SYNONYM','ENSEMBL'),
               SQL wildcard characters (_ and %) can be used to
               specify patterns.

  Example    :  my @dblinks = @{ $rnaproduct->get_all_DBLinks() };
                @dblinks = @{ $rnaproduct->get_all_DBLinks('mirbase%') };
                @dblinks = @{ $rnaproduct->get_all_DBLinks('%', 'ENSEMBL') };

  Description: This is here for consistancy with the Transcript
               and Gene classes.  It is a synonym for the
               get_all_DBEntries() method.

  Return type: Listref to Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_DBLinks {
  my $self = shift;
  return $self->get_all_DBEntries(@_);
}


=head2 get_all_object_xrefs

  Arg [1]    : (optional) String, external database name

  Arg [2]    : (optional) String, external_db type

  Example    : @oxrefs = @{ $rnaproduct->get_all_object_xrefs() };

  Description: Retrieves xrefs for this rnaproduct.

               This method will attempt to lazy-load xrefs from a
               database if an adaptor is available and no xrefs
               are present on the rnaproduct (i.e. they have not
               already been added or loaded).

                NB: This method is an alias for the
                    get_all_DBentries() method.

  Return type: Listref of Bio::EnsEMBL::DBEntry objects

  Status     : Stable

=cut

sub get_all_object_xrefs {
  my $self = shift;
  return $self->get_all_DBEntries(@_);
}


=head2 get_all_xrefs

  Arg [1]    : String database name (optional)
               SQL wildcard characters (_ and %) can be used to
               specify patterns.

  Example    : @xrefs = @{ $rnaproduct->get_all_xrefs() };
               @xrefs = @{ $rnaproduct->get_all_xrefs('mirbase%') };

  Description: This method is here for consistancy with the Gene
               and Transcript classes.  It is an alias for the
               get_all_DBLinks() method, which in turn directly
               calls get_all_DBEntries().

  Return type: Listref of Bio::EnsEMBL::DBEntry objects

  Status     : Stable

=cut

sub get_all_xrefs {
  my $self = shift;
  return $self->get_all_DBLinks(@_);
}


=head2 modified_date

  Arg [1]    : (optional) string $modified_date - modification date to set
  Example    : $rnaproduct->modified_date('2007-01-10 20:52:00');
  Description: Getter/setter for attribute modified_date
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub modified_date {
  my $self = shift;
  $self->{'modified_date'} = shift if (@_);
  return $self->{'modified_date'};
}


=head2 length

  Example    : print "RNA length =", $rnaproduct->length();
  Description: Retrieves the length of the nucleotide sequence represented
               by this RNAProduct object.
  Returntype : int
  Exceptions : none
  Caller     : webcode (protview etc.)
  Status     : Stable

=cut

sub length { ## no critic (Subroutines::ProhibitBuiltinHomonyms)
  my $self = shift;
  my $seq = $self->seq();

  return ($seq) ? CORE::length($seq) : 0;
}


=head2 seq

  Example    : print $rnaproduct->seq();
  Description: Retrieves a string representation of the nucleotide sequence
               of this RNAProduct.  This retrieves the transcript from the
               database and gets its sequence, or retrieves the sequence which
               was set via the constructor/setter.
  Returntype : string
  Exceptions : warning if the sequence is not set and cannot be retrieved from
               the database.
  Caller     : webcode (protview etc.)
  Status     : Stable

=cut

sub seq {
  my ($self, $sequence) = @_;

  if (defined($sequence)) {

    $self->{'seq'} = $sequence;

  } elsif (!defined($self->{'seq'})) {

    my $tr_seq = $self->transcript()->seq();
    if ($tr_seq->length() <= 0) {
      throw('Got no or empty sequence from the database');
    }
    $self->{'seq'} = $tr_seq->subseq($self->{'start'}, $self->{'end'});

  }

  return '' unless defined($self->{'seq'});
  return $self->{'seq'};
}


=head2 stable_id

  Arg [1]    : (optional) string $stable_id - stable ID to set
  Example    : $rnaproduct->stable_id('ENSM00090210');
  Description: Getter/setter for attribute stable_id.
               Unlike stable_id_version(), setting a new stable ID does NOT
               reset the version number.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub stable_id {
  my $self = shift;
  $self->{'stable_id'} = shift if (@_);
  return $self->{'stable_id'};
}


=head2 stable_id_version

  Arg [1]    : (optional) String - the stable ID with version to set
  Example    : $rnaproduct->stable_id("ENSM0059890.3");
  Description: Getter/setter for stable id with version for this rnaproduct.
               If the input string omits the version part, the version gets reset
               to undef; use stable_id() if you want to avoid this.
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub stable_id_version {
  my $self = shift;

  if (my $stable_id = shift) {
    # If there is at least one embedded period assume everything
    # beyond the last one is the version number. This may not work for
    # some species, if you are worried about ambiguity use stable_id() +
    # version() explicitly.
    my $vindex = rindex($stable_id, '.');
    ($self->{stable_id},
     $self->{version}) = ($vindex > 0 ?
			  (substr($stable_id, 0, $vindex),
			   substr($stable_id, $vindex + 1)) :
			  $stable_id, undef
			 );
  }

  return $self->{stable_id} . ($self->{version} ? ".$self->{version}" : '');
}


=head2 start

  Arg [1]    : (optional) int $start - start position to set
  Example    : $rnaproduct->start(17);
  Description: Getter/setter for the value of start, which is a position within
               the precursor Transcript.
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub start {
  my $self = shift;
  $self->{'start'} = shift if (@_);
  return $self->{'start'};
}


=head2 summary_as_hash

  Example       : $rnaproduct_summary = $rnaproduct->summary_as_hash();
  Description   : Retrieves a textual summary of this RNAProduct.
                  Not inherited from Feature.
  Returns       : hashref of arrays of descriptive strings
  Status        : Intended for internal use

=cut

sub summary_as_hash {
  my $self = shift;
  my %summary;
  my $id = $self->display_id;
  if ($self->version) {
    $id .= "." . $self->version;
  }
  $summary{'id'} = $id;
  $summary{'rnaproduct_id'} = $id;
  $summary{'genomic_start'} = $self->genomic_start;
  $summary{'genomic_end'} = $self->genomic_end;
  $summary{'length'} = $self->length;
  my $transcript = $self->transcript;
  $summary{'Parent'} = $transcript->display_id;
  return \%summary;
}


=head2 transcript

  Arg [1]       : Transcript object (optional)
  Description   : Sets or retrieves the transcript object associated
                  with this rnaproduct object.
  Exceptions    : Throws if there is no adaptor or no dbID defined for
                  the rnaproduct object.
  Returntype    : Bio::EnsEMBL::Transcript

=cut

sub transcript {
  my ($self, $transcript) = @_;

  if (defined($transcript)) {

    # Normal setter
    assert_ref($transcript, 'Bio::EnsEMBL::Transcript');
    $self->{'transcript'} = $transcript;
    weaken($self->{'transcript'});    # Avoid circular references.

  } elsif (@_ > 1) {

    # User has explicitly passed undef. Break connection to transcript.
    delete( $self->{'transcript'} );

  } elsif (!defined($self->{'transcript'})) {
    my $adaptor = $self->{'adaptor'};
    if (!defined($adaptor)) {
      throw("Adaptor not set for rnaproduct, cannot fetch its transcript");
    }

    my $dbID = $self->{'dbID'};
    if (!defined($dbID)) {
      throw("dbID not set for rnaproduct, cannot fetch its transcript.");
    }

    $self->{'transcript'} =
      $adaptor->db()->get_TranscriptAdaptor()
      ->fetch_by_rnaproduct_id($dbID);

    # Do not weaken the reference if we had to get the transcript from the
    # database. The user is probably working on rnaproducts directly,
    # not going through transcripts.
  }

  return $self->{'transcript'};
}


=head2 version

  Arg [1]    : (optional) string $version - version to set
  Example    : $rnaproduct->version(2);
  Description: Getter/setter for attribute version
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub version {
  my $self = shift;
  $self->{'version'} = shift if (@_);
  return $self->{'version'};
}

1;
