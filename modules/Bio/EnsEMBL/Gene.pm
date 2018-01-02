=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Gene - Object representing a genes

=head1 SYNOPSIS

  my $gene = Bio::EnsEMBL::Gene->new(
    -START  => 123,
    -END    => 1045,
    -STRAND => 1,
    -SLICE  => $slice
  );

  # print gene information
  print("gene start:end:strand is "
      . join( ":", map { $gene->$_ } qw(start end strand) )
      . "\n" );

  # set some additional attributes
  $gene->stable_id('ENSG000001');
  $gene->description('This is the gene description');

=head1 DESCRIPTION

A representation of a Gene within the Ensembl system. A gene is a set of one or
more alternative transcripts.

=head1 METHODS

=cut

package Bio::EnsEMBL::Gene;

use strict;

use POSIX;
use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Intron;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Feature);


=head2 new

  Arg [-START]  : 
       int - start postion of the gene
  Arg [-END]    : 
       int - end position of the gene
  Arg [-STRAND] : 
       int - 1,-1 tehe strand the gene is on
  Arg [-SLICE]  : 
       Bio::EnsEMBL::Slice - the slice the gene is on
  Arg [-STABLE_ID] :
        string - the stable identifier of this gene
  Arg [-VERSION] :
        int - the version of the stable identifier of this gene
  Arg [-EXTERNAL_NAME] :
        string - the external database name associated with this gene
  Arg [-EXTERNAL_DB] :
        string - the name of the database the external name is from
  Arg [-EXTERNAL_STATUS]:
        string - the status of the external identifier
  Arg [-DISPLAY_XREF]:
        Bio::EnsEMBL::DBEntry - The external database entry that is used
        to label this gene when it is displayed.
  Arg [-TRANSCRIPTS]:
        Listref of Bio::EnsEMBL::Transcripts - this gene's transcripts
  Arg [-CREATED_DATE]:
        string - the date the gene was created
  Arg [-MODIFIED_DATE]:
        string - the date the gene was last modified
  Arg [-DESCRIPTION]:
        string - the genes description
  Arg [-BIOTYPE]:
        string - the biotype e.g. "protein_coding"
  Arg [-SOURCE]:
        string - the genes source, e.g. "ensembl"
  Arg [-IS_CURRENT]:
        Boolean - specifies if this is the current version of the gene
  Arg [-CANONICAL_TRANSCRIPT]:
        Bio::EnsEMBL::Transcript - the canonical transcript of this gene
  Arg [-CANONICAL_TRANSCRIPT_ID]:
        integer - the canonical transcript dbID of this gene, if the
        transcript object itself is not available.

  Example    : $gene = Bio::EnsEMBL::Gene->new(...);
  Description: Creates a new gene object
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);
  my (
    $stable_id,               $version,
    $external_name,           $type,
    $external_db,             $external_status,
    $display_xref,            $description,
    $transcripts,             $created_date,
    $modified_date,           $confidence,
    $biotype,                 $source,
    $is_current,
    $canonical_transcript_id, $canonical_transcript
    )
    = rearrange( [
      'STABLE_ID',               'VERSION',
      'EXTERNAL_NAME',           'TYPE',
      'EXTERNAL_DB',             'EXTERNAL_STATUS',
      'DISPLAY_XREF',            'DESCRIPTION',
      'TRANSCRIPTS',             'CREATED_DATE',
      'MODIFIED_DATE',           'CONFIDENCE',
      'BIOTYPE',                 'SOURCE',
      'IS_CURRENT',
      'CANONICAL_TRANSCRIPT_ID', 'CANONICAL_TRANSCRIPT'
    ],
    @_
    );


  if ($transcripts) {
    $self->{'_transcript_array'} = $transcripts;
    $self->recalculate_coordinates();
  }

  $self->stable_id($stable_id);
  $self->{'created_date'}  = $created_date;
  $self->{'modified_date'} = $modified_date;

  $self->external_name($external_name) if ( defined $external_name );
  $self->external_db($external_db)     if ( defined $external_db );
  $self->external_status($external_status)
    if ( defined $external_status );
  $self->display_xref($display_xref) if ( defined $display_xref );
  $self->biotype($type)              if ( defined $type );
  $self->biotype($biotype)           if ( defined $biotype );
  $self->description($description);
  $self->source($source);

  # Default version
  if ( !defined($version) ) { $version = 1 }
  $self->{'version'} = $version;

  # default to is_current
  $is_current = 1 unless (defined($is_current));
  $self->{'is_current'} = $is_current;

  # Add the canonical transcript if we were given one, otherwise add the
  # canonical transcript internal ID if we were given one.
  if ( defined($canonical_transcript) ) {
    $self->canonical_transcript($canonical_transcript);
  } elsif ( defined($canonical_transcript_id) ) {
    $self->{'canonical_transcript_id'} = $canonical_transcript_id;
  }

  return $self;
}



=head2 external_name

  Arg [1]    : (optional) String - the external name to set
  Example    : $gene->external_name('BRCA2');
  Description: Getter/setter for attribute external_name.
  Returntype : String or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub external_name {
  my  $self  = shift;

  $self->{'external_name'} = shift if (@_);

  if (defined $self->{'external_name'}) {
    return $self->{'external_name'};
  }

  my $display_xref = $self->display_xref();

  if (defined $display_xref) {
    return $display_xref->display_id();
  } else {
    return undef;
  }
}



=head2 source

  Arg [1]    : (optional) String - the source to set
  Example    : $gene->source('ensembl');
  Description: Getter/setter for attribute source
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub source {
  my $self = shift;
  $self->{'source'} = shift if( @_ );
  return ( $self->{'source'} || "ensembl" );
}


=head2 external_db	

  Arg [1]    : (optional) String - name of external db to set
  Example    : $gene->external_db('HGNC');
  Description: Getter/setter for attribute external_db. The db is the one that 
               belongs to the external_name.  
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

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

  Arg [1]    : (optional) String - status of the external db
  Example    : $gene->external_status('KNOWNXREF');
  Description: Getter/setter for attribute external_status. The status of
               the external db of the one that belongs to the external_name.
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

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

  Arg [1]    : (optional) String - the description to set
  Example    : $gene->description('This is the gene\'s description');
  Description: Getter/setter for gene description
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub description {
    my $self = shift;
    $self->{'description'} = shift if( @_ );
    return $self->{'description'};
}


=head2 equals

  Arg [1]       : Bio::EnsEMBL::Gene gene
  Example       : if ($geneA->equals($geneB)) { ... }
  Description   : Compares two genes for equality.
                  The test for eqality goes through the following list
                  and terminates at the first true match:

                  1. If Bio::EnsEMBL::Feature::equals() returns false,
                     then the genes are *not* equal.
                  2. If the biotypes differ, then the genes are *not*
                     equal.
                  3. If both genes have stable IDs: if these are the
                     same, the genes are equal, otherwise not.
                  4. If both genes have the same number of transcripts
                     and if these are (when compared pair-wise sorted by
                     start-position and length) the same, then they are
                     equal, otherwise not.

  Return type   : Boolean (0, 1)

  Exceptions    : Thrown if a non-gene is passed as the argument.

=cut

sub equals {
  my ( $self, $gene ) = @_;

  if ( !defined($gene) ) { return 0 }
  if ( $self eq $gene ) { return 1 }

  assert_ref( $gene, 'Bio::EnsEMBL::Gene' );

  my $feature_equals = $self->SUPER::equals($gene);
  if ( defined($feature_equals) && $feature_equals == 0 ) {
    return 0;
  }

  if ( $self->biotype() ne $gene->biotype() ) {
    return 0;
  }

  if ( defined( $self->stable_id() ) && defined( $gene->stable_id() ) )
  {
    if   ( $self->stable_id() eq $gene->stable_id() ) { return 1 }
    else                                              { return 0 }
  }

  my @self_transcripts = sort {
    $a->start() <=> $b->start() ||
      $a->length() <=> $b->length()
  } @{ $self->get_all_Transcripts() };
  my @gene_transcripts = sort {
    $a->start() <=> $b->start() ||
      $a->length() <=> $b->length()
  } @{ $gene->get_all_Transcripts() };

  if ( scalar(@self_transcripts) != scalar(@gene_transcripts) ) {
    return 0;
  }

  while (@self_transcripts) {
    my $self_transcript = shift(@self_transcripts);
    my $gene_transcript = shift(@gene_transcripts);

    if ( !$self_transcript->equals($gene_transcript) ) {
      return 0;
    }
  }

  return 1;
} ## end sub equals

=head2 canonical_transcript

  Arg [1]    : (optional) Bio::EnsEMBL::Transcript - canonical_transcript object
  Example    : $gene->canonical_transcript($canonical_transcript);
  Description: Getter/setter for the canonical_transcript
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : Throws if argument is not a transcript object.
  Caller     : general
  Status     : Stable

=cut

sub canonical_transcript {
  my ( $self, $transcript ) = @_;

  if ( defined($transcript) ) {
    # We're attaching a new canonical transcript.

    assert_ref( $transcript, 'Bio::EnsEMBL::Transcript' );

    # If there's already a canonical transcript, make sure it doesn't
    # think it's still canonical.
    if ( defined( $self->{'canonical_transcript'} ) ) {
      $self->{'canonical_transcript'}->is_canonical(0);
    }

    $self->{'canonical_transcript'}    = $transcript;
    $self->{'canonical_transcript_id'} = $transcript->dbID();

    $transcript->is_canonical(1);

  } elsif (   !defined( $self->{'canonical_transcript'} )
            && defined( $self->{'canonical_transcript_id'} )
            && $self->{'canonical_transcript_id'} != 0 )
  {
    # We have not attached a canoncical transcript, but we have the dbID
    # of one.

    if ( defined( $self->adaptor() ) ) {
      my $transcript_adaptor =
        $self->adaptor()->db()->get_TranscriptAdaptor();

      my $canonical_transcript =
        $transcript_adaptor->fetch_by_dbID(
                                   $self->{'canonical_transcript_id'} );

      if ( defined($canonical_transcript) ) {
        # Recusive call...
        $self->canonical_transcript($canonical_transcript);
      }

    } else {
      warning(   "Gene has no adaptor "
               . "when trying to fetch canonical transcript." );
    }

  } ## end elsif ( !defined( $self->...))

  return $self->{'canonical_transcript'};
} ## end sub canonical_transcript


=head2 get_all_Attributes

  Arg [1]    : (optional) String $attrib_code
               The code of the attribute type to retrieve values for
  Example    : my ($author) = @{ $gene->get_all_Attributes('author') };
               my @gene_attributes = @{ $gene->get_all_Attributes };
  Description: Gets a list of Attributes of this gene.
               Optionally just get Attributes for given code.
  Returntype : Listref of Bio::EnsEMBL::Attribute
  Exceptions : warning if gene does not have attached adaptor and attempts lazy
               load.
  Caller     : general
  Status     : Stable

=cut

sub get_all_Attributes {
  my $self = shift;
  my $attrib_code = shift;

  if ( ! exists $self->{'attributes' } ) {
    if (!$self->adaptor() ) {
      return [];
    }

    my $attribute_adaptor = $self->adaptor->db->get_AttributeAdaptor();
    $self->{'attributes'} = $attribute_adaptor->fetch_all_by_Gene($self);
  }

  if ( defined $attrib_code ) {
    my @results = grep { uc($_->code()) eq uc($attrib_code) }
    @{$self->{'attributes'}};
    return \@results;
  } else {
    return $self->{'attributes'};
  }
}


=head2 add_Attributes

  Arg [1-N]  : list of Bio::EnsEMBL::Attribute's @attribs
               Attribute(s) to add
  Example    : my $attrib = Bio::EnsEMBL::Attribute->new(...);
               $gene->add_Attributes($attrib);
  Description: Adds an Attribute to the Gene. If you add an attribute before
               you retrieve any from database, lazy loading will be disabled.
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
     throw( "Argument to add_Attribute has to be an Bio::EnsEMBL::Attribute" );
    }
    push( @{$self->{'attributes'}}, $attrib );
  }

  return;
}


=head2 add_DBEntry

  Arg [1]    : Bio::EnsEMBL::DBEntry $dbe
               The dbEntry to be added
  Example    : my $dbe = Bio::EnsEMBL::DBEntery->new(...);
               $gene->add_DBEntry($dbe);
  Description: Associates a DBEntry with this gene. Note that adding DBEntries
               will prevent future lazy-loading of DBEntries for this gene
               (see get_all_DBEntries).
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


=head2 get_all_DBEntries

  Arg [1]    : (optional) String, external database name,
               SQL wildcard characters (_ and %) can be used to
               specify patterns.

  Arg [2]    : (optional) String, external_db type, can be one of
               ('ARRAY','ALT_TRANS','ALT_GENE','MISC','LIT','PRIMARY_DB_SYNONYM','ENSEMBL'),
               SQL wildcard characters (_ and %) can be used to
               specify patterns.

  Example    : my @dbentries = @{ $gene->get_all_DBEntries() };
               @dbentries = @{ $gene->get_all_DBEntries('Uniprot%') };
               @dbentries = @{ $gene->get_all_DBEntries('%', 'ENSEMBL') };}

  Description: Retrieves DBEntries (xrefs) for this gene.  This does
               *not* include DBEntries that are associated with the
               transcripts and corresponding translations of this
               gene (see get_all_DBLinks()).

               This method will attempt to lazy-load DBEntries
               from a database if an adaptor is available and no
               DBEntries are present on the gene (i.e. they have not
               already been added or loaded).

  Return type: Listref of Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : get_all_DBLinks, GeneAdaptor::store
  Status     : Stable

=cut

sub get_all_DBEntries {
  my ( $self, $db_name_exp, $ex_db_type ) = @_;

  my $cache_name = 'dbentries';

  if ( defined($db_name_exp) ) {
    $cache_name .= $db_name_exp;
  }

  if ( defined($ex_db_type) ) {
    $cache_name .= $ex_db_type;
  }

  # if not cached, retrieve all of the xrefs for this gene
  if ( !defined( $self->{$cache_name} ) && defined( $self->adaptor() ) )
  {
    $self->{$cache_name} =
      $self->adaptor()->db()->get_DBEntryAdaptor()
      ->fetch_all_by_Gene( $self, $db_name_exp, $ex_db_type );
  }

  $self->{$cache_name} ||= [];

  return $self->{$cache_name};
} ## end sub get_all_DBEntries

=head2 get_all_object_xrefs

  Arg [1]    : (optional) String, external database name

  Arg [2]    : (optional) String, external_db type

  Example    : @oxrefs = @{ $gene->get_all_object_xrefs() };

  Description: Retrieves xrefs for this gene.  This does *not*
               include xrefs that are associated with the
               transcripts or corresponding translations of this
               gene (see get_all_xrefs()).

               This method will attempt to lazy-load xrefs from a
               database if an adaptor is available and no xrefs are
               present on the gene (i.e. they have not already been
               added or loaded).

                NB: This method is an alias for the
                    get_all_DBentries() method.

  Return type: Listref of Bio::EnsEMBL::DBEntry objects

  Status     : Stable

=cut

sub get_all_object_xrefs {
  my $self = shift;
  return $self->get_all_DBEntries(@_);
}

=head2 get_all_DBLinks

  Arg [1]    : String database name (optional)
               SQL wildcard characters (_ and %) can be used to
               specify patterns.

  Arg [2]    : (optional) String, external database type, can be one of
               ('ARRAY','ALT_TRANS','ALT_GENE','MISC','LIT','PRIMARY_DB_SYNONYM','ENSEMBL'),
               SQL wildcard characters (_ and %) can be used to
               specify patterns.

  Example    : @dblinks = @{ $gene->get_all_DBLinks() };
               @dblinks = @{ $gene->get_all_DBLinks('Uniprot%') };
               @dblinks = @{ $gene->get_all_DBLinks('%', 'ENSEMBL') };}

  Description: Retrieves *all* related DBEntries for this gene. This
               includes all DBEntries that are associated with the
               transcripts and corresponding translations of this
               gene.

               If you only want to retrieve the DBEntries
               associated with the gene (and not the transcript
               and translations) then you should use the
               get_all_DBEntries() call instead.

               Note: Each entry may be listed more than once.  No
               uniqueness checks are done.  Also if you put in an
               incorrect external database name no checks are done
               to see if this exists, you will just get an empty
               list.

  Return type: Listref of Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_DBLinks {
  my ( $self, $db_name_exp, $ex_db_type ) = @_;

  my @links =
    @{ $self->get_all_DBEntries( $db_name_exp, $ex_db_type ) };

  # Add all of the transcript and translation xrefs to the return list.
  foreach my $transcript ( @{ $self->get_all_Transcripts() } ) {
    push( @links,
          @{$transcript->get_all_DBLinks( $db_name_exp, $ex_db_type ) }
    );
  }

  return \@links;
}

=head2 get_all_xrefs

  Arg [1]    : String database name (optional)
               SQL wildcard characters (_ and %) can be used to
               specify patterns.

  Example    : @xrefs = @{ $gene->get_all_xrefs() };
               @xrefs = @{ $gene->get_all_xrefs('Uniprot%') };

  Description: Retrieves *all* related xrefs for this gene.  This
               includes all xrefs that are associated with the
               transcripts and corresponding translations of this
               gene.

               If you want to retrieve the xrefs associated
               with only the gene (and not the transcript
               or translations) then you should use the
               get_all_object_xrefs() method instead.

               Note: Each entry may be listed more than once.  No
               uniqueness checks are done.  Also if you put in an
               incorrect external database name no checks are done
               to see if this exists, you will just get an empty
               list.

                NB: This method is an alias for the
                    get_all_DBLinks() method.

  Return type: Listref of Bio::EnsEMBL::DBEntry objects

  Status     : Stable

=cut

sub get_all_xrefs {
  my $self = shift;
  return $self->get_all_DBLinks(@_);
}

=head2 get_all_Exons

  Example    : my @exons = @{ $gene->get_all_Exons };
  Description: Returns a set of all the exons associated with this gene.
  Returntype : Listref of Bio::EnsEMBL::Exon objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub get_all_Exons {
  my $self = shift;

  my %h;
  my @out = ();

  foreach my $trans ( @{$self->get_all_Transcripts} ) {
    foreach my $e ( @{$trans->get_all_Exons} ) {
      $h{$e->start()."-".$e->end()."-".$e->strand()."-".$e->phase()."-".$e->end_phase()} = $e;
    }
  }

  push @out, values %h;

  return \@out;
}

=head2 get_all_Introns

  Arg [1]    : none
  Example    : my @introns = @{$gene->get_all_Introns()};
  Description: Returns an listref of the introns in this gene in order.
               i.e. the first intron in the listref is the 5prime most exon in
               the gene.
  Returntype : listref to Bio::EnsEMBL::Intron objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Introns {
  my $self = shift;

  my %h;
  my @out = ();
  my @introns;

  foreach my $trans ( @{$self->get_all_Transcripts} ) {
    my @exons = @{ $trans->get_all_Exons() };
    for (my $i = 0; $i < scalar(@exons) - 1; $i++) {
      my $intron = new Bio::EnsEMBL::Intron($exons[$i], $exons[$i+1]);
      push (@introns, $intron);
    }
  }

  return \@introns;
}


=head2 get_all_homologous_Genes

  Arg[1]     : String The compara synonym to use when looking for a database in the
               registry. If not provided we will use the very first compara database
               we find.
  Description: Queries the Ensembl Compara database and retrieves all
               Genes from other species that are orthologous.
               REQUIRES properly setup Registry conf file. Meaning that
               one of the aliases for each core db has to be "Genus species"
               e.g. "Homo sapiens" (as in the name column in genome_db table
               in the compara database).

               The data is cached in this Object for faster re-retreival.
  Returntype : listref [
                        Bio::EnsEMBL::Gene,
                        Bio::EnsEMBL::Compara::Homology,
                        string $species # needed as cannot get spp from Gene 
                       ]
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_homologous_Genes {
  my ($self, $db_synonym) = @_;

  #Look for DBAdaptors which have a group of compara; these are compara DBAs. 
  #If given a synonym 
  my %args = (-GROUP => 'compara');
  $args{-SPECIES} = $db_synonym if $db_synonym;
  my ($compara_dba) = @{Bio::EnsEMBL::Registry->get_all_DBAdaptors(%args)};
  unless( $compara_dba ) {
    throw("No compara found in Bio::EnsEMBL::Registry. Please fully populate the Registry or construct a Bio::EnsEMBL::Compara::DBSQL::DBAdaptor");
  }
  my $compara_species = $compara_dba->species();
  if( exists( $self->{'homologues'}->{$compara_species} ) ){
    return $self->{'homologues'}->{$compara_species};
  }
  $self->{'homologues'}->{$compara_species} = [];

  # Get the compara 'member' corresponding to self
  my $member_adaptor   = $compara_dba->get_adaptor('GeneMember');
  my $query_member = $member_adaptor->fetch_by_stable_id($self->stable_id);
  unless( $query_member ){ return $self->{'homologues'}->{$compara_species} };

  # Get the compara 'homologies' corresponding to 'member'
  my $homology_adaptor = $compara_dba->get_adaptor('Homology');
  my @homolos = @{$homology_adaptor->fetch_all_by_Member($query_member)};
  unless( scalar(@homolos) ){ return $self->{'homologues'}->{$compara_species} };

  # Get the ensembl 'genes' corresponding to 'homologies'
  foreach my $homolo( @homolos ){
    foreach my $member( @{$homolo->get_all_GeneMembers} ){
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
      push @{$self->{'homologues'}->{$compara_species}}, [$hgene,$homolo,$hspecies];
    }
  }
  return $self->{'homologues'}->{$compara_species};
}

=head2 _clear_homologues

  Description: Removes any cached homologues from the Gene which could have been
               fetched from the C<get_all_homologous_Genes()> call.
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub _clear_homologues {
  my ($self) = @_;
  delete $self->{homologues};
}


=head2 biotype

  Arg [1]    : (optional) String - the biotype to set
  Example    : $gene->biotype("protein_coding");
  Description: Getter/setter for the attribute biotype
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub biotype {
  my $self = shift;

  $self->{'biotype'} = shift if( @_ );
  return ( $self->{'biotype'} || "protein_coding" );
}


=head2 add_Transcript

  Arg [1]    : Bio::EnsEMBL::Transcript $trans
               The transcript to add to the gene
  Example    : my $transcript = Bio::EnsEMBL::Transcript->new(...);
               $gene->add_Transcript($transcript);
  Description: Adds another Transcript to the set of alternatively
               spliced Transcripts of this gene. If it shares exons 
               with another Transcript, these should be object-identical.
  Returntype : none
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub add_Transcript {
   my ($self, $trans) = @_;

   if( !ref $trans || ! $trans->isa("Bio::EnsEMBL::Transcript") ) {
       throw("$trans is not a Bio::EnsEMBL::Transcript!");
   }

   $self->{'_transcript_array'} ||= [];
   push(@{$self->{'_transcript_array'}},$trans);

   $self->recalculate_coordinates();
}

sub remove_Transcript {
  my ($self,$trans) = @_;
  if( !ref $trans || ! $trans->isa("Bio::EnsEMBL::Transcript") ) {
       throw("$trans is not a Bio::EnsEMBL::Transcript!");
  }
  # Clean transcript from live data
  $self->get_all_Transcripts; # force lazy load.
  my $array = $self->{_transcript_array};
  my $db_id = $trans->dbID;
  @$array = grep { $_->dbID != $db_id } @$array;
  # Recalculate and store new gene coordinates
  $self->adaptor->update_coords($self);
}


=head2 get_all_Transcripts

  Example    : my @transcripts = @{ $gene->get_all_Transcripts };
  Description: Returns the Transcripts in this gene.
  Returntype : Listref of Bio::EnsEMBL::Transcript objects
  Warning    : This method returns the internal transcript array 
               used by this object. Avoid any modification
               of this array. We class use of shift and 
               reassignment of the loop variable when iterating
               this array as modification.

               Dereferencing the structure as shown in the example is
               a safe way of using this data structure.
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Transcripts {
  my $self = shift;

  if( ! exists $self->{'_transcript_array'} ) {
    if( defined $self->adaptor() ) {
      my $ta = $self->adaptor()->db()->get_TranscriptAdaptor();
      my $transcripts = $ta->fetch_all_by_Gene( $self );
      $self->{'_transcript_array'} = $transcripts;
    }
  }
  my @array_copy;
  if (defined $self->{'_transcript_array'}) {
    @array_copy = @{ $self->{'_transcript_array'} } ;
    return \@array_copy;
  }
  return;
}


=head2 get_all_alt_alleles

  Example    : my @alt_genes = @{ $gene->get_all_alt_alleles };
               foreach my $alt_gene (@alt_genes) {
                 print "Alternate allele: " . $alt_gene->stable_id() . "\n";
               }
  Description: Returns a listref of Gene objects that represent this Gene on
               an alternative haplotype. Empty list if there is no such
               Gene (eg there is no overlapping haplotype).
  Returntype : listref of Bio::EnsEMBL::Gene objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_alt_alleles {
  my $self = shift;
  my $result = $self->adaptor()->fetch_all_alt_alleles( $self );
  return $result;
}


=head2 version

  Arg [1]    : (optional) Int
               A version number for the stable_id
  Example    : $gene->version(2);
  Description: Getter/setter for version number
  Returntype : Int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub version {
  my $self = shift;
  $self->{'version'} = shift if(@_);
  return $self->{'version'};
}


=head2 stable_id

  Arg [1]    : (optional) String - the stable ID to set
  Example    : $gene->stable_id("ENSG0000000001");
  Description: Getter/setter for stable id for this gene.
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub stable_id {
  my $self = shift;
  $self->{'stable_id'} = shift if(@_);
  return $self->{'stable_id'};
}

=head2 stable_id_version

  Arg [1]    : (optional) String - the stable ID with version to set
  Example    : $gene->stable_id("ENSG0000000001.3");
  Description: Getter/setter for stable id with version for this gene.
  Returntype : String
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub stable_id_version {
    my $self = shift;
    if(my $stable_id = shift) {
	# See if there's an embedded period, assume that's a
	# version, might not work for some species but you
	# should use ->stable_id() and version() if you're worried
	# about ambiguity
	my $vindex = rindex($stable_id, '.');
	# Set the stable_id and version pair depending on if
	# we found a version delimiter in the stable_id
	($self->{stable_id}, $self->{version}) = ($vindex > 0 ?
						  (substr($stable_id,0,$vindex), substr($stable_id,$vindex+1)) :
						  $stable_id, undef);
    }
    return $self->{stable_id} . ($self->{version} ? ".$self->{version}" : '');
}

=head2 is_current

  Arg [1]    : Boolean $is_current
  Example    : $gene->is_current(1)
  Description: Getter/setter for is_current state of this gene.
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

  Arg [1]    : (optional) String - created date to set (as a UNIX time int)
  Example    : $gene->created_date('1141948800');
  Description: Getter/setter for attribute created_date
  Returntype : String
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

  Arg [1]    : (optional) String - modified date to set (as a UNIX time int)
  Example    : $gene->modified_date('1141948800');
  Description: Getter/setter for attribute modified_date
  Returntype : String
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

  Arg [1]    : String - coordinate system name to transform to
  Arg [2]    : String - coordinate system version
  Example    : my $new_gene = $gene->transform('supercontig');
  Description: Moves this gene to the given coordinate system. If this gene has
               Transcripts attached, they move as well.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : throw on wrong parameters
  Caller     : general
  Status     : Stable

=cut

sub transform {
  my $self = shift;

  # catch for old style transform calls
  if( !@_  || ( ref $_[0] && ($_[0]->isa( "Bio::EnsEMBL::Slice" ) or $_[0]->isa( "Bio::EnsEMBL::LRGSlice" )) )) {
    deprecate('Calling transform without a coord system name is deprecated.');
    return $self->_deprecated_transform(@_);
  }

  my $new_gene = $self->SUPER::transform(@_);

  if ( !defined($new_gene) ) {
    # check if this gene projects at all to requested coord system,
    #  if not we are done.
    my @segments = @{ $self->project(@_) };
    if ( !@segments ) {
      return undef;
    }
  }

  #
  # If you are transforming the gene then make sure the transcripts and exons are loaded
  #

  foreach my $tran (@{$self->get_all_Transcripts}){
    $tran->get_all_Exons();
  }

  if( exists $self->{'_transcript_array'} ) {
    my @new_transcripts;
    my ( $strand, $slice );
    my $low_start = POSIX::INT_MAX;
    my $hi_end = POSIX::INT_MIN;
    for my $old_transcript ( @{$self->{'_transcript_array'}} ) {
      my $new_transcript = $old_transcript->transform( @_ );
      # this can fail if gene transform failed  
      
      return undef unless $new_transcript;

      if( ! defined $new_gene ) {
	if( $new_transcript->start() < $low_start ) {
	  $low_start = $new_transcript->start();
	}
	if( $new_transcript->end() > $hi_end ) {
	  $hi_end = $new_transcript->end();
	}
	$slice = $new_transcript->slice();
	$strand = $new_transcript->strand();
      }
      push( @new_transcripts, $new_transcript );
    }

    if( ! defined $new_gene ) {
      %$new_gene = %$self;
      bless $new_gene, ref( $self );

      $new_gene->start( $low_start );
      $new_gene->end( $hi_end );
      $new_gene->strand( $strand );
      $new_gene->slice( $slice );
    }

    $new_gene->{'_transcript_array'} = \@new_transcripts;
  }

  if(exists $self->{attributes}) {
    $new_gene->{attributes} = [@{$self->{attributes}}];
  }

  return $new_gene;
}


=head2 transfer

  Arg [1]    : Bio::EnsEMBL::Slice $destination_slice
  Example    : my $new_gene = $gene->transfer($slice);
  Description: Moves this Gene to given target slice coordinates. If Transcripts
               are attached they are moved as well. Returns a new gene.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general
  Status     : Stable

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

  if(exists $self->{attributes}) {
    $new_gene->{attributes} = [@{$self->{attributes}}];
  }
  
  return $new_gene;
}


=head2 display_xref

  Arg [1]    : (optional) Bio::EnsEMBL::DBEntry - the display xref to set
  Example    : $gene->display_xref($db_entry);
  Description: Getter/setter display_xref for this gene.
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


=head2 display_id

  Example    : print $gene->display_id();
  Description: This method returns a string that is considered to be
               the 'display' identifier. For genes this is (depending on
               availability and in this order) the stable Id, the dbID or an
               empty string.
  Returntype : String
  Exceptions : none
  Caller     : web drawing code
  Status     : Stable

=cut

sub display_id {
  my $self = shift;
  return $self->{'stable_id'} || $self->dbID || '';
}


=head2 recalculate_coordinates

  Example    : $gene->recalculate_coordinates;
  Description: Called when transcript added to the gene, tries to adapt the
               coords for the gene.
  Returntype : none
  Exceptions : none
  Caller     : internal
  Status     : Stable

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


=head2 get_all_DASFactories

  Example    : $dasref = $prot->get_all_DASFactories
  Description: Retrieves a listref of registered DAS objects
              TODO: Abstract to a DBLinkContainer obj
  Returntype : [ DAS_objects ]
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_DASFactories {
   my $self = shift;
   return [ $self->adaptor()->db()->_each_DASFeatureFactory ];
}


=head2 get_all_DAS_Features

  Example    : $features = $prot->get_all_DAS_Features;
  Description: Retrieves a hash reference to a hash of DAS feature
               sets, keyed by the DNS, NOTE the values of this hash
               are an anonymous array containing:
                (1) a pointer to an array of features
                (2) a pointer to the DAS stylesheet
  Returntype : hashref of Bio::SeqFeatures
  Exceptions : none
  Caller     : webcode
  Status     : Stable

=cut

sub get_all_DAS_Features{
  my ($self, @args) = @_;
  my $slice = $self->feature_Slice;
  return $self->SUPER::get_all_DAS_Features($slice);
}


=head2 load

  Arg [1]       : Boolean $load_xrefs
                  Load (or don't load) xrefs.  Default is to load xrefs.
  Example       : $gene->load();
  Description   : The Ensembl API makes extensive use of
                  lazy-loading.  Under some circumstances (e.g.,
                  when copying genes between databases), all data of
                  an object needs to be fully loaded.  This method
                  loads the parts of the object that are usually
                  lazy-loaded.  It will also call the equivalent
                  method on all the transcripts of the gene.
  Returns       : 

=cut

sub load {
  my ( $self, $load_xrefs ) = @_;

  if ( !defined($load_xrefs) ) { $load_xrefs = 1 }

  foreach my $transcript ( @{ $self->get_all_Transcripts() } ) {
    $transcript->load($load_xrefs);
  }

  $self->analysis();
  $self->get_all_Attributes();
  $self->stable_id();
  $self->canonical_transcript();

  if ($load_xrefs) {
    $self->get_all_DBEntries();
  }
}

=head2 flush_Transcripts

  Description : Empties out caches and unsets fields of this Gene.
                Beware of further actions without adding some new transcripts.
  Example     : $gene->flush_Transcripts();

=cut

sub flush_Transcripts {
    my $self = shift;
    $self->{'_transcript_array'} = [];
    $self->{'canonical_transcript_id'} = undef;
    $self->{'canonical_transcript'} = undef;
    return;
}

=head2 is_ref

  Description: getter setter for the gene attribute is_ref
  Arg [1]    : (optional) 1 or 0
  return     : boolean

=cut

sub is_reference{
  my ( $self, $is_ref) = @_;

  if(defined($is_ref)){
    $self->{'is_ref'} = $is_ref;
  }
  else{
    $self->{'is_ref'} = $self->adaptor->is_ref($self->dbID);
  }	
  return $self->{'is_ref'};
}

=head2 summary_as_hash

  Example       : $gene_summary = $gene->summary_as_hash();
  Description   : Extends Feature::summary_as_hash
                  Retrieves a summary of this Gene object.
	                  
  Returns       : hashref of arrays of descriptive strings
  Status        : Intended for internal use
=cut

sub summary_as_hash {
  my $self = shift;
  my $summary_ref = $self->SUPER::summary_as_hash;
  $summary_ref->{'description'} = $self->description;
  $summary_ref->{'biotype'} = $self->biotype;
  $summary_ref->{'Name'} = $self->external_name if $self->external_name;
  $summary_ref->{'logic_name'} = $self->analysis->logic_name() if defined $self->analysis();
  $summary_ref->{'source'} = $self->source();
  $summary_ref->{'gene_id'} = $summary_ref->{'id'};

  ## Will only work for for merged species
  my $havana_gene = $self->havana_gene();
  $summary_ref->{'havana_gene'} = $havana_gene->display_id() if defined $havana_gene;
  $summary_ref->{'havana_version'} = $havana_gene->version() if defined $havana_gene;

  ## Stable identifier of the parent gene this gene was projected from
   my $proj_parent_attributes = $self->get_all_Attributes("proj_parent_g");
    if (@{$proj_parent_attributes}) {
      $summary_ref->{'projection_parent_gene'} = $proj_parent_attributes->[0]->value;
    }
  return $summary_ref;
}

=head2 havana_gene

  Example       : $havana_gene = $transcript->havana_gene();
  Description   : Locates the corresponding havana gene
  Returns       : Bio::EnsEMBL::DBEntry
=cut

sub havana_gene {
  my $self = shift;
  my @otts = @{ $self->get_all_DBEntries('Vega_gene') };
  my $ott;
  foreach my $xref (@otts) {
    if ($xref->display_id() =~ /OTT/) {
      $ott = $xref;
      last;
    }
  }
  return $ott;
}




###########################
# DEPRECATED METHODS FOLLOW
###########################

=head2 fetch_coded_for_regulatory_factors

  Arg [1]    : none
  Example    : $gene->fetch_coded_for_regulatory_factors()
  Description: DEPRECATED: Fetches any regulatory_factors that are coded for by
               this gene.
  Returntype : Listref of Bio::Ensembl::RegulatoryFactor
  Exceptions :
  Caller     : ?
  Status     : At Risk
             : under development

=cut

sub fetch_coded_for_regulatory_factors {

  my ($self) = @_;

  my $rfa = $self->adaptor->db->get_RegulatoryFactorAdaptor();

  return $rfa->fetch_factors_coded_for_by_gene($self);

}


1;

