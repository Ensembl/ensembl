
#
# Ensembl module for Bio::EnsEMBL::DBSQL::SliceAdaptor
#
# Cared for by Ewan Birney <ensembl-dev@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::SlicetAdaptor - Adaptors for slices

=head1 SYNOPSIS



=head1 DESCRIPTION

Factory for getting out slices of assemblies.

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email ensembl-dev@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


package Bio::EnsEMBL::DBSQL::SliceAdaptor;
use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Mapper;

use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);

@ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');



sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  $self->{'_name_cache'} = {};
  $self->{'_id_cache'} = {};

  return $self;
}


=head2 fetch_by_region

  Arg [1]    : string $coord_system_name
               The name of the coordinate system of the slice to be created
               This may be a name of an acutal coordinate system or an alias
               to a coordinate system.  Valid aliases are 'seqlevel' or
               'toplevel'.
  Arg [2]    : string $seq_region_name
               The name of the sequence region that the slice will be
               created on
  Arg [3]    : int $start (optional, default = 1)
               The start of the slice on the sequence region
  Arg [4]    : int $end (optional, default = seq_region length)
               The end of the slice on the sequence region
  Arg [5]    : int $strand (optional, default = 1)
               The orientation of the slice on the sequence region
  Arg [6]    : string $version (optional, default = default version)
               The version of the coordinate system to use (e.g. NCBI33)
  Example    : $slice = $slice_adaptor->fetch_by_region('chromosome');
  Description: Creates a slice object on the given chromosome and coordinates.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_region {
  my ($self, $coord_system_name, $seq_region_name,
      $start, $end, $strand, $version) = @_;

  throw('seq_region_name argument is required') if(!$seq_region_name);

  my $csa = $self->db()->get_CoordSystemAdaptor();
  my $coord_system = $csa->fetch_by_name($coord_system_name,$version);

  #check the cache so we only go to the db if necessary
  my $name_cache = $self->{'_name_cache'};
  my $key = lc(join(':',$seq_region_name,
                    $coord_system->name(),
                    $coord_system->version));

  my $length;

  if(exists($name_cache->{$key})) {
    $length = $name_cache->{$key}->[1];
  } else {
    my $sth = $self->prepare("SELECT seq_region_id, length " .
                             "FROM seq_region " .
                             "WHERE name = ? AND coord_system_id = ?");

    #force seq_region_name cast to string so mysql cannot treat as int
    $sth->execute("$seq_region_name", $coord_system->dbID());

    if($sth->rows() != 1) {
      $version ||= '';
      $seq_region_name ||= '';
      $coord_system_name ||= '';
      throw("Cannot create slice on non-existant or ambigous seq_region:\n" .
            "  coord_system=[$coord_system_name],\n" .
            "  name=[$seq_region_name],\n" .
            "  version=[$version]");
    }

    my $id;
    ($id, $length) = $sth->fetchrow_array();
    $sth->finish();

    #cache results to speed up future queries
    $name_cache->{$key} = [$id,$length];
    $self->{'_id_cache'}->{$id} = [$seq_region_name, $length, $coord_system];
  }

  $start = 1 if(!defined($start));
  $strand = 1 if(!defined($strand));
  $end = $length if(!defined($end));

  if($end < $start) {
    throw('start [$start] must be less than or equal to end [$end]');
  }

  return Bio::EnsEMBL::Slice->new(-COORD_SYSTEM    => $coord_system,
                                  -SEQ_REGION_NAME => $seq_region_name,
                                  -START           => $start,
                                  -END             => $end,
                                  -STRAND          => $strand,
                                  -ADAPTOR         => $self);
}




=head2 fetch_by_seq_region_id

  Arg [1]    : string $seq_region_id
               The internal identifier of the seq_region to create this slice
               on
  Example    : $slice = $slice_adaptor->fetch_by_seq_region_id(34413);
  Description: Creates a slice object of an entire seq_region using the
               seq_region internal identifier to resolve the seq_region.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_seq_region_id {
  my ($self, $seq_region_id) = @_;

  my $id_cache = $self->{'_id_cache'};

  my ($name, $length, $cs);

  if(exists $id_cache->{$seq_region_id}) {
    ($name, $length, $cs) = @{$id_cache->{$seq_region_id}};
  } else {
    my $sth = $self->prepare("SELECT name, length, coord_system_id " .
                             "FROM seq_region " .
                             "WHERE seq_region_id = ?");

    $sth->execute($seq_region_id);

    if($sth->rows() != 1) {
      throw("Cannot create slice on non-existant or ambigous seq_region:" .
            "  seq_region_id=[$seq_region_id],\n");
    }

    my $cs_id;
    ($name, $length, $cs_id) = $sth->fetchrow_array();
    $sth->finish();

    $cs = $self->db->get_CoordSystemAdaptor->fetch_by_dbID($cs_id);

    #cache results to speed up repeated queries
    $id_cache->{$seq_region_id} = [$name, $length, $cs];
    my $key = lc(join(':', $name, $cs->name, $cs->version));
    $self->{'_name_cache'}->{$key} = [$seq_region_id, $length];
  }

  return Bio::EnsEMBL::Slice->new(-COORD_SYSTEM    => $cs,
                                  -SEQ_REGION_NAME => $name,
                                  -START           => 1,
                                  -END             => $length,
                                  -STRAND          => 1,
                                  -ADAPTOR         => $self);
}



=head2 get_seq_region_id

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to fetch a seq_region_id for
  Example    : $srid = $slice_adaptor->get_seq_region_id($slice);
  Description: Retrieves the seq_region id (in this database) given a slice
               Seq region ids are not stored on the slices themselves
               because they are intended to be somewhat database independant
               and seq_region_ids vary accross databases.
  Returntype : int
  Exceptions : throw if the seq_region of the slice is not in the db
               throw if incorrect arg provided
  Caller     : BaseFeatureAdaptor

=cut

sub get_seq_region_id {
  my $self = shift;
  my $slice = shift;

  if(!$slice || !ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Slice argument is required');
  }

  my $cs_name = $slice->coord_system->name();
  my $cs_version = $slice->coord_system->version();
  my $seq_region_name = $slice->seq_region_name();

  my $key = lc(join(':', $seq_region_name,$cs_name,$cs_version));

  my $name_cache = $self->{'_name_cache'};

  if(exists($name_cache->{$key})) {
    return $name_cache->{$key}->[0];
  }

  my $csa = $self->db()->get_CoordSystemAdaptor();
  my $coord_system = $csa->fetch_by_name($cs_name,$cs_version);

  my $sth = $self->prepare("SELECT seq_region_id, length " .
                           "FROM seq_region " .
                           "WHERE name = ? AND coord_system_id = ?");

  #force seq_region_name cast to string so mysql cannot treat as int
  $sth->execute("$seq_region_name", $coord_system->dbID());

  if($sth->rows() != 1) {
    throw("Non existant or ambigous seq_region:\n" .
          "  coord_system=[$cs_name],\n" .
          "  name=[$seq_region_name],\n" .
          "  version=[$cs_version]");
  }

  my($seq_region_id, $length) = $sth->fetchrow_array();
  $sth->finish();

  #cache information for future requests
  $name_cache->{$key} = [$seq_region_id, $length];
  $self->{'_id_cache'}->{$seq_region_id} =
    [$seq_region_name, $length, $coord_system];

  return $seq_region_id;
}


=head2 get_seq_region_attribs

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to fetch attributes for
  Example    : %attribs = %{$slice_adaptor->get_seq_region_attribs($slice)};
               ($htg_phase) = @{$attribs->{'htg_phase'} || []};
               @synonyms    = @{$attribs->{'synonym'} || []};
  Description: Retrieves a reference to a hash containing attrib code values
               and listref value keys.
  Returntype : hashref
  Exceptions : throw if the seq_region of the slice is not in the db
               throw if incorrect arg provided
  Caller     : Bio::EnsEMBL::Slice

=cut

sub get_seq_region_attribs {
  my $self = shift;
  my $slice = shift;

  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Slice argument is required');
  }

  my $srid = $self->get_seq_region_id($slice);

  if(!$srid) {
    throw('Slice is not on a seq_region stored in this database.');
  }

  $self->{'_attribs_cache'} ||= {};
  if($self->{'_attribs_cache'}->{$srid}) {
    return $self->{'_attribs_cache'}->{$srid};
  }

  my $sth = $self->prepare('SELECT at.code, sra.value ' .
                           'FROM   seq_region_attrib sra, attrib_type at ' .
                           'WHERE  sra.seq_region_id = ? ' .
                           'AND    at.attrib_type_id = sra.attrib_type_id');

  $sth->execute($srid);

  my($code, $attrib);
  $sth->bind_columns(\$code, \$attrib);

  my %attrib_hash;
  while($sth->fetch()) {
    $attrib_hash{$code} ||= [];
    push @{$attrib_hash{$code}}, $attrib;
  }

  $sth->finish();
  $self->{'_attribs_cache'} = \%attrib_hash;
  return \%attrib_hash;
}


=head2 fetch_all

  Arg [1]    : string $coord_system_name
               The name of the coordinate system to retrieve slices of.
               This may be a name of an acutal coordinate system or an alias
               to a coordinate system.  Valid aliases are 'seqlevel' or
               'toplevel'.
  Arg [2]    : string $coord_system_version (optional)
               The version of the coordinate system to retrieve slices of
  Example    : @chromos = @{$slice_adaptor->fetch_all('chromosome','NCBI33')};
               @contigs = @{$slice_adaptor->fetch_all('contig')};
  Description: Retrieves slices of all seq_regions for a given coordinate
               system.  This is analagous to the methods fetch_all which were
               formerly on the ChromosomeAdaptor, RawContigAdaptor and
               CloneAdaptor classes.  Slices fetched span the entire
               seq_regions and are on the forward strand.
  Returntype : listref of Bio::EnsEMBL::Slices
  Exceptions : thrown if invalid coord system is provided
  Caller     : general

=cut

sub fetch_all {
  my $self = shift;
  my $cs_name = shift;
  my $cs_version = shift || '';

  my $csa = $self->db->get_CoordSystemAdaptor();
  my $cs = $csa->fetch_by_name($cs_name, $cs_version);

  my $sth = $self->prepare('SELECT seq_region_id, name, length ' .
                        'FROM   seq_region ' .
                        'WHERE  coord_system_id =?');

  $sth->execute($cs->dbID);

  my ($seq_region_id, $name, $length);
  $sth->bind_columns(\$seq_region_id, \$name, \$length);

  my $name_cache = $self->{'_name_cache'};
  my $id_cache   = $self->{'_id_cache'};

  my $cs_key = lc($cs->name().':'.$cs_version);

  my @out;
  while($sth->fetch()) {
    #cache values for future reference
    my $key = lc($name) . ':'. $cs_key;
    $name_cache->{$key} = [$seq_region_id, $length];
    $id_cache->{$seq_region_id} = [$name, $length, $cs];
    push @out, Bio::EnsEMBL::Slice->new(-START  => 1,
                                        -END    => $length,
                                        -STRAND => 1,
                                        -SEQ_REGION_NAME => $name,
                                        -COORD_SYSTEM => $cs,
                                        -ADAPTOR => $self);
  }

  return \@out;
}


=head2 fetch_all_non_redundant

  Arg [1]    : none
  Example    : @all = @{$slice_adaptor->fetch_all_non_redundant()};
  Description: Retrieves all non-redundant slices, i.e. those which have
               the attribute 'nonredundant' set
  Returntype : listref of Bio::EnsEMBL::Slices
  Exceptions : none
  Caller     : general

=cut

sub fetch_all_non_redundant {

  my $self = shift;

  my $sth = $self->prepare("SELECT s.name, s.length, c.coord_system_id " .
			   "FROM seq_region s, coord_system c, seq_region_attrib sra, attrib_type at " .
			   "WHERE s.coord_system_id=c.coord_system_id " .
			   "AND at.code='nonredundant' " .
			   "AND at.attrib_type_id=sra.attrib_type_id " .
			   "AND sra.seq_region_id=s.seq_region_id");

  $sth->execute();

  my ($name, $length, $cs_id);
  $sth->bind_columns(\$name, \$length, \$cs_id);

  # Slice expects a CoordSystem object
  my $cs_adaptor = $self->db->get_CoordSystemAdaptor();


  my @out;
  while($sth->fetch()) {
    my $cs = $cs_adaptor->fetch_by_dbID($cs_id);
    push @out, Bio::EnsEMBL::Slice->new(-START  => 1,
                                        -END    => $length,
                                        -STRAND => 1,
                                        -SEQ_REGION_NAME => $name,
                                        -COORD_SYSTEM => $cs,
                                        -ADAPTOR => $self);
  }

  return \@out;

}


sub deleteObj {
  my $self = shift;

  $self->SUPER::deleteObj;

  $self->{'_id_cache'} = {};
  $self->{'_name_cache'} = {};
  $self->{'_exc_cache'} = {};
}



=head2 fetch_by_chr_band

 Title   : fetch_by_chr_band
 Usage   :
 Function: create a Slice representing a series of bands
 Example :
 Returns :
 Args    : the band name

=cut

sub fetch_by_chr_band {
  my ($self,$chr,$band) = @_;

  my $chr_slice = $self->fetch_by_region('chromosome', $chr);

  my $seq_region_id = $self->get_seq_region_id($chr_slice);

  my $sth = $self->db->prepare
        ("select min(k.chr_start), max(k.chr_end) " .
         "from karyotype as k " .
         "where k.seq_region_id = ? and k.band like ?");

  $sth->execute( $seq_region_id, "$band%" );
  my ( $slice_start, $slice_end) = $sth->fetchrow_array;

  if(defined $slice_start) {
    return $self->fetch_by_region('chromosome',$chr,$slice_start,$slice_end);
  }

  throw("Band not recognised in database");
}



=head2 fetch_by_transcript_stable_id

  Arg [1]    : string $transcriptid
               The stable id of the transcript around which the slice is 
               desired
  Arg [2]    : (optional) int $size
               The length of the flanking regions the slice should encompass 
               on either side of the transcript (0 by default)
  Example    : $slc = $sa->fetch_by_transcript_stable_id('ENST00000302930',10);
  Description: Creates a slice around the region of the specified transcript. 
               If a context size is given, the slice is extended by that 
               number of basepairs on either side of the
               transcript.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : Thrown if the transcript is not in the database.
  Caller     : general

=cut

sub fetch_by_transcript_stable_id{
  my ($self,$transcriptid,$size) = @_;

  throw('Transcript argument is required.') if(!$transcriptid);

  my $ta = $self->db->get_TranscriptAdaptor;
  my $transcript = $ta->fetch_by_stable_id($transcriptid);

  throw("Transcript [$transcriptid] does not exist in DB.") if(!$transcript);

  return $self->fetch_by_Feature($transcript, $size);
}



=head2 fetch_by_transcript_id

  Arg [1]    : int $transcriptid
               The unique database identifier of the transcript around which 
               the slice is desired
  Arg [2]    : (optional) int $size
               The length of the flanking regions the slice should encompass 
               on either side of the transcript (0 by default)
  Example    : $slc = $sa->fetch_by_transcript_id(24, 1000);
  Description: Creates a slice around the region of the specified transcript.
               If a context size is given, the slice is extended by that
               number of basepairs on either side of the
               transcript.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : throw on incorrect args
               throw if transcript is not in database
  Caller     : general

=cut

sub fetch_by_transcript_id {
  my ($self,$transcriptid,$size) = @_;

  throw('Transcript id argument is required.') if(!$transcriptid);

  my $transcript_adaptor = $self->db()->get_TranscriptAdaptor();
  my $transcript = $transcript_adaptor->fetch_by_dbID($transcriptid);

  throw("Transcript [$transcriptid] does not exist in DB.") if(!$transcript);

  return $self->fetch_by_Feature($transcript, $size);
}



=head2 fetch_by_gene_stable_id

  Arg [1]    : string $geneid
               The stable id of the gene around which the slice is
               desired
  Arg [2]    : (optional) int $size
               The length of the flanking regions the slice should encompass
               on either side of the gene (0 by default)
  Example    : $slc = $sa->fetch_by_transcript_stable_id('ENSG00000012123',10);
  Description: Creates a slice around the region of the specified gene.
               If a context size is given, the slice is extended by that
               number of basepairs on either side of the gene.
               The slice will be created in the genes native coordinate system.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : throw on incorrect args
               throw if transcript does not exist
  Caller     : general

=cut

sub fetch_by_gene_stable_id {
  my ($self,$geneid,$size) = @_;

  throw('Gene argument is required.') if(!$geneid);

  my $gene_adaptor = $self->db->get_GeneAdaptor();
  my $gene = $gene_adaptor->fetch_by_stable_id($geneid);

  throw("Gene [$geneid] does not exist in DB.") if(!$gene);

  return $self->fetch_by_Feature($gene, $size);
}



=head2 fetch_by_Feature

  Arg [1]    : Bio::EnsEMBL::Feature $feat
               The feature to fetch the slice around
  Arg [2]    : int size (optional)
               The desired number of flanking basepairs around the feature.
  Example    : $slice = $slice_adaptor->fetch_by_Feature($feat, 100);
  Description: Retrieves a slice around a specific feature.  All this really
               does is return a resized version of the slice that the feature
               is already on. Note that slices returned from this method
               are always on the forward strand of the seq_region regardless of
               the strandedness of the feature passed in.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : throw if the feature does not have an attached slice
               throw if feature argument is not provided
  Caller     : fetch_by_gene_stable_id, fetch_by_transcript_stable_id,
               fetch_by_gene_id, fetch_by_transcript_id

=cut

sub fetch_by_Feature{
  my ($self, $feature, $size) = @_;

  $size ||= 0;

  if(!ref($feature) || !$feature->isa('Bio::EnsEMBL::Feature')) {
    throw('Feature argument expected.');
  }

  my $slice = $feature->slice();
  if(!$slice || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Feature must be attached to a valid slice.');
  }

  my $fstart = $feature->start();
  my $fend   = $feature->end();
  if(!defined($fstart) || !defined($fend)) {
    throw('Feature must have defined start and end.');
  }

  #convert the feature slice coordinates to seq_region coordinates
  my $slice_start  = $slice->start();
  my $slice_end    = $slice->end();
  my $slice_strand = $slice->strand();
  if($slice_start != 1 || $slice_strand != 1) {
    if($slice_strand == 1) {
      $fstart = $fstart + $slice_start - 1;
      $fend   = $fend   + $slice_start - 1;
    } else {
      my $tmp_start = $fstart;
      $fstart = $slice_end - $fend      + 1;
      $fend   = $slice_end - $tmp_start + 1;
    }
  }

  #return a new slice covering the region of the feature
  return Bio::EnsEMBL::Slice->new(-seq_region_name => $slice->seq_region_name,
                                  -coord_system    => $slice->coord_system,
                                  -start           => $fstart - $size,
                                  -end             => $fend + $size,
                                  -strand          => 1,
                                  -adaptor         => $self);
}



=head2 fetch_by_misc_feature_attribute

  Arg [1]    : string $attribute_type
               The code of the attribute type
  Arg [2]    : (optional) string $attribute_value
               The value of the attribute to fetch by
  Arg [3]    : (optional) int $size
               The amount of flanking region around the misc feature desired.
  Example    : $slice = $sa->fetch_by_misc_feature_attribute('superctg',
                                                             'NT_030871');
               $slice = $sa->fetch_by_misc_feature_attribute('synonym',
                                                             'AL00012311',
                                                             $flanking);
  Description: Fetches a slice around a MiscFeature with a particular
               attribute type and value. If no value is specified then
               the feature with the particular attribute is used.
               If no size is specified then 0 is used.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : Throw if no feature with the specified attribute type and value
               exists in the database
               Warning if multiple features with the specified attribute type
               and value exist in the database.
  Caller     : webcode

=cut

sub fetch_by_misc_feature_attribute {
  my ($self, $attrib_type_code, $attrib_value, $size) = @_;

  my $mfa = $self->db()->get_MiscFeatureAdaptor();

  my $feats = $mfa->fetch_all_by_attribute_type_value($attrib_type_code,
                                                   $attrib_value);

  if(@$feats == 0) {
    throw("MiscFeature with $attrib_type_code=$attrib_value does " .
          "not exist in DB.");
  }

  if(@$feats > 1) {
    warning("MiscFeature with $attrib_type_code=$attrib_value is " .
            "ambiguous - using first one found.");
  }

  my ($feat) = @$feats;

  return $self->fetch_by_Feature($feat, $size);
}


=head2 get_exceptions_by_slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Example    : 
 Description: get Bio::EnsEMBL::Mapper::Coordinates where on given slice are
              PAR or overlapping haplotypes
  Returntype : listref of Bio::EnsEMBL::Mapper::Coordinate
  Exceptions : none
  Caller     : Slice

=cut


=head2 fetch_normalized_slice_projection

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Example    :  ( optional )
  Description: gives back a project style result. The returned slices 
               represent the areas to which there are symlinks for the 
               given slice. start, end show which area on given slice is 
               symlinked
  Returntype : [[start,end,$slice][]]
  Exceptions : none
  Caller     : BaseFeatureAdaptor

=cut


sub fetch_normalized_slice_projection {
  my $self = shift;
  my $slice = shift;
  
  if( $self->{'_exc_cache'}->{$slice->name()} ) {
    return $self->{'_exc_cache'}->{$slice->name()};
  }
  
  my $result = [];
  my $sql = "
    SELECT seq_region_id, seq_region_start, seq_region_end,
           exc_type, exc_seq_region_id, exc_seq_region_start,
           exc_seq_region_end
      FROM assembly_exception
     WHERE seq_region_id = ?";

  my $slice_seq_region_id = $self->get_seq_region_id( $slice );

  my $sth = $self->prepare( $sql );
  $sth->execute( $slice_seq_region_id );

  my (@haps, @pars);

  my ( $seq_region_id, $seq_region_start, $seq_region_end,
       $exc_type, $exc_seq_region_id, $exc_seq_region_start,
       $exc_seq_region_end );
  $sth->bind_columns( \$seq_region_id, \$seq_region_start, \$seq_region_end,
                      \$exc_type, \$exc_seq_region_id, \$exc_seq_region_start,
                      \$exc_seq_region_end );

  while( $sth->fetch() ) {
    # need overlapping PAR and all HAPs if any
    if( $exc_type eq "PAR" ) {
      if( $seq_region_start <= $slice->end() && 
          $seq_region_end >= $slice->start() ) {
        push( @pars, [ $seq_region_start, $seq_region_end, $exc_seq_region_id,
                       $exc_seq_region_start, $exc_seq_region_end ] );
      }
    } else {
      push( @haps, [ $seq_region_start, $seq_region_end, $exc_seq_region_id,
                     $exc_seq_region_start, $exc_seq_region_end ] );
    }
  }

  if( @pars || @haps ) {
    my @syms;
    
    if( @haps > 1 ) {
      my @sort_haps = sort { $a->[1] <=> $b->[1] } @haps;
      throw( "More than one HAP region not supported yet" );
    } elsif( @haps == 1 ) {
      my $hap = $haps[0];

      my $seq_reg_slice = $self->fetch_by_seq_region_id($slice_seq_region_id);
      my $exc_slice = $self->fetch_by_seq_region_id( $hap->[2] );

      #
      # lengths of haplotype and reference in db may be different
      # we want to use the maximum possible length for the mapping
      # between the two systems
      #
      my $len1 = $seq_reg_slice->length();
      my $len2 = $exc_slice->length();
      my $max_len = ($len1 > $len2) ? $len1 : $len2;

      #the inserted region can differ in length, but mapped sections
      #need to be same lengths
      my $diff = $hap->[4] - $hap->[1];
      
      # we want the region of the haplotype INVERTED
      push( @syms, [ 1, $hap->[0]-1, $hap->[2], 1, $hap->[3] - 1 ] );
      push( @syms, [ $hap->[1]+1, $max_len - $diff, 
                     $hap->[2], $hap->[4] + 1, $max_len ] );
      
    }
    
    # for now haps and pars should not be both there, but in theory we 
    # could handle it here by cleverly merging the pars into the existing syms,
    # for now just:
    push( @syms, @pars );

    my $mapper = Bio::EnsEMBL::Mapper->new( "sym", "org" );
    for my $sym ( @syms ) {
      $mapper->add_map_coordinates( $slice_seq_region_id, $sym->[0], $sym->[1],
                                    1, $sym->[2], $sym->[3], $sym->[4] );
    }

    my @linked = $mapper->map_coordinates( $slice_seq_region_id,
					   $slice->start(), $slice->end(), 
					   $slice->strand(), "sym" );
    
    # gaps are regions where there is no mapping to another region
    my $rel_start = 1;

    #if there was only one coord and it is a gap, we know it is just the
    #same slice with no overlapping symlinks
    if(@linked == 1 && $linked[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
      return [[1,$slice->length, $slice]];
    }

    for my $coord ( @linked ) {
      if( $coord->isa( "Bio::EnsEMBL::Mapper::Gap" )) {
        my $exc_slice = Bio::EnsEMBL::Slice->new
          (-START        => $coord->start(),
           -END          => $coord->end(),
           -STRAND       => $slice->strand(),
           -COORD_SYSTEM => $slice->coord_system(),
           -ADAPTOR      => $self,
           -SEQ_REGION_NAME => $slice->seq_region_name);
        push( @$result, [ $rel_start, $coord->length()+$rel_start-1, 
                          $exc_slice ] );
      } else {
        my $exc_slice = $self->fetch_by_seq_region_id( $coord->id() );
        my $exc2_slice = Bio::EnsEMBL::Slice->new
          (
           -START  => $coord->start(),
           -END    => $coord->end(),
           -STRAND => $coord->strand(),
           -SEQ_REGION_NAME => $exc_slice->seq_region_name(),
           -COORD_SYSTEM => $exc_slice->coord_system(),
           -ADAPTOR => $self
          );
	
        push( @$result, [ $rel_start, $coord->length() + $rel_start - 1, 
                          $exc2_slice ] );
      }
      $rel_start += $coord->length();
    }
  } else {
    #just return this slice, there were no haps or pars
    return  [[1,$slice->length, $slice]];
  }

  # cache and return
  $self->{'_exc_cache'}->{$slice->name()} = $result;
  return $result;
}



#####################################
# sub DEPRECATED METHODs
#####################################

=head2 fetch_by_mapfrag

 Function: DEPRECATED use fetch_by_misc_feature_attribute('synonym',$mapfrag)

=cut

sub fetch_by_mapfrag{
   my ($self,$mymapfrag,$flag,$size) = @_;
   deprecate('Use fetch_by_misc_feature_attribute instead');
   $flag ||= 'fixed-width'; # alt.. 'context'
   $size ||= $flag eq 'fixed-width' ? 100000 : 0;
   return $self->fetch_by_misc_feature_attribute('synonym',$mymapfrag,$size);
}



=head2 fetch_by_chr_start_end

  Description: DEPRECATED use fetch_by_region instead

=cut

sub fetch_by_chr_start_end {
  my ($self,$chr,$start,$end) = @_;
  deprecate('Use fetch_by_region() instead');

  #assume that by chromosome the user actually meant top-level coord
  #system since this is the old behaviour of this deprecated method
  my $csa = $self->db->get_CoordSystemAdaptor();
  my $cs = $csa->fetch_top_level();

  return $self->fetch_by_region($cs->name,$chr,$start,$end,1,$cs->version);
}



=head2 fetch_by_contig_name

  Description: Deprecated. Use fetch_by_region(), Slice::project(), 
               Slice::expand() instead

=cut

sub fetch_by_contig_name {
  my ($self, $name, $size) = @_;

  deprecate('Use fetch_by_region(), Slice::project() and Slice::expand().');

  #previously wanted chromosomal slice on a given contig.  Assume this means
  #a top-level slice on a given seq_region in the seq_level coord system
  my $csa = $self->db()->get_CoordSystemAdaptor();
  my $top_level = $csa->fetch_top_level();
  my $seq_level = $csa->fetch_sequence_level();

  my $seq_lvl_slice = $self->fetch_by_region($seq_level->name(), $name);

  my @projection = @{$seq_lvl_slice->project($top_level->name(),
                                             $top_level->version())};
  if(@projection == 0) {
    warning("contig $name is not used in ".$top_level->name().' assembly.');
    return undef;
  }

  if(@projection > 1) {
    warning("$name is mapped to multiple locations in " . $top_level->name());
  }

  return $projection[0]->[2]->expand($size, $size);
}


=head2 fetch_by_clone_accession

  Description: DEPRECATED.  Use fetch_by_region, Slice::project, Slice::expand
               instead.

=cut

sub fetch_by_clone_accession{
  my ($self,$name,$size) = @_;

  deprecate('Use fetch_by_region(), Slice::project() and Slice::expand().');

  my $csa = $self->db()->get_CoordSystemAdaptor();
  my $top_level = $csa->fetch_top_level();
  my $clone_cs = $csa->fetch_by_name('clone');

  if(!$clone_cs) {
    warning('Clone coordinate system does not exist for this species');
    return undef;
  }

  #this unfortunately needs a version on the end to work
  if(! ($name =~ /\./)) {
    my $sth = $self->prepare("SELECT sr.name " .
                             "FROM   seq_region sr, coord_system cs " .
                             "WHERE  cs.name = 'clone' " .
                             "AND    cs.coord_system_id = sr.coord_system_id ".
                             "AND    sr.name LIKE '$name.%'");
    $sth->execute();
    if(!$sth->rows()) {
      $sth->finish();
      throw("Clone $name not found in database");
    }

    ($name) = $sth->fetchrow_array();

    $sth->finish();
  }

  my $clone = $self->fetch_by_region($clone_cs->name(), $name);
  my @projection = @{$clone->project($top_level->name(),
                                     $top_level->version())};
  if(@projection == 0) {
    warning("clone $name is not used in " . $top_level->name() . ' assembly.');
    return undef;
  }

  if(@projection > 1) {
    warning("$name is mapped to multiple locations in " . $top_level->name());
  }

  return $projection[0]->[2]->expand($size, $size);
}


=head2 fetch_by_supercontig_name

  Description: DEPRECATED. Use fetch_by_region(), Slice::project() and
               Slice::expand() instead

=cut

sub fetch_by_supercontig_name {
  my ($self,$name, $size) = @_;

  deprecate('Use fetch_by_region(), Slice::project() and Slice::expand().');

  my $csa = $self->db()->get_CoordSystemAdaptor();
  my $top_level = $csa->fetch_top_level();
  my $sc_level = $csa->fetch_by_name('supercontig');

  if(!$sc_level) {
    warning('No supercontig coordinate system exists for this species.');
    return undef;
  }

  my $sc_slice = $self->fetch_by_region($sc_level->name(),$name);

  my @projection = @{$sc_slice->project($top_level->name(),
                                        $top_level->version())};
  if(@projection == 0) {
    warning("Supercontig $name is not used in " . $top_level->name() .
         ' assembly.');
    return undef;
  }

  if(@projection > 1) {
    warning("$name is mapped to multiple locations in " . $top_level->name());
  }

  return $projection[0]->[2]->expand($size, $size);
}




=head2 list_overlapping_supercontigs

  Description: DEPRECATED use Slice::project instead

=cut

sub list_overlapping_supercontigs {
   my ($self,$slice) = @_;

   deprecate('Use Slice::project() instead.');

   my $csa = $self->db()->get_CoordSystemAdaptor();
   my $top_level = $csa->fetch_top_level();
   my $sc_level = $csa->fetch_by_name('supercontig');

   if(!$sc_level) {
     warning('No supercontig coordinate system exists for this species.');
     return undef;
   }

   my @out;
   foreach my $seg (@{$slice->project($sc_level->name(), $sc_level->version)}){
     push @out, $seg->[2]->seq_region_name();
   }

   return \@out;
}



=head2 fetch_by_chr_name

  Description: DEPRECATED. Use fetch by region instead

=cut

sub fetch_by_chr_name{
   my ($self,$chr_name) = @_;
   deprecate('Use fetch_by_region() instead.');

   return $self->fetch_by_region('toplevel',$chr_name);
}


1;
