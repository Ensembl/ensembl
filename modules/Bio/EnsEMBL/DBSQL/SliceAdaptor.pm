
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

Bio::EnsEMBL::DBSQL::SliceAdaptor - Adaptors for slices

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

use Bio::EnsEMBL::Utils::Exception qw(throw deprecate);

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
  my $key = join(':',$seq_region_name,
                 $coord_system->name(),
                 $coord_system->version);

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
      throw("Cannot create slice on non-existant or ambigous seq_region:" .
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
    my $key = join(':', $name, $cs->name, $cs->version);
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

  my $key = join(':', $seq_region_name,$cs_name,$cs_version);

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
    throw("Non existant or ambigous seq_region:" .
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



sub deleteObj {
  my $self = shift;

  $self->SUPER::deleteObj;

  $self->{'_id_cache'} = {};
  $self->{'_name_cache'} = {};
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

    my $type = $self->db->assembly_type();

    warn( "XX>" ,$chr, "--", $band );
    my $sth = $self->db->prepare("
        select min(k.chr_start), max(k.chr_end)
          from chromosome as c, karyotype as k
         where c.chromosome_id = k.chromosome_id and c.name=? and k.band like ?
    ");
    $sth->execute( $chr, "$band%" );
    my ( $slice_start, $slice_end) = $sth->fetchrow_array;

    warn( $chr, "--", $band );
    unless( defined($slice_start) ) {
       my $sth = $self->db->prepare("
           select min(k.chr_start), max(k.chr_end)
             from chromosome as c, karyotype as k
            where c.chromosome_id = k.chromosome_id and k.band like ?
       ");
       $sth->execute( "$band%" );
       ( $slice_start, $slice_end) = $sth->fetchrow_array;
    }

   if(defined $slice_start) {
        return new Bio::EnsEMBL::Slice(
           -chr_name  => $chr,
           -chr_start => $slice_start,
           -chr_end   => $slice_end,
           -strand    => 1,
           -assembly_type => $type
        );
    }

    $self->throw("Band not recognised in database");
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
               transcript.  Throws if the transcript is not golden.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_transcript_stable_id{
  my ($self,$transcriptid,$size) = @_;

  # Just get the dbID, then fetch slice by that
  my $ta = $self->db->get_TranscriptAdaptor;
  my $transcript_obj = $ta->fetch_by_stable_id($transcriptid);
  my $dbID = $transcript_obj->dbID;
  
  return $self->fetch_by_transcript_id($dbID, $size);
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
  Exceptions : thrown on incorrect args
  Caller     : general

=cut

sub fetch_by_transcript_id {
  my ($self,$transcriptid,$size) = @_;

  unless( defined $transcriptid ) {
    $self->throw("Must have transcriptid id to fetch Slice of transcript");
  }

  $size = 0 unless(defined $size);
   
  my $ta = $self->db->get_TranscriptAdaptor;
  my $transcript_obj = $ta->fetch_by_dbID($transcriptid);
  
  my %exon_transforms;
  
  my $emptyslice;
  for my $exon ( @{$transcript_obj->get_all_Exons()} ) {
    $emptyslice = Bio::EnsEMBL::Slice->new( '-empty'   => 1,
					    '-adaptor' => $self,
					    '-ASSEMBLY_TYPE' =>
					    $self->db->assembly_type);     
    my $newExon = $exon->transform( $emptyslice );
    $exon_transforms{ $exon } = $newExon;
  }
  
  $transcript_obj->transform( \%exon_transforms );
  
  my $start = $transcript_obj->start() - $size;
  my $end = $transcript_obj->end() + $size;
  
  if($start < 1) {
    $start = 1;
  }
  
  my $slice = $self->fetch_by_chr_start_end($emptyslice->chr_name,
					    $start, $end);
  return $slice;
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
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_gene_stable_id{
   my ($self,$geneid,$size) = @_;

   if( !defined $geneid ) {
       $self->throw("Must have gene id to fetch Slice of gene");
   }
   if( !defined $size ) {$size=0;}

   my ($chr_name,$start,$end) = $self->_get_chr_start_end_of_gene($geneid);

   if( !defined $start ) {
     my $type = $self->db->assembly_type()
       or $self->throw("No assembly type defined");
     $self->throw("Gene [$geneid] is not on the golden path '$type'. " .
		  "Cannot build Slice.");
   }
     
   $start -= $size;
   $end += $size;
   
   if($start < 1) {
     $start = 1;
   }

   return $self->fetch_by_chr_start_end($chr_name, $start, $end);
}






=head2 fetch_by_mapfrag

 Title   : fetch_by_mapfrag
 Usage   : $slice = $slice_adaptor->fetch_by_mapfrag('20');
 Function: Creates a slice of a "mapfrag"
 Returns : Slice object
 Args    : chromosome name


=cut

sub fetch_by_mapfrag{
   my ($self,$mymapfrag,$flag,$size) = @_;

   $flag ||= 'fixed-width'; # alt.. 'context'
   $size ||= $flag eq 'fixed-width' ? 200000 : 0;
   unless( $mymapfrag ) {
       $self->throw("Mapfrag name argument required");
   }

   my( $chr_start,$chr_end);
  
   #set the end of the slice to the end of the chromosome
   my $ca = $self->db()->get_MapFragAdaptor();
   my $mapfrag = $ca->fetch_by_synonym($mymapfrag);
   return undef unless defined $mapfrag;

   if( $flag eq 'fixed-width' ) {
       my $halfsize = int( $size/2 );
       $chr_start = $mapfrag->seq_start - $halfsize;
       $chr_end   = $mapfrag->seq_start + $size - $halfsize;
   } else {
       $chr_start     = $mapfrag->seq_start - $size;
       $chr_end       = $mapfrag->seq_end   + $size;
   }
   my $type = $self->db->assembly_type();

   my $slice = Bio::EnsEMBL::Slice->new
     (
      -chr_name      => $mapfrag->seq,
      -chr_start     => $chr_start,
      -chr_end       => $chr_end,
      -assembly_type => $type,
      -adaptor       => $self
     );

   return $slice;
}






=head2 _get_chr_start_end_of_gene

 Title   : get_Gene_chr_bp
 Usage   : 
 Function: 
 Returns :  
 Args    :


=cut


sub _get_chr_start_end_of_gene {
  my ($self,$geneid) =  @_;
  
  my $type = $self->db->assembly_type()
    or $self->throw("No assembly type defined");
  
  my $sth = $self->db->prepare("SELECT  
   if(a.contig_ori=1,(e.contig_start-a.contig_start+a.chr_start),
                    (a.chr_start+a.contig_end-e.contig_end)),
   if(a.contig_ori=1,(e.contig_end-a.contig_start+a.chr_start),
                    (a.chr_start+a.contig_end-e.contig_start)),
     chr.name
  
                    FROM    exon e,
                        transcript tr,
                        exon_transcript et,
                        assembly a,
                        gene_stable_id gsi,
                        chromosome chr
                    WHERE e.exon_id=et.exon_id 
                    AND et.transcript_id =tr.transcript_id 
                    AND a.contig_id=e.contig_id 
                    AND a.type = '$type' 
                    AND tr.gene_id = gsi.gene_id
                    AND gsi.stable_id = '$geneid'
                    AND a.chromosome_id = chr.chromosome_id" 
                    );
   $sth->execute();

   my ($start,$end,$chr);
   my @start;
   while ( my @row=$sth->fetchrow_array){
      ($start,$end,$chr)=@row;
       push @start,$start;
       push @start,$end;
   }   
   
   my @start_sorted=sort { $a <=> $b } @start;

   $start=shift @start_sorted;
   $end=pop @start_sorted;

   return ($chr,$start,$end);      
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
  my $cs = $csa->get_top_coord_system();

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
    warn("$name is mapped to multiple locations in " . $top_level->name());
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

  my $clone = $self->fetch_by_region($clone_cs->name(), $name);
  my @projection = @{$clone->project($top_level->name(),
                                     $top_level->version())};
  if(@projection == 0) {
    warn("clone $name is not used in " . $top_level->name() . ' assembly.');
    return undef;
  }

  if(@projection > 1) {
    warn("$name is mapped to multiple locations in " . $top_level->name());
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
    warn('No supercontig coordinate system exists for this species.');
    return undef;
  }

  my $sc_slice = $self->fetch_by_region($sc_level->name(),$name);

  my @projection = @{$sc_slice->project($top_level->name(),
                                        $top_level->version())};
  if(@projection == 0) {
    warn("Supercontig $name is not used in " . $top_level->name() .
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
   foreach my $seg ($slice->project($sc_level->name(), $sc_level->version)){
     push @out, $seg->[2]->seq_region_name();
   }

   return \@out;
}



=head2 fetch_by_chr_name

  Description: DEPRECATED. Use fetch by region instead

=cut

sub fetch_by_chr_name{
   my ($self,$chr_name) = @_;
   deprecate('Use fetch_by_region() instead');

   my $csa = $self->db()->get_CoordSystemAdaptor();
   my $cs = $csa->fetch_top_level();
   return $self->fetch_by_region($cs->name());
}


1;






