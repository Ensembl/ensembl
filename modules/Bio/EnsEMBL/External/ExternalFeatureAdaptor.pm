#
# EnsEMBL module for Bio::EnsEMBL::External::ExternalFeatureAdaptor
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::External::ExternalFeatureAdaptor - Allows features created externally from EnsEMBL in a single coordinate system to be retrieved in several other (EnsEMBL-style) coordinate systems. This is intended to be a replacement for the old Bio::EnsEMBL::DB::ExternalFeatureFactoryI interface.

=head1 SYNOPSIS

  $database_adaptor = 
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => 'kaka.sanger.ac.uk',
                                        -dbname => 'homo_sapiens_core_9_30',
                                        -pass   => 'anonymous' );

  $xf_adaptor = new ExternalFeatureAdaptorSubClass;

  #Connect the EnsEMBL core database:
  $xf_adaptor->db($database_adaptor);

  #get some features in RawContig coords
  @feats = @{$xf_adaptor->fetch_all_by_contig_name('AC000087.2.1.42071')};

  #another way of doing the same thing
  $contig_adaptor = $database_adaptor->get_RawContigAdaptor;
  $contig = $contig_adaptor->fetch_by_name('AC000087.2.1.42071');
  @feats = @{$xf_adaptor->fetch_all_by_RawContig($contig);

  #get some features in assembly coords
  @feats = @{$xf_adaptor->fetch_all_by_chr_start_end('X', 100000, 200000)};

  #get some features in clone coords
  @feats = @{$xf_adaptor->fetch_all_by_clone_accession('AC000087')};

  #another way of doing the same thing
  $clone_adaptor = $db_adaptor->get_CloneAdaptor;
  $clone = $clone_adaptor->fetch_by_accession('AC000087');
  @feats = ${$xf_adaptor->fetch_all_by_Clone($clone);

  #Add the adaptor to the ensembl core dbadaptor (implicitly sets db attribute)
  $database_adaptor->add_ExternalFeatureAdaptor($xf_adaptor);

  #get some features in Slice coords
  $slice_adaptor = $database_adaptor->get_SliceAdaptor;
  $slice = $slice_adaptor->fetch_by_chr_start_end(1,100000,200000);
  @feats = @{$xf_adaptor->fetch_all_by_Slice($slice)};

  #now features can be retrieved directly from slice or RawContig
  @feats = @{$slice->get_all_ExternalFeatures};
  @feats = @{$contig->get_all_ExternalFeatures};
  

=head1 DESCRIPTION

This class is intended to be used as a method of getting external features into
EnsEMBL.  To work, this class must be extended and must implement the
the coordinate_systems method.  As well, the subclass is required to implement
a single fetch method so that the external features may be retrieved.  
By implementing a single fetch_method in a single coordinate system all
of the other ExternalFeatureAdaptor fetch methods become available for 
retrieving the data in several different coordinate systems.

The coordinate_systems method should return a list of strings indicating which
coordinate system(s) have been implemented.  If a given string is returned 
from the coordinate_systems method then the corresponding fetch method must be 
implemented.  The reverse is also true: if a fetch method is implemented then
coordinate_systems must return the appropriate string in its list of return 
values.  The following are the valid coordinate system values and the 
corresponding fetch methods that must be implemented:

  COORD SYSTEM STRING   FETCH_METHOD      
  -------------------   ------------
  'ASSEMBLY'            fetch_all_by_chr_start_end
  'CLONE'               fetch_all_by_clone_accession
  'CONTIG'              fetch_all_by_contig_name
  'SLICE'               fetch_all_by_Slice

For convenience fetch_all_by_RawContig and fetch_all_by_Clone methods are
also available.  These methods may be overridden if desired but due to
internal dependencies the fetch_all_by_contig_name and 
fetch_all_by_clone_accession must also be overridden if these methods are
altered.  See the method descriptions for more detail.

The objects returned by the fetch methods should be Bio::SeqFeature object,
though only the start, end, strand and attach_seq methods are actually used
by the ExternalFeatureAdaptor.  The objects which are returned by the 
ExternalFeature adaptor will be altered by the functions called.

Before the non-overridden ExternalFeatureAdaptor fetch methods may be called
an EnsEMBL core database adaptor must be attached to the ExternalFeatureAdaptor
.  This database adaptor is required to perform the remappings between various
coordinate system.  This may be done implicitly by adding the 
ExternalFeatureAdaptor to the database adaptor through a call to the 
DBAdaptor add_ExternalFeatureAdaptor method or explicitly by calling the 
ExternalFeatureAdaptor db method.


=head1 CONTACT

Post questions to the EnsEMBL developer list: <ensembl-dev@ebi.ac.uk>

=head1 AUTHOR

Graham McVicker

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

use strict;

package Bio::EnsEMBL::External::ExternalFeatureAdaptor;

use vars qw(@ISA);

use Bio::EnsEMBL::Root;


@ISA = qw(Bio::EnsEMBL::Root);


=head2 new

  Arg [1]    : none
  Example    : $xfa = new Bio::EnsEMBL::External::ExternalFeatureAdaptor;
  Description: Creates a new ExternalFeatureAdaptor object.  You may wish to
               extend this constructor and provide your own set of paremeters.
  Returntype : Bio::EnsEMBL::External::ExternalFeatureAdaptor
  Exceptions : none
  Caller     : general

=cut

sub new {
  my $class = shift;

  if(ref $class) {
    return bless {}, ref $class;
  }

  return bless {}, $class;
}



=head2 ensembl_db

  Arg [1]    : (optional) Bio::EnsEMBL::DBSQL::DBAdaptor
  Example    : $external_feature_adaptor->ensembl_db($new_val);
  Description: none
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : internal

=cut

sub ensembl_db {
  my ($self, $value) = @_;

  if($value) {
    #avoid potentially nasty memory leaks
    if(ref $value && $value->isa("Bio::EnsEMBL::Container")) {
      $self->{'ensembl_db'} = $value->_obj;
    } else {
      $self->{'ensembl_db'} = $value;
    }
  }     

  return $self->{'ensembl_db'};
}



=head2 coordinate_systems

  Arg [1]    : none
  Example    : @implemented_coord_systems = $ext_adaptor->coordinate_systems;
  Description: ABSTRACT method. Must be implemented by all 
               ExternalFeatureAdaptor subclasses.  This method returns a list
               of coordinate systems which are implemented by the subclass. 
               A minimum of on valid coordinate system must be implemented.
               Valid coordinate systems are: 'SLICE', 'ASSEMBLY', 'CONTIG',
               and 'CLONE'.
  Returntype : list of strings
  Exceptions : none
  Caller     : internal

=cut

sub coordinate_systems {
  my $self = shift;
  
  $self->throw("abstract method coordinate_systems not implemented\n");

  return '';
}


=head2 track_name

  Arg [1]    : none
  Example    : $track_name = $xf_adaptor->track_name;
  Description: Currently this is not really used.  In the future it may be 
               possible to have ExternalFeatures automatically displayed by
               the EnsEMBL web code.  By default this method returns 
               'External features' but you are encouraged to override this 
               method and provide your own meaningful name for the features
               your adaptor provides.  This also allows you to distinguish the
               type of features retrieved from RawContigs or Slices.  See
               the PODs for Bio::EnsEMBL::Slice::get_all_ExternalFeatures and 
               Bio::EnsEMBL::DBSQL::DBAdaptor::add_ExternalFeatureAdaptor 
               methods. 
  Returntype : string
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::DBAdaptor::add_ExternalFeatureAdaptor

=cut

sub track_name {
  my $self = shift;

  return 'External features';
}



=head2 feature_type

  Arg [1]    : none
  Example    : $feature_type = $xf_adaptor->track_name
  Description: Currently this is not used.  In the future it may be possible
               to have ExternalFeatures automatically displayed by the EnsEMBL
               web code.  This method would then be used do determine the 
               type of glyphs used to draw the features which are returned
               from this external adaptor.
  Returntype : string
  Exceptions : none
  Caller     : none

=cut

sub feature_type {
  my $self = shift;
  
  return qw(SIMPLE);
}



=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Example    : @features = @{$ext_adaptor->fetch_all_by_Slice($slice)};
  Description: Retrieves all features which lie in the region defined
               by $slice in slice coordinates.  
  
               If this method is overridden then the coordinate_systems method
               must return 'SLICE' as one of its values.  

               This method will work as is (i.e. without overriding it) 
               providing at least one of the other fetch methods is overridden.
  Returntype : reference to a list of Bio::SeqFeature objects in the Slice
               coordinate system
  Exceptions : Thrown on incorrect input arguments
  Caller     : general, fetch_all_by_chr_start_end

=cut

sub fetch_all_by_Slice {
  my ($self, $slice) = @_;

  unless($slice && ref $slice && $slice->isa('Bio::EnsEMBL::Slice')) {
    $self->throw("[$slice] is not a Bio::EnsEMBL::Slice");
  }

  my $out = [];

  my $slice_start  = $slice->chr_start;
  my $slice_end    = $slice->chr_end;
  my $slice_strand = $slice->strand;
  my $slice_chr    = $slice->chr_name;

  if($self->_supported('SLICE')) {
    $self->throw("ExternalFeatureAdaptor supports SLICE coordinate system" .
		 " but fetch_all_by_Slice not implemented");
  } 

  #fetch the features in assembly coords
  $out = $self->fetch_all_by_chr_start_end($slice_chr,$slice_start,$slice_end,$slice_strand);

  #convert from assembly coords to slice coords
  my($f_start, $f_end, $f_strand);
  foreach my $f (@$out) {
    if($slice_strand == 1) {
      $f_start  = $f->start - $slice_start + 1;
      $f_end    = $f->end   - $slice_start + 1;
      $f_strand = $f->strand;
    } else {
      $f_start  = $slice_end - $f->end   + 1;
      $f_end    = $slice_end - $f->start + 1;
      $f_strand = $f->strand * -1;
    }
    
    $f->start($f_start);
    $f->end($f_end);
    $f->strand($f_strand);
    $f->attach_seq($slice);
  }
  
  return $out;
}
  


=head2 fetch_all_by_RawContig

  Arg [1]    : Bio::EnsEMBL::RawContig $contig
  Arg [2]    : (optional) int $start
               The start coordinate of the RawContig to retrieve features from.
  Arg [3]    : (optional) int $end
               The end coordinate of the RawContig to retrieve features from.
  Example    : @features = @{$self->fetch_all_by_RawContig($contig)};
  Description: Retrieves features on the region defined by the $contig arg in 
               RawContig coordinates. 

               If this method is overridden then it is also necessary to 
               override the fetch_all_by_contig_name method due to 
               interdependencies.  As well, if this method is overridden then 
               the coordinate_systems method must return 'CONTIG' as one of 
               its values.  

               It is probably more useful to only override the 
               fetch_all_by_contig_name instead of both methods. This method 
               will work as is - providing at least one other fetch method has 
               been overridden.               
  Returntype : reference to a list of Bio::SeqFeature objects in the RawContig
               coordinate system.
  Exceptions : thrown if the input argument is incorrect
  Caller     : general, fetch_all_by_contig_name, 
               fetch_all_by_Clone, fetch_all_by_chr_start_end

=cut

sub fetch_all_by_RawContig {
  my ($self, $contig, $start, $end) = @_;

  unless($contig && ref $contig && $contig->isa('Bio::EnsEMBL::RawContig')) {
    $self->throw("[$contig] is not a Bio::EnsEMBL::RawContig");
  }

  my $out = [];

  if($self->_supported('CONTIG')) {
    return $self->fetch_all_by_contig_name($contig->name);
  }

  if($self->_supported('CLONE')) {
    #retrieve features in clone coordinates

    #if we didn't get start/end, fetch features from whole clone
    my $offset = $contig->embl_offset;
    my $length = $contig->length;
    $start ||= $offset;
    $end ||= $offset + $length - 1;
    
    my $feats = $self->fetch_all_by_Clone($contig->clone, $start, $end);
    
    my ($start, $end);

    #convert from clone coordinates to contig coordinates
    foreach my $f (@$feats) {
      $start = $f->start - $offset + 1;
      $end   = $f->end   - $offset + 1;
      
      #skip features which are not entirely on this contig
      next if($start < 1 || $end > $length);
      
      $f->start($start);
      $f->end($end);
      $f->attach_seq($contig);
      push(@$out, $f);
    }

    return $out;
  }

  unless($self->_supported('SLICE') || $self->_supported('ASSEMBLY')) {
    $self->_no_valid_coord_system;
  }

  unless($self->ensembl_db) {
    $self->throw('DB attribute not set.  This value must be set for the ' .
		 'ExternalFeatureAdaptor to function correctly');
  }

  my $asma = $self->ensembl_db->get_AssemblyMapperAdaptor;
  my $mapper = $asma->fetch_by_type($self->ensembl_db->assembly_type);
  
  #map the whole contig to the assembly
  my @mapped = $mapper->map_coordinates_to_assembly($contig->dbID, 1, 
						    $contig->length, 1);
  
  my $chr;

  foreach my $coord (@mapped) {
    #skip over parts of the contig that do not map to the assembly
    next if ($coord->isa('Bio::EnsEMBL::Mapper::Gap'));
    
    $chr = $coord->id;

    #retrieve the features in assembly coordinates
    push @$out, @{$self->fetch_all_by_chr_start_end($coord->id, 
						$coord->start, 
						$coord->end,
                                                $coord->strand)};
  }
    
  #map each feature from assembly coords back to raw contig coords
  foreach my $f (@$out) {
    @mapped =$mapper->map_coordinates_to_rawcontig($chr,
						   $f->start, 
						   $f->end,
						   $f->strand);
    
    #skip if feature does not map cleanly to our contig
    next if(@mapped != 1);
    next if($mapped[0]->isa('Bio::EnsEMBL::Mapper::Gap'));
    
    unless($mapped[0]->id == $contig->dbID) {
      $self->throw('Error mapping feature from assembly to contig');
    }
    
    $f->start ($mapped[0]->start);
    $f->end   ($mapped[0]->end);
    $f->strand($mapped[0]->strand);
    $f->attach_seq($contig);
  }

  return $out;
}


=head2 fetch_all_by_contig_name

  Arg [1]    : string $contig_name
  Example    : @fs = @{$self->fetch_all_by_contig_name('AB00879.1.1.39436')};
  Description: Retrieves features on the contig defined by the name 
               $contig_name in RawContig coordinates.

               If this method is overridden then the coordinate_systems 
               method must return 'CONTIG' as one of its values. If the 
               fetch_all_by_Contig method is overridden then this method
               should also be overridden to chain calls to that method (and
               avoid throwing an exception when this method is called).  

               This method will work as is (i.e. without being overridden) 
               providing at least one other fetch method has 
               been overridden and the fetch_all_by_Contig method has not been
               overridden.               
  Returntype : reference to a list of Bio::SeqFeature objects in the RawContig
               coordinate system.
  Exceptions : thrown if the input argument is incorrect
               thrown if the coordinate_systems method returns the value 
               'CONTIG' and this method has not been overridden.
  Caller     : general, fetch_all_by_Contig 

=cut

sub fetch_all_by_contig_name {
  my ($self, $contig_name) = @_;

  unless($contig_name) {
    $self->throw("contig_name argument not defined");
  }

  if($self->_supported('CONTIG')) {
    $self->throw("ExternalFeatureAdaptor supports CONTIG coordinate system" .
		 " but fetch_all_by_contig_name is not implemented");
  }

  unless($self->ensembl_db) {
    $self->throw('DB attribute not set.  This value must be set for the ' .
		 'ExternalFeatureAdaptor to function correctly');
  }

  my $contig_adaptor = $self->ensembl_db->get_RawContigAdaptor; 
  my $contig = $contig_adaptor->fetch_by_name($contig_name);

  unless($contig) {
    $self->warn("ExternalFeatureAdaptor::fetch_all_by_contig_name: Contig " .
		"[$contig_name] not found\n");
    return [];
  }
    

  return $self->fetch_all_by_RawContig($contig);
}




=head2 fetch_all_by_Clone

  Arg [1]    : Bio::EnsEMBL::Clone $clone
  Arg [2]    : (optional) int $clone_start the start of the clonal region
               interested in. This information may be used to speed up
               the query, or may be ignored.
  Arg [3]    : (optional) int $clone_end the end of the clonal region 
               interested in.  This information may be used to speed up 
               the query, or may be ignored.
  Example    : @features = @{$self->fetch_all_by_Clone($clone)};
  Description: Retrieves features on the region defined by the $clone arg in 
               Clone coordinates. 
               
               If this method is overridden then it is also necessary to 
               override the fetch_all_by_clone_accession method due to 
               interdependencies.  As well, if this method is overridden then 
               the coordinate_systems method must return 'CLONE' as one of its
               values.  
              
               It is probably more useful to override only the 
               fetch_all_by_clone_accession method rather than both of these 
               methods. This method will work as is - providing at least one 
               other fetch method has been overridden.               
  Returntype : reference to a list of Bio::SeqFeature objects in the Clone
               coordinate system
  Exceptions : thrown if the input argument is incorrect
               thrown if the coordinate systems method does not return any 
               valid values.
  Caller     : general, fetch_all_by_clone_accession, 
               fetch_all_by_RawContig

=cut

sub fetch_all_by_Clone {
  my ($self, $clone, $clone_start, $clone_end) = @_;

  unless($clone && ref $clone && $clone->isa('Bio::EnsEMBL::Clone')) {
    $self->throw("Clone must be a Bio::EnsEMBL::Clone");
  }

  if($self->_supported('CLONE')) {
    return $self->fetch_all_by_clone_accession($clone->id, $clone->embl_id,
					       $clone_start, $clone_end);
  }
  
  unless($self->_supported('CONTIG') || 
	 $self->_supported('SLICE')  ||
	 $self->_supported('ASSEMBLY')) {
    $self->_no_valid_coord_system;
  }

  my $out = [];

  #retrieve features from each contig in the clone
  my $contigs = $clone->get_all_Contigs;
  my $offset;
  foreach my $contig (@$contigs) {
    $offset = $contig->embl_offset;
    foreach my $f (@{$self->fetch_all_by_RawContig($contig)}) {
      #convert each feature to clone coordinates
      $f->start($f->start - $offset + 1);
      $f->end($f->end - $offset + 1);
      #$f->attach_seq($clone); #this might work in future...
      push @$out, $f;
    }
  }

  return $out;
}



=head2 fetch_all_by_clone_accession

  Arg [1]    : string $acc
               The EMBL accession number of the clone to fetch features from.
  Arg [2]    : (optional) string $ver
  Arg [3]    : (optional) int $start
  Arg [4]    : (optional) int $end
 
  Example    : @fs = @{$self->fetch_all_by_clone_accession('AC000093')};
  Description: Retrieves features on the clone defined by the $acc arg in 
               Clone coordinates. 
               
               If this method is overridden then the coordinate_systems method
               must return 'CLONE' as one of its values. The arguments 
               start, end, version are passed if this method is overridden and
               can optionally be used to reduce the scope of the query and 
               improve performance.  

               If the fetch_all_by_Clone method has been overridden then this
               method must also be overridden to chain calls to that method
               (and avoid throwing an exception).
               This method will work as is - providing at least one other 
               fetch method has been overridden and the fetch_all_by_Clone 
               method has not been overridden.
  Returntype : reference to a list of Bio::SeqFeature objects in the Clone
               coordinate system
  Exceptions : thrown if the input argument is incorrect
               thrown if the coordinate system method returns the value 'CLONE'
               and this method is not overridden.
               thrown if the coordinate systems method does not return any 
               valid values.
  Caller     : general, fetch_all_by_clone_accession, 
               fetch_all_by_RawContig

=cut

sub fetch_all_by_clone_accession {
  my ($self, $acc, $version, $start, $end) = @_;

  unless($acc) {
    $self->throw("clone accession argument not defined");
  }

  if($self->_supported('CLONE')) {
    $self->throw('ExternalFeatureAdaptor supports CLONE coordinate system ' .
		 'but does not implement fetch_all_by_clone_accession');
  }

  unless($self->ensembl_db) {
    $self->throw('DB attribute not set.  This value must be set for the ' .
		 'ExternalFeatureAdaptor to function correctly');
  }


  my $clone_adaptor = $self->ensembl_db->get_CloneAdaptor;
  my $clone = $clone_adaptor->fetch_by_accession($acc);


  unless($clone) {
    $self->warn("ExternalFeatureAdaptor::fetch_all_by_clone_accession: Clone "
		. "[$acc] not found\n");
    return [];
  }
    

  return $self->fetch_all_by_Clone($clone);
}



=head2 fetch_all_by_chr_start_end

  Arg [1]    : string $chr_name
               The name of the chromosome to retrieve features from
  Arg [2]    : int $start
               The start coordinate of the chromosomal region to retrieve
               features from.
  Arg [3]    : int $end
               The end coordinate of the chromosomal region to retrieve 
               features from.
  Arg [4]    : (optional) int $strand
               The strand of the chromosomal region to retrieve features from.
  Example    : @features
  Description: Retrieves features on the region defined by the $chr_name,
               $start, and $end args in assembly (chromosomal) coordinates. 

               If this method is overridden then the coordinate_systems method
               must return 'ASSEMBLY' as one of its values.  

               This method will work as is (i.e. without overriding it) 
               providing at least one of the other fetch methods is overridden.
  Returntype : reference to a list of Bio::SeqFeatures 
  Exceptions : Thrown if the coordinate_systems method returns ASSEMBLY as a 
               value and this method is not overridden.  
               Thrown if any of the input arguments are incorrect
  Caller     : general, fetch_all_by_Slice, fetch_all_by_RawContig

=cut

sub fetch_all_by_chr_start_end {
  my ($self, $chr_name, $start, $end, $strand) = @_;

  unless($chr_name && defined $start && defined $end && $start < $end) {
    $self->throw("Incorrect start [$start] end [$end] or chr [$chr_name] arg");
  }

  unless($self->ensembl_db) {
    $self->throw('DB attribute not set.  This value must be set for the ' .
		 'ExternalFeatureAdaptor to function correctly');
  }

  my $out = [];

  my $asma = $self->ensembl_db->get_AssemblyMapperAdaptor;
  my $mapper = $asma->fetch_by_type($self->ensembl_db->assembly_type);
  my $contig_adaptor = $self->ensembl_db->get_RawContigAdaptor;

  if($self->_supported('ASSEMBLY')) {
    $self->throw("ExternalFeatureAdaptor supports ASSEMBLY coordinate system".
		 " but fetch_all_by_chr_start_end is not implemented");

  } 

  my $slice_adaptor = $self->ensembl_db->get_SliceAdaptor;
  #get a slice of the whole chromosome
  my $chrom_slice = $slice_adaptor->fetch_by_chr_name($chr_name);

  if($self->_supported('SLICE')) {

    #fetch by slice and convert to assembly coords
    my $slice = $slice_adaptor->fetch_by_chr_start_end($chr_name,$start,$end);
    $out = $self->fetch_all_by_Slice($slice);
    
    foreach my $f (@$out) {
      $f->start($start + $f->start - 1);
      $f->end  ($start + $f->end   - 1);
      $f->attach_seq($chrom_slice);
    }

    return $out;
  }
  
  #fetch via rawcontig and convert to assembly coords

  #Figure out what contigs we are overlapping
  my @cids = $mapper->list_contig_ids($chr_name, $start, $end);
  #convert start/end from assembly to contig coords if we know the strand
  my %contig_coords;
  if ($strand) {
    my @mapped_contig = $mapper->map_coordinates_to_rawcontig($chr_name,
                                                              $start,
                                                              $end,
                                                              $strand);
    foreach my $mc (@mapped_contig) {
      next if $mc->isa('Bio::EnsEMBL::Mapper::Gap');
      $contig_coords{$mc->id} = [$mc->start, $mc->end];
    }
  }
    
  foreach my $cid (@cids) {
    my $contig = $contig_adaptor->fetch_by_dbID($cid);
    #retrieve features of each contig that we are overlapping
    my $feats = $self->fetch_all_by_RawContig($contig, @{$contig_coords{$cid}});
    
    foreach my $f (@$feats) {
      #convert each feature from contig coords to assembly coords
      my @mapped = $mapper->map_coordinates_to_assembly($cid,
							$f->start, 
							$f->end, 
							$f->strand);
      
      #if maps to multiple locations in assembly, skip feature
      next if(@mapped > 1);
	
      #if maps to a gap, skip
      next if($mapped[0]->isa('Bio::EnsEMBL::Mapper::Gap'));
      
      my $m_start  = $mapped[0]->start;
      my $m_end    = $mapped[0]->end;
      my $m_strand = $mapped[0]->strand;

      #skip features which do not overlap this assembly region
      next if($m_start > $end || $m_end < $start);

      $f->start($m_start);
      $f->end($m_end);
      $f->strand($m_strand);
      $f->attach_seq($chrom_slice);

      push @$out, $f;
    }
  }
    
  return $out;
}


=head2 _no_valid_coord_system

  Arg [1]    : none
  Example    : none
  Description: PRIVATE method - throws an error with a descriptive message
  Returntype : none
  Exceptions : always thrown
  Caller     : internal

=cut

sub _no_valid_coord_system {
  my $self = shift;

  $self->throw("This ExternalFeatureAdaptor does not support a known " .
		"coordinate system.\n Valid coordinate systems are: " .
		"[SLICE, ASSEMBLY, CONTIG, CLONE].\n This External Adaptor " . 
                "supports: [" . join(', ', $self->coordinate_systems) . "]");
}  




=head2 _supported

  Arg [1]    : string $system 
  Example    : print "CONTIG system supported" if($self->_supported('CONTIG'));
  Description: PRIVATE method. Tests if the coordinate system defined by
               the $system argument is implemented.
  Returntype : boolean
  Exceptions : none
  Caller     : internal

=cut

sub _supported {
  my ($self, $system) = @_;

  #construct the hash of supported features if it has not been already
  unless(exists $self->{_supported}) {
    $self->{_supported} = {};
    foreach my $coord_system ($self->coordinate_systems) {
      $self->{_supported}->{$coord_system} = 1;
    }
  }

  return $self->{_supported}->{$system};
}
  


1;
