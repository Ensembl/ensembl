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

Bio::EnsEMBL::External::ExternalFeatureAdaptor

=head 1 SUMMARY

Allows features created externally from Ensembl in a single
coordinate system to be retrieved in several other (Ensembl-style)
coordinate systems. This is intended to be a replacement for the old
Bio::EnsEMBL::DB::ExternalFeatureFactoryI interface.

=head1 SYNOPSIS

  $database_adaptor = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => 'kaka.sanger.ac.uk',
    -dbname => 'homo_sapiens_core_9_30',
    -pass   => 'anonymous'
  );

  $xf_adaptor = new ExternalFeatureAdaptorSubClass;

  # Connect the Ensembl core database:
  $xf_adaptor->db($database_adaptor);

  # get some features in vontig coords
  @feats =
    @{ $xf_adaptor->fetch_all_by_contig_name('AC000087.2.1.42071') };

  # get some features in assembly coords
  @feats =
    @{ $xf_adaptor->fetch_all_by_chr_start_end( 'X', 100000, 200000 ) };

  # get some features in clone coords
  @feats = @{ $xf_adaptor->fetch_all_by_clone_accession('AC000087') };

  # Add the adaptor to the ensembl core dbadaptor (implicitly sets db
  # attribute)
  $database_adaptor->add_ExternalFeatureAdaptor($xf_adaptor);

  # get some features in Slice coords
  $slice_adaptor = $database_adaptor->get_SliceAdaptor;
  $slice =
    $slice_adaptor->fetch_all_by_chr_start_end( 1, 100000, 200000 );
  @feats = @{ $xf_adaptor->fetch_all_by_Slice($slice) };

  # now features can be retrieved directly from Slice
  @feats = @{ $slice->get_all_ExternalFeatures };

=head1 DESCRIPTION

This class is intended to be used as a method of getting external
features into EnsEMBL.  To work, this class must be extended and must
implement the the coordinate_systems method.  As well, the subclass
is required to implement a single fetch method so that the external
features may be retrieved.  By implementing a single fetch_method in a
single coordinate system all of the other ExternalFeatureAdaptor fetch
methods become available for retrieving the data in several different
coordinate systems.

The coordinate_systems method should return a list of strings indicating
which coordinate system(s) have been implemented.  If a given string is
returned from the coordinate_systems method then the corresponding fetch
method must be implemented.  The reverse is also true: if a fetch method
is implemented then coordinate_systems must return the appropriate
string in its list of return values.  The following are the valid
coordinate system values and the corresponding fetch methods that must
be implemented:

  COORD SYSTEM STRING   FETCH_METHOD
  -------------------   ------------
  'ASSEMBLY'            fetch_all_by_chr_start_end
  'CLONE'               fetch_all_by_clone_accession
  'CONTIG'              fetch_all_by_contig_name
  'SUPERCONTIG'         fetch_all_by_supercontig_name
  'SLICE'               fetch_all_by_Slice

The objects returned by the fetch methods should be EnsEMBL or BioPerl
style Feature objects.  These objects MUST have start, end and strand
methods.

Before the non-overridden ExternalFeatureAdaptor fetch methods may
be called an EnsEMBL core database adaptor must be attached to the
ExternalFeatureAdaptor . This database adaptor is required to perform
the remappings between various coordinate system.  This may be done
implicitly by adding the ExternalFeatureAdaptor to the database adaptor
through a call to the DBAdaptor add_ExternalFeatureAdaptor method or
explicitly by calling the ExternalFeatureAdaptor ensembl_db method.

=head1 METHODS

=cut

package Bio::EnsEMBL::External::ExternalFeatureAdaptor;

use strict;

use Bio::EnsEMBL::Utils::Exception qw(warning throw);


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
    $self->{'ensembl_db'} = $value;
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

  throw("abstract method coordinate_systems not implemented\n");

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
               type of features retrieved from Slices.  See
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
    throw("[$slice] is not a Bio::EnsEMBL::Slice");
  }

  my $out = [];

  my $csa = $self->ensembl_db->get_CoordSystemAdaptor();

  my $slice_start  = $slice->start;
  my $slice_end    = $slice->end;
  my $slice_strand = $slice->strand;
  my $slice_seq_region  = $slice->seq_region_name;
  my $slice_seq_region_id  = $slice->get_seq_region_id;
  my $coord_system = $slice->coord_system;

  if($self->_supported('SLICE')) {
    throw("ExternalFeatureAdaptor supports SLICE coordinate system" .
		 " but fetch_all_by_Slice not implemented");
  }

  my %features;
  my $from_coord_system;

  my $fetch_method;

  #
  # Get all of the features from whatever coord system they are computed in
  #
  if($self->_supported('CLONE')) {
    $fetch_method = sub {
      my $self = shift;
      my $name = shift;
      my ($acc, $ver) = split(/\./, $name);
      $self->fetch_all_by_clone_accession($acc,$ver,@_);
    };
    $from_coord_system = $csa->fetch_by_name('clone');
  } elsif($self->_supported('ASSEMBLY')) {
    $from_coord_system = $csa->fetch_by_name('chromosome');
    $fetch_method = $self->can('fetch_all_by_chr_start_end');
  } elsif($self->_supported('CONTIG')) {
    $from_coord_system = $csa->fetch_by_name('contig');
    $fetch_method = $self->can('fetch_all_by_contig_name');
  } elsif($self->_supported('SUPERCONTIG')) {
    $from_coord_system = $csa->fetch_by_name('supercontig');
    $fetch_method = $self->can('fetch_all_by_supercontig_name');
  } else {
    $self->_no_valid_coord_systems();
  }

  if($from_coord_system->equals($coord_system)) {
    $features{$slice_seq_region} = &$fetch_method($self, $slice_seq_region,
                                                  $slice_start,$slice_end);
  } else {
    foreach my $segment (@{$slice->project($from_coord_system->name,
                                           $from_coord_system->version)}) {
      my ($start,$end,$pslice) = @$segment;
      $features{$pslice->seq_region_name } ||= [];
      push @{$features{$pslice->seq_region_name }},
           @{&$fetch_method($self, $pslice->seq_region_name,
                            $pslice->start(),
                            $pslice->end())};
    }
  }

  my @out;

  if(!$coord_system->equals($from_coord_system)) {             
    my $asma = $self->ensembl_db->get_AssemblyMapperAdaptor();
    my $mapper = $asma->fetch_by_CoordSystems($from_coord_system,
                                              $coord_system);
    my %slice_cache;
    my $slice_adaptor = $self->ensembl_db->get_SliceAdaptor();
    my $slice_setter;

    #convert the coordinates of each of the features retrieved
    foreach my $fseq_region (keys %features) {
      my $feats = $features{$fseq_region};
      next if(!$feats);
      $slice_setter = _guess_slice_setter($feats) if(!$slice_setter);

      foreach my $f (@$feats) {
        my($sr_id, $start, $end, $strand) = 
          $mapper->fastmap($fseq_region,$f->start,$f->end,$f->strand,
                           $from_coord_system);
        
        #maps to gap
        next if(!defined($sr_id));

        #maps to unexpected seq region, probably error in the externally
        if($sr_id ne $slice_seq_region_id) {
          warning("Externally created Feature mapped to [$sr_id] " .
                  "which is not on requested seq_region_id [$slice_seq_region_id]");
          next;
        }

        #update the coordinates of the feature
        &$slice_setter($f,$slice);
        $f->start($start);
        $f->end($end);
        $f->strand($strand);
        push @out, $f;
      }
    }
  } else {
    #we already know the seqregion the featues are on, we just have
    #to put them on the slice
    @out = @{$features{$slice_seq_region}};
    my $slice_setter = _guess_slice_setter(\@out);

    foreach my $f (@out) {
      &$slice_setter($f,$slice);
    }
  }

  # convert from assembly coords to slice coords
  # handle the circular slice case
  my $seq_region_len = $slice->seq_region_length();
  foreach my $f (@out) {
    my($f_start, $f_end, $f_strand);

    if ($slice->strand == 1) { # Positive strand
      $f_start = $f->start - $slice_start + 1;
      $f_end   = $f->end - $slice_start + 1;
      $f_strand = $f->strand;

      if ($slice->is_circular()) { # Handle cicular chromosomes
	if ($f_start > $f_end) { # Looking at a feature overlapping the chromsome origin.
	  if ($f_end > $slice_start) {
	    # Looking at the region in the beginning of the chromosome.
	    $f_start -= $seq_region_len;
	  }

	  if ($f_end < 0) {
	    $f_end += $seq_region_len;
	  }
	} else {
	  if ($slice_start > $slice_end && $f_end < 0) {
	    # Looking at the region overlapping the chromosome origin and 
	    # a feature which is at the beginning of the chromosome.
	    $f_start += $seq_region_len;
	    $f_end   += $seq_region_len;
	  }
	}
      }
    } else { # Negative strand
      my ($seq_region_start, $seq_region_end) = ($f->start, $f->end);
      $f_start = $slice_end - $seq_region_end + 1;
      $f_end = $slice_end - $seq_region_start + 1;
      $f_strand = $f->strand * -1;

      if ($slice->is_circular()) {
	if ($slice_start > $slice_end) { # slice spans origin or replication
	  if ($seq_region_start >= $slice_start) {
	    $f_end += $seq_region_len;
	    $f_start += $seq_region_len 
	      if $seq_region_end > $slice_start;

	  } elsif ($seq_region_start <= $slice_end) {
	    # do nothing
	  } elsif ($seq_region_end >= $slice_start) {
	    $f_start += $seq_region_len;
	    $f_end += $seq_region_len;
	  } elsif ($seq_region_end <= $slice_end) {
	    $f_end += $seq_region_len
	      if $f_end < 0;
	  } elsif ($seq_region_start > $seq_region_end) {
	    $f_end += $seq_region_len;
	  } else { }
	} else {
	  if ($seq_region_start <= $slice_end and $seq_region_end >= $slice_start) {
	    # do nothing
	  } elsif ($seq_region_start > $seq_region_end) {
	    if ($seq_region_start <= $slice_end) {
	      $f_start -= $seq_region_len;
	    } elsif ($seq_region_end >= $slice_start) {
	      $f_end += $seq_region_len;
	    } else { }
	  }
	}
      }
    }

    $f->start($f_start);
    $f->end($f_end);
    $f->strand($f_strand);    
  }
  
  return \@out;
}
  

sub _guess_slice_setter {
  my $features = shift;

  #we do not know what type of features these are.  They might
  #be bioperl features or old ensembl features, hopefully they are new
  #style features.  Try to come up with a setter method for the
  #slice.

  return undef if(!@$features);

  my ($f) = @$features;

  my $slice_setter;
  foreach my $method (qw(slice contig attach_seq)) {
    last if($slice_setter = $f->can($method));
  }
    
  if(!$slice_setter) {
    if($f->can('seqname')) {
      $slice_setter = sub { $_[0]->seqname($_[1]->seq_region_name()); };
    } else {
      $slice_setter = sub{} if(!$slice_setter);
    }
  }

  return $slice_setter;
}


=head2 fetch_all_by_contig_name

  Arg [1]    : string $contig_name
  Arg [2]    : int $start (optional)
               The start of the region on the contig to retrieve features on
               if not specified the whole of the contig is used.
  Arg [3]    : int $end (optional) 
               The end of the region on the contig to retrieve features on
               if not specified the whole of the contig is used.
  Example    : @fs = @{$self->fetch_all_by_contig_name('AB00879.1.1.39436')};
  Description: Retrieves features on the contig defined by the name 
               $contig_name in contig coordinates.

               If this method is overridden then the coordinate_systems 
               method must return 'CONTIG' as one of its values. 

               This method will work as is (i.e. without being overridden) 
               providing at least one other fetch method has 
               been overridden.               
  Returntype : reference to a list of Bio::SeqFeature objects in the contig
               coordinate system.
  Exceptions : thrown if the input argument is incorrect
               thrown if the coordinate_systems method returns the value 
               'CONTIG' and this method has not been overridden.
  Caller     : general, fetch_all_by_Slice

=cut

sub fetch_all_by_contig_name {
  my ($self, $contig_name, $start, $end) = @_;

  unless($contig_name) {
    throw("contig_name argument not defined");
  }

  if($self->_supported('CONTIG')) {
    throw("ExternalFeatureAdaptor supports CONTIG coordinate system" .
		 " but fetch_all_by_contig_name is not implemented");
  }

  unless($self->ensembl_db) {
    throw('DB attribute not set.  This value must be set for the ' .
		 'ExternalFeatureAdaptor to function correctly');
  }

  my $slice_adaptor = $self->ensembl_db->get_SliceAdaptor();
  my $slice = $slice_adaptor->fetch_by_region('contig', $contig_name,
                                             $start, $end);
  return $self->fetch_all_by_Slice($slice);
}



=head2 fetch_all_by_supercontig_name

  Arg [1]    : string $supercontig_name
  Arg [2]    : int $start (optional)
               The start of the region on the contig to retrieve features on
               if not specified the whole of the contig is used.
  Arg [3]    : int $end (optional) 
               The end of the region on the contig to retrieve features on
               if not specified the whole of the contig is used.
  Example    : @fs = @{$self->fetch_all_by_contig_name('NT_004321')};
  Description: Retrieves features on the contig defined by the name 
               $supercontigname in supercontig coordinates.

               If this method is overridden then the coordinate_systems 
               method must return 'SUPERCONTIG' as one of its values. 

               This method will work as is (i.e. without being overridden)
               providing at least one other fetch method has 
               been overridden.
  Returntype : reference to a list of Bio::SeqFeature objects in the contig
               coordinate system.
  Exceptions : thrown if the input argument is incorrect
               thrown if the coordinate_systems method returns the value
               'SUPERCONTIG' and this method has not been overridden.
  Caller     : general, fetch_all_by_Slice

=cut


sub fetch_all_by_supercontig_name {
  my ($self, $supercontig_name, $start, $end) = @_;

  unless($supercontig_name) {
    throw("supercontig_name argument not defined");
  }

  if($self->_supported('SUPERCONTIG')) {
    throw("ExternalFeatureAdaptor supports SUPERCONTIG coordinate system" .
		 " but fetch_all_by_supercontig_name is not implemented");
  }

  unless($self->ensembl_db) {
    throw('DB attribute not set.  This value must be set for the ' .
		 'ExternalFeatureAdaptor to function correctly');
  }

  my $slice_adaptor = $self->ensembl_db->get_SliceAdaptor();
  my $slice = $slice_adaptor->fetch_by_region('supercontig', $supercontig_name,
                                             $start, $end);
  return $self->fetch_all_by_Slice($slice);
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

               This method will work as is - providing at least one other
               fetch method has been overridden.
  Returntype : reference to a list of Bio::SeqFeature objects in the Clone
               coordinate system
  Exceptions : thrown if the input argument is incorrect
               thrown if the coordinate system method returns the value 'CLONE'
               and this method is not overridden.
               thrown if the coordinate systems method does not return any 
               valid values.
  Caller     : general, fetch_all_by_clone_accession

=cut

sub fetch_all_by_clone_accession {
  my ($self, $acc, $version, $start, $end) = @_;

  unless($acc) {
    throw("clone accession argument not defined");
  }

  if($self->_supported('CLONE')) {
    throw('ExternalFeatureAdaptor supports CLONE coordinate system ' .
		 'but does not implement fetch_all_by_clone_accession');
  }

  unless($self->ensembl_db) {
    throw('DB attribute not set.  This value must be set for the ' .
		 'ExternalFeatureAdaptor to function correctly');
  }

  if(defined($version)) {
    $acc = "$acc.$version";
  } elsif(!$acc =~ /\./) {
    $acc = "$acc.1";
  }

  my $slice_adaptor = $self->ensembl_db->get_SliceAdaptor;

  my $slice = $slice_adaptor->fetch_by_region('clone', $acc, $start, $end);

  return $self->fetch_all_by_Slice($slice);
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
  Caller     : general, fetch_all_by_Slice

=cut

sub fetch_all_by_chr_start_end {
  my ($self, $chr_name, $start, $end) = @_;

  unless($chr_name && defined $start && defined $end && $start < $end) {
    throw("Incorrect start [$start] end [$end] or chr [$chr_name] arg");
  }

  unless($self->ensembl_db) {
    throw('DB attribute not set.  This value must be set for the ' .
		 'ExternalFeatureAdaptor to function correctly');
  }

  my $slice_adaptor = $self->ensembl_db->get_SliceAdaptor();

  my $slice = $slice_adaptor->fetch_by_region('toplevel', $chr_name, $start,
                                              $end);

  return $self->fetch_all_by_Slice($slice);
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

  throw("This ExternalFeatureAdaptor does not support a known " .
		"coordinate system.\n Valid coordinate systems are: " .
		"[SLICE,ASSEMBLY,SUPERCONTIG,CONTIG,CLONE].\n This External Adaptor " . 
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
