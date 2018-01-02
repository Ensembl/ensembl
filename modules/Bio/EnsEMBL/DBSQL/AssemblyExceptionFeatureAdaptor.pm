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

Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor

=head1 SYNOPSIS

  my $assembly_exception_feature_adaptor =
    $database_adaptor->get_AssemblyExceptionFeatureAdaptor();

  @assembly_exception_features =
    $assembly_exception_feature_adaptor->fetch_all_by_Slice($slice);

=head1 DESCRIPTION

Assembly Exception Feature Adaptor - database access for assembly
exception features.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor;

use strict;
use warnings;
no warnings qw(uninitialized);

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::AssemblyExceptionFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Cache;

our @ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

# set the number of slices you'd like to cache
our $ASSEMBLY_EXCEPTION_FEATURE_CACHE_SIZE = 100;

=head2 new

  Arg [1]    : list of args @args
               Superclass constructor arguments
  Example    : none
  Description: Constructor which just initializes internal cache structures
  Returntype : Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor
  Exceptions : none
  Caller     : implementing subclass constructors
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  # initialize an LRU cache for slices
  my %cache;
  tie(%cache, 'Bio::EnsEMBL::Utils::Cache',
    $ASSEMBLY_EXCEPTION_FEATURE_CACHE_SIZE);

  $self->{'_aexc_slice_cache'} = \%cache;

  return $self;
}

=head2 fetch_all

  Arg [1]    : none
  Example    : my @axfs = @{$axfa->fetch_all()};
  Description: Retrieves all assembly exception features which are in the
               database and builds internal caches of the features.
  Returntype : reference to list of Bio::EnsEMBL::AssemblyExceptionFeatures
  Exceptions : none
  Caller     : fetch_by_dbID, fetch_by_Slice
  Status     : Stable

=cut

sub fetch_all {
  my $self = shift;

  # this is the "global" cache for all assembly exception features in the db
  if(defined($self->{'_aexc_cache'})) {
    return $self->{'_aexc_cache'};
  }

  my $statement = qq(
  SELECT    ae.assembly_exception_id,
            ae.seq_region_id,
            ae.seq_region_start,
            ae.seq_region_end,
            ae.exc_type,
            ae.exc_seq_region_id,
            ae.exc_seq_region_start,
            ae.exc_seq_region_end,
            ae.ori
  FROM      assembly_exception ae,
            coord_system cs,
            seq_region sr
  WHERE     cs.species_id = ?
    AND     sr.coord_system_id = cs.coord_system_id
    AND     sr.seq_region_id = ae.seq_region_id);

  my $sth = $self->prepare($statement);

  $sth->bind_param( 1, $self->species_id(), SQL_INTEGER );

  $sth->execute();

  my ($ax_id, $sr_id, $sr_start, $sr_end,
      $x_type, $x_sr_id, $x_sr_start, $x_sr_end, $ori);

  $sth->bind_columns(\$ax_id, \$sr_id, \$sr_start, \$sr_end,
                     \$x_type, \$x_sr_id, \$x_sr_start, \$x_sr_end, \$ori);

  my @features;
  my $sa = $self->db()->get_SliceAdaptor();

  $self->{'_aexc_dbID_cache'} = {};

  while($sth->fetch()) {
    my $slice   = $sa->fetch_by_seq_region_id($sr_id);
    my $x_slice = $sa->fetch_by_seq_region_id($x_sr_id);

    # each row creates TWO features, each of which has alternate_slice
    # pointing to the "other" one

   
    my $a = Bio::EnsEMBL::AssemblyExceptionFeature->new
          ('-dbID'            => $ax_id,
           '-start'           => $sr_start,
           '-end'             => $sr_end,
           '-strand'          => 1,
           '-adaptor'         => $self,
           '-slice'           => $slice,
           '-alternate_slice' => $x_slice->sub_Slice($x_sr_start, $x_sr_end),
           '-type'            => $x_type);
   
    push @features, $a;
    $self->{'_aexc_dbID_cache'}->{$ax_id} = $a;

    push @features, Bio::EnsEMBL::AssemblyExceptionFeature->new
          ('-dbID'            => $ax_id,
           '-start'           => $x_sr_start,
           '-end'             => $x_sr_end,
           '-strand'          => 1,
           '-adaptor'         => $self,
           '-slice'           => $x_slice,
           '-alternate_slice' => $slice->sub_Slice($sr_start, $sr_end),
           '-type'            => "$x_type REF" );
  }

  $sth->finish();

  $self->{'_aexc_cache'} = \@features;
  
  return \@features;
}


=head2 fetch_by_dbID

  Arg [1]    : int $dbID
  Example    : my $axf = $axfa->fetch_by_dbID(3);
  Description: Retrieves a single assembly exception feature via its internal
               identifier.  Note that this only retrieves one of the two
               assembly exception features which are represented by a single
               row in the assembly_exception table.
  Returntype : Bio::EnsEMBL::AssemblyExceptionFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  if(!exists($self->{'_aexc_dbID_cache'})) {
    # force loading of cache
    $self->fetch_all();
  }

  return $self->{'_aexc_dbID_cache'}->{$dbID};
}


=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Example    : my @axfs = @{$axfa->fetch_all_by_Slice($slice)};
  Description: Retrieves all assembly exception features which overlap the
               provided slice.  The returned features will be in coordinate
               system of the slice.
  Returntype : reference to list of Bio::EnsEMBL::AssemblyException features
  Exceptions : none
  Caller     : Feature::get_all_alt_locations, general
  Status     : Stable

=cut

sub fetch_all_by_Slice {
  my $self = shift;
  my $slice = shift;

  my $key= uc($slice->name());

  # return features from the slice cache if present
  if(exists($self->{'_aexc_slice_cache'}->{$key})) {
    return $self->{'_aexc_slice_cache'}->{$key};
  }

  my $all_features = $self->fetch_all();

  my $mcc = $self->db()->get_MetaCoordContainer();
  my $css = $mcc->fetch_all_CoordSystems_by_feature_type('assembly_exception');

  my @features;

  my $ma = $self->db()->get_AssemblyMapperAdaptor();

  foreach my $cs (@$css) {
    my $mapper;
    if($cs->equals($slice->coord_system)) {
      $mapper = undef;
    } else {
      $mapper = $ma->fetch_by_CoordSystems($cs,$slice->coord_system());
    }

    push @features, @{ $self->_remap($all_features, $mapper, $slice) };
  }

  $self->{'_aexc_slice_cache'}->{$key} = \@features;

  return \@features;
}


#
# Given a list of features checks if they are in the correct coord system
# by looking at the first features slice.  If they are not then they are
# converted and placed on the slice.
#
# Note that this is a re-implementation of a method with the same name in
# BaseFeatureAdaptor, and in contrast to the latter which maps features in
# place, this method returns a remapped copy of each feature. The reason for
# this is to get around conflicts with caching.
#
sub _remap {
  my ($self, $features, $mapper, $slice) = @_;

  # check if any remapping is actually needed
  if(@$features && (!$features->[0]->isa('Bio::EnsEMBL::Feature') ||
                    $features->[0]->slice == $slice)) {
    return $features;
  }
  # remapping has not been done, we have to do our own conversion from
  # to slice coords

  my @out;

  my $slice_start = $slice->start();
  my $slice_end   = $slice->end();
  my $slice_strand = $slice->strand();
  my $slice_cs    = $slice->coord_system();

  my ($seq_region, $start, $end, $strand, $seq_region_name);

  my $slice_seq_region = $slice->seq_region_name();

  foreach my $f (@$features) {
    # since feats were obtained in contig coords, attached seq is a contig
    my $fslice = $f->slice();
    if(!$fslice) {
      throw("Feature does not have attached slice.\n");
    }
    my $fseq_region = $fslice->seq_region_name();
    my $fcs = $fslice->coord_system();
    if(!$slice_cs->equals($fcs)) {
      # slice of feature in different coord system, mapping required
      ($seq_region, $start, $end, $strand) =
        $mapper->fastmap($fseq_region,$f->start(),$f->end(),$f->strand(),$fcs);
      # undefined start means gap
      next if(!defined $start);

      my $slice_adaptor = $self->db()->get_SliceAdaptor();
      $seq_region_name = $slice_adaptor->fetch_by_seq_region_id($seq_region)->seq_region_name;

    } else {
      $start           = $f->start();
      $end             = $f->end();
      $strand          = $f->strand();
      $seq_region_name = $f->slice->seq_region_name();
    }

    # maps to region outside desired area
    next if ($start > $slice_end) || ($end < $slice_start) || 
      ($slice_seq_region ne $seq_region_name);

    # create new copies of successfully mapped feaatures with shifted start,
    # end and strand
    my ($new_start, $new_end);
    my $seq_region_len = $slice->seq_region_length();

    if ($slice_strand == 1) { # Positive strand
		
      $new_start = $start - $slice_start + 1;
      $new_end   = $end - $slice_start + 1;

      if ($slice->is_circular()) {
        # Handle circular chromosomes.

        if ($new_start > $new_end) {
    	  # Looking at a feature overlapping the chromsome origin.
          if ($new_end > $slice_start) {
            # Looking at the region in the beginning of the chromosome.
            $new_start -= $seq_region_len;
          }
          if ($new_end < 0) {
            $new_end += $seq_region_len;
          }
        } else {
          if (   $slice_start > $slice_end && $new_end < 0) {
            # Looking at the region overlapping the chromosome
            # origin and a feature which is at the beginning of the
            # chromosome.
            $new_start += $seq_region_len;
            $new_end   += $seq_region_len;
          }
        }
      }   ## end if ($dest_slice->is_circular...)
    } else {  # Negative strand

      $new_start = $slice_end - $end + 1;
      $new_end = $slice_end - $start + 1;

      if ($slice->is_circular()) {

        if ($slice_start > $slice_end) { 
        # slice spans origin or replication

          if ($start >= $slice_start) {
            $new_end += $seq_region_len;
            $new_start += $seq_region_len if $end > $slice_start;

          } elsif ($start <= $slice_end) {
            # do nothing
          } elsif ($end >= $slice_start) {
            $new_start += $seq_region_len;
            $new_end += $seq_region_len;

          } elsif ($end <= $slice_end) {
            $new_end += $seq_region_len if $new_end < 0;
          } elsif ($start > $end) {
            $new_end += $seq_region_len;
          } else {

          }
              
        } else {

          if ($start <= $slice_end and $end >= $slice_start) {
            # do nothing
          } elsif ($start > $end) {
            if ($start <= $slice_end) {
          
              $new_start -= $seq_region_len;

            } elsif ($end >= $slice_start) {
              $new_end += $seq_region_len;
            } else {

            }
          }
        }

      }

    }    ## end else [ if ($dest_slice_strand...)]
    push @out, Bio::EnsEMBL::AssemblyExceptionFeature->new(
            '-dbID'            => $f->dbID,
            '-start'           => $new_start,
            '-end'             => $new_end,
            '-strand'          => $strand * $slice_strand,
            '-adaptor'         => $self,
            '-slice'           => $slice,
            '-alternate_slice' => $f->alternate_slice,
            '-type'            => $f->type,
    );
  }  # end foreach assembly exception

  return \@out;
}

=head2 store

    Arg[1]       : Bio::EnsEMBL::AssemblyException $asx
    Arg[2]       : Bio::EnsEMBL::AssemblyException $asx2

    Example      : $asx = Bio::EnsEMBL::AssemblyExceptionFeature->new(...)
                   $asx2 = Bio::EnsEMBL::AssemblyExceptionFeature->new(...)
                   $asx_seq_region_id = $asx_adaptor->store($asx);
    Description:  This stores a assembly exception feature in the 
                  assembly_exception table and returns the assembly_exception_id.
                  Needs 2 features: one pointing to the Assembly_exception, and the
                  other pointing to the region in the reference that is being mapped to
                  Will check that start, end and type are defined, and the alternate
                  slice is present as well.
    ReturnType:   int
    Exceptions:   throw if assembly exception not defined (needs start, end,
		  type and alternate_slice) of if $asx not a Bio::EnsEMBL::AssemblyException
    Caller:       general
    Status:       Stable

=cut

sub store{
  my $self = shift;
  my $asx = shift;
  my $asx2 = shift;


  if (! $asx->isa('Bio::EnsEMBL::AssemblyExceptionFeature')){
    throw("$asx is not a Ensembl assemlby exception -- not stored");
  }
  #if already present, return ID in the database
  my $db = $self->db();
  if ($asx->is_stored($db)){
    return $asx->dbID();
  }
  #do some checkings for the object
  #at the moment, the orientation is always 1
  if (! $asx->start || ! $asx->end ){
    throw("Assembly exception does not have coordinates");
  }
  if ($asx->type !~ /PAR|HAP|PATCH_NOVEL|PATCH_FIX/){
    throw("Only types of assembly exception features valid are PAR, HAP, PATCH_FIX or PATCH_NOVEL");
  }
  if ( !($asx->alternate_slice->isa('Bio::EnsEMBL::Slice')) ){
    throw("Alternate slice should be a Bio::EnsEMBL::Slice");
  }
  #now check the other Assembly exception feature, the one pointing to the REF
  # region
  if (!$asx2->isa('Bio::EnsEMBL::AssemblyExceptionFeature')){
    throw("$asx2 is not a Ensembl assemlby exception -- not stored");
  }
  if (! $asx2->start || ! $asx2->end ){
    throw("Assembly exception does not have coordinates");
  }
  if ($asx2->type !~ /HAP REF|PAR REF|PATCH_NOVEL REF|PATCH_FIX REF/){
    throw("$asx2 should have  type of assembly exception features HAP REF, PAR REF, PATCH_FIX REF or PATCH_NOVEL REF");
  }
  if (! ($asx2->alternate_slice->isa('Bio::EnsEMBL::Slice')) ){
    throw("Alternate slice should be a Bio::EnsEMBL::Slice");
  }
  #finally check that both features are pointing to each other slice
  if ($asx->slice != $asx2->alternate_slice || $asx->alternate_slice != $asx2->slice){
    throw("Slice and alternate slice in both features are not pointing to each other");
  }
  #prepare the SQL
  my $asx_sql = q{
	INSERT INTO assembly_exception( seq_region_id, seq_region_start,
					seq_region_end, 
					exc_type, exc_seq_region_id,
					exc_seq_region_start, exc_seq_region_end,
					ori)
	    VALUES (?, ?, ?, ?, ?, ?, ?, 1)
	};
  
  my $asx_st = $self->prepare($asx_sql);
  my $asx_id = undef;
  my $asx_seq_region_id;
  my $asx2_seq_region_id;
  my $original = $asx;
  my $original2 = $asx2;
  #check all feature information
  ($asx, $asx_seq_region_id) = $self->_pre_store($asx);
  ($asx2, $asx2_seq_region_id) = $self->_pre_store($asx2);
  
  #and store it
  $asx_st->bind_param(1, $asx_seq_region_id, SQL_INTEGER);
  $asx_st->bind_param(2, $asx->start(), SQL_INTEGER);
  $asx_st->bind_param(3, $asx->end(), SQL_INTEGER);
  $asx_st->bind_param(4, $asx->type(), SQL_VARCHAR);
  $asx_st->bind_param(5, $asx2_seq_region_id, SQL_INTEGER);
  $asx_st->bind_param(6, $asx2->start(), SQL_INTEGER);
  $asx_st->bind_param(7, $asx2->end(), SQL_INTEGER);
  
  $asx_st->execute();
  $asx_id = $self->last_insert_id('assembly_exception_id', undef, 'assembly_exception');
  
  #finally, update the dbID and adaptor of the asx and asx2
  $original->adaptor($self);
  $original->dbID($asx_id);
  $original2->adaptor($self);
  $original2->dbID($asx_id);
  #and finally update dbID cache with new assembly exception
  $self->{'_aexc_dbID_cache'}->{$asx_id} = $original;
  #and update the other caches as well
  push @{$self->{'_aexc_slice_cache'}->{uc($asx->slice->name)}},$original, $original2;
  push @{$self->{'_aexc_cache'}}, $original, $original2;
  
  return $asx_id;
}

#
# Helper function containing some common feature storing functionality
#
# Given a Feature this will return a copy (or the same feature if no changes 
# to the feature are needed) of the feature which is relative to the start
# of the seq_region it is on. The seq_region_id of the seq_region it is on
# is also returned.
#
# This method will also ensure that the database knows which coordinate
# systems that this feature is stored in.
# Since this adaptor doesn't inherit from BaseFeatureAdaptor, we need to copy
# the code
#

sub _pre_store {
  my $self    = shift;
  my $feature = shift;

  if(!ref($feature) || !$feature->isa('Bio::EnsEMBL::Feature')) {
    throw('Expected Feature argument.');
  }

  $self->_check_start_end_strand($feature->start(),$feature->end(),
                                 $feature->strand());


  my $db = $self->db();

  my $slice_adaptor = $db->get_SliceAdaptor();
  my $slice = $feature->slice();

  if(!ref($slice) || !($slice->isa('Bio::EnsEMBL::Slice') or $slice->isa('Bio::EnsEMBL::LRGSlice'))  ) {
    throw('Feature must be attached to Slice to be stored.');
  }

  # make sure feature coords are relative to start of entire seq_region

  if($slice->start != 1 || $slice->strand != 1) {
    #move feature onto a slice of the entire seq_region
    $slice = $slice_adaptor->fetch_by_region($slice->coord_system->name(),
                                             $slice->seq_region_name(),
                                             undef, #start
                                             undef, #end
                                             undef, #strand
                                             $slice->coord_system->version());

    $feature = $feature->transfer($slice);

    if(!$feature) {
      throw('Could not transfer Feature to slice of ' .
            'entire seq_region prior to storing');
    }
  }

  # Ensure this type of feature is known to be stored in this coord system.
  my $cs = $slice->coord_system;

   my $mcc = $db->get_MetaCoordContainer();

  $mcc->add_feature_type($cs, 'assembly_exception', $feature->length);

  my $seq_region_id = $slice_adaptor->get_seq_region_id($slice);

  if(!$seq_region_id) {
    throw('Feature is associated with seq_region which is not in this DB.');
  }

  return ($feature, $seq_region_id);
}

#
# helper function used to validate start/end/strand and 
# hstart/hend/hstrand etc.
#
sub _check_start_end_strand {
  my $self = shift;
  my $start = shift;
  my $end   = shift;
  my $strand = shift;

  #
  # Make sure that the start, end, strand are valid
  #
  if(int($start) != $start) {
    throw("Invalid Feature start [$start].  Must be integer.");
  }
  if(int($end) != $end) {
    throw("Invalid Feature end [$end]. Must be integer.");
  }
  if(int($strand) != $strand || $strand < -1 || $strand > 1) {
    throw("Invalid Feature strand [$strand]. Must be -1, 0 or 1.");
  }
  if($end < $start) {
    throw("Invalid Feature start/end [$start/$end]. Start must be less " .
          "than or equal to end.");
  }

  return 1;
}

=head2 remove

  Arg [1]    : $asx Bio::EnsEMBL::AssemblyFeatureException 
  Example    : $asx_adaptor->remove($asx);
  Description: This removes a assembly exception feature from the database.  
  Returntype : none
  Exceptions : thrown if $asx arg does not implement dbID(), or if
               $asx->dbID is not a true value
  Caller     : general
  Status     : Stable

=cut

#again, this method is generic in BaseFeatureAdaptor, but since this class
#is not inheriting, need to copy&paste
sub remove {
  my ($self, $feature) = @_;

  if(!$feature || !ref($feature) || !$feature->isa('Bio::EnsEMBL::AssemblyExceptionFeature')) {
    throw('AssemblyExceptionFeature argument is required');
  }

  if(!$feature->is_stored($self->db)) {
    throw("This feature is not stored in this database");
  }

  my $asx_id = $feature->dbID();
  my $key = uc($feature->slice->name);
  my $sth = $self->prepare("DELETE FROM assembly_exception WHERE assembly_exception_id = ?");
  $sth->bind_param(1,$feature->dbID,SQL_INTEGER);
  $sth->execute();
  
  #and clear cache
  #and finally update dbID cache
  delete $self->{'_aexc_dbID_cache'}->{$asx_id};
  #and remove from cache feature
  my @features;
  foreach my $asx (@{$self->{'_aexc_slice_cache'}->{$key}}){
      if ($asx->dbID != $asx_id){
	  push @features, $asx;
      }
  }
  $self->{'_aexc_slice_cache'}->{$key} = \@features;
  @features = ();
  foreach my $asx (@{$self->{'_aexc_cache'}}){
      if ($asx->dbID != $asx_id){
	  push @features, $asx;
      }
  }
  $self->{'_aexc_cache'} = \@features;

#unset the feature dbID ad adaptor
  $feature->dbID(undef);
  $feature->adaptor(undef);

  return;
}


1;
