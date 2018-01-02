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

Bio::EnsEMBL::Utils::Collector

=head1 SYNOPSIS

  # Inherit this base module in your feature specific Collector
  # instance:

  package Bio::EnsEMBL::Funcgen::Collector::ResultFeature;
  use base('Bio::EnsEMBL::Utils::Collector');

  # ... and define package config variables
  $Bio::EnsEMBL::Funcgen::Collector::bin_model = 'SIMPLE';
  $Bio::EnsEMBL::Funcgen::Collector::window_sizes =
    [ 30, 65, 130, 260, 450, 648, 950, 1296 ];
  # Could replace 30 with 0 here for low density data at natural resolution

  $Bio::EnsEMBL::Utils::Collector::bin_method =
    'count';    # only used by collector

  $Bio::EnsEMBL::Utils::Collector::packed_size = 2;

  # ... or simply use this module in a script either defining package
  # config variables, or passing as parameters to the constructor:

  my $collector =
    Bio::EnsEMBL::Utils::BaseCollector->new( -pack_template => 'v' );

  $Bio::EnsEMBL::Funcgen::Collector::pack_template = 'v';

  # Config variables can also be over-ridden by passing a config hash to
  # the store_window_bins_by_Slice() method:

  $collector->store_window_bins_by_Slice( $slice, (
                                            -pack_template => 'v',
                                            -packed_size   => 2 ) );

  # NOTE: Over-riding default config variables can cause problems when
  # storing or fetching data. e.g. Fetch may revert to using defaults or
  # table partitions may not match window sizes.

=head1 DESCRIPTION

This package is the base Collector class which contains generic
getter/setter methods along with the main 'collecting' methods which
perform the majority of the work in generating compressed data
collections optimised for web display.  The bins produced are aimed at
storage in a BLOB representing an entire seq_region i.e. even bins with
no features/null data are encoded as a 0 score.  Non-BLOB collections
are currently not supported.

If your Collection class defines a Bio::EnsEMBL::Feature, then its
adaptor should inherit from the relevant Collection class.

The minimum prerequisites of the input features/data are that they have
a start() and end() method.  For instance a Bio::EnsEMBL::Features
generated from a database or parsed from a flat file.

NOTE: This Collector does not have a lightweight mode previously used
for dynamic/on the fly collecting i.e. it does not take advantage of
bypassing object creation via the related BaseFeatureAdaptor method.

=cut

package Bio::EnsEMBL::Utils::Collector;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument  ('rearrange');
use Bio::EnsEMBL::Utils::Exception ('throw');

### Global package config vars

# Defaults
our $max_view_width = 1000000;   # Max bp width in location/detailed view
our $max_data_type_size = 16777216;    # Default is 16MB for long BLOB
# This is really a guide value as this should be set in the inheriting
# Collector class by deducting the rest of the row size from this value.
# Is is upto the inheritor to handle checking whether this size has been 
# exceeded.

# NOTE: Theoretically the min window size is: slice_length/(16777216/2)
# So for human chr1:  249,250,621/(16777216/2) = 29.7 => 30.  However,
# this size does not seem to directly translate to the MySQL
# max_allowed_packet_size.  Increasing max_allowed_packet_size to 64MB
# solves this issue, and substr operation doesn't appear to incur any of
# the potential memory(4*) usage issues.

# Others global package variables which are set in the inheriting
# Collector class.
our ( $bin_model,   $bin_method, $pack_template,
      $packed_size, $window_sizes );

=head2 new

  Args       : None
  Example    :

    my $collector = Bio::EnsEMBL::XXX::Collector::FEATURE->new();
    $collector->store_windows_by_Slice($slice);

    # Where XXX is, e.g. Compara, FuncGen etc.

  Description: Simple new method to enable use of collector
               when not inherited by a descendant of
               Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor

  Returntype : Bio::EnsEMBL::XXX::Collector
  Exceptions : None
  Caller     : Collector script
  Status     : At Risk

=cut

sub new {
  return bless {}, $_[0];   # Simply blesses this class as an empty hash.

  # Do not set anything here, as will not be first in ISA for
  # FeatureAdaptors.  Hence, not guaranteed to be called.
}

=head2 new_assembly

  Args       : optional - string assembly version e.g. GRCh37
  Example    : $collector->new_assembly('GRCh37');
  Description: Getter/Setter for new assembly version which should be
               used to project only 0 wsize Collections.
  Returntype : string
  Exceptions : None
  Caller     : store_window_bins_by_Slice() or
               write_collection() in inheriting Collector class.
  Status     : At Risk

=cut

sub new_assembly {
  my ( $self, $new_asm ) = @_;

  if ( defined($new_asm) ) {
    $self->{'new_assembly'} = $new_asm;
  }

  return $self->{'new_assembly'};
}

### Setter/Getter methods for basic/mandatory config
#   Can also be set using package variables in the inheriting
#   Collector/Adaptor or run script.  Allows over-riding of defaults set
#   in Adaptor/Collector.

# Package variables used here instead of attrs to enable easy
# default config in inheriting class/script method. Provided
# for easy/standardised fetch access outside of this package
# i.e. Collectors/Adaptors

=head2 max_data_type_size

  Args       : optional - int  Maximum size of collection in bytes
  Example    : $collector->max_data_type_size($new_max_size);
  Description: Getter/Setter for max_data_type_size, default is
               currently set at in this class as 16777216 (16MB), for
               long BLOB.  This is used by the write_collection()
               method to determine when to build and store a compressed
               collection.
  Returntype : int
  Exceptions : None
  Caller     : bins_per_record() and
               write_collection() in inheriting Collector class.
  Status     : At Risk

=cut

sub max_data_type_size {
  my ( $self, $size ) = @_;

  # Validate is sensible integer

  if ( defined($size) ) {
    if ( $size !~ /^\d+$/ ) {
      throw("max_data_type_size must be a integer of bytes, not $size");
    }
    $max_data_type_size = $size;
  } elsif ( !defined($max_data_type_size) ) {
    # This should never happen as we have defaults in this module.
    throw(   'You must define a '
           . '$Bio::EnsEMBL::Utils::Collector::max_data_type_size '
           . 'or pass -max_data_type_size config' );
  }

  return $max_data_type_size;
}

=head2 max_view_width

  Args       : optional - int  Maximum width of view
  Example    : $collector->max_view_width($new_max_width);
  Description: Getter/Setter for max_view_width, default is currently
               set at in this class as 500000bp, for maximum level of
               zoom permitted by location view.
  Returntype : int
  Exceptions : None
  Caller     : general
  Status     : At Risk

=cut

sub max_view_width {
  my ( $self, $size ) = @_;

  # Validate is sensible integer

  if ( defined($size) ) {
    if ( $size !~ /^\d+$/ ) {
      throw("max_view_width must be a integer, not $size");
    }
    $max_view_width = $size;
  } elsif ( !defined $max_view_width ) {
    # This should never happen as we have defaults in this module.
    throw(   'You must define a '
           . '$Bio::EnsEMBL::Utils::Collector::max_view_width '
           . 'or pass -max_view_width config' );
  }

  return $max_view_width;
}

=head2 bin_method

  Args[0]    : optional - string name of bin method e.g. 'max_magnitude'
  Args[1]    : optional - Bio::EnsEMBL::Funcgen::Parsers::InputSet
  Example    : my $bin_method = $self->bin_method();
  Description: Getter/Setter for bin_method, default is normally set in
               the inheriting Collector class either by package variable
               or by passing a config hash via the store methods.
  Returntype : string
  Exceptions : Throws if cannot set by package variable
  Caller     : general
  Status     : At Risk

=cut

sub bin_method {
  my ( $self, $bmethod, $config ) = @_;

  if ( defined($bmethod) ) {
    $bin_method = $bmethod;
  }

  if ( !defined($bin_method) ) {
    throw(   'You must define a '
           . '$Bio::EnsEMBL::Utils::Collector::bin_method '
           . 'or pass -bin_method config' );
  }

  if ( !$self->can( "_calculate_" . $bin_method ) ) {
    throw("$bin_method is not a valid/available binning method");
  }

  my $set_up_method = "_set_up_" . $bin_method;
  if ( $self->can($set_up_method) ) {
    $self->$set_up_method($config);
  }

  return $bin_method;
}

=head2 bin_model

  Args       : optional - string bin model e.g. SIMPLE or COMPLEX
  Example    : my $bin_model = $self->bin_model;
  Description: Getter/Setter for bin_model, default should be set in
               inheriting Collector class.  Currently only supports
               'SIMPLE' bin model.
  Returntype : string
  Exceptions : Throws if bin_model is not SIMPLE
  Caller     : general
  Status     : At Risk

=cut

sub bin_model {
  my ( $self, $bmodel ) = @_;

  if ( defined($bmodel) ) {
    $bin_model = $bmodel;
  }

  if ( !defined($bin_model) ) {
    throw(   'You must define a '
           . '$Bio::EnsEMBL::Utils::Collector::bin_model '
           . 'or pass -bin_model config' );
  }

  if ( $bin_model ne 'SIMPLE' ) {
    throw(   'Bio::EnsEMBL::Utils::Collector does not yet support '
           . 'non-SIMPLE bin models' );
  }

  return $bin_model;
}


=head2 window_sizes

  Args       : optional - arrayref of window sizes
  Example    :

    foreach my $wsize ( @{ $collector->window_sizes } )
    {    # Do some collecting
    }

  Description: Getter/Setter for window_sizes.  Default should be set
               in inheriting Collector (if the config is dynamic),
               FeatureAdaptor class or script using package variable or
               this method.
               NOTE: Redefining these may cause a mismatch with the
               table partition definition.
  Returntype : arrayref of ints
  Exceptions : Throws if cannot set a valid array of int window sizes
  Caller     : general
  Status     : At Risk - rename bin_sizes?

=cut

sub window_sizes {
  my ( $self, $sizes ) = @_;

  if ( defined($sizes) ) {
    $window_sizes = $sizes;
  }

  if ( !(    ref($window_sizes)
          && ( ref($window_sizes) eq 'ARRAY' )
          && ( scalar(@$window_sizes) > 0 ) ) )
  {
    throw('Must pass -windows_sizes in the config '
        . 'or define $Bio::EnsEMBL::Utils::Collector::window_sizes '
        . 'in your Collector as an array ref of integer window_sizes' );
  }

  return $window_sizes;
}




=head2 has_window_size

  Args       : int - window size to validate
  Example    : if( $collector->has_window_size('30') ){
                   #Do something wrt to 30bp window size
               }

  Description: Simple utility method to validate whether this Collector
               has a given window_size
  Returntype : Boolean
  Exceptions : Throws if window size not specified
  Caller     : general
  Status     : At Risk

=cut


sub has_window_size{
  my ( $self, $size ) = @_;

  if(! defined $size){
	throw('You must pass a window size to validate');
  }

  return grep(/$size/, @$window_sizes); 
}




### Getter/Setters for BLOB collection config
# NOTE: Overriding the defaults here may cause a mismatch when the data
# is retrieved.

=head2 pack_template

  Args       : optional - string perl 'pack' template 
  Example    : $self->pack_template('v');
  Description: Getter/Setter for pack_template. Default should be set
               in inheriting Collector (if the config is dynamic),
               FeatureAdaptor class or script using package variable or
               this method.
  Returntype : string
  Exceptions : Throws if cannot set pack_template from package variable
  Caller     : FeatureAdaptor::_obj_from_sth
  Status     : At Risk

=cut

sub pack_template {
  my ( $self, $template ) = @_;

  if ( defined($template) ) {
    $pack_template = $template;
  }

  if ( !defined($pack_template) ) {
    throw(   'Must pass a per score '
           . '-pack_template in the config '
           . 'or define $Bio::EnsEMBL::Utils::Collector::pack_template '
           . 'in your Collector' );
  }

  return $pack_template;
}

=head2 packed_size

  Args       : optional - int size of perl 'pack' template in bytes
  Example    : $self->packed_size(2);
  Description: Getter/Setter for packed_size.  Default should be set
               in inheriting Collector (if the config is dynamic),
               FeatureAdaptor class or script using package variable or
               this method.
  Returntype : string
  Exceptions : Throws if cannot set pack_template from pacakge variable
  Caller     : current_packed_size() and
               FeatureAdaptor::_obj_from_sth()
  Status     : At Risk

=cut

sub packed_size {
  my ( $self, $size ) = @_;

  if ( defined($size) ) {
    $packed_size = $size;
  }

  if ( !defined($packed_size) ) {
    throw(   'Must pass -packed_size(wrt to pack_template) config '
           . 'or define $Bio::EnsEMBL::Utils::Collector::packed_size '
           . 'in your Collector' );
  }

  if ( $packed_size !~ /^\d+$/ ) {
    throw(   "$packed_size is not an integer, "
           . "must pass a size integer for packed_size "
           . "which specifies size of pack_template:\t"
           . $pack_template );
  }

  return $packed_size;
}

=head2 bins_per_record

  Example    : my $bin_per_records = $self->bin_per_record
  Description: Simple method to calculate the max number of bins
               allowed per record given the current config.
  Returntype : int
  Exceptions : None
  Caller     :
  Status     : At Risk

=cut

sub bins_per_record {
  return int( $max_data_type_size/$packed_size );
}


=head2 current_packed_size

  Arg[0]     : int - window size
  Example    : my $cps = $self->current_packed_size($wsize);
  Description: Simple method to calculate the max number of bins
               allowed per record given the current config.
  Returntype : int
  Exceptions : None
  Caller     :
  Status     : At Risk

=cut

sub current_packed_size {
  my ( $self, $wsize ) = @_;
  return ( scalar( @{ $self->score_cache($wsize) } )*$packed_size );
}


=head2 score_cache

  Arg[0]     : int - window size
  Example    : my $cps = $self->current_packed_size($wsize);
  Description: Handles caching of bin scores for each window size
  Returntype : arrayref
  Exceptions : Throws if no window size defined
  Caller     : current_packed_size() and store_collection()
  methods
  Status     : At Risk

=cut

sub score_cache {
  my ( $self, $wsize, $scores ) = @_;

  if ( !defined($wsize) ) {
    throw('Must pass a window size argument');
  }

  $self->{'score_cache'}{$wsize} ||= [];

  if ( defined($scores) ) {
    push( @{ $self->{'score_cache'}{$wsize} }, @{$scores} );
  }

  return $self->{'score_cache'}{$wsize};
}


=head2 collection_start

  Arg[0]     : int - window_size
  Arg[1]     : optional int - seq_region_start
  Example    : my $coll_start->(150);
  Description: Getter/Setter collection seq_region_start
  Returntype : int
  Exceptions : Throws if no window size defined
  Caller     : store_window_bin_by_Slice() and write_collection()
  Status     : At Risk

=cut

sub collection_start {
  my ( $self, $wsize, $sr_start ) = @_;

  if ( !defined($wsize) ) {
    throw('Must pass a window size argument');
  }

  if ( defined($sr_start) ) {
    $self->{'collection_start'}{$wsize} = $sr_start;
  }

  return $self->{'collection_start'}{$wsize};
}


=head2 collection_end

  Arg[0]     : int - window_size
  Arg[1]     : optional int - seq_region_end
  Example    : my $coll_end->(150);
  Description: Getter/Setter collection seq_region_end
  Returntype : int
  Exceptions : Throws if no window size defined
  Caller     : inheriting Collector write_collection method
  Status     : At Risk

=cut

sub collection_end{
  my ($self, $wsize, $sr_end) = @_;
  throw('Must pass a window size argument') if ! defined $wsize;

  if(defined $sr_end){
	$self->{'collection_end'}{$wsize} = $sr_end;
  }
  else{
	return $self->{'collection_end'}{$wsize};
  }
}


=head2 collection_strand

  Arg[0]     : int - window_size
  Arg[1]     : optional int - seq_region_strand
  Example    : my $coll_start->(0);
  Description: Getter/Setter collection seq_region_strand
  Returntype : int
  Exceptions : Throws if no window size defined
  Caller     : inheriting Collector write_collection method
  Status     : At Risk - Collections are currently strandless

=cut

sub collection_strand {
  my ( $self, $wsize, $strand ) = @_;

  if ( !defined($wsize) ) {
    throw('Must pass a window size argument');
  }

  if ( defined $strand ) {
    $self->{'collection_strand'}{$wsize} = $strand;
  }

  return $self->{'collection_strand'}{$wsize};
}


### Here follows the actual working methods

=head2 _get_Slice_chunks

  Description: Defines the optimal set of slice chunks to use for
               generating collections such that redundant fetches
               are minimized.
  Returntype : hashref of window_size chunk size pairs
  Exceptions : Throws if no window sizes or max_view_width defined
  Caller     : store_window_bin_by_Slice()
  Status     : At Risk

=cut

sub _get_Slice_chunks {
  my $self = shift;

  if ( !defined($window_sizes) || !defined($max_view_width) ) {
    throw(   'You must pass both a window_size array ref '
           . 'and max_view_width arguments' );
  }

  if ( !defined( $self->{'_slice_chunks'} ) ) {
    # Calulate sensible slice length based on window sizes
    my @wsizes = sort { $a <=> $b } @$window_sizes;

    # Handle calculating only 0 wsize
    if ( scalar(@wsizes) == 1
         && $wsizes[0] == 0 )
    {
      return { $max_view_width => [0] };
    }

    my $multiplier = int( $max_view_width/$wsizes[$#wsizes] );

    my $chunk_length  = $multiplier*$wsizes[$#wsizes];
    my $not_divisible = 1;

    my %chunk_windows;   # Registry of chunk lengths to run with windows
    my %workable_chunks = map { $_ => {} } @wsizes;

    # get rid of natural resolution as this will always work
    delete $workable_chunks{'0'};

    while ( $not_divisible && $chunk_length != 0 ) {
      $not_divisible = 0;

      foreach my $wsize (@wsizes) {
        if ( $wsize == 0 ) {
          # Special wsize for normal data
          next;
        }

        # Set not divisible if modulus is true
        if ( $chunk_length % $wsize ) {
          $not_divisible = 1;
        } else {
          $workable_chunks{$wsize}{$chunk_length} = [];
        }
      }

      # Gradually shrink the length until we find a workable slice
      # length for all windows.
      if ($not_divisible) {
        $chunk_length -= $wsizes[$#wsizes];
      }
    }

    my %chunk_sets;

    if ( $chunk_length == 0 ) {
      print "Could not find chunk length "
        . "for all window sizes, "
        . "attempting to subset windows "
        . "using alternate slice length\n";

      foreach my $wsize ( keys(%workable_chunks) ) {
        # Loop through windows, seeing if they are workable in the other
        # windows.

        foreach my $chunk ( keys( %{ $workable_chunks{$wsize} } ) ) {

          foreach my $other_wsize ( keys %workable_chunks ) {
            next if $wsize == $other_wsize;

            if ( exists( $workable_chunks{$other_wsize}{$chunk} ) ) {
              # only push it onto the other wsize, as we will do the
              # reverse later
              $chunk_sets{$chunk}{$wsize} = undef;
            }
          }
        }
      }

      # %chunk_sets represents co-occurence of wsizes with repect to
      # chunks.  Take the set which has the most windows and the longest
      # chunk.  Then get the largest which handles the rest.

      # define possible set lengths
      my $i = 0;
      my %set_lengths;

      map { $set_lengths{$i} = []; $i++ } @wsizes;

      # get rid of natural resolution as this will always work
      delete $set_lengths{'0'};

      # Store chunks lengths for each set size
      foreach my $chunk ( keys(%chunk_sets) ) {
        my $set_size = scalar( values( %{ $chunk_sets{$chunk} } ) );
        push( @{ $set_lengths{$set_size} }, $chunk );
      }

      # Get the biggest set with the longest length;

      # Scalar here as we are disregarding natural resolution of 0 in
      # loop.
      my $largest_size      = scalar(@wsizes);
      my $found_largest_set = 0;

      while ( !$found_largest_set ) {
        $largest_size--;

        if ( scalar( @{ $set_lengths{$largest_size} } ) > 0 ) {
          $found_largest_set = 1;
        }
      }

      my ($largest_chunk) =
        sort { $b <=> $a } @{ $set_lengths{$largest_size} };

      my @largest_windows = keys %{ $chunk_sets{$largest_chunk} };
      @{ $chunk_windows{$largest_chunk} } = @largest_windows;

      print "Largest chunk $largest_chunk($largest_size) "
        . "contains windows: @largest_windows\n";

      my %remaining_windows = map { $_ => {} } @wsizes;

      # get rid of natural resolution as this will always work
      delete $remaining_windows{'0'};

      map { delete $remaining_windows{$_} } @largest_windows;
      my $remaining_set_size = scalar( keys(%remaining_windows) );

      # Use array here for practicality, would need to maintain hash if
      # we need to iterate.
      my @rwindows = keys(%remaining_windows);

      # Could be one window, but this will not be in the co-occurence
      # hash %chunk_sets.
      my $next_chunk;

      if ( scalar(@rwindows) == 1 ) {
        my ($last_window) = @rwindows;
        # Find a suitably large chunk for this one window.
        $multiplier = int( 500000/$last_window );
        $next_chunk = $multiplier*$last_window;
      } else {

        foreach my $chunk ( sort { $b <=> $a }
                            @{ $set_lengths{$remaining_set_size} } )
        {
          my $seen_count = 0;

          foreach my $rwindow (@rwindows) {
            if ( grep /$rwindow/,
                 ( values( %{ $chunk_sets{$chunk} } ) ) )
            {
              $seen_count++;
            }
          }

          if ( $seen_count == $remaining_set_size ) {
            $next_chunk = $chunk;
            last;
          }
        }

      }

      @{ $chunk_windows{$next_chunk} } = @rwindows;

      if ( defined($next_chunk) ) {
        print "Found next chunk length $next_chunk "
          . "contains remaining windows:\t@rwindows\n";
      } else {
        warn "Need to write iterative method for set definition";
        throw(   'Could not find workable slice length '
               . 'for remaining windows: '
               . join( ', ', @rwindows ) );
      }
    } else {
      @{ $chunk_windows{$chunk_length} } = keys(%workable_chunks);
      print "Found workable chunk length $chunk_length  "
        . "for all window sizes:\t"
        . join( ' ', @{ $chunk_windows{$chunk_length} } ) . "\n";
    }

    $self->{'_slice_chunks'} = \%chunk_windows;
  } ## end if ( !defined( $self->...))

  return $self->{'_slice_chunks'};
} ## end sub _get_Slice_chunks




=head2 set_config

  Arg[0]     : optional hash - parameter hash(see above methods for more info):

        WINDOW_SIZES        => array ref - subset of defined window
                               sizes
        BIN_METHOD          => string
        MAX_VIEW_WIDTH      => int
        MAX_DATA_TYPE_SIZE  => int
        PACK_TEMPLATE       => string
        PACKED_SIZE         => int
        BIN_MODEL           => string
        NEW_ASSEMBLY        => string
        METHOD_CONFIG       => hash of method specific config params
        SKIP_ZERO_WINDOW    => boolean - skips generation of 0 wsize
                               this is used if already generated
                               from an assembly projection.

               NOTE: Over-riding any of the default config may cause
               problems when storing or retrieving Collection data,
               except sub sets of default window sizes.

  Description: This method replaces the constructor as new will not be
               called for Adaptor based Collectors.
               Separating this from the store method is currently 
               redundant as jobs are normally submitetd in Slice based 
               jobs. However, this will be required if the store method 
               is further seaprated into fetch/generate and store methods
  Returntype : None
  Exceptions : Throws if no window sizes or max_view_width defined
  Caller     : Inheritor Collector e.g. Bio::EnsEMBL::Funcgen:Collector::ResultFeature
               or script.
  Status     : At Risk

=cut

sub set_config {
  my ( $self, %config ) = @_;

  my ( $wsizes,       $bmethod,  $mv_width,
       $md_type_size, $template, $psize,
       $bmodel,       $new_assm, $skip_zero_window,
       $method_config )
    = rearrange( [ 'WINDOW_SIZES',     'BIN_METHOD',
                   'MAX_VIEW_WIDTH',   'MAX_DATA_TYPE_SIZE',
                   'PACK_TEMPLATE',    'PACKED_SIZE',
                   'BIN_MODEL',        'NEW_ASSEMBLY',
                   'SKIP_ZERO_WINDOW', 'METHOD_CONFIG' ],
                 %config );

  ### VAILDATE/SET VARS/CONFIG

  # Attrs used in this method
  $self->bin_method( $bmethod, $method_config );
  $self->bin_model($bmodel);
  $self->window_sizes($wsizes);

  # Set to undef if we have empty array? To change this we need to
  # pass the config hash -window_sizes conditionally
  # This currently overwrite the defaults!
  #  if ( ref($window_sizes) eq 'ARRAY'
  #       && scalar( @{$window_sizes} ) == 0 )
  #  {
  #    $window_sizes = undef;
  #  }

  # Attrs used in other (store) methods
  $self->pack_template($template);
  $self->packed_size($psize);
  $self->max_data_type_size($md_type_size);
  $self->max_view_width($mv_width);

  # Other vars
  $self->new_assembly($new_assm);
  $self->{'_only_natural'} = 0;
  $self->{'_store_natural'} = grep /^0$/, @$window_sizes;

  ### Set window_sizes

  if ( $self->new_assembly() ) {
    print "Assembly projection may cause problems "
      . "for large Collections, "
      . "defaulting to window_sizes = (0)\n";


	if ( $skip_zero_window ) {
	  throw(   "You cannot -skip_zero_window or "
			   . "omit 0 from -window_sizes "
			   . "when projecting to a new assembly($new_assm) "
			   . "which should only be generated using window_size=0" );
	}
	



    # Then build the bins on the projected 0 level single Features

    # Test we haven't explicity set window_sizes to be something else
    if ( defined($wsizes)
         && !( scalar(@$wsizes) == 1 && $wsizes->[0] == 0 ) )
    {
      throw(   "You have set window_sizes config "
             . "which are not safe when projecting to "
             . "a new assembly($new_assm), "
             . "please omit window_sizes config or set to 0" );
    }

    $self->window_sizes( [0] );
  } else {

    if ( $wsizes && $skip_zero_window && 
		 ( grep /^0$/, @$wsizes )) {
	  #Only test passed params not default config

      throw(   "You have specied skip_zero_window "
             . "and window_size 0 in your parameters, "
             . "please remove one of these" );
    } 
	elsif ( defined($window_sizes) && !grep /^0$/, @$window_sizes ) {
      $skip_zero_window = 1;
      # re-add 0 window as we need this to build the collections
	  # see ...
      unshift( @{$window_sizes}, 0 );
    }
  }

 
  if ( $self->{'_store_natural'} && scalar( @{$window_sizes} ) == 1 ) {
    $self->{'_only_natural'} = 1;
  }
  if ($skip_zero_window) {
    $self->{'_store_natural'} = 0;
  }

  return;
} ## end sub set_config

=head2 store_window_bins_by_Slice

  Arg[0]     : Bio::EnsEMBL:Slice
  Example    : $collector->store_window_bins_by_Slice($slice);
  Description: This is the main run method, it loops through
               optimal slice chunks from _define_window_chunks,
               calls _bin_features_by_Slice as appropriate and
               calls write_collection in the inheriting Collector
               class/script.
  Returntype : None
  Exceptions : Throws if Bio::EnsEMBL::Slice is not defined
  Caller     : store methods in inheriting Collector class/script
  Status     : At Risk

=cut

sub store_window_bins_by_Slice {
  my ( $self, $slice ) = @_;

  warn "Need to be careful here "
    . "about cleaning start end strand caches between "
    . "serially run slices";

  if ( !(    defined($slice)
          && ref($slice)
          && $slice->isa('Bio::EnsEMBL::Slice') ) )
  {
    throw('You must pass a valid Bio::EnsEMBL::Slice');
  }

  # Rollback previously stored features.
  # Change 'can' to empty method stubb with pod ???
  if ( $self->can('rollback_Features_by_Slice') ) {
    $self->rollback_Features_by_Slice($slice);
  } else {
    warn ref($self)
      . " cannot rollback_Features_by_Slice. "
      . "This may result in storage failure "
      . "or duplicate Collections if there is pre-existing data";
  }

  ### PROCESS CHUNKS
  my %chunk_windows = %{ $self->_get_Slice_chunks };
  my (%counts);
  my $store_natural = $self->{'_store_natural'};
  my $only_natural  = $self->{'_only_natural'};
  $counts{0} = 0;    # Set natural res count to 0
  my $slice_end       = $slice->end;
  my $orig_start      = $slice->start;
  my $region          = $slice->coord_system_name;
  my $version         = $slice->coord_system->version;
  my $seq_region_name = $slice->seq_region_name;
  my $strand          = $slice->strand;

  # Warn if this is not a full slice.  Version needed in case we are
  # projecting from a non-default version slice
  my $full_slice =
    $slice->adaptor->fetch_by_region( $region, $seq_region_name, undef,
                                      undef, undef, $version );

  if (    ( $full_slice->start() != $orig_start )
       || ( $full_slice->end() != $slice_end ) )
  {
    warn "Generating collections using sub-Slices "
      . "can result in data issues/artifacts";
    # Last chunk might not be the correct window length.  Test
    # slices less than chunk length can cause failures in
    # _bin_features_by_window_sizes others?
  }

  # Set the initial collection_start to orig_start.  This is not the
  # case for 0 wsize where it must always be the true feature start.
  for my $wsize (@$window_sizes) {
    if ( $wsize == 0 ) { next }
    $self->collection_start( $wsize, $orig_start );

    # Also reset collection end and score cache in case we are running
    # serially.
    $self->{collection_end}{$wsize} = undef;
    $self->{'score_cache'}{$wsize} = [];
  }

  my $first_chunk_length = 1;

  foreach my $chunk_length ( sort keys %chunk_windows ) {
    print "Processing windows "
      . join( ', ', @{ $chunk_windows{$chunk_length} } )
      . " with chunk length $chunk_length\n";

    # Set window counts to 0
    map $counts{$_} = 0, @{ $chunk_windows{$chunk_length} };

    # May need to reset flat file parser handle or other caches via
    # inheriting Collector
    if ( !$first_chunk_length ) {
      # Change 'can' to empty method stubb with pod???
      if ( $self->can('reinitialise_input') ) {
        $self->reinitialise_input();
      }
    }

    $first_chunk_length = 0;

    # Now walk through slice using slice length chunks and build all
    # windows in each chunk.
    my $in_slice  = 1;
    my $start_adj = 0;
    my ( $sub_slice, $sub_end, $features, $bins );
    my $sub_start    = 1;
    my $slice_length = $slice->length();

    # Always create in local coords for fetch
    # Then change to seq_region coords for store if required

    while ($in_slice) {
      $sub_start += $start_adj;
      $sub_end = $sub_start + $chunk_length - 1;

      if ( $sub_end >= $slice_length ) {
        # Surplus bins are removed in store/write_collection in caller
        $in_slice = 0;
      }

      $sub_slice =
        $slice->adaptor->fetch_by_region( $region, $seq_region_name,
                                          $sub_start + $orig_start - 1,
                                          $sub_end + $orig_start - 1,
                                          $strand, $version );

      # Can't subslice as this will not clip if we go over the length of
      # the slice, unlike normal slice fetching.  Will clipping the end
      # to the slice end cause any problems here?  How will this affect
      # bin clipping?

      ### Grab features and shift chunk coords
      $features = $self->get_Features_by_Slice($sub_slice);

      # warn "Binning "
      #   . scalar(@$features)
      #   . " Features for chunk length $chunk_length, on Slice "
      #   . $sub_slice->name;

      if ( ( @{$features} )
         && ref( $features->[0] ) =~ /Bio::EnsEMBL::Utils::Collection/ )
      {
        # Would need to create base module with generic methods:
        # window_size, ...

      # Check that the returned feature/collections support window_size.
      # All Collections should be able to

        if ( $features->[0]->can('window_size') ) {
          if ( $features->[0]->window_size != 0 ) {
            throw(   "You are trying to generated Collections from "
                   . "a non-zero window sized Collection:\t"
                   . $features->[1]->{'window_size'} );
          }

        # This should never happen
        # if ( !$skip_zero_window ) {
        #   throw( 'You have retrieved data from a Collection '
        #        . 'which without using -skip_zero_window '
        #        . 'i.e. you are trying to generate overwrite '
        #        . 'the data you are generating the Collections from' );
        # }

        } else {
          throw(   'Something is wrong, '
                 . 'the Collection you have retrieved '
                 . 'does not support the method window_size' );
        }
      } ## end if ( ( @{$features} ) ...)

      # Set collection start here for 0 window_size
      if (    @{$features}
           && $store_natural
           && !defined( $self->collection_start(0) ) )
      {
        $self->collection_start( 0,
                                 $features->[0]->start + $sub_start );
      }

      if ($in_slice) {
        $start_adj = $chunk_length;
      }

      # Collect features into wsize bins
      if ( !$only_natural ) {
        # Get hashref of wsize=>bin array pairs
        $bins =
          $self->_bin_features_by_Slice_window_sizes(
                         -slice        => $sub_slice,
                         -window_sizes => $chunk_windows{$chunk_length},
                         -features     => $features, );
      }

      # Handle 0 wsize
      if ($store_natural) {
        foreach my $feature ( @{$features} ) {
          $counts{0}++;

          if ( $bin_model eq 'SIMPLE' ) {
            $self->collection_start( 0, $feature->start + $sub_start );

            $self->write_collection(
              0,
              $slice,    # Pass Slice to sub-slice when storing
              $feature->end + $sub_start,
              $feature->strand,   # Need to pass strand for 0 resolution
              $feature->scores, );
          }
        }

        print "Window size 0 (natural resolution) has "
          . scalar( @{$features} )
          . " feature bins for:\t"
          . $sub_slice->name . "\n";
      }

      # Now store collections for wsizes >0
      my $num_bins;

      foreach my $wsize ( sort keys( %{$bins} ) ) {
        $num_bins = scalar( @{ $bins->{$wsize} } );
        $counts{$wsize} += $num_bins;

        if ( $bin_model eq 'SIMPLE' ) {
          $self->write_collection(
            $wsize,
            $slice,
            #$sub_start,
            $sub_end,
            $slice->strand,    # This is most likely 1!
             # Override this woth 0 in descendant Collector if required.
            $bins->{$wsize}, );

        } else {
          throw(   'Bio::EnsEMBL::Utils::Collector '
                 . 'does not yet support non-SIMPLE bin models' );
          # i.e. More than one score
        }
      }
    } ## end while ($in_slice)

    # Turn off storing of natural resolution for next chunk length sets
    $store_natural = 0;
  } ## end foreach my $chunk_length ( ...)

  # Write last collections for each wsize

  foreach my $wsize (@$window_sizes) {

    if (    ( $wsize == 0 && !$store_natural )
         || ( $wsize != 0 && $only_natural ) )
    {
      next;
    }

    print "Writing final $wsize window_size collection, "
      . "this may result in slightly different "
      . "bin numbers from counts due to removing "
      . "overhanging bins past end of slice\n";
    $self->write_collection( $wsize, $slice );
  }

  # Print some counts
  foreach my $wsize ( sort ( keys %counts ) ) {
    print "Generated "
      . $counts{$wsize}
      . " bins for window size $wsize for "
      . $slice->name . "\n";
    # Some may have failed to store if we are projecting to a new
    # assembly.
  }

  return;
} ## end sub store_window_bins_by_Slice

=head2 _bin_features_by_Slice_window_sizes

  Args[0]    : Bio::EnsEMBL::Slice
  Args[1]    : ARRAYREF of window sizes
  Args[2]    : ARRAYREF of features with start and end method
               e.g. Bio::EnsEMBL::Features
  Example    :

      $bins =
        $self->_bin_features_by_window_sizes(
                         -slice        => $slice,
                         -window_sizes => $chunk_windows{$chunk_length},
                         -features     => $features, );

  Description: Bins feature scores for a given list of window sizes and
               predefined method.
  Returntype : HASHREF of scores per bin per window size
  Exceptions : None
  Caller     : store_window_bins_by_Slice
  Status     : At Risk

=cut

sub _bin_features_by_Slice_window_sizes {
  my ( $self, @args ) = @_;

  my ( $slice, $wsizes, $features ) =
    rearrange( [ 'SLICE', 'WINDOW_SIZES', 'FEATURES' ], @args );

  # Generate these once in caller?
  my $calc_method = '_calculate_' . $bin_method;
  my $post_method = '_post_process_' . $bin_method;

  # Do this conditional on the Collection type i.e. is
  # collection seq_region blob then no else yes Would need
  # $Bio::EnsEMBL::Utils::Collector::collection_format=BLOB|STANDARD
  # if ( !defined($features) || !@{$features} ) { return {} }

  # Set up some hashes to store data by window_size
  my ( %bins, %nbins, %bin_counts );
  my $slice_start  = $slice->start();
  my $slice_length = $slice->length();

  # Set up some bin data for the windows
  foreach my $wsize (@$wsizes) {
    $nbins{$wsize} = int( $slice_length/$wsize );    # int rounds down
         # nbins is index of the bin not the 'number'
         # Unless $slice_length is a multiple!
    if ( !( $slice_length % $wsize ) ) { $nbins{$wsize}-- }

    # Create default bins with 0
    $bins{$wsize} = [];
    map { $bins{$wsize}->[$_] = 0 } ( 0 .. $nbins{$wsize} );

    # Set bin counts to 0 for each bin
    $bin_counts{$wsize} = [];

    # This is adding an undef to the start of the array!?
    map { $bin_counts{$wsize}->[ ($_) ] = 0 } @{ $bins{$wsize} };

    foreach my $bin ( @{ $bins{$wsize} } ) {
      $bin_counts{$wsize}->[$bin] = 0;
    }
  }

  my $feature_index = 0;
  my ( $bin_index, @bin_masks );

  foreach my $feature ( @{$features} ) {
    # Set up the bins for each window size

    foreach my $wsize (@$wsizes) {
      my $start_bin = int( ( $feature->start )/$wsize );
      my $end_bin   = int( ( $feature->end )/$wsize );

      if ( $end_bin > $nbins{$wsize} ) {
        $end_bin = $nbins{$wsize};
      }

      $self->$calc_method( $feature, $start_bin, $end_bin,
                           $wsize,   \%bins,     \%bin_counts );
    }

  }

  # Now do post processing of bins if required
  if ( $self->can($post_method) ) {
    $self->$post_method( \%bins, \%bin_counts );
  }

  return \%bins;
} ## end sub _bin_features_by_Slice_window_sizes
# end sub _bin_features_by_Slice


### Here follows the bin methods
# These may also be defined in the inheriting Collector class.  No tests
# as these are internal and require speed.


=head2 _calculate_count

  Args[0]    : feature e.g. Bio::EnsEMBL::Feature
  Args[1]    : int - start bin
  Args[2]    : int - end bin
  Args[3]    : int - window_size
  Args[4]    : hashref - score bins
  Example    : $self->$calc_method
  Description: Adds count to bins which this feature overlaps
  Returntype : None
  Exceptions : None
  Caller     : _bin_features_by_window_sizes
  Status     : At Risk

=cut

sub _calculate_count {
  my ( $self, $feature, $start_bin, $end_bin, $wsize, $bins_ref ) = @_;

  my $bin_index;

  for ( $bin_index = $start_bin; $bin_index <= $end_bin; ++$bin_index )
  {
    $bins_ref->{$wsize}->[$bin_index]++;
  }

  return;
}


=head2 _calculate_average_score

  Args[0]    : feature e.g. Bio::EnsEMBL::Feature
  Args[1]    : int - start bin
  Args[2]    : int - end bin
  Args[3]    : int - window_size
  Args[4]    : hashref - score bins
  Example    : $self->$calc_method
  Description: Adds score to bins which this feature overlaps
  Returntype : None
  Exceptions : None
  Caller     : _bin_features_by_window_sizes
  Status     : At Risk

=cut


sub _calculate_average_score {
  my ( $self, $feature, $start_bin, $end_bin, $wsize, $bins_ref,
       $bin_counts_ref )
    = @_;

  # This is simple an average of all the scores for features which
  # overlap this bin.  No weighting with respect to the bin or the
  # feature.

  my $score = $self->get_score_by_Feature($feature);

  for ( my $bin_index = $start_bin;
        $bin_index <= $end_bin;
        ++$bin_index )
  {
    # We should really push onto array here so we can have median or
    # mean.

    $bins_ref->{$wsize}->[$bin_index] += $score;
    $bin_counts_ref->{$wsize}->[$bin_index]++;
  }

  return;
}


=head2 _post_process_average_score

  Args[0]    : hashref - score bins
  Args[1]    : hashref - count bins
  Example    : $self->$post_method
  Description: Post processes bins to calculate average score
  Returntype : None
  Exceptions : None
  Caller     : _bin_features_by_window_sizes
  Status     : At Risk

=cut	

sub _post_process_average_score {
  my ( $self, $bins_ref, $bin_counts_ref ) = @_;

  foreach my $wsize ( keys %{$bins_ref} ) {
    foreach my $bin_index ( 0 .. $#{ $bins_ref->{$wsize} } ) {

      if ( $bin_counts_ref->{$wsize}->[$bin_index] ) {
        $bins_ref->{$wsize}->[$bin_index] /=
          $bin_counts_ref->{$wsize}->[$bin_index];
      }

    }
  }

  return;
}


=head2 _calculate_max_magnitude

  Args[0]    : feature e.g. Bio::EnsEMBL::Feature
  Args[1]    : int - start bin
  Args[2]    : int - end bin
  Args[3]    : int - window_size
  Args[4]    : hashref - score bins
  Example    : $self->$calc_method
  Description: Sets max +/-ve scores for bins which this feature overlaps
  Returntype : None
  Exceptions : None
  Caller     : _bin_features_by_window_sizes
  Status     : At Risk

=cut

sub _calculate_max_magnitude {
  my ( $self, $feature, $start_bin, $end_bin, $wsize, $bins_ref ) = @_;

  my $score = $self->get_score_by_Feature($feature);

  # Max magnitude
  # Take the highest value +ve or -ve score
  for ( my $bin_index = $start_bin;
        $bin_index <= $end_bin;
        ++$bin_index )
  {

    # We really need to capture the lowest -ve and higest +ve scores
    # here and post process to pick between them.

    $bins_ref->{$wsize}->[$bin_index] ||= [ 0, 0 ];    #-ve, +ve

    if ( $score < $bins_ref->{$wsize}->[$bin_index]->[0] ) {
      $bins_ref->{$wsize}->[$bin_index]->[0] = $score;
    } elsif ( $score > $bins_ref->{$wsize}->[$bin_index][1] ) {
      $bins_ref->{$wsize}->[$bin_index]->[1] = $score;
    }
  }

  return;
} ## end sub _calculate_max_magnitude


=head2 _post_process_max_magnitude

  Args[0]    : hashref - score bins
  Args[1]    : hashref - count bins
  Example    : $self->$post_method
  Description: Post processes bins to pick largest +ve or -ve score
  Returntype : None
  Exceptions : None
  Caller     : _bin_features_by_window_sizes
  Status     : At Risk

=cut	

sub _post_process_max_magnitude {
  my ( $self, $bins_ref ) = @_;

  # Take the highest value +ve or -ve score

  foreach my $wsize ( keys %{$bins_ref} ) {
    foreach my $bin_index ( 0 .. $#{ $bins_ref->{$wsize} } ) {

      # Have potential for no listref in a given bin

      # default value if we haven't seen anything is 0
      # Actually want an array of -ve +ve values

      if ( $bins_ref->{$wsize}->[$bin_index] ) {
        my $tmp_minus = -$bins_ref->{$wsize}->[$bin_index]->[0];

        if ( $tmp_minus > $bins_ref->{$wsize}->[$bin_index]->[1] ) {
          $bins_ref->{$wsize}->[$bin_index] =
            $bins_ref->{$wsize}->[$bin_index]->[0];
        } else {
          $bins_ref->{$wsize}->[$bin_index] =
            $bins_ref->{$wsize}->[$bin_index]->[1];
        }

      }

    }
  }
  return;
} ## end sub _post_process_max_magnitude


=head2 _calculate_RPKM

  Args[0]    : feature e.g. Bio::EnsEMBL::Feature
  Args[1]    : int - start bin
  Args[2]    : int - end bin
  Args[3]    : int - window_size
  Args[4]    : hashref - score bins
  Example    : $self->$calc_method
  Description: Stores counts to calculate Read Per Kb per Million(RPKM)
  Returntype : None
  Exceptions : None
  Caller     : _bin_features_by_window_sizes
  Status     : At Risk

=cut

sub _calculate_RPKM {
  my ( $self, $feature, $start_bin, $end_bin, $wsize, $bins_ref ) = @_;

  $self->_calculate_count( $feature, $start_bin, $end_bin,
                           $wsize,   $bins_ref );

  return;
}


=head2 _post_process_RPKM

  Args[0]    : hashref - score bins
  Args[1]    : hashref - count bins
  Example    : $self->$post_method
  Description: Post processes bins to calculate average score
  Returntype : None
  Exceptions : None
  Caller     : _bin_features_by_window_sizes
  Status     : At Risk

=cut	

sub _post_process_RPKM {
  my ( $self, $bins_ref ) = @_;

  #10^9 x C / NGB
  #C = Reads overlapping bin
  #N = Total reads in the experiment
  #G = Length of bin in bps 
  #(don't really have to account for non-ref/HAPs or gender here 
  #as should be close enough, CellTypes/gender differences will be miniscule)
  #B = length of each bin

  foreach my $wsize(keys %{$bins_ref}){
	
	foreach my $bin_index(0..$#{$bins_ref->{$wsize}}){
	  $bins_ref->{$wsize}->[$bin_index] =  
		((10**9) *  
		 $bins_ref->{$wsize}->[$bin_index])/(($self->_RPKM_factor($wsize)) * $wsize);
	}
  }

  return;

}


=head2 _set_up_RPKM

  Args[0]    : hashref - method config e.g
                 {
                  DNADB         => Bio::EnsEMBL::DBSQL::DBAdaptor,
                  TOTAL_FEATURE => $total_feature_count,
                 }

  Example    : $self->$set_up_method($config);
  Description: Sets the RPKM factor
  Returntype : None
  Exceptions : Throws is required config params are not set
  Caller     : bin_method
  Status     : At Risk

=cut	


sub _set_up_RPKM{
  my ($self, $config) = @_;

  my ($dnadb, $total_features) = rearrange([ 'DNADB', 'TOTAL_FEATURES'], %{$config});
  
  #Test specifically here to notify about config hash
  if(! $total_features){
	throw("For RPKM you must pass a valid 'total_features' ".
		  "as part of the method config hash.");
  }
  
  if(! $dnadb){
	throw("For RPKM you must pass 'dnadb' as part of the method config hash.");
  }

  foreach my $wsize(@{$self->window_sizes}){
	#Should never have 0 here
	$self->_RPKM_factor($wsize, ($wsize * $total_features));	#N*G

	warn "setting $wsize RPKM factor($wsize * $total_features) to ".	
	  $self->_RPKM_factor($wsize);
  }
  
  return;
} ## end sub _set_up_RPKM


=head2 _RPKM_factor

  Args[0]    : int - RPKM factor i.e.  (Total reads in the experiment *
               Genome length)
  Example    : $self->_RPKM_factor($wsize, $factor);
  Description: Gets/Sets the RPKM factor
  Returntype : int
  Exceptions : None
  Caller     : _set_up_RPKM, _post_process_RPKM
  Status     : At Risk

=cut	

sub _RPKM_factor{
  my ($self, $wsize, $factor) = @_;

  if (! defined $wsize){
	throw("You must pass at least window_size to get or set the RPKM factor");
  }

  if(defined $factor){
	$self->{'RPKM_factor'}{$wsize} = $factor;
  }
  elsif(! exists $self->{'RPKM_factor'}{$wsize}){
	#This should never happen unless the window sizes 
	#are redefined after initialisation
	throw("You have requested an RPKM factor for a window_size".
		  " which has not been set:\t$wsize");
  }

  return $self->{'RPKM_factor'}{$wsize};
}

=head2 get_diploid_genome_length_by_gender

  Args[0]    : string - RPKM factor i.e.  (Total reads in the experiment *
               Genome length)
  Args[1]    : string - gender e.g. male or female
  Example    :

      my $glength =
        $self->get_diploid_genome_length_by_gender( $dnadb, $gender );

  Description: Gets the gender specific diploid genome length,
               including non-ref but not including haplotypes.  Only
               handles species with X/Y sex chromosomes.
  Returntype : int
  Exceptions : None
  Caller     : _set_up_RPKM, _post_process_RPKM
  Status     : At Risk - Move to and export from generic Utils Slice module???

=cut	

sub get_diploid_genome_length_by_gender {
  my ( $dnadb, $gender ) = @_;

  my %sex_chrs = ( 'Y' => 'male',
                   'X' => 'female', );

  my $dip_length = 0;

  if (!(( ref($dnadb) && $dnadb->isa('Bio::EnsEMBL::DBSQL::DBAdaptor') )
        && $dnadb->grou() eq 'core'
        && ( defined $gender && $gender =~ /(male|female)/ ) ) )
  {
    throw(   "Must provide valid "
           . "Bio::EnsEMBL::DBSQL::DBAdaptor($dnadb) and "
           . "gender ($gender) arguments" );
  }

  my @ref_slices = $dnadb->get_SliceAdaptor->fetch_all('toplevel');

  # Include non-ref(unassembled), but omit haps/lrgs(i.e. redundant)

  foreach my $slice (
     @{ $dnadb->get_SliceAdaptor->fetch_all( 'toplevel', undef, 1, 1 ) }
    )
  {
    # Include duplicated region for true diploid length

    # Skip haps/lrgs
    if ( ( $slice->coord_system->name() eq 'chromosome'
           && !$slice->is_reference() )
         || $slice->coord_system->name() eq 'lrg' )
    {
      next;
    }

    if ( exists( $sex_chrs{ $slice->seq_region_name() } ) ) {
      if ( $gender eq 'male' ) {
        $dip_length += $slice->length;
      } elsif ( $sex_chrs{ $slice->seq_region_name } eq 'male' ) {
        next;
      }
    }

    $dip_length += 2*$slice->length;
  }

  return $dip_length;
} ## end sub get_diploid_genome_length_by_gender


1;
