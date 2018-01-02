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

Bio::EnsEMBL::DBSQL::SliceAdaptor - A database aware adaptor responsible for
the creation of Slice objects.

=head1 SYNOPSIS

  use Bio::EnsEMBL::Utils::Slice qw(split_Slices);
  use Bio::EnsEMBL::Registry;

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
  );

  $slice_adaptor =
    Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "slice" );

  # get a slice on the entire chromosome X
  $chr_slice = $slice_adaptor->fetch_by_region( 'chromosome', 'X' );

  # get a slice for each clone in the database
  foreach $cln_slice ( @{ $slice_adaptor->fetch_all('clone') } ) {
    # do something with clone
  }

  # get a slice which is part of NT_004321
  $spctg_slice =
    $slice_adaptor->fetch_by_region( 'supercontig', 'NT_004321',
    200_000, 600_000 );

  # get all non-redundant slices from the highest possible coordinate
  # systems
  $slices = $slice_adaptor->fetch_all('toplevel');

  # include non-reference regions
  $slices = $slice_adaptor->fetch_all( 'toplevel', undef, 1 );

  # include duplicate regions
  $slices = $slice_adaptor->fetch_all( 'toplevel', undef, 0, 1 );

  # split up a list of slices into smaller slices
  $overlap    = 1000;
  $max_length = 1e6;
  $slices     = split_Slices( $slices, $max_length, $overlap );

  # store a list of slice names in a file
  open( FILE, ">$filename" ) or die("Could not open file $filename");
  foreach my $slice (@$slices) {
    print FILE $slice->name(), "\n";
  }
  close FILE;

  # retrieve a list of slices from a file
  open( FILE, $filename ) or die("Could not open file $filename");
  while ( $name = <FILE> ) {
    chomp($name);
    $slice = $slice_adaptor->fetch_by_name($name);
    # do something with slice
  }

=head1 DESCRIPTION

This module is responsible for fetching Slices representing genomic
regions from a database.  A Details on how slices can be used are in the
Bio::EnsEMBL::Slice module.

=head1 METHODS

=cut


package Bio::EnsEMBL::DBSQL::SliceAdaptor;
use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::CircularSlice;
use Bio::EnsEMBL::Mapper;
use Bio::EnsEMBL::LRGSlice;
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning stack_trace_dump);
use Bio::EnsEMBL::ProjectionSegment;
use Scalar::Util qw/looks_like_number/;
use Bio::EnsEMBL::Utils::Scalar qw/assert_integer/;

@ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  # use a cache which is shared and also used by the assembly
  # mapper adaptor

  my $seq_region_cache = $self->db->get_SeqRegionCache();

  $self->{'sr_name_cache'} = $seq_region_cache->{'name_cache'};
  $self->{'sr_id_cache'}   = $seq_region_cache->{'id_cache'};

  $self->{'lrg_region_test'} = undef;
  my $meta_container = $self->db->get_MetaContainer();
  my @values = $meta_container->list_value_by_key("LRG");
  if(scalar(@values) and $values[0]->[0]){
    $self->{'lrg_region_test'} = $values[0]->[0];
  }
  return $self;
}


=head2 fetch_by_region

  Arg [1]    : string $coord_system_name (optional)
               The name of the coordinate system of the slice to be created
               This may be a name of an actual coordinate system or an alias
               to a coordinate system.  Valid aliases are 'seqlevel' or
               'toplevel'.
  Arg [2]    : string $seq_region_name
               The name of the sequence region that the slice will be
               created on.
  Arg [3]    : int $start (optional, default = 1)
               The start of the slice on the sequence region
  Arg [4]    : int $end (optional, default = seq_region length)
               The end of the slice on the sequence region
  Arg [5]    : int $strand (optional, default = 1)
               The orientation of the slice on the sequence region
  Arg [6]    : string $version (optional, default = default version)
               The version of the coordinate system to use (e.g. NCBI33)
  Arg [7]    : boolean $no_fuzz (optional, default = undef (false))
               If true (non-zero), do not use "fuzzy matching" (see below).
  Example    : $slice = $slice_adaptor->fetch_by_region('chromosome', 'X');
               $slice = $slice_adaptor->fetch_by_region('clone', 'AC008066.4');
  Description: Retrieves a slice on the requested region.  At a minimum the
               name the name of the seq_region to fetch must be provided.

               If no coordinate system name is provided than a slice on the
               highest ranked coordinate system with a matching
               seq_region_name will be returned.  If a version but no
               coordinate system name is provided, the same behaviour will
               apply, but only coordinate systems of the appropriate version
               are considered.  The same applies if the 'toplevel' coordinate
               system is specified, however in this case the version is
               ignored.  The coordinate system should always be specified if
               it is known, since this is unambiguous and faster.

               Some fuzzy matching is performed if no exact match for
               the provided name is found.  This allows clones to be
               fetched even when their version is not known.  For
               example fetch_by_region('clone', 'AC008066') will
               retrieve the sequence_region with name 'AC008066.4'.

               The fuzzy matching can be turned off by setting the
               $no_fuzz argument to a true value.

               If the requested seq_region is not found in the database undef
               is returned.

  Returntype : Bio::EnsEMBL::Slice or undef
  Exceptions : throw if no seq_region_name is provided
               throw if invalid coord_system_name is provided
               throw if start > end is provided
  Caller     : general
  Status     : Stable

=cut


#
# ARNE: This subroutine needs simplification!! 
#
sub fetch_by_region {
  my ( $self, $coord_system_name, $seq_region_name, $start, $end,
       $strand, $version, $no_fuzz )
    = @_;

  assert_integer($start, 'start') if defined $start;
  assert_integer($end, 'end') if defined $end;

  if ( !defined($start) )  { $start  = 1 }
  if ( !defined($strand) ) { $strand = 1 }

  if ( !defined($seq_region_name) ) {
    throw('seq_region_name argument is required');
  }

  my $cs;
  my $csa = $self->db->get_CoordSystemAdaptor();

  if ( defined($coord_system_name) ) {
    $cs = $csa->fetch_by_name( $coord_system_name, $version );

    if ( !defined($cs) ) {
      throw( sprintf( "Unknown coordinate system:\n"
                        . "name='%s' version='%s'\n",
                      $coord_system_name, $version ) );
    }

    # fetching by toplevel is same as fetching w/o name or version
    if ( $cs->is_top_level() ) {
      $cs      = undef;
      $version = undef;
    }

  } ## end if ( defined($coord_system_name...))

  my $constraint;
  my $sql;
  my @bind_params;
  my $key;

  if ( defined($cs) ) {
    $sql = sprintf( "SELECT sr.name, sr.seq_region_id, sr.length, %d "
                      . "FROM seq_region sr ",
                    $cs->dbID() );

    $constraint = "AND sr.coord_system_id = ?";
    push( @bind_params, [ $cs->dbID(), SQL_INTEGER ] );

    $key = "$seq_region_name:" . $cs->dbID();
  } else {
    $sql =
      "SELECT sr.name, sr.seq_region_id, sr.length, cs.coord_system_id "
      . "FROM seq_region sr, coord_system cs ";

    $constraint = "AND sr.coord_system_id = cs.coord_system_id "
      . "AND cs.species_id = ? ";
    push( @bind_params, [ $self->species_id(), SQL_INTEGER ] );

    if ( defined($version) ) {
      $constraint .= "AND cs.version = ? ";
      push( @bind_params, [ $version, SQL_VARCHAR ] );
    }

    $constraint .= "ORDER BY cs.rank ASC";
  }

  # check the cache so we only go to the db if necessary
  my $length;
  my $arr;

  if ( defined($key) ) { $arr = $self->{'sr_name_cache'}->{$key} }

  if ( defined($arr) ) {
    $length = $arr->[3];
  } else {
    my $sth =
      $self->prepare( $sql . "WHERE sr.name = ? " . $constraint );

    unshift( @bind_params, [ $seq_region_name, SQL_VARCHAR ] );

    my $pos = 0;
    foreach my $param (@bind_params) {
      $sth->bind_param( ++$pos, $param->[0], $param->[1] );
    }

    $sth->execute();
    my @row = $sth->fetchrow_array();
    $sth->finish();

    unless ( @row ) {

      # try synonyms
      my $syn_sql = "select s.name, cs.name, cs.version from seq_region s join seq_region_synonym ss using (seq_region_id) join coord_system cs using (coord_system_id) where ss.synonym like ? and cs.species_id =? ";
      if (defined $coord_system_name && defined $cs) {
        $syn_sql .= "AND cs.name = '" . $coord_system_name . "' ";
      }
      if (defined $version) {
        $syn_sql .= "AND cs.version = '" . $version . "' ";
      }
      my $syn_sql_sth = $self->prepare($syn_sql);
      my $escaped_seq_region_name = $seq_region_name;
      my $escape_char = $self->dbc->db_handle->get_info(14);
      $escaped_seq_region_name =~ s/([_%])/$escape_char$1/g;
      $syn_sql_sth->bind_param(1, "$escaped_seq_region_name%", SQL_VARCHAR);
      $syn_sql_sth->bind_param(2, $self->species_id(), SQL_INTEGER);
      $syn_sql_sth->execute();
      my ($new_name, $new_coord_system, $new_version);
      $syn_sql_sth->bind_columns( \$new_name, \$new_coord_system, \$new_version);
            
      if($syn_sql_sth->fetch){
        $syn_sql_sth->finish;
        if ((not defined($cs)) || ($cs->name eq $new_coord_system && $cs->version eq $new_version)) {
            return $self->fetch_by_region($new_coord_system, $new_name, $start, $end, $strand, $new_version, $no_fuzz);
        } elsif ($cs->name ne $new_coord_system) {
            warning("Searched for a known feature on coordinate system: ".$cs->dbID." but found it on: ".$new_coord_system.
            "\n No result returned, consider searching without coordinate system or use toplevel.");
            return;
        }
        
      }
      $syn_sql_sth->finish;


      if ($no_fuzz) { return; }

      # Do fuzzy matching, assuming that we are just missing a version
      # on the end of the seq_region name.

      $sth =
        $self->prepare( $sql . " WHERE sr.name LIKE ? " . $constraint );

      $bind_params[0] =
        [ sprintf( '%s.%%', $seq_region_name ), SQL_VARCHAR ];

      $pos = 0;
      foreach my $param (@bind_params) {
        $sth->bind_param( ++$pos, $param->[0], $param->[1] );
      }

      $sth->execute();

      my $prefix_len = length($seq_region_name) + 1;
      my $high_ver   = undef;
      my $high_cs    = $cs;

      # Find the fuzzy-matched seq_region with the highest postfix
      # (which ought to be a version).

      my ( $tmp_name, $id, $tmp_length, $cs_id );
      $sth->bind_columns( \( $tmp_name, $id, $tmp_length, $cs_id ) );

      my $i = 0;

      while ( $sth->fetch ) {
        my $tmp_cs =
          ( defined($cs) ? $cs : $csa->fetch_by_dbID($cs_id) );

        # cache values for future reference
        my $arr = [ $id, $tmp_name, $cs_id, $tmp_length ];
        $self->{'sr_name_cache'}->{"$tmp_name:$cs_id"} = $arr;
        $self->{'sr_id_cache'}->{"$id"}                = $arr;

        my $tmp_ver = substr( $tmp_name, $prefix_len );

        # skip versions which are non-numeric and apparently not
        # versions
        if ( $tmp_ver !~ /^\d+$/ ) { next }

        # take version with highest num, if two versions match take one
        # with highest ranked coord system (lowest num)
        if ( !defined($high_ver)
          || $tmp_ver > $high_ver
          || ( $tmp_ver == $high_ver && $tmp_cs->rank < $high_cs->rank )
          )
        {
          $seq_region_name = $tmp_name;
          $length          = $tmp_length;
          $high_ver        = $tmp_ver;
          $high_cs         = $tmp_cs;
        }

        $i++;
      } ## end while ( $sth->fetch )
      $sth->finish();

      # warn if fuzzy matching found more than one result
      if ( $i > 1 ) {
        warning(
          sprintf(
            "Fuzzy matching of seq_region_name "
              . "returned more than one result.\n"
              . "You might want to check whether the returned seq_region\n"
              . "(%s:%s) is the one you intended to fetch.\n",
            $high_cs->name(), $seq_region_name ) );
      }

      $cs = $high_cs;

      # return if we did not find any appropriate match:
      if ( !defined($high_ver) ) { return; }

    } else {

      my ( $id, $cs_id );
      ( $seq_region_name, $id, $length, $cs_id ) = @row;

      # cache to speed up for future queries
      my $arr = [ $id, $seq_region_name, $cs_id, $length ];
      $self->{'sr_name_cache'}->{"$seq_region_name:$cs_id"} = $arr;
      $self->{'sr_id_cache'}->{"$id"}                       = $arr;
      $cs = $csa->fetch_by_dbID($cs_id);
    }
  } ## end else [ if ( defined($arr) ) ]

  if ( !defined($end) ) { $end = $length }

  #If this was given then check if we've got a circular seq region otherwise
  #let it fall through to the normal Slice method
  if ( $end + 1 < $start ) {
    my $cs_id = $cs->dbID();
    my $seq_region_id = $self->{'sr_name_cache'}->{"$seq_region_name:$cs_id"}->[0];
    if($self->is_circular($seq_region_id)) {
      my $new_sl =
        Bio::EnsEMBL::CircularSlice->new(
                                     -COORD_SYSTEM    => $cs,
                                     -SEQ_REGION_NAME => $seq_region_name,
                                     -SEQ_REGION_LENGTH => $length,
                                     -START             => $start,
                                     -END               => $end,
                                     -STRAND            => 1,
                                     -ADAPTOR           => $self );
  
      return $new_sl;
    }
  }

  if ( defined( $self->{'lrg_region_test'} )
       and substr( $cs->name, 0, 3 ) eq $self->{'lrg_region_test'} )
  {
    return
      Bio::EnsEMBL::LRGSlice->new( -COORD_SYSTEM    => $cs,
                                   -SEQ_REGION_NAME => $seq_region_name,
                                   -SEQ_REGION_LENGTH => $length,
                                   -START             => $start,
                                   -END               => $end,
                                   -STRAND            => $strand,
                                   -ADAPTOR           => $self );
  } else {
    return
      Bio::EnsEMBL::Slice->new_fast( {
                                  'coord_system'    => $cs,
                                  'seq_region_name' => $seq_region_name,
                                  'seq_region_length' => $length,
                                  'start'             => $start,
                                  'end'               => $end,
                                  'strand'            => $strand,
                                  'adaptor'           => $self } );
  }
} ## end sub fetch_by_region

=head2 fetch_by_toplevel_location

  Arg [1]     : string $location
                Ensembl formatted location. Can be a format like 
                C<name:start-end>, C<name:start..end>, C<name:start:end>, 
                C<name:start>, C<name>. We can also support strand 
                specification as a +/- or 1/-1. 
                
                Location names must be separated by a C<:>. All others can be
                separated by C<..>, C<:> or C<->.
  Arg[2]      : boolean $no_warnings
                Suppress warnings from this method
  Arg[3]      : boolean $no_fuzz
                Stop fuzzy matching of sequence regions from occuring
  Arg[4]      : boolean $ucsc
                If we are unsuccessful at retriving a location retry taking any 
                possible chr prefix into account e.g. chrX and X are treated as
                equivalents
  Example     : my $slice = $sa->fetch_by_toplevel_location('X:1-10000')
                my $slice = $sa->fetch_by_toplevel_location('X:1-10000:-1')
  Description : Converts an Ensembl location/region into the sequence region
                name, start and end and passes them onto C<fetch_by_region()>. 
                The code assumes that the required slice is on the top level
                coordinate system. The code assumes that location formatting
                is not perfect and will perform basic cleanup before parsing.
  Returntype  : Bio::EnsEMBL::Slice
  Exceptions  : If $location is false otherwise see C<fetch_by_location()>
                or C<fetch_by_region()>
  Caller      : General
  Status      : Beta

=cut

sub fetch_by_toplevel_location {
  my ($self, $location, $no_warnings, $no_fuzz, $ucsc) = @_;
  return $self->fetch_by_location($location, 'toplevel', undef, $no_warnings, $no_fuzz, $ucsc);
}

=head2 fetch_by_location

  Arg [1]     : string $location
                Ensembl formatted location. Can be a format like 
                C<name:start-end>, C<name:start..end>, C<name:start:end>, 
                C<name:start>, C<name>. We can also support strand 
                specification as a +/- or 1/-1. 
                
                Location names must be separated by a C<:>. All others can be
                separated by C<..>, C<:>, C<_> or C<->.
  Arg[2]      : String $coord_system_name
                The coordinate system to retrieve
  Arg[3]      : String $coord_system_version
                Optional parameter. Version of the coordinate system to fetch
  Arg[4]      : boolean $no_warnings
                Suppress warnings from this method
  Arg[5]      : boolean $no_fuzz
                Stop fuzzy matching of sequence regions from occuring
  Arg[6]      : boolean $ucsc
                If we are unsuccessful at retriving a location retry taking any 
                possible chr prefix into account e.g. chrX and X are treated as
                equivalents
  Example     : my $slice = $sa->fetch_by_location('X:1-10000','chromosome')
                my $slice = $sa->fetch_by_location('X:1-10000:-1','toplevel')
  Description : Converts an Ensembl location/region into the sequence region
                name, start and end and passes them onto C<fetch_by_region()>. 
                The code assumes that location formatting is not perfect and 
                will perform basic cleanup before parsing.
  Returntype  : Bio::EnsEMBL::Slice
  Exceptions  : If $location or coordinate system is false otherwise 
                see C<fetch_by_region()>
  Caller      : General
  Status      : Beta

=cut

sub fetch_by_location {
  my ($self, $location, $coord_system_name, $coord_system_version, $no_warnings, $no_fuzz, $ucsc) = @_;
  
  throw "No coordinate system name specified" unless $coord_system_name;
  
  my ($seq_region_name, $start, $end, $strand) = $self->parse_location_to_values($location, $no_warnings);

  if(! $seq_region_name) {
    return;
  }
    
  if(defined $start && defined $end && $start > $end) {
    throw "Cannot request a slice whose start is greater than its end. Start: $start. End: $end";
  }
  
  my $slice = $self->fetch_by_region($coord_system_name, $seq_region_name, $start, $end, $strand, $coord_system_version, $no_fuzz);
  if(! defined $slice) {
    if($ucsc) {
      my $ucsc_seq_region_name = $seq_region_name;
      $ucsc_seq_region_name =~ s/^chr//;
      if($ucsc_seq_region_name ne $seq_region_name) {
        $slice = $self->fetch_by_region($coord_system_name, $ucsc_seq_region_name,  $start, $end, $strand, $coord_system_version, $no_fuzz);
        return if ! defined $slice; #if we had no slice still then bail
      }
      else {
        return; #If it was not different then we didn't have the prefix so just return (same bail as before)
      }
    }
    else {
      return; #We didn't have a slice and no UCSC specifics are being triggered
    }
  }
  
  my $srl = $slice->seq_region_length();
  my $name = $slice->seq_region_name();
  if(defined $start && $start > $srl) {
    throw "Cannot request a slice whose start ($start) is greater than $srl for $name.";
  }
  if(defined $end && $end > $srl) {
    warning "Requested end ($end) is greater than $srl for $name. Resetting to $srl" if ! $no_warnings;
    $slice->{end} = $srl;
  }
  
  return $slice;
}

=head2 parse_location_to_values

  Arg [1]     : string $location
                Ensembl formatted location. Can be a format like 
                C<name:start-end>, C<name:start..end>, C<name:start:end>, 
                C<name:start>, C<name>. We can also support strand 
                specification as a +/- or 1/-1. 
                
                Location names must be separated by a C<:>. All others can be
                separated by C<..>, C<:> C<_>, or C<->.

                If the start is negative, start will be reset to 1 (e.g.: 1: -10-1,000')
                If both start and end are negative, returns undef (e.g.: 1: -10--1,000')
  Arg[2]      : boolean $no_warnings
                Suppress warnings from this method
  Arg[3]      : boolean $no_errors
                Supress errors being thrown from this method
  Example			: my ($name, $start, $end, $strand) = $sa->parse_location_to_values('X:1..100:1);
  Description	: Takes in an Ensembl location String and returns the parsed
                values
  Returntype 	: List. Contains name, start, end and strand 

=cut


sub parse_location_to_values {
  my ($self, $location, $no_warnings, $no_errors) = @_;
  
  throw 'You must specify a location' if ! $location;
  
  #cleanup any nomenclature like 1 000 or 1,000
  my $number_seps_regex = qr/\s+|,/;
  my $separator_regex = qr/(?:-|[.]{2}|\:|_)?/; # support -, .., : and _ as separators
  my $hgvs_nomenclature_regex = qr/(?:g\.)?/; # check for HGVS looking locations e.g. X:g.1-100
  my $number_regex = qr/[0-9, EMKG]+/xmsi;
  my $number_regex_signed = qr/-?[0-9, EMKG]+/xmsi; # to capture negative locations as sometimes we end up in negative location if the location is padded
  my $strand_regex = qr/[+-1]|-1/xms;
  
  my $regex = qr/^((?:\w|\.|_|-)+) \s* :? \s* $hgvs_nomenclature_regex ($number_regex_signed)? $separator_regex ($number_regex)? $separator_regex ($strand_regex)? $/xms;
  my ($seq_region_name, $start, $end, $strand);
  if(($seq_region_name, $start, $end, $strand) = $location =~ $regex) {
    
    if(defined $strand) {
      if(!looks_like_number($strand)) {
        $strand = ($strand eq '+') ? 1 : -1;
      }
    }
    
    if(defined $start) {
      $start =~ s/$number_seps_regex//g; 
      if($start < 1) {
        warning "Start was less than 1 (${start}) which is not allowed. Resetting to 1"  if ! $no_warnings;

        unless(defined $end) {
          # We will reach here only when the location is given without start and '-' is used as seperator eg: 1:-10 (expected to return 1:1-10)
          $end = abs($start);   	
        }
          $start = 1;
      }
    }
    if(defined $end) {
      $end =~ s/$number_seps_regex//g;
      if($end < 1) {
        throw "Cannot request negative or 0 end indexes through this interface. Given $end but expected something greater than 0" unless $no_errors;
      }
    }
    
  }
  
  return ($seq_region_name, $start, $end, $strand);
}

=head2 fetch_by_region_unique

  Arg [1]    : string $coord_system_name (optional)
               The name of the coordinate system of the slice to be created
               This may be a name of an actual coordinate system or an alias
               to a coordinate system.  Valid aliases are 'seqlevel' or
               'toplevel'.
  Arg [2]    : string $seq_region_name
               The name of the sequence region that the slice will be
               created on.
  Arg [3]    : int $start (optional, default = 1)
               The start of the slice on the sequence region
  Arg [4]    : int $end (optional, default = seq_region length)
               The end of the slice on the sequence region
  Arg [5]    : int $strand (optional, default = 1)
               The orientation of the slice on the sequence region
  Arg [6]    : string $version (optional, default = default version)
               The version of the coordinate system to use (e.g. NCBI33)
  Arg [7]    : boolean $no_fuzz (optional, default = undef (false))
               If true (non-zero), do not use "fuzzy matching" (see below).
  Example    : $slice = $slice_adaptor->fetch_by_region_unique('chromosome', 'HSCHR6_MHC_COX');
  Description: Retrieves a slice on the requested region but returns only the unique
               parts of the slice.  At a minimum the
               name the name of the seq_region to fetch must be provided.

               If no coordinate system name is provided than a slice on the
               highest ranked coordinate system with a matching
               seq_region_name will be returned.  If a version but no
               coordinate system name is provided, the same behaviour will
               apply, but only coordinate systems of the appropriate version
               are considered.  The same applies if the 'toplevel' coordinate
               system is specified, however in this case the version is
               ignored.  The coordinate system should always be specified if
               it is known, since this is unambiguous and faster.

               Some fuzzy matching is performed if no exact match for
               the provided name is found.  This allows clones to be
               fetched even when their version is not known.  For
               example fetch_by_region('clone', 'AC008066') will
               retrieve the sequence_region with name 'AC008066.4'.

               The fuzzy matching can be turned off by setting the
               $no_fuzz argument to a true value.

               If the requested seq_region is not found in the database undef
               is returned.

  Returntype : listref Bio::EnsEMBL::Slice
  Exceptions : throw if no seq_region_name is provided
               throw if invalid coord_system_name is provided
               throw if start > end is provided
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_region_unique {
  my $self = shift;

  my @out   = ();
  my $slice = $self->fetch_by_region(@_);


  if ( !exists( $self->{'asm_exc_cache'} ) ) {
    $self->_build_exception_cache();
  }

  if ( exists(
          $self->{'asm_exc_cache'}->{ $self->get_seq_region_id($slice) }
       ) )
  {
    # Dereference symlinked assembly regions.  Take out any regions
    # which are symlinked because these are duplicates.
    my @projection =
      @{ $self->fetch_normalized_slice_projection($slice) };

    foreach my $segment (@projection) {
      if ( $segment->[2]->seq_region_name() eq $slice->seq_region_name()
        && $segment->[2]->coord_system->equals( $slice->coord_system ) )
      {
        push( @out, $segment->[2] );
      }
    }
  } else {
    @out = ($slice);
  }

  return \@out;
} ## end sub fetch_by_region_unique

=head2 fetch_by_name

  Arg [1]    : string $name
  Example    : $name  = 'chromosome:NCBI34:X:1000000:2000000:1';
               $slice = $slice_adaptor->fetch_by_name($name);
               $slice2 = $slice_adaptor->fetch_by_name($slice3->name());
  Description: Fetches a slice using a slice name (i.e. the value returned by
               the Slice::name method).  This is useful if you wish to 
               store a unique identifier for a slice in a file or database or
               pass a slice over a network.
               Slice::name allows you to serialise/marshall a slice and this
               method allows you to deserialise/unmarshal it.

               Returns undef if no seq_region with the provided name exists in
               the database.

  Returntype : Bio::EnsEMBL::Slice or undef
  Exceptions : throw if incorrent arg provided
  Caller     : Pipeline
  Status     : Stable

=cut

sub fetch_by_name {
  my $self = shift;
  my $name = shift;

  if(!$name) {
    throw("name argument is required");
  }

  my @array = split(/:/,$name);

  if(scalar(@array) < 3 || scalar(@array) > 6) {
    throw("Malformed slice name [$name].  Format is " .
        "coord_system:version:name:start:end:strand");
  }

  # Rearrange arguments to suit fetch_by_region

  my @targetarray;

  $targetarray[0]=$array[0];
  $targetarray[5]=(($array[1]&&$array[1] ne "")?$array[1]:undef);
  $targetarray[1]=(($array[2]&&$array[2] ne "")?$array[2]:undef);
  $targetarray[2]=(($array[3]&&$array[3] ne "")?$array[3]:undef);
  $targetarray[3]=(($array[4]&&$array[4] ne "")?$array[4]:undef);
  $targetarray[4]=(($array[5]&&$array[5] ne "")?$array[5]:undef);
  return $self->fetch_by_region(@targetarray);
}



=head2 fetch_by_seq_region_id

  Arg [1]    : string $seq_region_id
               The internal identifier of the seq_region to create this slice
               on
  Arg [2]    : optional start
  Arg [3]    : optional end
  Arg [4]    : optional strand
  Example    : $slice = $slice_adaptor->fetch_by_seq_region_id(34413);
  Description: Creates a slice object of an entire seq_region using the
               seq_region internal identifier to resolve the seq_region.
               Returns undef if no such slice exists.
  Returntype : Bio::EnsEMBL::Slice or undef
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_seq_region_id {
  my ( $self, $seq_region_id, $start, $end, $strand, $check_prior_ids ) = @_;

  my $csa = $self->db->get_CoordSystemAdaptor();
  my $arr = $self->{'sr_id_cache'}->{$seq_region_id};
  my ( $name, $length, $cs, $cs_id );


  if ( $arr && defined( $arr->[2] ) ) {
    ( $name, $cs_id, $length ) = ( $arr->[1], $arr->[2], $arr->[3] );
    $cs = $self->db->get_CoordSystemAdaptor->fetch_by_dbID($cs_id);
  } else {
    my $sth =
      $self->prepare(   "SELECT sr.name, sr.coord_system_id, sr.length "
                      . "FROM seq_region sr "
                      . "WHERE sr.seq_region_id = ? " );

    $sth->bind_param( 1, $seq_region_id, SQL_INTEGER );
    $sth->execute();

    my @row = $sth->fetchrow_array();
    unless ( @row ) {
      # This could have been an old seq region id so see if we can
      # translate it into a more recent version.
      if($check_prior_ids) {
        if(exists $csa->{_external_seq_region_mapping}->{$seq_region_id}) {
          my $new_seq_region_id = $csa->{_external_seq_region_mapping}->{$seq_region_id};
          # No need to pass check prior ids flag because it's a 1 step relationship
          return $self->fetch_by_seq_region_id($new_seq_region_id, $start, $end, $strand);
        }
      }
      return undef;
    }

    ( $name, $cs_id, $length ) = @row;
    $sth->finish();

    $cs = $self->db->get_CoordSystemAdaptor->fetch_by_dbID($cs_id);

    #cache results to speed up repeated queries
    my $arr = [ $seq_region_id, $name, $cs_id, $length ];

    $self->{'sr_name_cache'}->{"$name:$cs_id"} = $arr;
    $self->{'sr_id_cache'}->{"$seq_region_id"} = $arr;
  }

  return
    Bio::EnsEMBL::Slice->new_fast({ 
	                      'coord_system'     => $cs,
                              'seq_region_name'  => $name,
                              'seq_region_length'=> $length,
                              'start'            => $start || 1,
                              'end'              => $end || $length,
                              'strand'           => $strand || 1,
                              'adaptor'           => $self} );
} ## end sub fetch_by_seq_region_id



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
  Status     : Stable

=cut

sub get_seq_region_id {
  my $self = shift;
  my $slice = shift;

  if(!$slice || !ref($slice) || !($slice->isa('Bio::EnsEMBL::Slice') or $slice->isa('Bio::EnsEMBL::LRGSlice'))) {
    throw('Slice argument is required');
  }
  
  my $seq_region_name = $slice->seq_region_name();
  my $key = $seq_region_name.":".$slice->coord_system->dbID();
  my $arr = $self->{'sr_name_cache'}->{"$key"};

  if( $arr ) {
    return $arr->[0];
  }

  my $cs_id = $slice->coord_system->dbID();

  my $sth = $self->prepare("SELECT seq_region_id, length " .
                           "FROM seq_region " .
                           "WHERE name = ? AND coord_system_id = ?");

  #force seq_region_name cast to string so mysql cannot treat as int
  $sth->bind_param(1,"$seq_region_name",SQL_VARCHAR);
  $sth->bind_param(2,$cs_id,SQL_INTEGER);
  $sth->execute();

  my @row = $sth->fetchrow_array();
  unless ( @row ) {
    throw("No-existent seq_region [$seq_region_name] in coord system [$cs_id]");
  }
  my @more = $sth->fetchrow_array();
  if ( @more ) {
    throw("Ambiguous seq_region [$seq_region_name] in coord system [$cs_id]");
  }

  my($seq_region_id, $length) = @row;
  $sth->finish();

  #cache information for future requests
  $arr = [ $seq_region_id, $seq_region_name, $cs_id, $length ];

  $self->{'sr_name_cache'}->{"$seq_region_name:$cs_id"} = $arr;
  $self->{'sr_id_cache'}->{"$seq_region_id"} = $arr;

  return $seq_region_id;
}



=head2 fetch_all

  Arg [1]    : string $coord_system_name
               The name of the coordinate system to retrieve slices of.
               This may be a name of an acutal coordinate system or an alias
               to a coordinate system.  Valid aliases are 'seqlevel' or
               'toplevel'.
  Arg [2]    : string $coord_system_version (optional)
               The version of the coordinate system to retrieve slices of
  Arg [3]    : bool $include_non_reference (optional)
               If this argument is not provided then only reference slices
               will be returned. If set, both reference and non reference
               slices will be returned.
  Arg [4]    : int $include_duplicates (optional)
               If set duplicate regions will be returned.
               
               NOTE: if you do not use this option and you have a PAR
               (pseudo-autosomal region) at the beginning of your seq_region
               then your slice will not start at position 1, so coordinates
               retrieved from this slice might not be what you expected.

  Arg[5]     : bool $include_lrg (optional)  (default 0)
               If set lrg regions will be returned aswell.


  Example    : @chromos = @{$slice_adaptor->fetch_all('chromosome','NCBI33')};
               @contigs = @{$slice_adaptor->fetch_all('contig')};

               # get even non-reference regions
               @slices = @{$slice_adaptor->fetch_all('toplevel',undef,1)};

               # include duplicate regions (such as pseudo autosomal regions)
               @slices = @{$slice_adaptor->fetch_all('toplevel', undef,0,1)};

  Description: Retrieves slices of all seq_regions for a given coordinate
               system.  This is analagous to the methods fetch_all which were
               formerly on the ChromosomeAdaptor, RawContigAdaptor and
               CloneAdaptor classes.  Slices fetched span the entire
               seq_regions and are on the forward strand.
               If the coordinate system with the provided name and version
               does not exist an empty list is returned.
               If the coordinate system name provided is 'toplevel', all
               non-redundant toplevel slices are returned (note that any
               coord_system_version argument is ignored in that case).

               Retrieved slices can be broken into smaller slices using the
               Bio::EnsEMBL::Utils::Slice module.

  Returntype : listref of Bio::EnsEMBL::Slices
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all {
  my $self = shift;
  my $cs_name = shift;
  my $cs_version = shift || '';

  my ($include_non_reference, $include_duplicates, $include_lrg) = @_;

  #
  # verify existance of requested coord system and get its id
  #
  my $csa       = $self->db->get_CoordSystemAdaptor();
  my $orig_cs   = $csa->fetch_by_name($cs_name, $cs_version);

  return [] if ( !$orig_cs );

  my %bad_vals=();


  #
  # Get a hash of non reference seq regions
  #
  if ( !$include_non_reference ) {
    my $sth =
      $self->prepare(   'SELECT sr.seq_region_id '
                      . 'FROM seq_region sr, seq_region_attrib sra, '
                      . 'attrib_type at, coord_system cs '
                      . 'WHERE at.code = "non_ref" '
                      . 'AND sra.seq_region_id = sr.seq_region_id '
                      . 'AND at.attrib_type_id = sra.attrib_type_id '
                      . 'AND sr.coord_system_id = cs.coord_system_id '
                      . 'AND cs.species_id = ?' );

    $sth->bind_param( 1, $self->species_id(), SQL_INTEGER );
    $sth->execute();

    my ($seq_region_id);
    $sth->bind_columns( \$seq_region_id );

    while ( $sth->fetch() ) {
      $bad_vals{$seq_region_id} = 1;
    }
  }

  #
  # if we do not want lrg's then add them to the bad list;
  #
  if ( !$include_lrg ) {
    my $sth =
      $self->prepare(   'SELECT sr.seq_region_id '
                      . 'FROM seq_region sr, seq_region_attrib sra, '
                      . 'attrib_type at, coord_system cs '
                      . 'WHERE at.code = "LRG" '
                      . 'AND sra.seq_region_id = sr.seq_region_id '
                      . 'AND at.attrib_type_id = sra.attrib_type_id '
                      . 'AND sr.coord_system_id = cs.coord_system_id '
                      . 'AND cs.species_id = ?' );

    $sth->bind_param( 1, $self->species_id(), SQL_INTEGER );
    $sth->execute();

    my ($seq_region_id);
    $sth->bind_columns( \$seq_region_id );

    while ( $sth->fetch() ) {
      $bad_vals{$seq_region_id} = 1;
    }
  }

  #
  # Retrieve the seq_regions from the database
  #

  my $sth;
  if ( $orig_cs->is_top_level() ) {
    $sth =
      $self->prepare(   'SELECT sr.seq_region_id, sr.name, '
                      . 'sr.length, sr.coord_system_id '
                      . 'FROM seq_region sr, seq_region_attrib sra, '
                      . 'attrib_type at, coord_system cs '
                      . 'WHERE at.code = "toplevel" '
                      . 'AND at.attrib_type_id = sra.attrib_type_id '
                      . 'AND sra.seq_region_id = sr.seq_region_id '
                      . 'AND sr.coord_system_id = cs.coord_system_id '
                      . 'AND cs.species_id = ?' );

    $sth->bind_param( 1, $self->species_id(), SQL_INTEGER );
    $sth->execute();
  } else {
    $sth =
      $self->prepare(   'SELECT sr.seq_region_id, sr.name, '
                      . 'sr.length, sr.coord_system_id '
                      . 'FROM seq_region sr '
                      . 'WHERE sr.coord_system_id = ?' );

    $sth->bind_param( 1, $orig_cs->dbID, SQL_INTEGER );
    $sth->execute();
  }

  my ( $seq_region_id, $name, $length, $cs_id );
  $sth->bind_columns( \( $seq_region_id, $name, $length, $cs_id ) );

  my $cache_count = 0;

  my @out;
  while($sth->fetch()) {
    if(!defined($bad_vals{$seq_region_id})){
      my $cs = $csa->fetch_by_dbID($cs_id);

      if(!$cs) {
        throw("seq_region $name references non-existent coord_system $cs_id.");
      }

      #cache values for future reference, but stop adding to the cache once we
      #we know we have filled it up
      if($cache_count < $Bio::EnsEMBL::Utils::SeqRegionCache::SEQ_REGION_CACHE_SIZE) {
        my $arr = [ $seq_region_id, $name, $cs_id, $length ];

        $self->{'sr_name_cache'}->{"$name:$cs_id"} = $arr;
        $self->{'sr_id_cache'}->{"$seq_region_id"} = $arr;

        $cache_count++;
      }

      my $slice = Bio::EnsEMBL::Slice->new_fast({
	  'start'           => 1,
          'end'             => $length,
          'strand'          => 1,
         'seq_region_name'  => $name,
         'seq_region_length'=> $length,
         'coord_system'     => $cs,
         'adaptor'          => $self});

      if(!defined($include_duplicates) or !$include_duplicates){
        # test if this slice *could* have a duplicate (exception) region
        $self->_build_exception_cache() if(!exists $self->{'asm_exc_cache'});
        if(exists $self->{asm_exc_cache}->{$seq_region_id}) {

          # Dereference symlinked assembly regions.  Take out
          # any regions which are symlinked because these are duplicates
          my @projection = @{$self->fetch_normalized_slice_projection($slice)};
          foreach my $segment ( @projection) {
            if($segment->[2]->seq_region_name() eq $slice->seq_region_name() &&
               $segment->[2]->coord_system->equals($slice->coord_system)) {
              push @out, $segment->[2];
            }
          }
        } else {
          # no duplicate regions
          push @out, $slice;
        }
      } else {
        # we want duplicates anyway so do not do any checks
        push @out, $slice;
      }
    }
  }

  return \@out;
}

=head2 fetch_all_by_genome_component

  Arg [1]    : string $genome_component_name
               The name of the genome component to retrieve slices of.
  Example    : @slices = @{$slice_adaptor->fetch_all_by_genome_component('A')};
  Description: Returns the list of all top level slices for a a given 
               genome component
  Returntype : listref of Bio::EnsEMBL::Slices
  Exceptions : If argument is not provided or is not a valid genome
               component
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_genome_component {
  my $self = shift;
  my $genome_component = shift;
  defined $genome_component or
    throw "Undefined genome component";

  # check the provided genome component is valid
  my $gc = $self->db->get_adaptor('GenomeContainer');
  my $is_valid_component = grep { $_ eq $genome_component } 
    @{$gc->get_genome_components};
  throw "Invalid genome component"
    unless $is_valid_component;
  
  #
  # Retrieve the toplevel seq_regions from the database
  #
  my $sth =
    $self->prepare(   "SELECT sr.seq_region_id, sr.name, sr.length, sr.coord_system_id "
		      . "FROM seq_region sr "
		      . "JOIN seq_region_attrib sa1 USING (seq_region_id) "
		      . "JOIN attrib_type a1 ON sa1.attrib_type_id = a1.attrib_type_id "
		      . "JOIN seq_region_attrib sa2 USING (seq_region_id) "
		      . "JOIN attrib_type a2 ON sa2.attrib_type_id = a2.attrib_type_id "
		      . "WHERE sa2.value=? and a1.code='toplevel' and a2.code='genome_component'"
		  );
  
  if (looks_like_number($genome_component)) {
    $sth->bind_param( 1, $genome_component, SQL_INTEGER );
  } else {
    $sth->bind_param( 1, $genome_component, SQL_VARCHAR );
  }
  
  $sth->execute();

  my ( $seq_region_id, $name, $length, $cs_id );
  $sth->bind_columns( \( $seq_region_id, $name, $length, $cs_id ) );

  my $cache_count = 0;

  my @out;
  my $csa = $self->db->get_CoordSystemAdaptor();

  while($sth->fetch()) {
    my $cs = $csa->fetch_by_dbID($cs_id);
    throw("seq_region $name references non-existent coord_system $cs_id.")
	  unless $cs;
    
    # cache values for future reference, but stop adding to the cache once we
    # we know we have filled it up
    if($cache_count < $Bio::EnsEMBL::Utils::SeqRegionCache::SEQ_REGION_CACHE_SIZE) {
      my $arr = [ $seq_region_id, $name, $cs_id, $length ];

      $self->{'sr_name_cache'}->{"$name:$cs_id"} = $arr;
      $self->{'sr_id_cache'}->{"$seq_region_id"} = $arr;

      $cache_count++;
    }

    push @out, Bio::EnsEMBL::Slice->new_fast({
					      'start'           => 1,
					      'end'             => $length,
					      'strand'          => 1,
					      'seq_region_name'  => $name,
					      'seq_region_length'=> $length,
					      'coord_system'     => $cs,
					      'adaptor'          => $self});
  }

  return \@out;
}

=head2 get_genome_component_for_slice

  Arg [1]    : An object of type Bio::EnsEMBL::Slice
  Example    : my $component = $slice->get_genome_component();
  Description: Returns the genome component of a slice
  Returntype : Scalar; the identifier of the genome component of the slice
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_genome_component_for_slice {
  my ($self, $slice) = @_;

  throw "Undefined slice" unless defined $slice;
  throw "Argument is not a slice"
    unless $slice->isa("Bio::EnsEMBL::Slice");

  my $seq_region_id = $self->get_seq_region_id($slice);

  my $sth =
    $self->prepare(   "SELECT sa.value "
		      . "FROM seq_region_attrib sa "
		      . "JOIN seq_region sr USING (seq_region_id) "
		      . "JOIN attrib_type at ON sa.attrib_type_id = at.attrib_type_id "
		      . "WHERE sr.seq_region_id=? AND at.code='genome_component'"
		  );
  $sth->bind_param( 1, $seq_region_id, SQL_INTEGER );
  $sth->execute();

  my $genome_component;
  $sth->bind_columns( \( $genome_component ) );
  $sth->fetch();

  return $genome_component;
}

=head2 fetch_all_karyotype

  Example    : my $top = $slice_adptor->fetch_all_karyotype()
  Description: returns the list of all slices which are part of the karyotype
  Returntype : listref of Bio::EnsEMBL::Slices
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_karyotype {
  my $self = shift;

  my $csa       = $self->db->get_CoordSystemAdaptor();

  my $sth = 
    $self->prepare( 'SELECT sr.seq_region_id, sr.name, '
                      . 'sr.length, sr.coord_system_id, sra.value '
                      . 'FROM seq_region sr, seq_region_attrib sra, '
                      . 'attrib_type at, coord_system cs '
                      . 'WHERE at.code = "karyotype_rank" '
                      . 'AND at.attrib_type_id = sra.attrib_type_id '
                      . 'AND sra.seq_region_id = sr.seq_region_id '
                      . 'AND sr.coord_system_id = cs.coord_system_id '
                      . 'AND cs.species_id = ?');
  $sth->bind_param( 1, $self->species_id(), SQL_INTEGER );
  $sth->execute();
  my ( $seq_region_id, $name, $length, $cs_id, $rank );
  $sth->bind_columns( \( $seq_region_id, $name, $length, $cs_id, $rank ) );

  my @out;
  while($sth->fetch()) {
    my $cs = $csa->fetch_by_dbID($cs_id);

    my $slice = Bio::EnsEMBL::Slice->new_fast({
          'start'           => 1,
          'end'             => ($length+0),
          'strand'          => 1,
         'seq_region_name'  => $name,
         'seq_region_length'=> ($length+0),
         'coord_system'     => $cs,
         'adaptor'          => $self,
         'karyotype'        => 1,
         'karyotype_rank'   => ($rank+0),
         });

    push @out, $slice;
  }

  #Sort using Perl as value in MySQL is a text field and not portable
  @out = sort { $a->{karyotype_rank} <=> $b->{karyotype_rank} } @out;

  return \@out;
}

=head2 is_toplevel
  Arg        : int seq_region_id 
  Example    : my $top = $slice_adptor->is_toplevel($seq_region_id)
  Description: Returns 1 if slice is a toplevel slice else 0
  Returntype : int
  Caller     : Slice method is_toplevel
  Status     : At Risk

=cut

sub is_toplevel {
  my $self = shift;
  my $id   = shift;

  my $sth = $self->prepare(
            "SELECT at.code from seq_region_attrib sra, attrib_type at "
              . "WHERE sra.seq_region_id = ? "
              . "AND at.attrib_type_id = sra.attrib_type_id "
              . "AND at.code = 'toplevel'" );

  $sth->bind_param( 1, $id, SQL_INTEGER );
  $sth->execute();

  my $code;
  $sth->bind_columns( \$code );

  while ( $sth->fetch ) {
    $sth->finish;
    return 1;
  }

  $sth->finish;
  return 0;
}


=head2 has_karyotype
  Arg        : int seq_region_id 
  Example    : my $karyotype = $slice_adptor->has_karyotype($seq_region_id)
  Description: Returns 1 if slice is a part of a karyotype else 0
  Returntype : int
  Caller     : Slice method has_karyotype
  Status     : At Risk

=cut

sub has_karyotype {
  my $self = shift;
  my $id   = shift;

  my $sth = $self->prepare(
            "SELECT at.code from seq_region_attrib sra, attrib_type at "
              . "WHERE sra.seq_region_id = ? "
              . "AND at.attrib_type_id = sra.attrib_type_id "
              . "AND at.code = 'karyotype_rank'" );

  $sth->bind_param( 1, $id, SQL_INTEGER );
  $sth->execute();

  my $code;
  $sth->bind_columns( \$code );

  while ( $sth->fetch ) {
    $sth->finish;
    return 1;
  }

  $sth->finish;
  return 0;
}

=head2 get_karyotype_rank
  Arg        : int seq_region_id 
  Example    : my $rank = $slice_adptor->get_karyotype_rank($seq_region_id)
  Description: Returns the rank of a slice if it is part of the karyotype else 0
  Returntype : int
  Caller     : Slice method get_karyotype_rank
  Status     : At Risk

=cut

sub get_karyotype_rank {
  my $self = shift;
  my $id   = shift;

  my $sth = $self->prepare(
            "SELECT sra.value from seq_region_attrib sra, attrib_type at "
              . "WHERE sra.seq_region_id = ? "
              . "AND at.attrib_type_id = sra.attrib_type_id "
              . "AND at.code = 'karyotype_rank'" );

  $sth->bind_param( 1, $id, SQL_INTEGER );
  $sth->execute();

  my $code;
  $sth->bind_columns( \$code );

  my $rank = $sth->fetchrow_array();
  $sth->finish();

  return $rank;
}



=head2 is_reference
  Arg        : int seq_region_id 
  Example    : my $reference = $slice_adptor->is_reference($seq_region_id)
  Description: Returns 1 if slice is a reference slice else 0
  Returntype : int
  Caller     : Slice method is_reference
  Status     : At Risk

=cut

sub is_reference {
  my $self = shift;
  my $id   = shift;

  my $sth = $self->prepare(
            "SELECT at.code from seq_region_attrib sra, attrib_type at "
              . "WHERE sra.seq_region_id = ? "
              . "AND at.attrib_type_id = sra.attrib_type_id "
              . "AND at.code = 'non_ref'" );

  $sth->bind_param( 1, $id, SQL_INTEGER );
  $sth->execute();

  my $code;
  $sth->bind_columns( \$code );

  while ( $sth->fetch ) {
    $sth->finish;
    return 0;
  }

  $sth->finish;
  return 1;
}

=head2 is_circular

  Arg[1]      : int seq_region_id
  Example     : my $circular = $slice_adptor->is_circular($seq_region_id);
  Description : Indicates if the sequence region was circular or not
  Returntype  : Boolean
  
=cut

sub is_circular {
  my ($self, $id) = @_;
  
  if (! defined $self->{is_circular}) {
    $self->_build_circular_slice_cache();
  }
  
  return 0 if $self->{is_circular} == 0;
  return (exists $self->{circular_sr_id_cache}->{$id}) ? 1 : 0;
}


=head2 fetch_by_chr_band

 Title   : fetch_by_chr_band
 Usage   :
 Function: create a Slice representing a series of bands
 Example :
 Returns : Bio::EnsEMBL::Slice
 Args    : the band name
 Status     : Stable

=cut

sub fetch_by_chr_band {
  my ( $self, $chr, $band ) = @_;

  my $chr_slice = $self->fetch_by_region( 'toplevel', $chr );
  my $seq_region_id = $self->get_seq_region_id($chr_slice);

  my $sth =
    $self->prepare(   'SELECT MIN(k.seq_region_start), '
                    . 'MAX(k.seq_region_end) '
                    . 'FROM karyotype k '
                    . 'WHERE k.seq_region_id = ? '
                    . 'AND k.band LIKE ?' );

  $sth->bind_param( 1, $seq_region_id, SQL_INTEGER );
  $sth->bind_param( 2, "$band%",       SQL_VARCHAR );
  $sth->execute();

  my ( $slice_start, $slice_end ) = $sth->fetchrow_array;

  if ( defined $slice_start ) {
    return
      $self->fetch_by_region( 'toplevel',   $chr,
                              $slice_start, $slice_end );
  }

  throw("Band not recognised in database");
} ## end sub fetch_by_chr_band



=head2 fetch_by_exon_stable_id

  Arg [1]    : string $exonid
               The stable id of the exon around which the slice is 
               desired
  Arg [2]    : (optional) int $size
               The length of the flanking regions the slice should encompass 
               on either side of the exon (0 by default)
  Example    : $slc = $sa->fetch_by_exon_stable_id('ENSE00000302930',10);
  Description: Creates a slice around the region of the specified exon. 
               If a context size is given, the slice is extended by that 
               number of basepairs on either side of the exon.
               
               The slice will be created in the exon's native coordinate system
               and in the forward orientation.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : Thrown if the exon is not in the database.
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_exon_stable_id{
  my ($self,$exonid,$size) = @_;

  throw('Exon argument is required.') if(!$exonid);

  my $ea = $self->db->get_ExonAdaptor;
  my $exon = $ea->fetch_by_stable_id($exonid);

  throw("Exon [$exonid] does not exist in DB.") if(!$exon);

  return $self->fetch_by_Feature($exon, $size);
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
               
               The slice will be created in the transcript's native coordinate
               system and in the forward orientation.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : Thrown if the transcript is not in the database.
  Caller     : general
  Status     : Stable

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
               
               The slice will be created in the transcript's native coordinate
               system and in the forward orientation.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : throw on incorrect args
               throw if transcript is not in database
  Caller     : general
  Status     : Stable

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
  Example    : $slc = $sa->fetch_by_gene_stable_id('ENSG00000012123',10);
  Description: Creates a slice around the region of the specified gene.
               If a context size is given, the slice is extended by that
               number of basepairs on either side of the gene.
               
               The slice will be created in the gene's native coordinate system
               and in the forward orientation.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : throw on incorrect args
               throw if transcript does not exist
  Caller     : general
  Status     : Stable

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
               The size may also be provided as a percentage of the feature 
               size such as 200% or 80.5%.
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
  Status     : Stable

=cut

sub fetch_by_Feature{
  my ($self, $feature, $size) = @_;

  $size ||= 0;

  if(!ref($feature) || !$feature->isa('Bio::EnsEMBL::Feature')) {
    throw('Feature argument expected.');
  }

  my $slice = $feature->slice();
  if(!$slice || !($slice->isa('Bio::EnsEMBL::Slice') or $slice->isa('Bio::EnsEMBL::LRGSlice') )) {
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

  ## Size may be stored as a %age of the length of the feature
  ## Size = 100% gives no context
  ## Size = 200% gives context - 50% the size of the feature either side of 
  ## feature

  $size = int( ($1-100)/200 * ($fend-$fstart+1) ) if( $size =~/([\d+\.]+)%/ );

  #return a new slice covering the region of the feature
  my $S = Bio::EnsEMBL::Slice->new_fast({
    'seq_region_name'   => $slice->seq_region_name,
    'seq_region_length' => $slice->seq_region_length,
    'coord_system'      => $slice->coord_system,
    'start'             => $fstart - $size,
    'end'               => $fend + $size,
    'strand'            => 1,
     'adaptor'          => $self});
  $S->{'_raw_feature_strand'}  = $feature->strand * $slice_strand if $feature->can('strand');
  return $S;
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
  Status     : Stable

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

=head2 fetch_by_misc_feature_set

  Arg [1]    : string $attribute_type
               The code of the attribute type
  Arg [2]    : (optional) string $attribute_value
               The value of the attribute to fetch by
  Arg [3]    : (optional) the name of the set
  Arg [4]    : (optional) int $size
               The amount of flanking region around the misc feature desired.
  Example    : $slice = $sa->fetch_by_misc_feature_set('clone',
                                                        'RP11-411G9'
                                                        'tilepath');
  Description: Fetches a slice around a MiscFeature with a particular
               attribute type, value and set. If no value is specified then
               the feature with the particular attribute is used.
               A size can be specified to include flanking region
               If no size is specified then 0 is used.
  Returntype : Bio::EnsEMBL::Slice
  Exceptions : Throw if no feature with the specified attribute type, value and set
               exists in the database
               Warning if multiple features with the specified attribute type, set
               and value exist in the database.
  Caller     : webcode
  Status     : Stable

=cut

sub fetch_by_misc_feature_set {
  my ($self, $attrib_type_code, $attrib_value, $misc_set, $size) = @_;

  my $mfa = $self->db()->get_MiscFeatureAdaptor();

  my $feat = $mfa->fetch_by_attribute_set_value($attrib_type_code,
                                                $attrib_value,
                                                $misc_set);

  return $self->fetch_by_Feature($feat, $size);
}

=head2 fetch_normalized_slice_projection

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Arg [2]    : boolean $filter_projections 
               Optionally filter the projections to remove anything 
               which is the same sequence region as the given slice
  Example    :  ( optional )
  Description: gives back a project style result. The returned slices 
               represent the areas to which there are symlinks for the 
               given slice. start, end show which area on given slice is 
               symlinked
  Returntype : [[start,end,$slice][]]
  Exceptions : none
  Caller     : BaseFeatureAdaptor
  Status     : Stable

=cut


sub fetch_normalized_slice_projection {
  my $self = shift;
  my $slice = shift;
  my $filter_projections = shift;

  my $slice_seq_region_id = $self->get_seq_region_id( $slice );

  $self->_build_exception_cache() if(!exists($self->{'asm_exc_cache'}));

  my $result = $self->{'asm_exc_cache'}->{$slice_seq_region_id};

  $result ||= [];

  my (@haps, @pars);

  foreach my $row (@$result) {
    my ( $seq_region_id, $seq_region_start, $seq_region_end,
         $exc_type, $exc_seq_region_id, $exc_seq_region_start,
         $exc_seq_region_end ) = @$row;

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

  if(!@pars && !@haps) {
    #just return this slice, there were no haps or pars
    return  [bless ( [1,$slice->length, $slice], "Bio::EnsEMBL::ProjectionSegment")];
  }

  my @syms;
  if( @haps >= 1 ) {
    my @sort_haps = sort { $a->[1] <=> $b->[1] } @haps;

    my $count =0;
    my $chr_start = 1;
    my $hap_start = 1;
    my $last = 0;

    my $seq_region_slice = $self->fetch_by_seq_region_id($slice_seq_region_id);
    my $exc_slice = $self->fetch_by_seq_region_id( $sort_haps[0][2] );
    my $len1 = $seq_region_slice->length();
    my $len2 = $exc_slice->length();
    my $max_len = ($len1 > $len2) ? $len1 : $len2;

    while($count <= scalar(@sort_haps)  and !$last){
      my $chr_end;
      my $hap_end;
      if(defined($sort_haps[$count]) and defined($sort_haps[$count][0]) ){
        $hap_end = $sort_haps[$count][0]-1;
        $chr_end = $sort_haps[$count][3]-1
      }
      else{
  $last = 1;
  $hap_end = $len1;
  $chr_end = $len2;
  my $diff = ($hap_end-$hap_start)-($chr_end-$chr_start);
  if($diff > 0){
    push( @syms, [ $hap_start, $hap_end, $sort_haps[0][2], $chr_start, $chr_end+$diff] );  
  }
  elsif($diff < 0){
    push( @syms, [ $hap_start, $hap_end - $diff, $sort_haps[0][2], $chr_start, $chr_end] );  
  }
  else{
    push( @syms, [ $hap_start, $hap_end, $sort_haps[0][2], $chr_start, $chr_end] );  
  }	
  next;
      }
      if($hap_end and $hap_start < $len1){ # if hap at start or end of chromosome
        push( @syms, [ $hap_start, $hap_end, $sort_haps[0][2], $chr_start, $chr_end] );  
      }
      $chr_start = $chr_end + ($sort_haps[$count][4]-$sort_haps[$count][3]) + 2;
      $hap_start = $hap_end + ($sort_haps[$count][1]-$sort_haps[$count][0]) + 2;
      $count++;
    }


  }


  # for now haps and pars should not be both there, but in theory we 
  # could handle it here by cleverly merging the pars into the existing syms,
  # for now just:
  push( @syms, @pars );

  my $mapper = Bio::EnsEMBL::Mapper->new( "sym", "org" );
  my $count = 0;
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
    return [bless( [1,$slice->length, $slice], "Bio::EnsEMBL::ProjectionSegment" )];
  }

  my @out;
  for my $coord ( @linked ) {
    if( $coord->isa( "Bio::EnsEMBL::Mapper::Gap" )) {
      my $exc_slice = Bio::EnsEMBL::Slice->new_fast({
        'start'             => $coord->start(),
        'end'               => $coord->end(),
        'strand'            => $slice->strand(),
        'coord_system'      => $slice->coord_system(),
         'adaptor'          => $self,
         'seq_region_name'  => $slice->seq_region_name(),
         'seq_region_length' => $slice->seq_region_length()});
      push( @out, bless ( [ $rel_start, $coord->length()+$rel_start-1,
                        $exc_slice ], "Bio::EnsEMBL::ProjectionSegment") );
    } else {
      my $exc_slice = $self->fetch_by_seq_region_id( $coord->id() );
      my $exc2_slice = Bio::EnsEMBL::Slice->new_fast({
        
         'start'             => $coord->start(),
         'end'               => $coord->end(),
         'strand'            => $coord->strand(),
         'seq_region_name'   => $exc_slice->seq_region_name(),
         'seq_region_length' => $exc_slice->seq_region_length(),
         'coord_system'      => $exc_slice->coord_system(),
         'adaptor'           => $self
        });
	
      push( @out, bless( [ $rel_start, $coord->length() + $rel_start - 1,
                    $exc2_slice ], "Bio::EnsEMBL::ProjectionSegment") );
    }
    $rel_start += $coord->length();
  }
  
  if($filter_projections) {
    return $self->_filter_Slice_projections($slice, \@out);
  }
  return \@out;
}

=head2 _filter_Slice_projections

    Arg [1]     : Bio::EnsEMBL::Slice The slice the projections were made from
    Arg [2]     : Array The projections which were fetched from the previous slice
    Description : Removes any projections which occur within the same sequence 
                  region as the given Slice object
    Returntype  : ArrayRef Bio::EnsEMBL::ProjectionSegment; Returns an array
                  of projected segments
=cut

sub _filter_Slice_projections {
  my ($self, $slice, $projections) = @_;
  my @proj = @{ $projections };
  if ( !@proj ) {
    throw('Was not given any projections to filter. Database may have incorrect assembly_exception information loaded');
  }
  
  # Want to get features on the FULL original slice as well as any
  # symlinked slices.
  
  # Filter out partial slices from projection that are on same
  # seq_region as original slice.

  my $sr_id = $slice->get_seq_region_id();

  @proj = grep { $_->to_Slice->get_seq_region_id() != $sr_id } @proj;

  my $segment = bless( [ 1, $slice->length(), $slice ],
                       'Bio::EnsEMBL::ProjectionSegment' );
  push( @proj, $segment );
  return \@proj;
}


=head2 store

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Arg [2]    : (optional) $seqref reference to a string
               The sequence associated with the slice to be stored.
  Example    : $slice = Bio::EnsEMBL::Slice->new(...);
               $seq_region_id = $slice_adaptor->store($slice, \$sequence);
  Description: This stores a slice as a sequence region in the database
               and returns the seq region id. The passed in slice must
               start at 1, and must have a valid seq_region name and coordinate
               system. The attached coordinate system must already be stored in
               the database.  The sequence region is assumed to start at 1 and
               to have a length equalling the length of the slice.  The end of
               the slice must equal the seq_region_length.
               If the slice coordinate system is the sequence level coordinate
               system then the seqref argument must also be passed.  If the
               slice coordinate system is NOT a sequence level coordinate
               system then the sequence argument cannot be passed.
  Returntype : int 
  Exceptions : throw if slice has no coord system.
               throw if slice coord system is not already stored.
               throw if slice coord system is seqlevel and no sequence is 
                     provided.
               throw if slice coord system is not seqlevel and sequence is
                     provided.
               throw if slice does not start at 1
               throw if sequence is provided and the sequence length does not
                     match the slice length.
               throw if the SQL insert fails (e.g. on duplicate seq region)
               throw if slice argument is not passed
               throw if the slice end is not equal to seq_region_length
  Caller     : database loading scripts
  Status     : Stable

=cut



sub store {
  my $self = shift;
  my $slice = shift;
  my $seqref = shift;
  my $not_dna = shift;

  #
  # Get all of the sanity checks out of the way before storing anything
  #

  if(!ref($slice) || !($slice->isa('Bio::EnsEMBL::Slice') or $slice->isa('Bio::EnsEMBL::LRGSlice'))) {
    throw('Slice argument is required');
  }

  my $cs = $slice->coord_system();
  throw("Slice must have attached CoordSystem.") if(!$cs);

  my $db = $self->db();
  if(!$cs->is_stored($db)) {
    throw("Slice CoordSystem must already be stored in DB.") 
  }

  if($slice->start != 1 || $slice->strand != 1) {
    throw("Slice must have start==1 and strand==1.");
  }

  if($slice->end() != $slice->seq_region_length()) {
    throw("Slice must have end==seq_region_length");
  }

  my $sr_len = $slice->length();
  my $sr_name  = $slice->seq_region_name();

  if($sr_name eq '') {
    throw("Slice must have valid seq region name.");
  }

  if($cs->is_sequence_level() && !$not_dna) {
    if(!$seqref) {
      throw("Must provide sequence for sequence level coord system.");
    }
    if(ref($seqref) ne 'SCALAR') {
      throw("Sequence must be a scalar reference.");
    }
    my $seq_len = length($$seqref);

    if($seq_len != $sr_len) {
      throw("Sequence length ($seq_len) must match slice length ($sr_len).");
    }
  } else {
    if($seqref) {
      throw("Cannot provide sequence for non-sequence level seq regions.");
    }
  }

  #store the seq_region

  my $sth = $db->dbc->prepare("INSERT INTO seq_region " .
                              "            ( name, length, coord_system_id ) " .
                              "     VALUES ( ?, ?, ? )"
      );

  $sth->bind_param(1,$sr_name,SQL_VARCHAR);
  $sth->bind_param(2,$sr_len,SQL_INTEGER);
  $sth->bind_param(3,$cs->dbID,SQL_INTEGER);

  $sth->execute();

  my $seq_region_id = $self->last_insert_id('seq_region_id', undef, 'seq_region');

  if(!$seq_region_id) {
    throw("Database seq_region insertion failed.");
  }

  if($cs->is_sequence_level() && !$not_dna) {
    #store sequence if it was provided
    my $seq_adaptor = $db->get_SequenceAdaptor();
    $seq_adaptor->store($seq_region_id, $$seqref);
  }

  #synonyms
  if(defined($slice->{'synonym'})){
    foreach my $syn (@{$slice->{'synonym'}} ){
      $syn->seq_region_id($seq_region_id); # set the seq_region_id
      $syn->adaptor->store($syn);
    }
  }
  
  
  $slice->adaptor($self);

  return $seq_region_id;
}


sub update {
  my $self = shift;
  my $slice = shift;

  #
  # Get all of the sanity checks out of the way before storing anything
  #

  if(!ref($slice) || !($slice->isa('Bio::EnsEMBL::Slice') or $slice->isa('Bio::EnsEMBL::LRGSlice'))) {
    throw('Slice argument is required');
  }

  my $cs = $slice->coord_system();
  throw("Slice must have attached CoordSystem.") if(!$cs);

  my $db = $self->db();

  if($slice->start != 1 || $slice->strand != 1) {
    throw("Slice must have start==1 and strand==1.");
  }

  if($slice->end() != $slice->seq_region_length()) {
    throw("Slice must have end==seq_region_length");
  }

  my $sr_len = $slice->length();
  my $sr_name  = $slice->seq_region_name();

  if(!$sr_name) {
    throw("Slice must have valid seq region name.");
  }

  #update the seq_region

  my $seq_region_id = $slice->get_seq_region_id();
  my $update_sql = qq(
     UPDATE seq_region
        SET name = ?,
            length = ?,
            coord_system_id = ?
      WHERE seq_region_id = ?
  );

  my $sth = $db->dbc->prepare($update_sql);

  $sth->bind_param(1,$sr_name,SQL_VARCHAR);
  $sth->bind_param(2,$sr_len,SQL_INTEGER);
  $sth->bind_param(3,$cs->dbID,SQL_INTEGER);

  $sth->bind_param(4, $seq_region_id, SQL_INTEGER);

  $sth->execute();

  #synonyms
  if(defined($slice->{'synonym'})){
    foreach my $syn (@{$slice->{'synonym'}} ){
      $syn->seq_region_id($seq_region_id); # set the seq_region_id
      my $syn_adaptor = $db->get_SeqRegionSynonymAdaptor();
      $syn_adaptor->store($syn);
    }
  }
}

=head2 remove

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               The slice to remove from the database
  Example    : $slice_adaptor->remove($slice);
  Description: Removes a slice completely from the database.
               All associated seq_region_attrib are removed as well.
               If dna is attached to the slice, it is also removed.
  Returntype : none
  Exceptions : throw if slice has no coord system.
               throw if slice argument is not passed
               warning if slice is not stored in this database
  Caller     : general
  Status     : Stable

=cut


sub remove {
  my $self = shift;
  my $slice = shift;

  #
  # Get all of the sanity checks out of the way before storing anything
  #

  if(!ref($slice) || !($slice->isa('Bio::EnsEMBL::Slice') or $slice->isa('Bio::EnsEMBL::LRGSlice'))) {
    throw('Slice argument is required');
  }

  my $cs = $slice->coord_system();
  throw("Slice must have attached CoordSystem.") if(!$cs);

  my $db = $self->db();
  my $seq_region_id = $slice->get_seq_region_id();

  if ($cs->is_sequence_level()) {
    my $seq_adaptor = $db->get_SequenceAdaptor();
    $seq_adaptor->remove($seq_region_id);
  }

  if(defined($slice->{'synonym'})){
    foreach my $syn (@{$slice->{'synonym'}} ){
      $syn->seq_region_id($seq_region_id); # set the seq_region_id
      $syn->adaptor->remove($syn);
    }
  }

  my $attrib_adaptor = $self->db->get_AttributeAdaptor();
  $attrib_adaptor->remove_from_Slice($slice);
  $self->remove_assembly($slice);

  my $sr_name  = $slice->seq_region_name();

  #remove the seq_region

  my $sth = $db->dbc->prepare("DELETE FROM seq_region " .
                         "WHERE name = ? AND " .
                         "      coord_system_id = ?" );

  $sth->bind_param(1,$sr_name,SQL_VARCHAR);
  $sth->bind_param(2,$cs->dbID,SQL_INTEGER);

  $sth->execute();

  return;
}


=head2 remove_assembly

  Arg [1]    : Bio::EnsEMBL::Slice $asm_slice or $cmp_slice
  Example    : $slice_adaptor->remove_assembly( $slice );
  Description: Deletes from the assembly table 
               where asm or cmp corresponds to slice
               Do not call this method unless you really know what you are doing
  Returntype : none
  Exceptions : throw if slice has no coord system (cs).
               throw if no slice provided, or argument is not a slice
  Caller     : Internal
  Status     : Experimental

=cut


sub remove_assembly {
  my $self      = shift;
  my $slice = shift;

  if(!ref($slice) || !($slice->isa('Bio::EnsEMBL::Slice') or $slice->isa('Bio::EnsEMBL::LRGSlice'))) {
    throw('Slice argument is required');
  }

  my $cs = $slice->coord_system();
  throw("Slice must have attached CoordSystem.") if(!$cs);

  #
  # Checks complete. Delete the data
  #
  my $sth = $self->db->dbc->prepare
      ("DELETE FROM assembly " .
       "WHERE   asm_seq_region_id = ? OR " .
       "        cmp_seq_region_id = ? ");

  my $asm_seq_region_id = $self->get_seq_region_id( $slice );
  my $cmp_seq_region_id = $self->get_seq_region_id( $slice );

  $sth->bind_param(1,$asm_seq_region_id,SQL_INTEGER);
  $sth->bind_param(2,$cmp_seq_region_id,SQL_INTEGER);

  $sth->execute();
  $sth->finish();

  return;

}

=head2 fetch_assembly

  Arg [1]    : Bio::EnsEMBL::Slice $asm_slice
  Arg [2]    : Bio::EnsEMBL::Slice $cmp_slice
  Example    : $asm = $slice_adaptor->fetch_assembly( $slice1, $slice2 );
  Description: Fetches from the assembly table based on the 
               coordinates of the two slices supplied. 
               Returns a mapper object mapping the two slices
               Do not call this method unless you really know what you are doing
  Returntype : Bio::EnsEMBL::Mapper
  Exceptions : throw if either slice has no coord system (cs).
               throw if there is no mapping path between coord systems
               throw if there are existing mappings between either slice
               and the oposite cs
  Caller     : Internal
  Status     : Experimental

=cut


sub fetch_assembly {
  my $self      = shift;
  my $asm_slice = shift;
  my $cmp_slice = shift;

  if(!ref($asm_slice) || !($asm_slice->isa('Bio::EnsEMBL::Slice') or $asm_slice->isa('Bio::EnsEMBL::LRGSlice'))) {
    throw('Assembled Slice argument is required');
  }
  if(!ref($cmp_slice) || !($cmp_slice->isa('Bio::EnsEMBL::Slice') or $cmp_slice->isa('Bio::EnsEMBL::LRGSlice')) ) {
    throw('Assembled Slice argument is required');
  }

  my $asm_cs = $asm_slice->coord_system();
  throw("Slice must have attached CoordSystem.") if(!$asm_cs);
  my $cmp_cs = $cmp_slice->coord_system();
  throw("Slice must have attached CoordSystem.") if(!$cmp_cs);

  my @path =
    @{ $asm_cs->adaptor()->get_mapping_path( $asm_cs, $cmp_cs ) };

  if ( !@path ) {
    throw("No mapping path defined between "
        . $asm_cs->name() . " and "
        . $cmp_cs->name() );
  }

  #
  # Checks complete. Fetch the data
  #
  my $sth = $self->db->dbc->prepare
      ("SELECT * FROM assembly " .
       "WHERE   asm_seq_region_id = ? AND " .
       "        cmp_seq_region_id = ? AND " .
       "        asm_start = ? AND " .
       "        asm_end   = ? AND " .
       "        cmp_start = ? AND " .
       "        cmp_end   = ? AND " .
       "        ori       = ?" );

  my $asm_seq_region_id = $self->get_seq_region_id( $asm_slice );
  my $cmp_seq_region_id = $self->get_seq_region_id( $cmp_slice );
  my $ori = $asm_slice->strand * $cmp_slice->strand;

  $sth->bind_param(1,$asm_seq_region_id,SQL_INTEGER);
  $sth->bind_param(2,$cmp_seq_region_id,SQL_INTEGER);
  $sth->bind_param(3,$asm_slice->start,SQL_INTEGER);
  $sth->bind_param(4,$asm_slice->end,SQL_INTEGER);
  $sth->bind_param(5,$cmp_slice->start,SQL_INTEGER);
  $sth->bind_param(6,$cmp_slice->end,SQL_INTEGER);
  $sth->bind_param(7,$ori,SQL_INTEGER);

  $sth->execute();

  my $results = $sth->fetchrow_array();
  $sth->finish();

  my $mapper;
  if ($results) {
    $mapper = Bio::EnsEMBL::Mapper->new($asm_slice, $cmp_slice, $asm_cs, $cmp_cs);
  }

  return $mapper;

}

=head2 store_assembly

  Arg [1]    : Bio::EnsEMBL::Slice $asm_slice
  Arg [2]    : Bio::EnsEMBL::Slice $cmp_slice
  Example    : $asm = $slice_adaptor->store_assembly( $slice1, $slice2 );
  Description: Creates an entry in the assembly table based on the 
               coordinates of the two slices supplied. Returns a string 
               representation of the assembly that gets created.
  Returntype : string
  Exceptions : throw if either slice has no coord system (cs).
               throw unless the cs rank of the asm_slice is lower than the 
               cmp_slice.
               throw if there is no mapping path between coord systems
               throw if the lengths of each slice are not equal
               throw if there are existing mappings between either slice
               and the oposite cs
  Caller     : database loading scripts
  Status     : Experimental

=cut

sub store_assembly{
  my $self = shift;
  my $asm_slice = shift;
  my $cmp_slice = shift;

  #
  # Get all of the sanity checks out of the way before storing anything
  #

  if(!ref($asm_slice) || !($asm_slice->isa('Bio::EnsEMBL::Slice') or $asm_slice->isa('Bio::EnsEMBL::LRGSlice'))) {
    throw('Assembled Slice argument is required');
  }
  if(!ref($cmp_slice) || !($cmp_slice->isa('Bio::EnsEMBL::Slice') or $cmp_slice->isa('Bio::EnsEMBL::LRGSlice')) ) {
    throw('Assembled Slice argument is required');
  }

  my $asm_cs = $asm_slice->coord_system();
  throw("Slice must have attached CoordSystem.") if(!$asm_cs);
  my $cmp_cs = $cmp_slice->coord_system();
  throw("Slice must have attached CoordSystem.") if(!$cmp_cs);

  my @path =
    @{ $asm_cs->adaptor()->get_mapping_path( $asm_cs, $cmp_cs ) };

  if ( !@path ) {
    throw("No mapping path defined between "
        . $asm_cs->name() . " and "
        . $cmp_cs->name() );
  }

  if( $asm_slice->length != $cmp_slice->length ){
    throw("The lengths of the assembled and component slices are not equal" );
  }

  # For now we disallow any existing mappings between the asm slice and cmp
  # CoordSystem and vice-versa. 
  # Some cases of multiple mappings may be allowable by the API, but their 
  # logic needs to be coded below.

  my $asm_proj = $asm_slice->project( $cmp_cs->name, $cmp_cs->version );
  if( @$asm_proj && $cmp_cs->name ne 'lrg' && $asm_cs->name ne 'lrg'){
    throw("Regions of the assembled slice are already assembled ".
          "into the component CoordSystem" ); 
  }
  my $cmp_proj = $cmp_slice->project( $asm_cs->name, $asm_cs->version );
  if( @$cmp_proj && $cmp_cs->name ne 'lrg' && $asm_cs->name ne 'lrg'){
    throw("Regions of the component slice are already assembled ".
          "into the assembled CoordSystem" ); 
  }

  #
  # Checks complete. Store the data
  #
  my $sth = $self->db->dbc->prepare
      ("INSERT INTO assembly " .
       "      ( asm_seq_region_id, " .
       "        cmp_seq_region_id, " .
       "        asm_start, " .
       "        asm_end  , " .
       "        cmp_start, " .
       "        cmp_end  , " .
       "        ori       )" .
       "VALUES ( ?, ?, ?, ?, ?, ?, ? )");

  my $asm_seq_region_id = $self->get_seq_region_id( $asm_slice );
  my $cmp_seq_region_id = $self->get_seq_region_id( $cmp_slice );
  my $ori = $asm_slice->strand * $cmp_slice->strand;

  $sth->bind_param(1,$asm_seq_region_id,SQL_INTEGER);
  $sth->bind_param(2,$cmp_seq_region_id,SQL_INTEGER);
  $sth->bind_param(3,$asm_slice->start,SQL_INTEGER);
  $sth->bind_param(4,$asm_slice->end,SQL_INTEGER);
  $sth->bind_param(5,$cmp_slice->start,SQL_INTEGER);
  $sth->bind_param(6,$cmp_slice->end,SQL_INTEGER);
  $sth->bind_param(7,$ori,SQL_INTEGER);

  $sth->execute();

  #use Data::Dumper qw( Dumper );
  #warn Dumper( $self->db->{seq_region_cache} );
  #$self->db->{seq_region_cache} = undef;
  #$self->_cache_seq_regions();

  my $ama = $self->db->get_AssemblyMapperAdaptor();
  $ama->delete_cache();


  return $asm_slice->name . "<>" . $cmp_slice->name;

}


=head2 prepare

  Arg [1]    : String $sql
  Example    :  ( optional )
  Description: overrides the default adaptor prepare method.
               All slice sql will usually use the dna_db.
  Returntype : DBD::sth 
  Exceptions : none
  Caller     : internal, convenience method
  Status     : Stable

=cut

sub prepare {
  my ( $self, $sql ) = @_;
  return $self->db()->dnadb()->dbc->prepare($sql);
}

sub _build_exception_cache {
  my $self = shift;

  # build up a cache of the entire assembly exception table
  # it should be small anyway
  my $sth =
    $self->prepare( 'SELECT ae.seq_region_id, ae.seq_region_start, '
              . 'ae.seq_region_end, ae.exc_type, ae.exc_seq_region_id, '
              . 'ae.exc_seq_region_start, ae.exc_seq_region_end '
              . 'FROM assembly_exception ae, '
              . 'seq_region sr, coord_system cs '
              . 'WHERE sr.seq_region_id = ae.seq_region_id '
              . 'AND sr.coord_system_id = cs.coord_system_id '
              . 'AND cs.species_id = ?' );

  $sth->bind_param( 1, $self->species_id(), SQL_INTEGER );
  $sth->execute();

  my %hash;
  $self->{'asm_exc_cache'} = \%hash;

  my $row;
  while ( $row = $sth->fetchrow_arrayref() ) {
    my @result = @$row;
    $hash{ $result[0] } ||= [];
    push( @{ $hash{ $result[0] } }, \@result );
  }
  $sth->finish();
} ## end sub _build_exception_cache

=head2 cache_toplevel_seq_mappings

  Args       : none
  Example    : $slice_adaptor->cache_toplevel_seq_mappings();
  Description: caches all the assembly mappings needed for genes
  Returntype : None
  Exceptions : None
  Caller     : general
  Status     : At Risk
             : New experimental code

=cut

sub cache_toplevel_seq_mappings {
  my ($self) = @_;

  # Get the sequence level to map too

  my $sql = (<<SSQL);
  SELECT    name
  FROM  coord_system
  WHERE attrib like "%sequence_level%"
  AND   species_id = ?
SSQL

  my $sth = $self->prepare($sql);
  $sth->bind_param( 1, $self->species_id(), SQL_INTEGER );
  $sth->execute();

  my $sequence_level = $sth->fetchrow_array();

  $sth->finish();

  my $csa = $self->db->get_CoordSystemAdaptor();
  my $ama = $self->db->get_AssemblyMapperAdaptor();

  my $cs1 = $csa->fetch_by_name($sequence_level);

  #get level to map too.

  $sql = (<<LSQL);
  SELECT DISTINCT(cs.name)
  FROM  seq_region sr,
        seq_region_attrib sra,
        attrib_type at,
        coord_system cs
  WHERE sra.seq_region_id = sr.seq_region_id
  AND   sra.attrib_type_id = at.attrib_type_id
  AND   at.code = "toplevel"
  AND   cs.coord_system_id = sr.coord_system_id
  AND   cs.species_id = ?
LSQL

  $sth = $self->prepare($sql);
  $sth->bind_param( 1, $self->species_id(), SQL_INTEGER );
  $sth->execute();

  while ( my $csn = $sth->fetchrow_array() ) {
    if ( $csn eq $sequence_level ) { next }
    my $cs2 = $csa->fetch_by_name($csn);
    my $am = $ama->fetch_by_CoordSystems( $cs1, $cs2 );
    $am->register_all();
  }

} ## end sub cache_toplevel_seq_mappings


sub _build_circular_slice_cache {
  my $self = shift;

  # build up a cache of circular sequence region ids
  my $sth =
    	$self->prepare( "SELECT sra.seq_region_id FROM seq_region_attrib sra "
		  	. "INNER JOIN attrib_type at ON sra.attrib_type_id = at.attrib_type_id "
			. "INNER JOIN seq_region sr ON sra.seq_region_id = sr.seq_region_id "
			. "INNER JOIN coord_system cs ON sr.coord_system_id = cs.coord_system_id "
			. "WHERE code = 'circular_seq' and cs.species_id = ?");

  $sth->bind_param( 1, $self->species_id(), SQL_INTEGER );
  $sth->execute();

  my $id;
  my %hash;
  if ( ($id) = $sth->fetchrow_array() ) {
  	$self->{'circular_sr_id_cache'} = \%hash;
        $self->{'is_circular'} = 1;
	$hash{ $id } = $id;
 	while ( ($id) = $sth->fetchrow_array() ) {
    		$hash{ $id } = $id;
  	}
  } else {
	$self->{'is_circular'} = 0;
  }
  $sth->finish();
} ## end _build_circular_slice_cache


1;
