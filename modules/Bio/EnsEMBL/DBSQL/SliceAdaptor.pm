=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

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

  # include non-duplicate regions
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

  # retreive a list of slices from a file
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
use Scalar::Util qw/looks_like_number/;

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

  if ( !defined($start) )  { $start  = 1 }
  if ( !defined($strand) ) { $strand = 1 }

  if ( !defined($seq_region_name) ) {
    throw('seq_region_name argument is required');
  }

  my $cs;
  my $csa = $self->db->get_CoordSystemAdaptor();

  if ( defined($coord_system_name) ) {
    $cs = $csa->fetch_by_name( $coord_system_name, $version );

    ## REMOVE THESE THREE LINES WHEN STICKLEBACK DB IS FIXED!
    ## Anne/ap5 (2007-10-09):
    # The problem was that the stickleback genebuild called the
    # chromosomes 'groups', which meant they weren't being picked out by
    # the karyotype drawing code.  Apparently they are usually called
    # 'groups' in the stickleback community, even though they really are
    # chromosomes!

    if ( !defined($cs) && $coord_system_name eq 'chromosome' ) {
      $cs = $csa->fetch_by_name( 'group', $version );
    }

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

    if ( $sth->rows() == 0 ) {
      $sth->finish();


      # try synonyms
      my $syn_sql_sth = $self->prepare("select s.name, cs.name, cs.version from seq_region s join seq_region_synonym ss using (seq_region_id) join coord_system cs using (coord_system_id) where ss.synonym = ? and cs.species_id =?");
      $syn_sql_sth->bind_param(1, $seq_region_name, SQL_VARCHAR);
      $syn_sql_sth->bind_param(2, $self->species_id(), SQL_INTEGER);
      $syn_sql_sth->execute();
      my ($new_name, $new_coord_system, $new_version);
      $syn_sql_sth->bind_columns( \$new_name, \$new_coord_system, \$new_version);
            
      if($syn_sql_sth->fetch){
        $syn_sql_sth->finish;
        if (not defined($cs)) {
            return $self->fetch_by_region($new_coord_system, $new_name, $start, $end, $strand, $new_version, $no_fuzz);
        } elsif ($cs->dbID != $new_coord_system) {
            warning("Searched for a known feature on coordinate system: ".$cs->dbID." but found it on: ".$new_coord_system.
            "\n No result returned, consider searching without coordinate system or use toplevel.");
            return;
        }
        
      }
      $syn_sql_sth->finish;


      if ($no_fuzz) { return undef }

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
      if ( !defined($high_ver) ) { return undef }

    } else {

      my ( $id, $cs_id );
      ( $seq_region_name, $id, $length, $cs_id ) =
        $sth->fetchrow_array();
      $sth->finish();

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
  my ($self, $location, $no_warnings, $no_fuzz) = @_;
  return $self->fetch_by_location($location, 'toplevel', undef, $no_warnings, $no_fuzz);
}

=head2 fetch_by_location

  Arg [1]     : string $location
                Ensembl formatted location. Can be a format like 
                C<name:start-end>, C<name:start..end>, C<name:start:end>, 
                C<name:start>, C<name>. We can also support strand 
                specification as a +/- or 1/-1. 
                
                Location names must be separated by a C<:>. All others can be
                separated by C<..>, C<:> or C<->.
  Arg[2]      : String $coord_system_name
                The coordinate system to retrieve
  Arg[3]      : String $coord_system_version
                Optional parameter. Version of the coordinate system to fetch
  Arg[4]      : boolean $no_warnings
                Suppress warnings from this method
  Arg[5]      : boolean $no_fuzz
                Stop fuzzy matching of sequence regions from occuring
  Example     : my $slice = $sa->fetch_by_toplevel_location('X:1-10000')
                my $slice = $sa->fetch_by_toplevel_location('X:1-10000:-1')
  Description : Converts an Ensembl location/region into the sequence region
                name, start and end and passes them onto C<fetch_by_region()>. 
                The code assumes that the required slice is on the top level
                coordinate system. The code assumes that location formatting
                is not perfect and will perform basic cleanup before parsing.
  Returntype  : Bio::EnsEMBL::Slice
  Exceptions  : If $location or coordinate system is false otherwise 
                see C<fetch_by_region()>
  Caller      : General
  Status      : Beta

=cut

sub fetch_by_location {
  my ($self, $location, $coord_system_name, $coord_system_version, $no_warnings, $no_fuzz) = @_;
  
  throw "No coordinate system name specified" unless $coord_system_name;
  
  my ($seq_region_name, $start, $end, $strand) = $self->parse_location_to_values($location, $no_warnings);

  if(! $seq_region_name) {
    return;
  }
    
  if(defined $start && defined $end && $start > $end) {
    throw "Cannot request a slice whose start is greater than its end. Start: $start. End: $end";
  }
  
  my $slice = $self->fetch_by_region($coord_system_name, $seq_region_name, $start, $end, $strand, $coord_system_version, $no_fuzz);
  return unless $slice;
  
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
                separated by C<..>, C<:> or C<->.
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
  
  #cleanup any nomenclature like 1_000 or 1 000 or 1,000
  my $number_seps_regex = qr/\s+|,|_/;
  my $separator_regex = qr/(?:-|[.]{2}|\:)?/;
  my $number_regex = qr/[0-9,_ E]+/xms;
  my $strand_regex = qr/[+-1]|-1/xms;
  
  my $regex = qr/^((?:\w|\.|_|-)+) \s* :? \s* ($number_regex)? $separator_regex ($number_regex)? $separator_regex ($strand_regex)? $/xms;
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
        $start = 1;
      }
    }
    if(defined $end) {
      $end =~ s/$number_seps_regex//g;
      if($end < 1) {
        throw "Cannot request negative or 0 end indexes through this interface. Given $end but expected something greater than 0" unless $no_errors;
      }
    }
    
    if(defined $start && defined $end && $start > $end) {
      throw "Cannot request a slice whose start is greater than its end. Start: $start. End: $end" unless $no_errors;
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
  my ( $self, $seq_region_id, $start, $end, $strand ) = @_;

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

    if ( $sth->rows() == 0 ) { return undef }

    ( $name, $cs_id, $length ) = $sth->fetchrow_array();
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

  if($sth->rows() != 1) {
    throw("Non existant or ambigous seq_region:\n" .
          "  coord_system=[$cs_id],\n" .
          "  name=[$seq_region_name],\n");

  }

  my($seq_region_id, $length) = $sth->fetchrow_array();
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
               will be returned. If set, both reference and non refeference
               slices will be rerurned.
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

=head2 fetch_by_band

 Title   : fetch_by_band
 Usage   :
 Function: Does not work please use fetch_by_chr_band
 Example :
 Returns : Bio::EnsEMBL::Slice
 Args    : the band name
 Status     : AT RISK

=cut

sub fetch_by_band {
  my ($self,$band) = @_;

  my $sth = $self->dbc->prepare
        ("select s.name,max(k.seq_region_id)-min(k.seq_region_id, min(k.seq_region_start), max(k.seq_region_id) " .
         "from karyotype as k " .
         "where k.band like ? and k.seq_region_id = s.seq_region_id");

  $sth->bind_param(1,"$band%",SQL_VARCHAR);
  $sth->execute();
  my ( $seq_region_name, $discrepancy, $seq_region_start, $seq_region_end) = $sth->fetchrow_array;

  if($seq_region_name && $discrepancy>0) {
    throw("Band maps to multiple seq_regions");
  } else {
    return $self->fetch_by_region('toplevel',$seq_region_name,$seq_region_start,$seq_region_end);
  }
  throw("Band not recognised in database");
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
  Status     : Stable

=cut


sub fetch_normalized_slice_projection {
  my $self = shift;
  my $slice = shift;

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

    my $seq_reg_slice = $self->fetch_by_seq_region_id($slice_seq_region_id);
    my $exc_slice = $self->fetch_by_seq_region_id( $sort_haps[0][2] );
    my $len1 = $seq_reg_slice->length();
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

  return \@out;
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

  if(!$sr_name) {
    throw("Slice must have valid seq region name.");
  }

  if($cs->is_sequence_level()) {
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
                         "SET    name = ?, " .
                         "       length = ?, " .
                         "       coord_system_id = ?" );

  $sth->bind_param(1,$sr_name,SQL_VARCHAR);
  $sth->bind_param(2,$sr_len,SQL_INTEGER);
  $sth->bind_param(3,$cs->dbID,SQL_INTEGER);

  $sth->execute();

  my $seq_region_id = $sth->{'mysql_insertid'};

  if(!$seq_region_id) {
    throw("Database seq_region insertion failed.");
  }

  if($cs->is_sequence_level()) {
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


=head2 store_assembly

  Arg [1]    : Bio::EnsEMBL::Slice $asm_slice
  Arg [2]    : Bio::EnsEMBL::Slice $cmp_slice
  Example    : $asm = $slice_adaptor->store_assembly( $slice1, $slice2 );
  Description: Creates an entry in the analysis table based on the 
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

  unless( $asm_cs->rank < $cmp_cs->rank ){
    throw("Assembled Slice CoordSystem->rank must be lower than ".
          "the component Slice Coord_system" );
  }

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
  if( @$asm_proj ){
    throw("Regions of the assembled slice are already assembled ".
          "into the component CoordSystem" ); 
  }
  my $cmp_proj = $cmp_slice->project( $asm_cs->name, $asm_cs->version );
  if( @$cmp_proj ){
    throw("Regions of the component slice are already assembled ".
          "into the assembled CoordSystem" ); 
  }

  #
  # Checks complete. Store the data
  #
  my $sth = $self->db->dbc->prepare
      ("INSERT INTO assembly " .
       "SET     asm_seq_region_id = ?, " .
       "        cmp_seq_region_id = ?, " .
       "        asm_start = ?, " .
       "        asm_end   = ?, " .
       "        cmp_start = ?, " .
       "        cmp_end   = ?, " .
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
  my ($cs) = @{$csa->fetch_all()}; # get the highest coord system

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
  my $seq_level = $csa->fetch_sequence_level();

  my $seq_lvl_slice = $self->fetch_by_region($seq_level->name(), $name);

  if(!$seq_lvl_slice) {
    return undef;
  }

  my @projection = @{$seq_lvl_slice->project('toplevel')};

  if(@projection != 1) {
    warning("$name is mapped to multiple toplevel locations.");
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
  my $clone_cs = $csa->fetch_by_name('clone');

  if(!$clone_cs) {
    warning('Clone coordinate system does not exist for this species');
    return undef;
  }

  #this unfortunately needs a version on the end to work
  if(! ($name =~ /\./)) {
    my $sth =
      $self->prepare(  "SELECT sr.name "
                     . "FROM   seq_region sr, coord_system cs "
                     . "WHERE  cs.name = 'clone' "
                     . "AND    cs.coord_system_id = sr.coord_system_id "
                     . "AND    sr.name LIKE '$name.%'"
                     . "AND    cs.species_id = ?" );

    $sth->bind_param( 1, $self->species_id(), SQL_INTEGER );
    $sth->execute();

    if(!$sth->rows()) {
      $sth->finish();
      throw("Clone $name not found in database");
    }

    ($name) = $sth->fetchrow_array();

    $sth->finish();
  }

  my $clone = $self->fetch_by_region($clone_cs->name(), $name);
  return undef if(!$clone);

  my @projection = @{$clone->project('toplevel')};

  if(@projection != 1) {
    warning("$name is mapped to multiple locations.");
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
  my $sc_level = $csa->fetch_by_name('supercontig');

  if(!$sc_level) {
    warning('No supercontig coordinate system exists for this species.');
    return undef;
  }

  my $sc_slice = $self->fetch_by_region($sc_level->name(),$name);

  return undef if(!$sc_slice);

  my @projection = @{$sc_slice->project('toplevel')};

  if(@projection > 1) {
    warning("$name is mapped to multiple locations in toplevel");
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

   my $csa = $self->db->get_CoordSystemAdaptor();

   my $top_cs = @{$csa->fetch_all()};

   return $self->fetch_by_region($top_cs->name(),$chr_name,
                                 undef,undef,undef,$top_cs->version);
}






1;
