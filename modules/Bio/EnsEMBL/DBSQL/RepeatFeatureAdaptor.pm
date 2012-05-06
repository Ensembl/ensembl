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

Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor

=head1 SYNOPSIS

  $rfa = $database_adaptor->get_RepeatFeatureAdaptor();

  my $repeat  = $rfa->fetch_by_dbID(1234);
  my @repeats = @{ $rfa->fetch_all_by_Slice($slice) };

=head1 DESCRIPTION

This is an adaptor for the retrieval and storage of RepeatFeature
objects from the database.  Most of the implementation is in the
superclass BaseFeatureAdaptor.

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor;

use strict;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::RepeatFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Scalar qw/wrap_array/;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);


=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Arg [2]    : (optional) string $logic_name
               Limits RepeatFeatures obtained to those having an Analysis with
               of the specified logic_name.  If no logic name is specified
               Repeats of all analysis types are retrieved.
  Arg [3]    : (optional) string/array $repeat_type
               Limits RepeatFeatures obtained to those of specified 
               repeat_type
  Example    : @rfeats = @{$rfa->fetch_all_by_Slice($slice, undef, 'Type II Transposons')};
               @rfeats = @{$rfa->fetch_all_by_Slice($slice, undef, ['Type II Transposons', 'RNA repeats'])};
  Description: Retrieves repeat features overlapping the area designated by
               the provided slice argument.  Returned features will be in
               in the same coordinate system as the provided slice and will
               have coordinates relative to the slice start.
  Returntype : reference to a list of Bio::EnsEMBL::RepeatFeatures.
  Exceptions : throw on bad argument
  Caller     : Slice::get_all_RepeatFeatures
  Status     : Stable

=cut

sub fetch_all_by_Slice {
  my $self = shift;
  my $slice = shift;
  my $logic_name = shift;
  my $repeat_type = shift;

  my $constraint = '';

  # MySQL was optimising the query the incorrect way when joining to
  # the repeat_consensus table on type
  $self->_straight_join(1);

  if($repeat_type) {
    my $rta = wrap_array($repeat_type);
    if(scalar(@{$rta}) > 1) {
      $constraint .= sprintf('rc.repeat_type IN (%s)', join(q{,}, map {"'${_}'"} @{$rta}));
    }
    else {
      $constraint .= "rc.repeat_type = '${repeat_type}'";      
    }
  }

  my $result =
    $self->fetch_all_by_Slice_constraint($slice,$constraint,$logic_name);


  $self->_straight_join(0);

  return $result;
}


#  _tablename
#
#   Arg [1]    : none
#   Example    : none
#   Description: PROTECTED Implementation of abstract superclass method to 
#                provide the name of the tables to query 
#   Returntype : string
#   Exceptions : none
#   Caller     : internal


sub _tables {
  my $self = shift;

  return (['repeat_feature', 'r'], ['repeat_consensus', 'rc']);
}


# _columns
#
#   Arg [1]    : none
#   Example    : none
#   Description: PROTECTED Implementation of abstract superclass method to 
#                provide the name of the columns to query 
#   Returntype : list of strings
#   Exceptions : none
#   Caller     : internal

sub _columns {
  my $self = shift;

  return qw (r.repeat_feature_id
	     r.seq_region_id
	     r.seq_region_start
	     r.seq_region_end
	     r.seq_region_strand
	     r.repeat_consensus_id
	     r.repeat_start
	     r.repeat_end
	     r.analysis_id
	     r.score
	     rc.repeat_name
	     rc.repeat_class
	     rc.repeat_type
	     rc.repeat_consensus);
}


# _default_where_clause
#  Arg [1]    : none
#  Example    : none
#  Description: Overrides superclass method to provide an additional 
#               table joining constraint before the SQL query is performed.
#  Returntype : string
#  Exceptions : none
#  Caller     : generic_fetch
#

sub _default_where_clause {
  my $self = shift;

  return 'r.repeat_consensus_id = rc.repeat_consensus_id';
}



#  Description: PROTECTED implementation of abstract superclass method.
#               responsible for the creation of RepeatFeatures from a
#               hashref generated from an SQL query

sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;

  #
  # This code is ugly because an attempt has been made to remove as many
  # function calls as possible for speed purposes.  Thus many caches and
  # a fair bit of gymnastics is used.
  #

  my $rca = $self->db()->get_RepeatConsensusAdaptor();
  my $sa = $self->db()->get_SliceAdaptor();
  my $aa = $self->db->get_AnalysisAdaptor();

  my @features;
  my %rc_hash;
  my %analysis_hash;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;

  my($repeat_feature_id, $seq_region_id, $seq_region_start, $seq_region_end,
     $seq_region_strand, $repeat_consensus_id, $repeat_start, $repeat_end,
     $analysis_id, $score, $repeat_name, $repeat_class, $repeat_type,
     $repeat_consensus);

  $sth->bind_columns( \$repeat_feature_id, \$seq_region_id, \$seq_region_start,
                      \$seq_region_end, \$seq_region_strand,
                      \$repeat_consensus_id, \$repeat_start,\$repeat_end,
                      \$analysis_id, \$score, \$repeat_name, \$repeat_class,
                      \$repeat_type, \$repeat_consensus );

  my $asm_cs;
  my $cmp_cs;
  my $asm_cs_vers;
  my $asm_cs_name;
  my $cmp_cs_vers;
  my $cmp_cs_name;
  if($mapper) {
    $asm_cs = $mapper->assembled_CoordSystem();
    $cmp_cs = $mapper->component_CoordSystem();
    $asm_cs_name = $asm_cs->name();
    $asm_cs_vers = $asm_cs->version();
    $cmp_cs_name = $cmp_cs->name();
    $cmp_cs_vers = $cmp_cs->version();
  }

  my $dest_slice_start;
  my $dest_slice_end;
  my $dest_slice_strand;
  my $dest_slice_length;
  my $dest_slice_sr_name;
  my $dest_slice_sr_id;
  if($dest_slice) {
    $dest_slice_start  = $dest_slice->start();
    $dest_slice_end    = $dest_slice->end();
    $dest_slice_strand = $dest_slice->strand();
    $dest_slice_length = $dest_slice->length();
    $dest_slice_sr_name = $dest_slice->seq_region_name();
    $dest_slice_sr_id   = $dest_slice->get_seq_region_id();
  }

  FEATURE: while($sth->fetch()) {
    #create a repeat consensus object

    my $rc = $rc_hash{$repeat_consensus_id} ||=
      Bio::EnsEMBL::RepeatConsensus->new_fast
          ({'dbID'             => $repeat_consensus_id,
            'adaptor'          => $rca,
            'name'             => $repeat_name,
            'repeat_class'     => $repeat_class,
	          'repeat_type'      => $repeat_type,
            'repeat_consensus' => $repeat_consensus,
            'length'           => length($repeat_consensus)});

    #get the analysis object
    my $analysis = $analysis_hash{$analysis_id} ||=
      $aa->fetch_by_dbID($analysis_id);
    #need to get the internal_seq_region, if present
    $seq_region_id = $self->get_seq_region_id_internal($seq_region_id);
    my $slice = $slice_hash{"ID:".$seq_region_id};

    if(!$slice) {
      $slice = $sa->fetch_by_seq_region_id($seq_region_id);
      $slice_hash{"ID:".$seq_region_id} = $slice;
      $sr_name_hash{$seq_region_id} = $slice->seq_region_name();
      $sr_cs_hash{$seq_region_id} = $slice->coord_system();
    }

    my $sr_name = $sr_name_hash{$seq_region_id};
    my $sr_cs   = $sr_cs_hash{$seq_region_id};
    #
    # remap the feature coordinates to another coord system 
    # if a mapper was provided
    #
    if($mapper) {

      if (defined $dest_slice && $mapper->isa('Bio::EnsEMBL::ChainedAssemblyMapper')  ) {
	    ( $seq_region_id,  $seq_region_start,
	      $seq_region_end, $seq_region_strand )
		=
		$mapper->map( $sr_name, $seq_region_start, $seq_region_end,
                          $seq_region_strand, $sr_cs, 1, $dest_slice);

      } else {

	    ( $seq_region_id,  $seq_region_start,
	      $seq_region_end, $seq_region_strand )
		=
		$mapper->fastmap( $sr_name, $seq_region_start, $seq_region_end,
                          $seq_region_strand, $sr_cs );
      }

      #skip features that map to gaps or coord system boundaries
      next FEATURE if(!defined($seq_region_id));

      #get a slice in the coord system we just mapped to
#      if($asm_cs == $sr_cs || ($cmp_cs != $sr_cs && $asm_cs->equals($sr_cs))) {
        $slice = $slice_hash{"ID:".$seq_region_id} ||=
          $sa->fetch_by_seq_region_id($seq_region_id);
#      } else {
#        $slice = $slice_hash{"NAME:$sr_name:$asm_cs_name:$asm_cs_vers"} ||=
#          $sa->fetch_by_region($asm_cs_name, $sr_name, undef, undef, undef,
#                               $asm_cs_vers);
#      }
    }

    #
    # If a destination slice was provided convert the coords
    # If the dest_slice starts at 1 and is foward strand, nothing needs doing
    #
    if($dest_slice) {
      if($dest_slice_start != 1 || $dest_slice_strand != 1) {
        if($dest_slice_strand == 1) {
          $seq_region_start = $seq_region_start - $dest_slice_start + 1;
          $seq_region_end   = $seq_region_end   - $dest_slice_start + 1;
        } else {
          my $tmp_seq_region_start = $seq_region_start;
          $seq_region_start = $dest_slice_end - $seq_region_end + 1;
          $seq_region_end   = $dest_slice_end - $tmp_seq_region_start + 1;
          $seq_region_strand *= -1;
        }
      }

      #throw away features off the end of the requested slice
      if($seq_region_end < 1 || $seq_region_start > $dest_slice_length ||
	( $dest_slice_sr_id ne $seq_region_id )) {
	next FEATURE;
      }
      $slice = $dest_slice;
    }

    # Finally, create the new RepeatFeature.
    push( @features,
          $self->_create_feature_fast( 'Bio::EnsEMBL::RepeatFeature', {
                                         'dbID' => $repeat_feature_id,
                                         'analysis' => $analysis,
                                         'start'  => $seq_region_start,
                                         'end'    => $seq_region_end,
                                         'strand' => $seq_region_strand,
                                         'score'  => $score,
                                         'hstart' => $repeat_start,
                                         'hend'   => $repeat_end,
                                         'repeat_consensus' => $rc,
                                         'adaptor'          => $self,
                                         'slice'            => $slice
                                       } ) );

  }

  return \@features;
}


=head2 store

  Arg [1]    : list of Bio::EnsEMBL::RepeatFeatures $repeat_feature_id
               the list of repeat features to store in the database
  Example    : $repeat_feature_adaptor->store(@repeat_features);
  Description: stores a repeat feature in the database
  Returntype : none
  Exceptions : if the repeat features do not have attached sequences 
               or if repeat_consensus are not present 
  Caller     : general
  Status     : Stable

=cut

sub store {
  my( $self, @repeats ) = @_;

  my $db = $self->db();
  my $rca = $db->get_RepeatConsensusAdaptor();
  my $sa  = $db->get_SliceAdaptor();
  my ($cons, $db_id);

  my $sth = $self->prepare(qq{
    INSERT into repeat_feature( repeat_feature_id
				, seq_region_id
				, seq_region_start
				, seq_region_end
				, seq_region_strand
				, repeat_consensus_id
				, repeat_start
				, repeat_end
				, score
				, analysis_id )
      VALUES(NULL, ?,?,?,?,?,?,?,?,?)
    });

 FEATURE: foreach my $rf (@repeats) {
    if(!ref($rf) || !$rf->isa('Bio::EnsEMBL::RepeatFeature')) {
      throw('Expected RepeatFeature argument not [' . ref($rf) .'].');
    }

    if($rf->is_stored($db)) {
      warning("RepeatFeature [".$rf->dbID."] is already stored in this DB.");
      next FEATURE;
    }

    my $cons = $rf->repeat_consensus();
    throw("Must have a RepeatConsensus attached") if(!defined($cons));

    # for tandem repeats - simply store consensus and repeat
    # one pair per hit. don't need to check consensi stored
    # already. consensus has name and class set to 'trf'

    if ($cons->repeat_class eq 'trf') {

      # Look for matches already stored
      my @match = @{$rca->fetch_all_by_class_seq('trf', $cons->repeat_consensus)};
      if (@match) {
        $cons->dbID($match[0]->dbID());
      }
      else {
        $rca->store($cons);
      }

    } elsif ($cons->repeat_class eq 'Simple_repeat') {

      my $rcon = $cons->name;
      $rcon =~ s/\((\S+)\)n/$1/;   # get repeat element
      $cons->repeat_consensus($rcon);

      # Look for matches already stored
      my $match = $rca->fetch_by_name_class($cons->name, 'Simple_repeat');
      if ($match) {
        $cons->dbID($match->dbID());
      }
      else {
        $rca->store($cons);
      }

    } else {

      # for other repeats - need to see if a consensus is stored already
      if(!$cons->dbID) {
        my $match = ($rca->fetch_by_name($cons->name));
	
        if($match) {
          #set the consensus dbID to be the same as the database one
          $cons->dbID($match->dbID());
        } else {
          # if we don't match a consensus already stored create a fake one 
          # and set consensus to 'N' as null seq not allowed
          # FIXME: not happy with this, but ho hum ...
          warning("Can't find " . $cons->name . "\n");
          $cons->repeat_consensus("N");
          $rca->store($cons);
        }
      }
	
      #if (@match > 1) {
      #multiple consensi were matched
      #  $self->warn(@match . " consensi for " . $cons->name . "\n");
      #}
    }

    my $slice = $rf->slice();
    if(!ref($slice) || !($slice->isa("Bio::EnsEMBL::Slice") or $slice->isa('Bio::EnsEMBL::LRGSlice'))) {
      throw("RepeatFeature cannot be stored without an associated slice.");
    }

    my $original = $rf;
    my $seq_region_id;
    ($rf, $seq_region_id) = $self->_pre_store($rf);

    $sth->bind_param(1,$seq_region_id,SQL_INTEGER);
    $sth->bind_param(2,$rf->start,SQL_INTEGER);
    $sth->bind_param(3,$rf->end,SQL_INTEGER);
    $sth->bind_param(4,$rf->strand,SQL_TINYINT);
    $sth->bind_param(5,$rf->repeat_consensus->dbID,SQL_INTEGER);
    $sth->bind_param(6,$rf->hstart,SQL_INTEGER);
    $sth->bind_param(7,$rf->hend,SQL_INTEGER);
    $sth->bind_param(8,$rf->score,SQL_DOUBLE);
    $sth->bind_param(9,$rf->analysis->dbID,SQL_INTEGER);

    $sth->execute();

    my $db_id = $sth->{'mysql_insertid'}
      or throw("Didn't get an insertid from the INSERT statement");

    $original->dbID($db_id);
    $original->adaptor($self);
  }
}


=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$repeat_feature_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all repeat features in the current db
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub list_dbIDs {
   my ($self, $ordered) = @_;

   return $self->_list_dbIDs("repeat_feature", undef, $ordered);
}

1;





