#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor
#
# Copyright EMBL/EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor

=head1 SYNOPSIS

$rfa = $database_adaptor->get_RepeatFeatureAdaptor();

my $repeat = $rfa->fetch_by_dbID(1234);
my @repeats = @{$rfa->fetch_all_by_Slice($slice)};

=head1 DESCRIPTION

This is an adaptor for the retrieval and storage of RepeatFeature objects
from the RepeatFeature database.  Most of the implementation is in the
superclass BaseFeatureAdaptor.

=head1 AUTHOR - James Gilbert

=head1 CONTACT

Post questions to the EnsEMBL developer list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor;

use strict;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::RepeatFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);


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
     $analysis_id, $score, $repeat_name, $repeat_class,
     $repeat_consensus);

  $sth->bind_columns( \$repeat_feature_id, \$seq_region_id, \$seq_region_start,
                      \$seq_region_end, \$seq_region_strand,
                      \$repeat_consensus_id, \$repeat_start,\$repeat_end,
                      \$analysis_id, \$score, \$repeat_name, \$repeat_class,
                      \$repeat_consensus );

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
  if($dest_slice) {
    $dest_slice_start  = $dest_slice->start();
    $dest_slice_end    = $dest_slice->end();
    $dest_slice_strand = $dest_slice->strand();
    $dest_slice_length = $dest_slice->length();
  }

  FEATURE: while($sth->fetch()) {
    #create a repeat consensus object
    my $rc = $rc_hash{$repeat_consensus_id} ||=
      Bio::EnsEMBL::RepeatConsensus->new
          (-DBID             => $repeat_consensus_id,
           -ADAPTOR          => $rca,
           -NAME             => $repeat_name,
           -REPEAT_CLASS     => $repeat_class,
           -REPEAT_CONSENSUS => $repeat_consensus,
           -LENGTH           => length($repeat_consensus));

    #get the analysis object
    my $analysis = $analysis_hash{$analysis_id} ||=
      $aa->fetch_by_dbID($analysis_id);

    my $slice = $slice_hash{"ID:".$seq_region_id};

    if(!$slice) {
      $slice = $sa->fetch_by_seq_region_id($seq_region_id);
      $slice_hash{"ID:".$seq_region_id} = $slice;
      $sr_name_hash{$seq_region_id} = $slice->seq_region_name();
      $sr_cs_hash{$seq_region_id} = $slice->coord_system();
    }

    #
    # remap the feature coordinates to another coord system 
    # if a mapper was provided
    #
    if($mapper) {
      my $sr_name = $sr_name_hash{$seq_region_id};
      my $sr_cs   = $sr_cs_hash{$seq_region_id};

      ($sr_name,$seq_region_start,$seq_region_end,$seq_region_strand) =
        $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end,
                          $seq_region_strand, $sr_cs);

      #skip features that map to gaps or coord system boundaries
      next FEATURE if(!defined($sr_name));

      #get a slice in the coord system we just mapped to
      if($asm_cs == $sr_cs || ($cmp_cs != $sr_cs && $asm_cs->equals($sr_cs))) {
        $slice = $slice_hash{"NAME:$sr_name:$cmp_cs_name:$cmp_cs_vers"} ||=
          $sa->fetch_by_region($cmp_cs_name, $sr_name,undef, undef, undef,
                               $cmp_cs_vers);
      } else {
        $slice = $slice_hash{"NAME:$sr_name:$asm_cs_name:$asm_cs_vers"} ||=
          $sa->fetch_by_region($asm_cs_name, $sr_name, undef, undef, undef,
                               $asm_cs_vers);
      }
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

        #throw away features off the end of the requested slice
        if($seq_region_end < 1 || $seq_region_start > $dest_slice_length) {
          next FEATURE;
        }
      }
      $slice = $dest_slice;
    }

    #finally, create the new repeat feature
    push @features, Bio::EnsEMBL::RepeatFeature->new_fast
      ( { 'analysis'      =>  $analysis,
          'start'         =>  $seq_region_start,
          'end'           =>  $seq_region_end,
          'strand'        =>  $seq_region_strand,
          'score'         =>  $score,
          'hstart'        =>  $repeat_start,
          'hend'          =>  $repeat_end,
          'repeat_consensus' => $rc,
          'adaptor'       =>  $self,
          'slice'         =>  $slice,
          'dbID'          =>  $repeat_feature_id } );
  }

  return \@features;
}


=head2 store

  Arg [1]    : list of Bio::EnsEMBL::RepeatFeatures $repeat_feature_id
               the list of repeat features to store in the database
  Example    : $repeat_feature_adaptor->store(1234, @repeat_features);
  Description: stores a repeat feature in the database
  Returntype : none
  Exceptions : if the repeat features do not have attached sequences 
               or if repeat_consensus are not present 
  Caller     : general

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
      my @match = @{$rca->fetch_by_class_seq('trf', $cons->repeat_consensus)};
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
    if(!ref($slice) || !$slice->isa("Bio::EnsEMBL::Slice")) {
      throw("RepeatFeature cannot be stored without an associated slice.");
    }

    my $original = $rf;
    my $seq_region_id;
    ($rf, $seq_region_id) = $self->_pre_store($rf);

    $sth->execute($seq_region_id,
                  $rf->start,
                  $rf->end,
                  $rf->strand,
                  $rf->repeat_consensus->dbID(),
                  $rf->hstart,
                  $rf->hend,
                  $rf->score,
                  $rf->analysis->dbID);

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
  Returntype : list of ints
  Exceptions : none
  Caller     : ?

=cut

sub list_dbIDs {
   my ($self) = @_;

   return $self->_list_dbIDs("repeat_feature");
}

1;





