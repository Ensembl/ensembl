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

Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor

=head1 SYNOPSIS

  my $reg = 'Bio::EnsEMBL::Registry';

  $reg->
    load_registry_from_db( ...

  my $sfa =
    $reg->get_adaptor('homo sapiens', 'core', 'SimpleFeature');

  print ref($sfa), "\n";

  my $sf_aref =
    $sfa->fetch_all;

  print scalar @$sf_aref, "\n";

=head1 DESCRIPTION

Simple Feature Adaptor - database access for simple features

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);


=head2 store

  Arg [1]    : list of Bio::EnsEMBL::SimpleFeatures @sf
               the simple features to store in the database
  Example    : $simple_feature_adaptor->store(@simple_feats);
  Description: Stores a list of simple feature objects in the database
  Returntype : none
  Exceptions : thrown if @sf is not defined, if any of the features do not
               have an attached slice.
               or if any elements of @sf are not Bio::EnsEMBL::SimpleFeatures 
  Caller     : general
  Status     : Stable

=cut

sub store{
  my ($self,@sf) = @_;

  if( scalar(@sf) == 0 ) {
    throw("Must call store with list of SimpleFeatures");
  }

  my $sth = $self->prepare
    ("INSERT INTO simple_feature (seq_region_id, seq_region_start, " .
                                 "seq_region_end, seq_region_strand, " .
                                 "display_label, analysis_id, score) " .
     "VALUES (?,?,?,?,?,?,?)");

  my $db = $self->db();
  my $analysis_adaptor = $db->get_AnalysisAdaptor();

 FEATURE: foreach my $sf ( @sf ) {

    if( !ref $sf || !$sf->isa("Bio::EnsEMBL::SimpleFeature") ) {
      throw("SimpleFeature must be an Ensembl SimpleFeature, " .
            "not a [".ref($sf)."]");
    }

    if($sf->is_stored($db)) {
      warning("SimpleFeature [".$sf->dbID."] is already stored" .
              " in this database.");
      next FEATURE;
    }

    if(!defined($sf->analysis)) {
      throw("An analysis must be attached to the features to be stored.");
    }

    #store the analysis if it has not been stored yet
    if(!$sf->analysis->is_stored($db)) {
      $analysis_adaptor->store($sf->analysis());
    }

    my $original = $sf;
    my $seq_region_id;
    ($sf, $seq_region_id) = $self->_pre_store($sf);

    $sth->bind_param(1,$seq_region_id,SQL_INTEGER);
    $sth->bind_param(2,$sf->start,SQL_INTEGER);
    $sth->bind_param(3,$sf->end,SQL_INTEGER);
    $sth->bind_param(4,$sf->strand,SQL_TINYINT);
    $sth->bind_param(5,$sf->display_label,SQL_VARCHAR);
    $sth->bind_param(6,$sf->analysis->dbID,SQL_INTEGER);
    $sth->bind_param(7,$sf->score,SQL_DOUBLE);

    $sth->execute();

    $original->dbID($sth->{'mysql_insertid'});
    $original->adaptor($self);
  }
}


=head2 _tables

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               returns the names, aliases of the tables to use for queries
  Returntype : list of listrefs of strings
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub _tables {
  my $self = shift;
  
  return ['simple_feature', 'sf'];
}


=head2 _columns

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               returns a list of columns to use for queries
  Returntype : list of strings
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub _columns {
  my $self = shift;

  return qw( sf.simple_feature_id
             sf.seq_region_id sf.seq_region_start sf.seq_region_end
             sf.seq_region_strand sf.display_label sf.analysis_id sf.score );
}


=head2 _objs_from_sth

  Arg [1]    : hash reference $hashref
  Example    : none
  Description: PROTECTED implementation of superclass abstract method.
               creates SimpleFeatures from an executed DBI statement handle.
  Returntype : list reference to Bio::EnsEMBL::SimpleFeature objects
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub _objs_from_sth {
  my ($self, $sth, $mapper, $dest_slice) = @_;

  #
  # This code is ugly because an attempt has been made to remove as many
  # function calls as possible for speed purposes.  Thus many caches and
  # a fair bit of gymnastics is used.
  #

  my $sa = $self->db()->get_SliceAdaptor();
  my $aa = $self->db->get_AnalysisAdaptor();

  my @features;
  my %analysis_hash;
  my %slice_hash;
  my %sr_name_hash;
  my %sr_cs_hash;
  
  
  my($simple_feature_id,$seq_region_id, $seq_region_start, $seq_region_end,
     $seq_region_strand, $display_label, $analysis_id, $score);

  $sth->bind_columns(\$simple_feature_id,\$seq_region_id, \$seq_region_start,
                     \$seq_region_end, \$seq_region_strand, \$display_label,
                     \$analysis_id, \$score);
  
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
  my $dest_slice_seq_region_id;
  if($dest_slice) {
    $dest_slice_start  = $dest_slice->start();
    $dest_slice_end    = $dest_slice->end();
    $dest_slice_strand = $dest_slice->strand();
    $dest_slice_length = $dest_slice->length();
    $dest_slice_sr_name = $dest_slice->seq_region_name();
    $dest_slice_seq_region_id =$dest_slice->get_seq_region_id();
  }

  my $count =0;
  FEATURE: while($sth->fetch()) {
      $count++;
    #get the analysis object
    my $analysis = $analysis_hash{$analysis_id} ||=
      $aa->fetch_by_dbID($analysis_id);

    #need to get the internal_seq_region, if present
    $seq_region_id = $self->get_seq_region_id_internal($seq_region_id);
    #get the slice object
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
      if($asm_cs == $sr_cs || ($cmp_cs != $sr_cs && $asm_cs->equals($sr_cs))) {
        $slice = $slice_hash{"ID:".$seq_region_id} ||=
          $sa->fetch_by_seq_region_id($seq_region_id);
      } else {
        $slice = $slice_hash{"ID:".$seq_region_id} ||=
          $sa->fetch_by_seq_region_id($seq_region_id);
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
      }
       
      #throw away features off the end of the requested slice
      if($seq_region_end < 1 || $seq_region_start > $dest_slice_length ||
	( $dest_slice_seq_region_id != $seq_region_id )) {
#	print STDERR "IGNORED DUE TO CUTOFF  $dest_slice_seq_region_id ne $seq_region_id . $sr_name\n";
	next FEATURE;
      }
      $slice = $dest_slice;
    }

    push( @features,
          $self->_create_feature_fast(
                                    'Bio::EnsEMBL::SimpleFeature', {
                                      'start'    => $seq_region_start,
                                      'end'      => $seq_region_end,
                                      'strand'   => $seq_region_strand,
                                      'slice'    => $slice,
                                      'analysis' => $analysis,
                                      'adaptor'  => $self,
                                      'dbID'     => $simple_feature_id,
                                      'display_label' => $display_label,
                                      'score'         => $score
                                    } ) );

    }

  return \@features;
}


=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$simple_feature_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all simple features in the current db
  Arg[1]     : <optional> int. not 0 for the ids to be sorted by the seq_region.
  Returntype : list of ints
  Exceptions : none
  Caller     : ?
  Status     : Stable

=cut

sub list_dbIDs {
   my ($self, $ordered) = @_;

   return $self->_list_dbIDs("simple_feature", undef, $ordered);
}

1;
