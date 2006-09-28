# EnsEMBL module for DitagFeatureAdaptor
#
# Copyright EMBL-EBI/Wellcome Trust Sanger Center 2006
#
# You may distribute this module under the same terms as perl itself
#
# Cared for by EnsEMBL (ensembl-dev@ebi.ac.uk)

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Map::DBSQL::DitagFeatureAdaptor

=head1 SYNOPSIS

my $dfa = $db->get_DitagFeatureAdaptor;
my $ditagFeatures = $dfa->fetch_all_by_Slice($slice, "SME005");

foreach my $ditagFeature (@$ditagFeatures){
  print $ditagFeature->ditag_id . " " .
        $ditagFeature->slice    . " " . $ditagFeature->start . "-" .
        $ditagFeature->end      . " " . $ditagFeature->strand;
}

=head1 DESCRIPTION

Provides database interaction for the Bio::EnsEMBL::Map::DitagFeature object

=cut

package Bio::EnsEMBL::Map::DBSQL::DitagFeatureAdaptor;

use strict;
use vars ('@ISA');

use Bio::EnsEMBL::Map::Ditag;
use Bio::EnsEMBL::Map::DitagFeature;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 fetch_all

  Arg [1]    : none
  Example    : @all_tags = @{$ditagfeature_adaptor->fetch_all};
  Description: Retrieves all ditagFeatures from the database
  Returntype : listref of Bio::EnsEMBL::Map::DitagFeature
  Caller     : general

=cut

sub fetch_all {
  my $self = shift;

  my $sth = $self->prepare("SELECT ditag_feature_id, ditag_id, seq_region_id, seq_region_start, 
                            seq_region_end, seq_region_strand, analysis_id, hit_start, hit_end, 
                            hit_strand, cigar_line, ditag_side, ditag_pair_id 
                            FROM   ditag_feature" );
  $sth->execute;

  my $result = $self->_fetch($sth);

  return $result;
}


=head2 fetch_by_dbID

  Arg [1]    : ditagFeature dbID
  Example    : @my_tags = @{$ditagfeature_adaptor->fetch_by_dbID($my_id)};
  Description: Retrieves a ditagFeature from the database.
  Returntype : Bio::EnsEMBL::Map::DitagFeature
  Caller     : general

=cut

sub fetch_by_dbID {
  my ($self, $dbid) = @_;

  my $sth = $self->prepare("SELECT ditag_feature_id, ditag_id, seq_region_id, seq_region_start, 
                            seq_region_end, seq_region_strand, analysis_id, hit_start, hit_end, 
                            hit_strand, cigar_line, ditag_side, ditag_pair_id 
                            FROM   ditag_feature
                            WHERE  ditag_feature_id = ?" );
  $sth->execute($dbid);

  my $result = $self->_fetch($sth);

  return $result->[0];
}



=head2 fetch_all_by_ditagID

  Arg [1]    : ditag dbID
  Arg [2]    : (optional) ditag-pair dbID
  Arg [3]    : (optional) analysis ID
  Example    : @my_tags = @{$ditagfeature_adaptor->fetch_all_by_ditag_id($my_id)};
  Description: Retrieves all ditagFeatures from the database linking to a specific ditag-id
  Returntype : listref of Bio::EnsEMBL::Map::DitagFeature
  Caller     : general

=cut

sub fetch_all_by_ditagID {
  my ($self, $ditag_id, $ditag_pair_id, $analysis_id) = @_;

  my $arg = $ditag_id;
  my $sql = "SELECT ditag_feature_id, ditag_id, seq_region_id, seq_region_start, 
             seq_region_end, seq_region_strand, analysis_id, hit_start, hit_end, 
             hit_strand, cigar_line, ditag_side, ditag_pair_id 
             FROM   ditag_feature 
             WHERE  ditag_id = ? ";
  if($ditag_pair_id){
    $sql .= "AND ditag_pair_id = ? ";
    $arg .= ", $ditag_pair_id";
  }
  if($analysis_id){
    $sql .= "AND analysis_id = ? ";
    $arg .= ", $analysis_id";
  }
  $sql   .= "ORDER BY ditag_pair_id";
  my $sth = $self->prepare($sql);
  $sth->execute(split(",",$arg));

  my $result = $self->_fetch($sth);

  return $result;
}


=head2 fetch_all_by_type

  Arg [1]    : ditag type
  Example    : @my_tags = @{$ditagfeature_adaptor->fetch_all_by_type($type)};
  Description: Retrieves all ditagFeatures from the database linking to a specific ditag-type
  Returntype : listref of Bio::EnsEMBL::Map::DitagFeature
  Caller     : general

=cut

sub fetch_all_by_type {
  my ($self, $ditag_type) = @_;

  my $sth = $self->prepare("SELECT df.ditag_feature_id, df.ditag_id, df.seq_region_id, 
                            df.seq_region_start, df.seq_region_end, df.seq_region_strand, 
                            df.analysis_id, df.hit_start, df.hit_end, df.hit_strand, 
                            df.cigar_line, df.ditag_side, ditag_pair_id 
                            FROM   ditag_feature df, ditag d 
                            WHERE  df.ditag_id=d.ditag_id and d.type = ? 
                            ORDER BY df.ditag_id, df.ditag_pair_id" );
  $sth->execute($ditag_type);

  my $result = $self->_fetch($sth);

  return $result;
}



=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : (optional) ditag type (specific library)
  Arg [3]    : (optional) analysis logic_name
  Example    : $tags = $ditagfeature_adaptor->fetch_all_by_Slice($slice, "SME005");
  Description: Retrieves ditagFeatures from the database overlapping a specific region
               and (optional) of a specific ditag type or analysis.
  Returntype : listref of Bio::EnsEMBL::Map::DitagFeatures
  Caller     : general

=cut

sub fetch_all_by_Slice {
  my ($self, $slice, $tagtype, $logic_name) = @_;

  my @result;
  my $moresql;

  if(!ref($slice) || !$slice->isa("Bio::EnsEMBL::Slice")) {
    throw("Bio::EnsEMBL::Slice argument expected not $slice.");
  }

  #get affected ditag_feature_ids
  my $sql = "SELECT df.ditag_feature_id ".
            "FROM ditag_feature df ";
  if($tagtype){
    $sql .= ", ditag d ".
            "WHERE df.ditag_id=d.ditag_id ".
            "AND d.type = \"".$tagtype."\" AND ";
  }
  else{
    $sql .= "WHERE ";
  }
  if($logic_name){
    my $analysis = $self->db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
    if(!$analysis) {
      return undef;
    }
    $sql .= "df.analysis_id = ".$analysis->dbID()." AND ";
  }
  $sql .= "df.seq_region_id = ".$slice->get_seq_region_id.
          " AND df.seq_region_start <= ".$slice->end.
	  " AND df.seq_region_end   >= ".$slice->start;

  my $sth = $self->prepare($sql);
  $sth->execute();
  my @id_list = map {$_->[0]} @{$sth->fetchall_arrayref([0],undef)};

  #fetch ditagFeatures for these ids
  #whith splitting large queries into smaller batches
  my $max_size     = 1000;
  my $ids_to_fetch = "";

  while(@id_list) {
    my @ids;
    if(@id_list > $max_size) {
      @ids = splice(@id_list, 0, $max_size);
    } else {
      @ids = splice(@id_list, 0);
    }
    $ids_to_fetch = join(', ', @ids);

    my $sth = $self->prepare("SELECT ditag_feature_id, ditag_id, seq_region_id, seq_region_start, 
                              seq_region_end, seq_region_strand, analysis_id, hit_start, hit_end, 
                              hit_strand, cigar_line, ditag_side, ditag_pair_id 
                              FROM   ditag_feature
                              WHERE  ditag_feature_id IN(".$ids_to_fetch.")" );
    $sth->execute();

    my $result = $self->_fetch($sth);
    push(@result, @$result);
  }

  return \@result;
}


=head2 _fetch

  Arg [1]    : statement handler
  Description: generic sql-fetch function for the DitagFeature fetch methods
  Returntype : listref of Bio::EnsEMBL::Map::DitagFeatures
  Caller     : private

=cut

sub _fetch {
  my ($self, $sth) = @_;

  my ( $tag_id, $mothertag_id, $seqreg, $seqstart, $seqend, $strand, $analysis_id, $hit_start,
       $hit_end, $hit_strand, $cigar_line, $ditag_side, $ditag_pair_id );
  $sth->bind_columns( \$tag_id,        \$mothertag_id, \$seqreg,
                      \$seqstart,      \$seqend,       \$strand,
                      \$analysis_id,   \$hit_start,    \$hit_end,
                      \$hit_strand,    \$cigar_line,   \$ditag_side,
                      \$ditag_pair_id );

  my @ditags;

  while ( $sth->fetch ) {
    my $analysis_obj = $self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id);
    my $slice        = $self->db->get_SliceAdaptor->fetch_by_seq_region_id($seqreg);

    push @ditags,
      Bio::EnsEMBL::Map::DitagFeature->new( -dbid          => $tag_id,
                                            -slice         => $slice,
                                            -start         => $seqstart,
                                            -end           => $seqend,
                                            -strand        => $strand, 
                                            -analysis      => $analysis_obj,
                                            -hit_start     => $hit_start,
                                            -hit_end       => $hit_end,
                                            -hit_strand    => $hit_strand,
                                            -ditag_id      => $mothertag_id,
                                            -cigar_line    => $cigar_line,
                                            -ditag_side    => $ditag_side,
					    -ditag_pair_id => $ditag_pair_id,
                                            -adaptor       => $self,
                                            );
  }

  return \@ditags;
}


=head2 sequence

  Arg [1]    : dbID of DitagFeature
  Example    : $ditagfeature_adaptor->get_sequence($ditagFeature->dbID)
  Description: get the part of the sequence of a ditag,
               that is actully aligned to the genome.
  Returntype : string
  Exceptions : thrown if not all data needed for storing is populated in the
               ditag features
  Caller     : Bio::EnsEMBL::Map::DitagFeature

=cut

sub sequence {
  my ($self, $dbID) = @_;

  my $sequence = undef;
  my $db  = $self->db() or throw "Couldn t get database connection.";
  my $sql = "SELECT d.sequence, df.hit_start, df.hit_end, df.hit_strand ".
            "FROM ditag d, ditag_feature df ".
	    "WHERE df.ditag_id=d.ditag_id and df.ditag_feature_id = ?";
  my $sth = $db->dbc->prepare($sql);
  $sth->execute( $dbID );
  my ($seq, $start, $end, $strand) = $sth->fetchrow_array();
  if($seq and $start and $end and $strand){
    $sequence = substr($seq, ($start-1), ($end-$strand));
    if($strand == -1) {
      $sequence =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
    }
  }

  return $sequence;
}


=head2 store

  Arg [1]    : (Array ref of) Bio::EnsEMBL::Map::DitagFeature
  Example    : $ditagfeature_adaptor->store(@ditag_features);
  Description: Stores a single ditagFeature or
               a list of ditagFeatures in this database.
  Returntype : none
  Exceptions : thrown if not all data needed for storing is populated in the
               ditag features
  Caller     : general

=cut

sub store {
  my ( $self, $ditags ) = @_;

  if ( ref $ditags eq 'ARRAY' ) {
    if ( scalar(@$ditags) == 0 ) {
      throw( "Must call store with ditagFeature or list ref of ditagsFeature" );
    }
  } elsif ($ditags) {
    my @ditags;
    push @ditags, $ditags;
    $ditags = \@ditags;
  } else {
    throw( "Must call store with ditagFeature or list ref of ditagsFeature." );
  }

  my $db = $self->db() or throw "Couldn t get database connection.";

  my $sth1 = $self->prepare( "INSERT INTO ditag_feature( ditag_id, seq_region_id, seq_region_start, 
                              seq_region_end, seq_region_strand, analysis_id, hit_start, hit_end, 
                              hit_strand, cigar_line, ditag_side, ditag_pair_id ) 
                              VALUES( ?,?,?,?,?,?,?,?,?,?,?,? )" );
  my $sth2 = $self->prepare( "INSERT INTO ditag_feature( ditag_feature_ID, ditag_id, seq_region_id, 
                              seq_region_start, seq_region_end, seq_region_strand, analysis_id, hit_start, 
                              hit_end, hit_strand, cigar_line, ditag_side, ditag_pair_id ) 
                              VALUES( ?,?,?,?,?,?,?,?,?,?,?,?,? )" );
  my $sth3 = $self->prepare( "SELECT COUNT(*) FROM ditag_feature 
                              WHERE ditag_id = ?" );

TAG:
  foreach my $ditag (@$ditags) {

    if ( !ref $ditag || !$ditag->isa("Bio::EnsEMBL::Map::DitagFeature") ) {
      throw(   "Object must be an Ensembl DitagFeature, "
             . "not a " . ref($ditag) );
    }

    if ( $ditag->is_stored($db) ) {
      warning(   "DitagFeature " . $ditag->dbID .
                 " is already stored in this database,".
                 " maybe you ned to use the update() method." );
      next TAG;
    }

#    #check if more than x tags with this ditag id exist
#    $sth3->execute( $ditag->ditag_id );
#    my ($num) = $sth3->fetchrow_array();
#    if ( ($num) and ($num > 1) ) {
#      warning( "There are already at least 2 DitagFeatures relating to Ditag ".
#                 $ditag->ditag_id." stored in this database." );
#      if ( $num > 4 ) {
#        warning( "not storing" );
#        next TAG;
#      }
#    }

    if ( $ditag->dbID ) {
      $sth2->bind_param( 1,  $ditag->dbID,                      SQL_INTEGER );
      $sth2->bind_param( 2,  $ditag->ditag_id,                  SQL_INTEGER );
      $sth2->bind_param( 3, ($ditag->slice->get_seq_region_id), SQL_INTEGER );
      $sth2->bind_param( 4,  $ditag->start,                     SQL_INTEGER );
      $sth2->bind_param( 5,  $ditag->end,                       SQL_INTEGER );
      $sth2->bind_param( 6,  $ditag->strand,                    SQL_VARCHAR );
      $sth2->bind_param( 7,  $ditag->analysis->dbID,            SQL_INTEGER );
      $sth2->bind_param( 8,  $ditag->hit_start,                 SQL_INTEGER );
      $sth2->bind_param( 9,  $ditag->hit_end,                   SQL_INTEGER );
      $sth2->bind_param( 10, $ditag->hit_strand,                SQL_VARCHAR );
      $sth2->bind_param( 11, $ditag->cigar_line,                SQL_VARCHAR );
      $sth2->bind_param( 12, $ditag->ditag_side,                SQL_VARCHAR );
      $sth2->bind_param( 13, $ditag->ditag_pair_id,             SQL_VARCHAR );
      $sth2->execute();
    }
    else{
      $sth1->bind_param( 1,  $ditag->ditag_id,                  SQL_INTEGER );
      $sth1->bind_param( 2, ($ditag->slice->get_seq_region_id), SQL_INTEGER );
      $sth1->bind_param( 3,  $ditag->start,                     SQL_INTEGER );
      $sth1->bind_param( 4,  $ditag->end,                       SQL_INTEGER );
      $sth1->bind_param( 5,  $ditag->strand,                    SQL_VARCHAR );
      $sth1->bind_param( 6,  $ditag->analysis->dbID,            SQL_INTEGER );
      $sth1->bind_param( 7,  $ditag->hit_start,                 SQL_INTEGER );
      $sth1->bind_param( 8,  $ditag->hit_end,                   SQL_INTEGER );
      $sth1->bind_param( 9,  $ditag->hit_strand,                SQL_VARCHAR );
      $sth1->bind_param( 10, $ditag->cigar_line,                SQL_VARCHAR );
      $sth1->bind_param( 11, $ditag->ditag_side,                SQL_VARCHAR );
      $sth1->bind_param( 12, $ditag->ditag_pair_id,             SQL_VARCHAR );
      $sth1->execute();
      my $dbID = $sth1->{'mysql_insertid'};
      $ditag->dbID($dbID);
      $ditag->adaptor($self);
    }

  }
}


=head2 batch_store

  Arg [1]    : (Array ref of) Bio::EnsEMBL::Map::DitagFeatures
  Arg [2]    : bool have_dbIDs
  Example    : $ditagfeature_adaptor->batch_store(\@ditag_features);
  Description: Stores a list of ditagFeatures in this database.
               DitagFeatures are expected to have no dbID yet unless flag "have_dbIDs" is true.
               They are inserted in one combined INSERT for better performance.
  Returntype : none
  Exceptions : thrown if not all data needed for storing is given for the
               ditag features
  Caller     : general

=cut

sub batch_store {
  my ( $self, $ditags, $have_dbIDs ) = @_;

  my %tag_ids = ();
  my @good_ditags;
  my ($sql, $sqladd);
  my $inserts = 0;

  if ( ref $ditags eq 'ARRAY' ) {
    if ( scalar(@$ditags) == 0 ) {
      throw( "Must call store with ditagFeature or list ref of ditagsFeature" );
    }
  } elsif ($ditags) {
    my @ditags;
    push @ditags, $ditags;
    $ditags = \@ditags;
  } else {
    throw( "Must call store with ditagFeature or list ref of ditagsFeature." );
  }

  my $db = $self->db() or throw "Couldn t get database connection.";

  #check whether it s a DitagFeature object and is not stored already
  foreach my $ditag (@$ditags) {

    if ( !ref $ditag || !$ditag->isa("Bio::EnsEMBL::Map::DitagFeature") ) {
      throw(   "Object must be an Ensembl DitagFeature, "
             . "not a " . ref($ditag) );
    }
    if ( $ditag->is_stored($db) ) {
      warning(   "DitagFeature " . $ditag->dbID
                 . " is already stored in this database." );
      next;
    }
    $tag_ids{$ditag->ditag_id} = $ditag;
    push(@good_ditags, $ditag);
  }
  $ditags = undef;

  #create batch INSERT
  if($have_dbIDs){
    $sql = "INSERT INTO ditag_feature ( ditag_feature_id, ditag_id, seq_region_id, seq_region_start, ".
           "seq_region_end, seq_region_strand, analysis_id, hit_start, hit_end, ".
           "hit_strand, cigar_line, ditag_side, ditag_pair_id ) VALUES ";
    foreach my $ditag (@good_ditags) {
      $sqladd = "";
      if($inserts){ $sqladd = ", " }
      $sqladd .= "(". $ditag->ditag_feature_id.", ".$ditag->ditag_id.", ".($ditag->slice->get_seq_region_id).
	         ", ". $ditag->start.", ".$ditag->end.", '".$ditag->strand."', ".$ditag->analysis->dbID.", ".
                 $ditag->hit_start.", ".$ditag->hit_end.", '".$ditag->hit_strand."', '".$ditag->cigar_line."', '".
                 $ditag->ditag_side."', ".$ditag->ditag_pair_id.")";
      $sql .= $sqladd;
      $inserts++;
    }
  }
  else{
    $sql = "INSERT INTO ditag_feature ( ditag_id, seq_region_id, seq_region_start, ".
           "seq_region_end, seq_region_strand, analysis_id, hit_start, hit_end, ".
           "hit_strand, cigar_line, ditag_side, ditag_pair_id ) VALUES ";
    foreach my $ditag (@good_ditags) {
      $sqladd = "";
      if($inserts){ $sqladd = ", " }
      $sqladd .= "(". $ditag->ditag_id.", ".($ditag->slice->get_seq_region_id).", ". $ditag->start.", ".
                 $ditag->end.", '".$ditag->strand."', ".$ditag->analysis->dbID.", ".$ditag->hit_start.
                 ", ".$ditag->hit_end.", '".$ditag->hit_strand."', '".$ditag->cigar_line."', '".
                 $ditag->ditag_side."', ".$ditag->ditag_pair_id.")";
      $sql .= $sqladd;
      $inserts++;
    }
  }

  #STORE
  if($inserts){
    print STDERR "\nHave $inserts Features.\n";
    eval{
      $db->dbc->do($sql);
    };
    if($@){
      warning("Problem inserting ditag batch!".$@."\n");
    }
  }
  else{
    warn "Nothing stored!";
  }

}


=head2 update_ditagFeature

  Arg [1]    : ditagFeature to update
  Description: update an existing ditagFeature with new values
  Returntype : 1 on success

=cut

sub update_ditagFeature {
  my ($self, $ditagFeature) = @_;

  my $sth = $self->prepare( "UPDATE ditag_feature
                             SET ditag_id=?, seq_region_id=?, seq_region_start=?, seq_region_end=?,
                             seq_region_strand=?, analysis_id=?, hit_start=?, hit_end=?, hit_strand=?,
                             cigar_line=?, ditag_side=?, ditag_pair_id=?
                             where ditag_feature_id=?;" );
  my $result =$sth->execute(
                            $ditagFeature->ditag_id,
                            $ditagFeature->seq_region_id,
                            $ditagFeature->seq_region_start,
                            $ditagFeature->seq_region_end,
                            $ditagFeature->seq_region_strand,
                            $ditagFeature->hit_start,
                            $ditagFeature->hit_end,
                            $ditagFeature->hit_strand,
                            $ditagFeature->cigar_line,
                            $ditagFeature->ditag_side,
                            $ditagFeature->ditag_pair_id,
                            $ditagFeature->dbID,
                           );

  return $result;
}


=head2 list_dbIDs

  Args       : None
  Example    : my @feature_ids = @{$dfa->list_dbIDs()};
  Description: Gets an array of internal IDs for all DitagFeature objects in
               the current database.
  Returntype : List of ints
  Exceptions : None

=cut

sub list_dbIDs {
	my $self = shift;
	
	return $self->_list_dbIDs('ditag_feature');
}

1;
