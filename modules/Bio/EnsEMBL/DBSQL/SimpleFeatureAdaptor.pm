#
# EnsEMBL module for Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright EMBL/EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor 

=head1 SYNOPSIS

my $simple_feature_adaptor = $database_adaptor->get_SimpleFeatureAdaptor();
@simple_features = $simple_feature_adaptor->fetch_by_Slice($slice);

=head1 DESCRIPTION

Simple Feature Adaptor - database access for simple features 

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

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
  Example    : $simple_feature_adaptor->store(1234, @simple_feats);
  Description: Stores a list of simple feature objects in the database
  Returntype : none
  Exceptions : thrown if @sf is not defined, if any of the features do not
               have an attached slice.
               or if any elements of @sf are not Bio::EnsEMBL::SeqFeatures 
  Caller     : general

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
  my $slice_adaptor = $db->get_SliceAdaptor();
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

    my $slice = $sf->slice();
    if(!ref($slice) || !$slice->isa("Bio::EnsEMBL::Slice")) {
      throw("A slice must be attached to the features to be stored.");
    }

    # make sure that the feature coordinates are relative to
    # the start of the seq_region that the prediction transcript is on
    if($slice->start != 1 || $slice->strand != 1) {
      #move the feature onto a slice of the entire seq_region
      $slice = $slice_adaptor->fetch_by_region($slice->coord_system->name(),
                                               $slice->seq_region_name(),
                                               undef, #start
                                               undef, #end
                                               undef, #strand
                                              $slice->coord_system->version());

      $sf = $sf->transfer($slice);

      if(!$sf) {
        throw('Could not transfer SimpleFeature to slice of ' .
              'entire seq_region prior to storing');
      }
    }

    my $seq_region_id = $slice_adaptor->get_seq_region_id($slice);

    $sth->execute($seq_region_id, $sf->start, $sf->end, $sf->strand,
                  $sf->display_label, $sf->analysis->dbID, $sf->score);

    $sf->dbID($sth->{'mysql_insertid'});
    $sf->adaptor($self);
  }
}


=head2 _tablename

  Arg [1]    : none
  Example    : none
  Description: PROTECTED implementation of superclass abstract method
               returns the names, aliases of the tables to use for queries
  Returntype : list of listrefs of strings
  Exceptions : none
  Caller     : internal

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

=cut

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my $db = $self->db();
  my $aa = $db->get_AnalysisAdaptor();
  my $slice_adaptor = $db->get_SliceAdaptor();

  my @features;
  my %slice_cache;
  my %analysis_cache;

  my($simple_feature_id,$seq_region_id, $seq_region_start, $seq_region_end,
     $seq_region_strand, $display_label, $analysis_id, $score);

  $sth->bind_columns(\$simple_feature_id,\$seq_region_id, \$seq_region_start,
                     \$seq_region_end, \$seq_region_strand, \$display_label,
                     \$analysis_id, \$score);

  while($sth->fetch()) {
    my $slice = $slice_cache{$seq_region_id} ||=
      $slice_adaptor->fetch_by_seq_region_id($seq_region_id);

    my $analysis = $analysis_cache{$analysis_id} ||=
      $aa->fetch_by_dbID($analysis_id);

    push @features, Bio::EnsEMBL::SimpleFeature->new
      (-START => $seq_region_start,
       -END   => $seq_region_end,
       -STRAND => $seq_region_strand,
       -SLICE => $slice,
       -ANALYSIS => $analysis,
       -ADAPTOR => $self,
       -DBID => $simple_feature_id,
       -DISPLAY_LABEL => $display_label,
       -SCORE => $score);
  }

  return \@features;
}


=head2 list_dbIDs

  Arg [1]    : none
  Example    : @feature_ids = @{$simple_feature_adaptor->list_dbIDs()};
  Description: Gets an array of internal ids for all simple features in the current db
  Returntype : list of ints
  Exceptions : none
  Caller     : ?

=cut

sub list_dbIDs {
   my ($self) = @_;

   return $self->_list_dbIDs("simple_feature");
}

1;
