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


# Let the code begin...


package Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::SimpleFeature;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);


=head2 store

  Arg [1]    : list of Bio::EnsEMBL::SimpleFeatures @sf
               the simple features to store in the database
  Example    : $simple_feature_adaptor->store(1234, @simple_feats);
  Description: Stores a list of simple feature objects in the database
  Returntype : none
  Exceptions : thrown if @sf is not defined, if any of the features do not
               have an attached contig object, 
               or if any elements of @sf are not Bio::EnsEMBL::SeqFeatures 
  Caller     : general

=cut

sub store{
  my ($self,@sf) = @_;
  
  if( scalar(@sf) == 0 ) {
    $self->throw("Must call store with list of sequence features");
  }
  
  my $sth = 
    $self->prepare("INSERT INTO simple_feature (contig_id, contig_start,
                                                contig_end, contig_strand,
                                                display_label, analysis_id,
                                                score) 
                    VALUES (?,?,?,?,?,?,?)");

  foreach my $sf ( @sf ) {
    if( !ref $sf || !$sf->isa("Bio::EnsEMBL::SimpleFeature") ) {
      $self->throw("Simple feature must be an Ensembl SimpleFeature, " .
		   "not a [$sf]");
    }
    
    if( !defined $sf->analysis ) {
      $self->throw("Cannot store sequence features without analysis");
    }
    if( !defined $sf->analysis->dbID ) {
      $self->throw("I think we should always have an analysis object " .
		   "which has originated from the database. No dbID, " .
		   "not putting in!");
    }
    
    my $contig = $sf->entire_seq();
    unless(defined $contig && $contig->isa("Bio::EnsEMBL::RawContig")) {
      $self->throw("Cannot store feature without a Contig object attached via "
		   . "attach_seq\n");
    }

    # store analysis if not there already
    $self->db->get_AnalysisAdaptor->store($sf->analysis);

    $sth->execute($contig->dbID(), $sf->start, $sf->end, $sf->strand,
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
	     sf.contig_id sf.contig_start sf.contig_end sf.contig_strand
	     sf.display_label sf.analysis_id score );
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

  my $aa = $self->db()->get_AnalysisAdaptor();  
  my $rca = $self->db()->get_RawContigAdaptor();

  my @features = ();
  
  my $hashref;
  while($hashref = $sth->fetchrow_hashref()) {
    my $contig = $rca->fetch_by_dbID($hashref->{'contig_id'});
    my $analysis = $aa->fetch_by_dbID($hashref->{'analysis_id'});

    my $out = Bio::EnsEMBL::SimpleFeature->new();
    $out->start($hashref->{'contig_start'});
    $out->end($hashref->{'contig_end'});
    $out->strand($hashref->{'contig_strand'});
    $out->analysis($analysis);
    $out->display_label($hashref->{'display_label'});
    $out->attach_seq($contig); 

    if($hashref->{'score'}) {
      $out->score($hashref->{'score'});
    }
    
    $out->dbID($hashref->{'simple_feature_id'});

    push @features, $out;
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
