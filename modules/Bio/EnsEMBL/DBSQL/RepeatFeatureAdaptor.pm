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

$repeat_feature_adaptor = $database_adaptor->get_RepeatFeatureAdaptor();

=head1 DESCRIPTION

This is an adaptor for the retrieval and storage of RepeatFeature objects
from the RepeatFeature database.  Most of the implementation is in the
superclass BaseFeatureAdaptor. 

=head1 AUTHOR - James Gilbert

Email jgrg@ebi.ac.uk

=head1 CONTACT

Arne Stabenau - stabenau@ebi.ac.uk
Graham McVicker - mcvicker@ebi.ac.uk
Ewan Birney - birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor;

use strict;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::RepeatFeature;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);


=head2 _tablename

  Arg [1]    : none
  Example    : none
  Description: PROTECTED Implementation of abstract superclass method to 
               provide the name of the tables to query 
  Returntype : string
  Exceptions : none
  Caller     : internal

=cut

sub _tables {
  my $self = shift;

  return (['repeat_feature', 'r'], ['repeat_consensus', 'rc']);
}


=head2 _columns

  Arg [1]    : none
  Example    : none
  Description: PROTECTED Implementation of abstract superclass method to 
               provide the name of the columns to query 
  Returntype : list of strings
  Exceptions : none
  Caller     : internal

=cut

sub _columns {
  my $self = shift;
  
  return qw (r.repeat_feature_id
	     r.contig_id
	     r.contig_start
	     r.contig_end
	     r.contig_strand
	     r.repeat_consensus_id
	     r.repeat_start
	     r.repeat_end
	     r.analysis_id
	     r.score
	     rc.repeat_name
	     rc.repeat_class
	     rc.repeat_consensus);
}



=head2 _default_where_clause

  Arg [1]    : none
  Example    : none
  Description: Overrides superclass method to provide an additional 
               table joining constraint before the SQL query is performed. 
  Returntype : string
  Exceptions : none
  Caller     : generic_fetch

=cut

sub _default_where_clause {
  my $self = shift;

  return 'r.repeat_consensus_id = rc.repeat_consensus_id';
}


=head2 _obj_from_hashref

  Arg [1]    : Hashreference $hashref
  Example    : none 
  Description: PROTECTED implementation of abstract superclass method.  
               responsible for the creation of RepeatFeatures from a 
               hashref generated from an SQL query
  Returntype : listref of Bio::EnsEMBL::RepeatFeatures in contig coordinates
  Exceptions : none
  Caller     : internal

=cut

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my $rca = $self->db()->get_RepeatConsensusAdaptor();
  my $ca = $self->db()->get_RawContigAdaptor();
  my $aa = $self->db->get_AnalysisAdaptor();

  my @features;
  my %rc_hash;
  my %analysis_hash;
  my %contig_hash;

  my($repeat_feature_id, $contig_id, $contig_start, $contig_end, 
     $contig_strand, $repeat_consensus_id, $repeat_start, $repeat_end,
     $analysis_id, $score, $repeat_name, $repeat_class,
     $repeat_consensus);
  
  $sth->bind_columns( \$repeat_feature_id, \$contig_id, \$contig_start, 
                      \$contig_end, \$contig_strand, \$repeat_consensus_id, 
                      \$repeat_start,\$repeat_end, \$analysis_id, \$score, 
                      \$repeat_name, \$repeat_class,
                      \$repeat_consensus );

  my $rc;
  my $contig;
  my $analysis;

  while($sth->fetch()) {
    #create a repeat consensus object
    unless($rc = $rc_hash{$repeat_consensus_id}) {
      $rc = new Bio::EnsEMBL::RepeatConsensus;
      $rc->dbID($repeat_consensus_id);
      $rc->repeat_class($repeat_class);
      $rc->name($repeat_name);
      $rc->repeat_consensus($repeat_consensus);
      $rc->adaptor($rca);

      $rc_hash{$repeat_consensus_id} = $rc;
    }
    
    unless($analysis = $analysis_hash{$analysis_id}) {
      $analysis = $aa->fetch_by_dbID($analysis_id);
      $analysis_hash{$analysis_id} = $analysis;
    }

    unless($contig = $contig_hash{$contig_id}) {
      $contig = $ca->fetch_by_dbID($contig_id);
      $contig_hash{$contig_id} = $contig;
    }

    #create the new repeat feature
    push @features, Bio::EnsEMBL::RepeatFeature->new_fast(
			    { '_gsf_tag_hash'  =>  {},
			      '_gsf_sub_array' =>  [],
                              '_parse_h'       =>  {},
                              '_analysis'      =>  $analysis,
                              '_gsf_start'         =>  $contig_start,
                              '_gsf_end'           =>  $contig_end,
                              '_gsf_strand'        =>  $contig_strand,
                              '_gsf_score'         =>  $score,
                              '_hstart'        =>  $repeat_start,
                              '_hend'          =>  $repeat_end,
                              '_repeat_consensus' => $rc,
			      '_adaptor'       =>  $self,
			      '_gsf_seq'       =>  $contig,
                              '_db_id'         =>  $repeat_feature_id } );
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
  
  my $rca = $self->db->get_RepeatConsensusAdaptor;
  my ($cons, $db_id);

  my $sth = $self->prepare(qq{
    INSERT into repeat_feature( repeat_feature_id
				, contig_id
				, contig_start
				, contig_end
				, contig_strand
				, repeat_consensus_id
				, repeat_start
				, repeat_end
				, score
				, analysis_id )
      VALUES(NULL, ?,?,?,?,?,?,?,?,?)
    });

  foreach my $rf (@repeats) {
    $self->throw("Must have a RepeatConsensus attached")
      unless defined ($cons = $rf->repeat_consensus);
      
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
      unless ($cons->dbID) {
	my $match = ($rca->fetch_by_name($cons->name));
	
	if($match) {
	  #set the consensus dbID to be the same as the database one
	  $cons->dbID($match->dbID());
	} else {
	  # if we don't match a consensus already stored create a fake one 
	  # and set consensus to 'N' as null seq not allowed
	  # FIXME: not happy with this, but ho hum ...
	  $self->warn("Can't find " . $cons->name . "\n");
	  $cons->repeat_consensus("N");
	  $rca->store($cons);
	}
	
	#if (@match > 1) {
	  #multiple consensi were matched
	#  $self->warn(@match . " consensi for " . $cons->name . "\n");
	#}
      }
    }
    
    my $contig = $rf->entire_seq();

    unless(defined $contig && $contig->isa("Bio::EnsEMBL::RawContig")) {
      $self->throw("RepeatFeature cannot be stored without a contig " .
		   "attached via the attach_seq method");
    } unless($contig->dbID()) {
      $self->throw("RepeatFeature cannot be stored because attached contig " .
		   "does not have a dbID");
    }
    
    $sth->execute(
		  $contig->dbID(),
		  $rf->start,
		  $rf->end,
		  $rf->strand,
		  $rf->repeat_consensus->dbID(),
		  $rf->hstart,
		  $rf->hend,
		  $rf->score,
		  $rf->analysis->dbID,
		 );

    my $db_id = $sth->{'mysql_insertid'}
    or $self->throw("Didn't get an insertid from the INSERT statement");
    $rf->dbID($db_id);
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





