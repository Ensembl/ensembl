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

sub _tablename {
  my $self = shift;

  return 'repeat_feature r, repeat_consensus rc';
}


=head2 _tablename

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
	     r.repeat_id
	     r.repeat_start
	     r.repeat_end
	     r.analysis_id
	     r.score
	     rc.repeat_id
	     rc.repeat_name
	     rc.repeat_class
	     rc.repeat_consensus);
}



=head2 generic_fetch

  Arg [1]    : string $constraint
  Arg [2]    : string $logic_name
  Example    : @repeats = $repeat_feature_adaptor->generic_fetch('','');
  Description: Overrides superclass method to provide an additional 
               table joining constraint before the SQL query is performed. 
  Returntype : list of Bio::EnsEMBL::RepeatFeatures in contig coordinates
  Exceptions : none
  Caller     : internal

=cut

sub generic_fetch {
  my ($self, $constraint, $logic_name) = @_;

  #modify the constraint to join the repeat and repeat_consensus tables
  if($constraint) {
    $constraint .= ' AND r.repeat_id = rc.repeat_id';
  } else {
    $constraint = 'r.repeat_id = rc.repeat_id';
  }

  #invoke the super class method
  return $self->SUPER::generic_fetch($constraint, $logic_name);
}


=head2 _obj_from_hashref

  Arg [1]    : Hashreference $hashref
  Example    : none 
  Description: PROTECTED implementation of abstract superclass method.  
               responsible for the creation of RepeatFeatures from a 
               hashref generated from an SQL query
  Returntype : list of Bio::EnsEMBL::RepeatFeatures in contig coordinates
  Exceptions : none
  Caller     : internal

=cut

sub _obj_from_hashref {
  my ($self, $hashref) = @_;

  my $rca = $self->db()->get_RepeatConsensusAdaptor();

  #create a repeat consensus object
  my $rc = new Bio::EnsEMBL::RepeatConsensus;
  $rc->dbID($hashref->{'repeat_id'});
  $rc->repeat_class($hashref->{'repeat_class'});
  $rc->name($hashref->{'repeat_name'});
  $rc->repeat_consensus($hashref->{'repeat_consensus'});
  $rc->adaptor($rca);

  #get the analysis object for this repeat
  my $aa = $self->db->get_AnalysisAdaptor();
  my $analysis = $aa->fetch_by_dbID($hashref->{'analysis_id'});

  #create the new repeat feature
  my $r = new Bio::EnsEMBL::RepeatFeature;
  $r->dbID($hashref->{'repeat_feature_id'});
  $r->repeat_id($hashref->{'repeat_id'});
  $r->contig_id($hashref->{'contig_id'});

  $r->start($hashref->{'contig_start'});
  $r->end($hashref->{'contig_end'});

  $r->score($hashref->{'score'});
  $r->strand( $hashref->{'contig_strand'} );  
  $r->hstart( $hashref->{'repeat_start'} );
  $r->hend( $hashref->{'repeat_end'} );

  $r->analysis($analysis);
  $r->repeat_consensus($rc);
  $r->adaptor($self);

  #attach the appropriate contig to this sequence
  my $ca = $self->db()->get_RawContigAdaptor();
  my $contig = $ca->fetch_by_dbID($hashref->{'contig_id'});
  $r->attach_seq($contig);
 
  return $r;
}
  


=head2 store

  Arg [1]    : int $contig_id
               the database id of the contig to store this repeat_feature on
  Arg [2]    : list of Bio::EnsEMBL::RepeatFeatures $repeat_feature_id
               the list of repeat features to store in the database
  Example    : $repeat_feature_adaptor->store(1234, @repeat_features);
  Description: stores a repeat feature in the database
  Returntype : none
  Exceptions : thrown if $contig_id is not an int, or if repeat_consensus is
               not used
  Caller     : genreal

=cut

sub store {
  my( $self, $contig_id, @repeats ) = @_;
  
  my $rca = $self->db->get_RepeatConsensusAdaptor;
  my ($cons, $db_id);

  $self->throw("Can't store repeats without a contig_id (got '$contig_id')")
    unless $contig_id =~ /^\d+$/;

  my $sth = $self->prepare(qq{
    INSERT into repeat_feature( repeat_feature_id
				, contig_id
				, contig_start
				, contig_end
				, contig_strand
				, repeat_id
				, repeat_start
				, repeat_end
				, score
				, analysis_id )
      VALUES(NULL, ?,?,?,?,?,?,?,?,?)
    });
  foreach my $rf (@repeats) {
            
        unless ($rf->repeat_id){

   $self->throw("Must have a RepeatConsensus attached")
	unless defined ($cons = $rf->repeat_consensus);
      
      # for tandem repeats - simply store consensus and repeat
      # one pair per hit. don't need to check consensi stored
      # already. consensus has name and class set to 'trf'
 }
      if ($cons->repeat_class eq 'trf') {

	$rca->store($cons);
	$rf->repeat_id($cons->dbID);

      } elsif ($cons->repeat_class eq 'Simple_repeat') {

        my $rcon = $cons->name;
        $rcon =~ s/\((\S+)\)n/$1/;   # get repeat element
	$cons->repeat_consensus($rcon);
	$rca->store($cons);
	$rf->repeat_id($cons->dbID);

      } else {
	
	# for other repeats - need to see if a consensus is stored already
	
	unless ($rf->repeat_id) {

	  # need to get the consensus seq object for this repeat

	  unless ($cons->dbID) {
	    my @match = ($rca->fetch_by_name($cons->name));
	    
	    if (@match > 1) {
	      $self->warn(@match . " consensi for " . $cons->name . "\n");
	    } elsif (@match == 0) {
	      # if we don't match a consensus already stored
	      # create a fake one
	      # set consensus to 'N' as null seq not allowed
	      # FIXME: not happy with this, but ho hum ...
	      $self->warn("Can't find " . $cons->name . "\n");
	      $cons->repeat_consensus("N");
	      $rca->store($cons);
	    }
	    $db_id = ($rca->fetch_by_name($cons->name))[0]->dbID;
	    
	    $cons->dbID($db_id);
	  }
	  $rf->repeat_id($cons->dbID);
	}
      }
	#print STDERR "repeat = ".$rf->analysis." ".$rf->analysis->dbID."\n";
      $sth->execute(
		    $contig_id,
		    $rf->start,
		    $rf->end,
		    $rf->strand,
		    $rf->repeat_id,
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




1;





