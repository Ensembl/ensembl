
### Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor

package Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor;

use strict;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::RepeatFeature;
use Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptorI;
use Bio::LocationI;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor
	  Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptorI);

use implements qw(Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptorI);

sub fetch_by_RawContig {
    my( $self, $contig, $logic_name ) = @_;

    my @repeats = $self->fetch_by_contig_id($contig->dbID, $logic_name);
    foreach my $r (@repeats) {
        $r->attach_seq($contig);
    }
    return @repeats;
}


sub fetch_by_contig_id {
    my( $self, $contig_id, $logic_name ) = @_;
    my $constraint = "contig_id = $contig_id";

    if($logic_name){
      my $analysis = 
	$self->db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
      if($analysis->dbID()) {
	$self->warn("No analysis for logic name $logic_name\n");
	return ();
      }

      $constraint .= " AND analysis_id = ".$analysis->dbID;
    }

    return $self->_generic_fetch($constraint);
}


sub fetch_by_dbID {
    my( $self, $db_id ) = @_;

    my ($rf) = $self->_generic_fetch(
        qq{ repeat_feature_id = $db_id }
        );
    return $rf;
}



sub _generic_fetch {
    my( $self, $where_clause ) = @_;

    my( $repeat_feature_id,
        $contig_id,
        $contig_start,
        $contig_end,
        $contig_strand,
        $repeat_id,
        $repeat_start,
        $repeat_end,
        $analysis_id,
	$score,
	$repeat_name,
	$repeat_class,
	$repeat_consensus
        );

    my $sth = $self->prepare(qq{
      SELECT f.repeat_feature_id
	, f.contig_id
	, f.contig_start
        , f.contig_end
        , f.contig_strand
 	, f.repeat_id
        , f.repeat_start
        , f.repeat_end
        , f.analysis_id
	, f.score
	, c.repeat_name
	, c.repeat_class
	, c.repeat_consensus
        FROM repeat_feature f, repeat_consensus c
        WHERE }. $where_clause . " AND f.repeat_id = c.repeat_id");

    $sth->execute;
    $sth->bind_columns(
        \$repeat_feature_id,
        \$contig_id,
        \$contig_start,
        \$contig_end,
        \$contig_strand,
        \$repeat_id,
        \$repeat_start,
        \$repeat_end,
        \$analysis_id,
	\$score,
	\$repeat_name,
        \$repeat_class,
	\$repeat_consensus
        );

    my $rca = $self->db->get_RepeatConsensusAdaptor();
    my $aa  = $self->db->get_AnalysisAdaptor();
    my $rfa = $self->db->get_RepeatFeatureAdaptor();

    my @repeats;
    while ($sth->fetch) {        
      my $r = $self->_new_repeat($contig_start, $contig_end, $contig_strand, 
				 $repeat_start, $repeat_end, $score, 
				 $contig_id, $repeat_id, $repeat_name,
				 $repeat_class, $repeat_consensus,
				 $analysis_id, $rfa, $rca, 
				 $aa, $repeat_feature_id);
      push(@repeats, $r);
    }
    return( @repeats );
}




sub fetch_by_assembly_location{
  my ($self,$start,$end,$chr,$type, $logic_name) = @_;
  
  my $constraint;

  if($logic_name){
    my $analysis  = 
      $self->db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
    unless($analysis->dbID()) {
      $self->warn("No analysis for logic name $logic_name\n");
      return ();
    }
    $constraint = " analysis_id = ".$analysis->dbID;
  }
  
  my @repeats = 
    $self->fetch_by_assembly_location_constraint($start, $end, $chr, 
						 $type, $constraint);

  return @repeats;
}



sub fetch_by_Slice{
  my ($self, $slice, $logic_name) = @_;

  my $constraint;

  if($logic_name){
      my $analysis  = 
	$self->db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name);
      unless($analysis->dbID()) {
	$self->warn("No analysis for logic name $logic_name\n");
	return ();
      }

      $constraint = " analysis_id = ".$analysis->dbID;
    }
  

  my @repeats = 
    $self->fetch_by_assembly_location_constraint($slice->chr_start, 
						 $slice->chr_end, 
						 $slice->chr_name, 
						 $slice->assembly_type, 
						 $constraint);
 
  foreach my $r(@repeats){
    #convert from chromosomal to slice coordiantes
    my $start = ($r->start - ($slice->chr_start - 1));
    my $end = ($r->end - ($slice->chr_start - 1));
    $r->start($start);
    $r->end($end);
  }
  return(@repeats);
}



sub fetch_by_assembly_location_constraint{
  my ($self,$chr_start,$chr_end,$chr,$type,$constraint) = @_;
  
  if( !defined $type ) {
    $self->throw("Assembly location must be start,end,chr,type");
  }
  
  if( $chr_start !~ /^\d/ || $chr_end !~ /^\d/ ) {
    $self->throw("start/end must be numbers not $chr_start,$chr_end " .
		 "(have you typed the location in the right way around " .
		 "- start,end,chromosome,type)?");
  }
  
  my $mapper = $self->db->get_AssemblyMapperAdaptor->fetch_by_type($type);
  
  $mapper->register_region($chr,$chr_start,$chr_end);

  # determine the raw contigs that these repeats are on
  my @cids = $mapper->list_contig_ids($chr, $chr_start ,$chr_end);
  my $cid_list = join(',',@cids);
  my $sql = "contig_id in($cid_list) ";
  if($constraint){
    $sql .= "AND $constraint";
  }

  #fetch the repeats from the database
  my @repeats = $self->_generic_fetch(qq{$sql});
  my @result;

  #
  # Adjust each repeat so start and end are in chromosomal coordinates
  # 
  foreach my $r(@repeats){
    my @coord_list = 
      $mapper->map_coordinates_to_assembly($r->contig_id, $r->start, 
					   $r->end, $r->strand, "rawcontig");

    if(scalar(@coord_list) > 1){
      #$self->warn("this feature doesn't cleanly map skipping\n");
      next;
    }
    if($coord_list[0]->isa("Bio::EnsEMBL::Mapper::Gap")){
      #$self->warn("this feature is on a part of ".$r->contig_id.
      #            " which isn't on the golden path skipping");
      next;
    }
    if(!($coord_list[0]->start >= $chr_start) ||
       !($coord_list[0]->end <= $chr_end)){
      next;
    }

    $r->start($coord_list[0]->start());
    $r->end($coord_list[0]->end());
    push( @result, $r );
  }

  return @result;
}




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



sub _new_repeat{
  my($self, $start, $end, $strand, $hstart, $hend, $score, $contig_id, 
     $repeat_id, $repeat_name, $repeat_class, $repeat_consensus, 
     $analysis_id, $rfa, $rca, $aa, $dbID) = @_;

  #create a repeat consensus object
  my $rc = new Bio::EnsEMBL::RepeatConsensus;
  $rc->dbID($repeat_id);
  $rc->repeat_class($repeat_class);
  $rc->name($repeat_name);
  $rc->repeat_consensus($repeat_consensus);
  $rc->adaptor($rca);

  #get the analysis object for this repeat
  my $analysis = $aa->fetch_by_dbID($analysis_id);

  #create the new repeat feature
  my $r = new Bio::EnsEMBL::RepeatFeature;
  $r->dbID($dbID);
  $r->repeat_consensus($rc);
  $r->repeat_id($repeat_id);
  $r->contig_id( $contig_id);
  $r->start    ( $start  );
  $r->end      ( $end    );
  if($strand == -0){
    $strand = 0;
  }
  $r->score($score);
  $r->strand   ( $strand );  
  $r->hstart   ( $hstart  );
  $r->hend     ( $hend    );
  $r->analysis($analysis);

  return $r;
}

1;





__END__

=head1 NAME - Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

