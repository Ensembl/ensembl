#
# BioPerl module for Bio::EnsEMBL::DBSQL::BaseAlignFeatureAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor - Abstract Base class for 
                                          FeatureAdaptors

=head1 SYNOPSIS

Abstract class should not be instantiated.  Implementation of
abstract methods must be performed by subclasses.

=head1 DESCRIPTION

This is a base adaptor for feature adaptors. This base class is simply a way
of eliminating code duplication through the implementation of methods 
common to all feature adaptors.

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root

use Bio::EnsEMBL::DBSQL::BaseAdaptor;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

  
sub generic_fetch {
  my ($self, $constraint, $logic_name) = @_;
  
  my $tablename = $self->_tablename();
  my $columns = join(', ', $self->_columns());
  
  if($logic_name) {
    #determine the analysis id via the logic_name
    my $aa = $self->db->get_AnalysisAdaptor();
    my $analysis = $aa->fetch_by_logic_name($logic_name);
    unless(defined $analysis && $analysis->dbID() ) {
      $self->warn("No analysis for logic name $logic_name exists\n");
      return ();
    }
    
    my $analysis_id = $analysis->dbID();
    
    if($constraint) {
      $constraint .= " AND analysis_id = $analysis_id";
    } else {
      $constraint = " analysis_id = $analysis_id";
    }
  } 
      
  my $sql = 
    "SELECT $columns 
     FROM $tablename";

  if($constraint) {
     $sql .= " WHERE $constraint";
  }
  
  my $sth = $self->prepare($sql);
  $sth->execute();
  
  my $hashref;
  my @out;

  while($hashref = $sth->fetchrow_hashref()) {
    push @out, $self->_obj_from_hashref($hashref);
  }
  
  return @out;
}


=head2 fetch_by_dbID
  Args      : none
  Function  : Retrieves AlignFeature from database
  Returntype: BaseAlignFeature
  Exceptions: thrown: if _columns() or _table() not implemented by subclass
                      if $id arg is not defined
  Caller    : Slice

=cut

sub fetch_by_dbID{
  my ($self,$id) = @_;
  
  unless(defined $id) {
    $self->throw("fetch_by_dbID must have an id");
  }

  my $tablename = $self->_tablename();
  my $constraint = "${tablename}_id = $id";

  return $self->generic_fetch($constraint);
}



=head2 fetch_by_contig_id_constraint

 Title   : fetch_by_contig_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_by_contig_id_constraint {
  my ($self, $cid, $constraint, $logic_name) = @_;
  
  if( !defined $cid ) {
    $self->throw("fetch_by_contig_id_constraint must have an contig id");
  }

  if($constraint) {
    $constraint .= " AND contig_id = $cid";
  } else {
    $constraint = "contig_id = $cid";
  }

  return $self->generic_fetch($constraint, $logic_name);
}

   
sub fetch_by_contig_id{
  my ($self, $cid, $logic_name) = @_;

  #fetch by contig id constraint with empty constraint
  return $self->fetch_by_contig_id_constraint($cid, '',$logic_name);
}


sub fetch_by_contig_id_and_score{
  my($self, $cid, $score, $logic_name) = @_;

  my $constraint;

  if(!defined $score){
    $self->throw("need a score even if its 0\n");
  } else{
    $constraint = "score > $score";
  }
    
  my @features = 
    $self->fetch_by_contig_id_constraint($cid, $constraint, $logic_name);
  
  return @features;
}


sub fetch_by_Slice_constraint {
  my($self, $slice, $constraint, $logic_name) = @_;

  if(!$slice){
    $self->throw("need a slice to work\n");
  }
  unless($slice->isa("Bio::EnsEMBL::Slice")) {
    $self->throw("$slice isn't a slice");
  }

  my @features = 
    $self->fetch_by_assembly_location_constraint($slice->chr_start,
						 $slice->chr_end,
						 $slice->chr_name,
						 $slice->assembly_type,
						 $constraint,
						 $logic_name);

  #convert from chromosomal coordinates to slice coordinates
  foreach my $f (@features){
    my $start = ($f->start - ($slice->chr_start - 1));
    my $end = ($f->end - ($slice->chr_start - 1));
   
    $f->start($start);
    $f->end($end);
    $f->attach_seq($slice);
  }

  return @features;
}


sub fetch_by_Slice {
  my ($self, $slice, $logic_name) = @_;
  
  #fetch by constraint with empty constraint
  return $self->fetch_by_Slice_constraint($slice, '', $logic_name);
}


sub fetch_by_Slice_and_score {
  my ($self, $slice, $score, $logic_name) = @_;
  my $constraint;

  if(!defined $score) {
    $self->throw("need a score even if its 0\n");
  } else {
    $constraint = "score > $score";
  }

  return $self->fetch_by_Slice_constraint($slice, $constraint, $logic_name);
}  


=head2 fetch_by_assembly_location

 Title   : fetch_by_assembly_location
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_by_assembly_location{
  my ($self,$start,$end,$chr,$type, $logic_name) = @_;

  #fetch by assembly location constraint w/ empty constraint
  return $self->fetch_by_assembly_location_constraint($start, $end, $chr,
						      $type, '', $logic_name);
}


sub fetch_by_assembly_location_and_score{
  my ($self, $start, $end, $chr, $type, $score, $logic_name) = @_;
  my $constraint;

  if(!defined $score) {
    $self->throw("need a score even if its 0\n");
  } else {
    $constraint = "score > $score";
  }

  return $self->fetch_by_assembly_location_constraint($start, $end, $chr,$type,
						     $constraint, $logic_name);
}


=head2 fetch_by_assembly_location_constraint

 Title   : fetch_by_assembly_location_constraint
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_by_assembly_location_constraint {
  my ($self, $chr_start, $chr_end, $chr, $type, $constraint, $logic_name) = @_;

  if( !defined $type ) {
    $self->throw("Assembly location must be start,end,chr,type");
  }

  if( $chr_start !~ /^\d/ || $chr_end !~ /^\d/ ) {
    $self->throw("start/end must be numbers not $chr_start,$chr_end " .
		 "(have you typed the location in the right way around" .
		 "- start,end,chromosome,type)?");
  }
  
  my $mapper = $self->db->get_AssemblyMapperAdaptor->fetch_by_type($type);  
  my @cids = $mapper->list_contig_ids($chr, $chr_start ,$chr_end);
  
  if( scalar(@cids) == 0 ) {
    return ();
  }

  my $cid_list = join(',',@cids);

  #construct the constraint for the contig ids 
  if($constraint) {
    $constraint .= " AND contig_id IN ($cid_list)";
  } else {
    $constraint = "contig_id IN ($cid_list)";
  }
  
  my @features = $self->generic_fetch($constraint, $logic_name); 

  my @out;

  #convert the features to assembly coordinates from raw contig coordinates
  foreach my $f (@features) {
    #since feats were obtained in contig coords, attached seq is a contig
    my $contig_id = $f->entire_seq->dbID();
    my @coord_list = 
      $mapper->map_coordinates_to_assembly($contig_id, $f->start(),
					   $f->end(),$f->strand(),"rawcontig");
       
    # coord list > 1 - means does not cleanly map At the moment, skip
    if( scalar(@coord_list) > 1 ) {
      next;
    }
     
    my ($coord) = @coord_list;

    #maps to a gap in assembly?
    if($coord->isa("Bio::EnsEMBL::Mapper::Gap")){
      next;
    }

    #maps to region outside desired area
    if(!($coord->start() >= $chr_start) ||
       !($coord->end() <= $chr_end)) {
      next;
    }
        
    $f->start($coord->start());
    $f->end($coord->end());
    $f->strand($coord->strand());
    $f->seqname($coord->id());

    push(@out,$f);
  }
  
  return @out;
}

=head2 store

 Title   : store
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub store{
  my $self = @_;

  $self->throw("Abstract method store not defined by implementing subclass\n");
}


=head2 _tablename

  Args      : none
  Function  : abstract protected method implemented by subclass to provide 
              name of table for sql queries
  Returntype: string
  Exceptions: thrown if not implemented by subclass
  Caller    : Implementing AlignFeatureAdaptor subclasses

=cut

sub _tablename {
  my $self = shift;

  $self->throw("abstract method _tablename not defined by implementing" .
               " subclass of AlignFeatureAdaptor");
  return undef;
}

=head2 _columns

  Args      : none
  Function  : abstract protected method implemented by subclass to provide 
              column names for sql queries
  Returntype: string list
  Exceptions: thrown if not implemented by subclass
  Caller    : Implementing AlignFeatureAdaptor subclasses

=cut

sub _columns {
  my $self = shift;

  $self->throw("abstract method _columns not defined by implementing" .
               " subclass of AlignFeatureAdaptor");
}


=head2 _obj_from_hashref

  Args      : a  DBI hashref
  Function  : abstract protected method implemented by subclass to provide 
              object creation from a sql DBI hashref
  Returntype: BaseAlignFeature
  Exceptions: thrown if not implemented by subclass
  Caller    : Implementing AlignFeatureAdaptor subclasses

=cut

sub _obj_from_hashref {
  my $self = shift;

  $self->throw("abstract method _obj_from_hashref not defined by implementing"
             . " subclass of AlignFeatureAdaptor");
} 

1;


