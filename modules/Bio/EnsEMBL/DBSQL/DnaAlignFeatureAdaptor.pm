

#
# BioPerl module for Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor - Adaptor for DnaAlignFeatures

=head1 SYNOPSIS

    $pfadp = $dbadaptor->get_DnaAlignFeatureAdaptor();

    my @feature_array = $pfadp->fetch_by_contig_id($contig_numeric_id);

    my @feature_array = $pfadp->fetch_by_assembly_location($start,$end,$chr,'UCSC');
 
    $pfadp->store($contig_numeric_id,@feature_array);


=head1 DESCRIPTION


This is an adaptor for DNA features on DNA sequence. Like other
feature getting adaptors it has a number of fetch_ functions and a
store function.


=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::FeatureFactory;
use Bio::EnsEMBL::DnaDnaAlignFeature;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);
# new() can be inherited from Bio::EnsEMBL::Root


=head2 fetch_by_dbID

 Title   : fetch_by_dbID
 Function:
 Returns : 
 Args    :


=cut

sub fetch_by_dbID{
   my ($self,$id) = @_;

   if( !defined $id ) {
       $self->throw("fetch_by_dbID must have an id");
   }

   my $sth = $self->prepare("select d.contig_id,d.contig_start,d.contig_end,d.contig_strand,d.hit_start,d.hit_end,d.hit_strand,d.hit_name,d.cigar_line,d.analysis_id, d.score from dna_align_feature d where d.dna_align_feature_id = $id");
   $sth->execute();

   my ($contig_id,$start,$end,$strand,$hstart,$hend,$hstrand,$hname,$cigar,$analysis_id, $score) = $sth->fetchrow_array();

   if( !defined $contig_id ) {
       $self->throw("No simple feature with id $id");
   }

   my $contig = $self->db->get_RawContigAdaptor->fetch_by_dbID($contig_id);
   my $analysis = $self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id);
   my $out= $self->_new_feature($start,$end,$strand,$score,$hstart,$hend,$hstrand,$hname,$cigar,$analysis,$contig->name,$contig->seq);

   return $out;

}

=head2 fetch_by_contig_id_constraint

 Title   : fetch_by_contig_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_by_contig_id_constraint{
   my ($self,$cid, $constraint) = @_;

   if( !defined $cid ) {
       $self->throw("fetch_by_contig_id must have an contig id");
   }
   
   my $sql = "select d.contig_id,d.contig_start,d.contig_end,d.contig_strand,d.hit_start,d.hit_end,d.hit_strand,d.hit_name,d.cigar_line,d.analysis_id, d.score,d.perc_ident,d.evalue from dna_align_feature d where d.contig_id = $cid";
   if($constraint){
     $sql .= " AND $constraint";
   }
   #print $sql."\n";
   my $sth = $self->prepare($sql);
   $sth->execute();

   my ($contig_id,$start,$end,$strand,$hstart,$hend,$hstrand,$hname,$cigar,$analysis_id, $score,$perc_ident,$evalue);

   $sth->bind_columns(undef,\$contig_id,\$start,\$end,\$strand,\$hstart,\$hend,\$hstrand,\$hname,\$cigar,\$analysis_id, \$score,\$perc_ident,\$evalue);

   my @f;
   my $contig = $self->db->get_RawContigAdaptor->fetch_by_dbID($cid);
   my %ana;

   while( $sth->fetch ) {
       if( !defined $ana{$analysis_id} ) {
	   $ana{$analysis_id} = $self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id);
       }


       my $out= $self->_new_feature($start,$end,$strand,$score,$hstart,$hend,$hstrand,$hname,$cigar,$ana{$analysis_id},$perc_ident,$evalue,$contig->name,$contig->seq);

       push(@f,$out);
   }
   return @f;
}

sub fetch_by_contig_id{
  my ($self, $cid, $logic_name) = @_;
  
  my $analysis;
  my $constraint = undef;
 
  if($logic_name){
    my $aa = $self->db->get_AnalysisAdaptor($logic_name);
    $analysis = $aa->fetch_by_logic_name($logic_name);
    $constraint = " d.analysis_id = ".$analysis->dbID();
  }
  my @features = $self->fetch_by_contig_id_constraint($cid, $constraint);

  return @features;
}

sub fetch_by_contig_id_and_score{
  my($self, $cid, $score, $logic_name) = @_;

  my $analysis;
  my $constraint;
  if(!$cid){
    $self->throw("need a contig id or this won't work\n");
  }  
  if(!$score){
    $self->throw("need a score even if its 0\n");
  }else{
    $constraint = "d.score > $score";
  }
  if($logic_name){
   my $aa = $self->db->get_AnalysisAdaptor();
   $analysis = $aa->fetch_by_logic_name($logic_name);
   $constraint .= " and d.analysis_id = ".$analysis->dbID(); 
  }

  
  my @features = $self->fetch_by_contig_id_constraint($cid, $constraint);
  
  return @features;

}


sub fetch_by_contig_id_and_pid{
  my($self, $cid, $pid, $logic_name) = @_;

  my $analysis;
  my $constraint;
  if(!$cid){
    $self->throw("need a contig id or this won't work\n");
  }  
  if(!$pid){
    $self->throw("need a pid even if its 0\n");
  }else{
    $constraint = "d.perc_ident > $pid";
  }
  if($logic_name){
   my $aa = $self->db->get_AnalysisAdaptor();
   $analysis = $aa->fetch_by_logic_name($logic_name);
   $constraint .= " and d.analysis_id = ".$analysis->dbID(); 
  }

  
  my @features = $self->fetch_by_contig_id_constraint($cid, $constraint);
  
  return @features;

}


sub fetch_by_Slice{
  my($self, $slice, $logic_name) = @_;
  my $constraint;
  my $analysis;
  if(!$slice){
    $self->throw("need a slice to work\n");
  }
  unless ($slice->isa("Bio::EnsEMBL::Slice")) 
        {
            $self->throw("$slice isn't a slice");
        }
  if($logic_name){
   my $aa = $self->db->get_AnalysisAdaptor();
   $analysis = $aa->fetch_by_logic_name($logic_name);
   $constraint .= " d.analysis_id = ".$analysis->dbID(); 
  }
  
  my @features = $self->fetch_by_assembly_location_constraint($slice->chr_start,$slice->chr_end,$slice->chr_name,$slice->assembly_type, $constraint);
  
  my @out_f;
  
  foreach my $f(@features){
   
    my $start = ($f->start - ($slice->chr_start - 1));
    my $end = ($f->end - ($slice->chr_start - 1));
   
    my $out = $self->_new_feature($start,$end,$f->strand,$f->score,$f->hstart,$f->hend,$f->hstrand,$f->hseqname,$f->cigar_string,$f->analysis,$f->percent_id,$f->p_value,$f->seqname,undef);

    push(@out_f, $out);
  }

  return @out_f;
}

sub fetch_by_Slice_and_score {
  my ($self,$slice,$score, $logic_name) = @_;
  my $constraint;
  my $analysis;
  if(!$slice){
    $self->throw("need a slice to work\n");
  }
  unless ($slice->isa("Bio::EnsEMBL::Slice")) 
        {
            $self->throw("$slice isn't a slice");
        }
  if(!$score){
    $self->throw("need a score even if its 0\n");
  }else{
    $constraint .= "d.score > $score";
  }
  if($logic_name){
   my $aa = $self->db->get_AnalysisAdaptor();
   $analysis = $aa->fetch_by_logic_name($logic_name);
   $constraint .= " and d.analysis_id = ".$analysis->dbID(); 
  }
  
  my @features = $self->fetch_by_assembly_location_constraint($slice->chr_start,$slice->chr_end,$slice->chr_name,$slice->assembly_type, $constraint);

  my @out_f;

  foreach my $f(@features){
    my $start = ($f->start - ($slice->chr_start - 1));
    my $end = ($f->end - ($slice->chr_start - 1));

    my $out = $self->_new_feature($start,$end,$f->strand,$f->score,$f->hstart,$f->hend,$f->hstrand,$f->hseqname,$f->cigar_string,$f->analysis,$f->percent_id,$f->p_value,$f->seqname,undef);

    push(@out_f, $out);
  }

  return @out_f;

}  

sub fetch_by_Slice_and_pid {
  my ($self,$slice,$pid, $logic_name) = @_;
  my $constraint;
  my $analysis;
  if(!$slice){
    $self->throw("need a slice to work\n");
  }
  unless ($slice->isa("Bio::EnsEMBL::Slice")) 
        {
            $self->throw("$slice isn't a slice");
        }
  if(!$pid){
    $self->throw("need a pid even if its 0\n");
  }else{
    $constraint .= "d.perc_ident > $pid";
  }
  if($logic_name){
   my $aa = $self->db->get_AnalysisAdaptor();
   $analysis = $aa->fetch_by_logic_name($logic_name);
   $constraint .= " and d.analysis_id = ".$analysis->dbID(); 
  }
  my @features = $self->fetch_by_assembly_location_constraint($slice->chr_start,$slice->chr_end,$slice->chr_name,$slice->assembly_type, $constraint);

  my @out_f;

  foreach my $f(@features){
    my $start = ($f->start - ($slice->chr_start - 1));
    my $end = ($f->end - ($slice->chr_start - 1));

    my $out = $self->_new_feature($start,$end,$f->strand,$f->score,$f->hstart,$f->hend,$f->hstrand,$f->hseqname,$f->cigar_string,$f->analysis,$f->percent_id,$f->p_value,$f->seqname,undef);

    push(@out_f, $out);
  }

  return @out_f;

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
  
  if( !defined $type ) {
    $self->throw("Assembly location must be start,end,chr,type");
  }
  my $constraint = undef;
  my $analysis;
  if($logic_name){
    my $aa = $self->db->get_AnalysisAdaptor();
    $analysis = $aa->fetch_by_logic_name($logic_name);
    $constraint = " d.analysis_id = ".$analysis->dbID();  
  }
  return $self->fetch_by_assembly_location_constraint($start,$end,$chr,$type,$constraint);

}

sub fetch_by_assembly_location_and_score{
  my ($self,$start,$end,$chr,$type, $score, $logic_name) = @_;
  my $constraint;
  if( !defined $type ) {
    $self->throw("Assembly location must be start,end,chr,type");
  }
  if(!$score){
    $self->throw("need a score even if its 0\n");
  }else{
    $constraint .= "d.score > $score";
  }

  my $analysis;
  if($logic_name){
    my $aa = $self->db->get_AnalysisAdaptor();
    $analysis = $aa->fetch_by_logic_name($logic_name);
    $constraint .= " and d.analysis_id = ".$analysis->dbID();  
  }
  return $self->fetch_by_assembly_location_constraint($start,$end,$chr,$type,$constraint);

}


sub fetch_by_assembly_location_and_pid{
  my ($self,$start,$end,$chr,$type, $pid, $logic_name) = @_;
  my $constraint;
  if( !defined $type ) {
    $self->throw("Assembly location must be start,end,chr,type");
  }
  if(!$pid){
    $self->throw("need a pid even if its 0\n");
  }else{
    $constraint .= "d.perc_ident > $pid";
  }

  my $analysis;
  if($logic_name){
    my $aa = $self->db->get_AnalysisAdaptor();
    $analysis = $aa->fetch_by_logic_name($logic_name);
    $constraint .= " and d.analysis_id = ".$analysis->dbID();  
  }
  return $self->fetch_by_assembly_location_constraint($start,$end,$chr,$type,$constraint);

}

=head2 fetch_by_assembly_location_constraint

 Title   : fetch_by_assembly_location_constraint
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_by_assembly_location_constraint{
  my ($self,$chr_start,$chr_end,$chr,$type,$constraint) = @_;
  
  if( !defined $type ) {
    $self->throw("Assembly location must be start,end,chr,type");
  }

  if( $chr_start !~ /^\d/ || $chr_end !~ /^\d/ ) {
    $self->throw("start/end must be numbers not $chr_start,$chr_end (have you typed the location in the right way around - start,end,chromosome,type)?");
  }
  
  my $mapper = $self->db->get_AssemblyMapperAdaptor->fetch_by_type($type);
  
  $mapper->register_region($chr,$chr_start,$chr_end);

  my @cids = $mapper->list_contig_ids($chr, $chr_start ,$chr_end);
  
  # build the SQL
  

 

  if( scalar(@cids) == 0 ) {
    return ();
  }

  my $cid_list = join(',',@cids);

  my $sql = "select d.contig_id,d.contig_start,d.contig_end,d.contig_strand,d.hit_start,d.hit_end,d.hit_strand,d.hit_name,d.cigar_line,d.analysis_id,d.score,d.perc_ident,d.evalue from dna_align_feature d where d.contig_id in ($cid_list)";
  
  if($constraint) {
    $sql .=  " AND $constraint";
  }
  #print STDERR "SQL $sql\n";

  my $sth = $self->prepare($sql);

  $sth->execute();
  
  
  my ($contig_id,$start,$end,$strand,$hstart,$hend,$hstrand,$hname,$cigar,$analysis_id, $score,$perc_ident,$evalue);
  
  $sth->bind_columns(undef,\$contig_id,\$start,\$end,\$strand,\$hstart,\$hend,\$hstrand, \$hname,\$cigar,\$analysis_id, \$score,\$perc_ident,\$evalue);
  

  my @f;
  my %ana;
  my $counter = 0;
  while( $sth->fetch ) {
    # we whether this is sensible to use or not
    
    my @coord_list = $mapper->map_coordinates_to_assembly($contig_id, $start,$end,$strand,"rawcontig");
       
    # coord list > 1 - means does not cleanly map At the moment, skip
    if( scalar(@coord_list) > 1 ) {
      #$self->warn("maps to ".scalar(@coord_list)." coordinate objs not all of feature will be on golden path skipping\n");
      next;
      }
     
    if($coord_list[0]->isa("Bio::EnsEMBL::Mapper::Gap")){
      #$self->warn("this feature is on a part of $contig_id which isn't on the golden path skipping");
      next;
    }
    if(!($coord_list[0]->start >= $chr_start) ||
       !($coord_list[0]->end <= $chr_end)) {
      next;
    }
    if( !defined $ana{$analysis_id} ) {
      $ana{$analysis_id} = $self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id);
    }

    # ok, ready to build a sequence feature: do we want this relative or not?
    
  
    my $out= $self->_new_feature($coord_list[0]->start,$coord_list[0]->end,$coord_list[0]->strand,$score,$hstart,$hend,$hstrand,$hname,$cigar,$ana{$analysis_id},$perc_ident,$evalue,$coord_list[0]->id,undef);

    push(@f,$out);
  }
  
  return @f;

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
   my ($self,$contig_id,@sf) = @_;

   if( scalar(@sf) == 0 ) {
       $self->throw("Must call store with contig_id then sequence features");
   }

   if( $contig_id !~ /^\d+$/ ) {
       $self->throw("Contig_id must be a number, not [$contig_id]");
   }

   my $sth = $self->prepare("insert into dna_align_feature (contig_id,contig_start,contig_end,contig_strand,hit_start,hit_end,hit_strand,hit_name,cigar_line,analysis_id,score,evalue, perc_ident) values (?,?,?,?,?,?,?,?,?,?,?, ?, ?)");

   foreach my $sf ( @sf ) {
       if( !ref $sf || !$sf->isa("Bio::EnsEMBL::DnaDnaAlignFeature") ) {
	   $self->throw("feature must be an Ensembl DnaDnaAlignFeature, not a [$sf]");
       }

       if( !defined $sf->analysis ) {
	   $self->throw("Cannot store sequence features without analysis");
       }
       if( !defined $sf->analysis->dbID ) {
	   # maybe we should throw here. Shouldn't we always have an analysis from the database?
	   $self->throw("I think we should always have an analysis object which has originated from the database. No dbID, not putting in!");
       }
       #print STDERR "storing ".$sf->gffstring."\n";
       $sth->execute($contig_id,$sf->start,$sf->end,$sf->strand,$sf->hstart,$sf->hend,$sf->hstrand,$sf->hseqname,$sf->cigar_string,$sf->analysis->dbID,$sf->score, $sf->p_value, $sf->percent_id);
   }


}

=head2 Internal functions

Internal functions to the adaptor which you never need to call

=cut


sub _new_feature {
  my ($self,$start,$end,$strand,$score,$hstart,$hend,$hstrand,$hseqname,$cigar,$analysis,$perc_ident,$evalue,$seqname,$seq) = @_;

  if( !defined $seqname ) {
    $self->throw("Internal error - wrong number of arguments to new_feature");
  }
 
  my $f1 = Bio::EnsEMBL::SeqFeature->new();
  my $f2 = Bio::EnsEMBL::SeqFeature->new();

  $f1->start($start);
  $f1->end($end);
  $f1->strand($strand);
  $f1->score($score);
  $f1->percent_id($perc_ident);
  $f1->p_value($evalue);
  $f1->seqname($seqname);
  if( defined $seq ) {
    $f1->attach_seq($seq);
  }

  $f2->start($hstart);
  $f2->end($hend);
  $f2->strand($hstrand);
  $f2->percent_id($perc_ident);
  $f2->p_value($evalue);
  $f2->seqname($hseqname);

  $f1->analysis($analysis);
  $f2->analysis($analysis);


  my $out = Bio::EnsEMBL::DnaDnaAlignFeature->new( -cigar_string => $cigar, -feature1 => $f1, -feature2 => $f2);


  return $out;
}
    
1;


