

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

# Object preamble - inherits from Bio::Root::RootI

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::FeatureFactory;
use Bio::EnsEMBL::DnaDnaAlignFeature;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);
# new() can be inherited from Bio::Root::RootI


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

   my $sth = $self->prepare("select p.contig_id,p.contig_start,p.contig_end,p.contig_strand,p.hit_start,p.hit_end,p.hit_strand,p.hit_name,p.cigar_line,p.analysis_id from dna_align_feature p where p.dna_align_feature_id = $id");
   $sth->execute();

   my ($contig_id,$start,$end,$strand,$hstart,$hend,$hstrand,$hname,$cigar,$analysis_id) = $sth->fetchrow_array();

   if( !defined $contig_id ) {
       $self->throw("No simple feature with id $id");
   }

   my $contig = $self->db->get_RawContigAdaptor->fetch_by_dbID($contig_id);

   my $out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();;
   $out->start($start);
   $out->end($end);
   $out->strand($strand);

   $out->hstart($hstart);
   $out->hend($hend);
   $out->hstrand($hstrand);
   $out->hseqname($hname);
   $out->cigar($cigar);
   
   $out->analysis($self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id));

   $out->seqname($contig->id);
   $out->attach_seq($contig->seq);

   return $out;

}

=head2 fetch_by_contig_id

 Title   : fetch_by_contig_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_by_contig_id{
   my ($self,$cid) = @_;

   if( !defined $cid ) {
       $self->throw("fetch_by_contig_id must have an contig id");
   }

   my $sth = $self->prepare("select p.contig_id,p.contig_start,p.contig_end,p.contig_strand,p.hit_start,p.hit_end,p.hit_strand,p.hit_name,p.cigar_line,p.analysis_id, p.score from dna_align_feature p where p.contig_id = $cid");
   $sth->execute();

   my ($contig_id,$start,$end,$strand,$hstart,$hend,$hstrand,$hname,$cigar,$analysis_id, $score);

   $sth->bind_columns(undef,\$contig_id,\$start,\$end,\$strand,\$hstart,\$hend,\$hstrand,\$hname,\$cigar,\$analysis_id, \$score);

   my @f;
   my $contig = $self->db->get_RawContigAdaptor->fetch_by_dbID($cid);
   my %ana;

   while( $sth->fetch ) {
       if( !defined $ana{$analysis_id} ) {
	   $ana{$analysis_id} = $self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id);
       }


       my $out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();
       $out->start($start);
       $out->end($end);
       $out->strand($strand);
       $out->score($score);
       $out->hstart($hstart);
       $out->hend($hend);
       $out->hstrand($hstrand);
       $out->hseqname($hname);
       $out->cigar($cigar);

       $out->analysis($ana{$analysis_id});
       $out->seqname($contig->name);
       $out->attach_seq($contig->seq);
       push(@f,$out);
   }
   
   return @f;
}
sub fetch_by_contig_id_and_logic_name{

  my($self, $cid, $logic_name) = @_;

  
   if( !defined $cid ) {
       $self->throw("fetch_by_contig_id_and_logic_name must have an contig id");
   }

   if(!defined $logic_name){
     $self->throw("must provide a logic_name to fetch by logic name: $!\n");
   } 

   my $sth = $self->prepare("select p.contig_id,p.contig_start,p.contig_end,p.contig_strand,p.hit_start,p.hit_end,p.hit_name, p.cigar_line,p.analysis_id from dna_align_feature p, analysis a where p.contig_id = $cid and a.analysis_id = p.analysis_id and a.logic_name = '$logic_name'");
   $sth->execute();

   my ($contig_id,$start,$end,$strand,$hstart,$hend,$hname, $cigar,$analysis_id);

   $sth->bind_columns(undef,\$contig_id,\$start,\$end,\$strand,\$hstart,\$hend,\$hname,\$cigar,\$analysis_id);

   my @f;
   my $contig = $self->db->get_RawContigAdaptor->fetch_by_dbID($cid);
   my %ana;

   while( $sth->fetch ) {
       if( !defined $ana{$analysis_id} ) {
	   $ana{$analysis_id} = $self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id);
       }


       my $out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();;
       $out->start($start);
       $out->end($end);
       $out->strand($strand);

       $out->hstart($hstart);
       $out->hend($hend);
       $out->hseqname($hname);
       
       $out->cigar($cigar);

       $out->analysis($ana{$analysis_id});
       $out->seqname($contig->name);
       $out->attach_seq($contig->seq);
       push(@f,$out);
   }
   
   return @f;

}

sub fetch_by_contig_id_and_dbname{

  my($self, $cid, $db_name) = @_;

  
   if( !defined $cid ) {
       $self->throw("fetch_by_contig_id_and_logic_name must have an contig id");
   }

   if(!defined $db_name){
     $self->throw("must provide a db_name to fetch by db name: $!\n");
   }

   my $sth1 = $self->prepare("select analysis_id from analysis where db = '$db_name'");

   $sth1->execute();
   my @analysis_ids;
   my $array_ref;
   while($array_ref = $sth1->fetchrow_arrayref){
        my $analysis_id = $array_ref->[0];
        push(@analysis_ids, $analysis_id);
   }
   my $analysis_idlist = join(',', @analysis_ids);
   my $sth = $self->prepare("select p.contig_id,p.contig_start,p.contig_end,p.contig_strand,p.hit_start,p.hit_end,p.hit_name, p.cigar_line,p.analysis_id from dna_align_feature p where p.contig_id = $cid and p.analysis_id in($analysis_idlist)");
   $sth->execute();

   my ($contig_id,$start,$end,$strand,$hstart,$hend,$hname, $cigar,$analysis_id);

   $sth->bind_columns(undef,\$contig_id,\$start,\$end,\$strand,\$hstart,\$hend,\$hname,\$cigar,\$analysis_id);

   my @f;
   my $contig = $self->db->get_RawContigAdaptor->fetch_by_dbID($cid);
   my %ana;

   while( $sth->fetch ) {
       if( !defined $ana{$analysis_id} ) {
	   $ana{$analysis_id} = $self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id);
       }


       my $out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();;
       $out->start($start);
       $out->end($end);
       $out->strand($strand);

       $out->hstart($hstart);
       $out->hend($hend);
       $out->hseqname($hname);
       
       $out->cigar($cigar);

       $out->analysis($ana{$analysis_id});
       $out->seqname($contig->name);
       $out->attach_seq($contig->seq);
       push(@f,$out);
   }
   
   return @f;

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
  my ($self,$start,$end,$chr,$type) = @_;
  
  if( !defined $type ) {
    $self->throw("Assembly location must be start,end,chr,type");
  }
  
  if( $start !~ /^\d/ || $end !~ /^\d/ ) {
    $self->throw("start/end must be numbers not $start,$end (have you typed the location in the right way around - start,end,chromosome,type");
  }
  
  my $mapper = $self->db->get_AssemblyMapperAdaptor->fetch_by_type($type);
  
  $mapper->register_region($chr,$start,$end);

  my @cids = $mapper->list_contig_ids($chr, $start ,$end);
  
  # build the SQL
  print STDERR "have @cids contig ids\n";
  my $cid_list = join(',',@cids);
  my $sql = "select p.contig_id,p.contig_start,p.contig_end,p.contig_strand,p.hit_start,p.hit_end,p.hit_strand,p.hit_name,p.cigar_line,a.gff_source,a.gff_feature from dna_align_feature p where p.contig_id in ($cid_list)";
  
  my $sth = $self->prepare("select p.contig_id,p.contig_start,p.contig_end,p.contig_strand,p.hit_start,p.hit_end,p.hit_strand,p.hit_name,p.cigar_line, p.analysis_id from dna_align_feature p where p.contig_id in ($cid_list)");
  $sth->execute();
  
  
  my ($contig_id,$start,$end,$strand,$hstart,$hend,$hstrand,$hname,$cigar,$analysis_id);
  
  $sth->bind_columns(undef,\$contig_id,\$start,\$end,\$strand,\$hstart,\$hend,\$hstrand, \$hname,\$cigar,\$analysis_id);
  

  my @f;
  my %ana;
  my $counter = 0;
  while( $sth->fetch ) {
    # we whether this is sensible to use or not
    
    my @coord_list = $mapper->map_coordinates_to_assembly($contig_id, $start,$end,$strand,"rawcontig");
       
    # coord list > 1 - means does not cleanly map. At the moment, skip
    if( scalar(@coord_list) > 1 ) {
      next;
    }
    
    if( !defined $ana{$analysis_id} ) {
      $ana{$analysis_id} = $self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id);
    }

    # ok, ready to build a sequence feature: do we want this relative or not?
    
    my $out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();;
    if($coord_list[0]->isa("Bio::EnsEMBL::Mapper::Gap")){
      #print STDERR "the gap is from ".$coord_list[0]->start." to ".$coord_list[0]->end." on contig ".$contig_id."\n";
      $self->warn("feature is on a piece of contig not on golden path or in a gap skipping as not needed\n");
      $counter++;
      next;
    }
    
    $out->start($coord_list[0]->start);
    $out->end($coord_list[0]->end);
    $out->strand($coord_list[0]->strand);
    $out->seqname($coord_list[0]->id);
    
    $out->hstart($hstart);
    $out->hend($hend);
    $out->hstrand($hstrand);
    $out->hseqname($hname);
    $out->cigar($cigar);
    
    $out->analysis($ana{$analysis_id});
       
    push(@f,$out);
  }
  #print STDERR "have ".$counter." gaps\n";
  return @f;

}
sub fetch_by_assembly_location_and_dbname{

    my ($self,$start,$end,$chr,$type,$db_name) = @_;

   if( !defined $type ) {
       $self->throw("Assembly location must be start,end,chr,type");
   }

   if( $start !~ /^\d/ || $end !~ /^\d/ ) {
       $self->throw("start/end must be numbers not $start,$end (have you typed the location in the right way around - start,end,chromosome,type");
   }

   if(!defined $db_name){
  
      $self->throw($db_name." not defined must have a db name to get features\n");

   } 

   my $mapper = $self->db->get_AssemblyMapperAdaptor->fetch_by_type($type);
   $mapper->register_region($start,$end,$chr);

   my @cids = $mapper->list_contig_ids($start,$end,$chr);

   # build the SQL

   my $cid_list = join(',',@cids);
    

    my $sth1 = $self->prepare("select analysis_id from analysis where db = '$db_name'");
   

   $sth1->execute();
   my @analysis_ids;
   my $array_ref;
   while($array_ref = $sth1->fetchrow_arrayref){
        my $analysis_id = $array_ref->[0];
        push(@analysis_ids, $analysis_id);
   }
   my $analysis_idlist = join(',', @analysis_ids);
   my $sth = $self->prepare("select p.contig_id,p.contig_start,p.contig_end,p.contig_strand,p.hit_start,p.hit_end,p.hit_name, p.cigar_line,p.analysis_id from dna_align_feature p where p.contig_id in($cid_list) and a.analysis_id = p.analysis_id and p.analysis_id in($analysis_idlist)");
   $sth->execute();

   my ($contig_id,$start,$end,$strand,$hstart,$hend,$hname, $cigar,$analysis_id);

   $sth->bind_columns(undef,\$contig_id,\$start,\$end,\$strand,\$hstart,\$hend,\$hname,\$cigar,\$analysis_id);

   my @f;
    my %ana;

   while( $sth->fetch ) {
       # we whether this is sensible to use or not
       my @coord_list = $mapper->map_coordinates_to_assembly($start,$end,$strand,$contig_id,"rawcontig");
       
       # coord list > 1 - means does not cleanly map. At the moment, skip
       if( scalar(@coord_list) > 1 ) {
	   next;
       }

       if( !defined $ana{$analysis_id} ) {
	   $ana{$analysis_id} = $self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id);
       }

       # ok, ready to build a sequence feature: do we want this relative or not?

       my $out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();;
       $out->start($coord_list[0]->start);
       $out->end($coord_list[0]->end);
       $out->strand($coord_list[0]->strand);
       $out->seqname($coord_list[0]->seqname);

       $out->hstart($hstart);
       $out->hend($hend);
       $out->hseqname($hname);
       $out->cigar($cigar);

       $out->analysis($ana{$analysis_id});
       
       push(@f,$out); 
   }
   
   return @f;

}

sub fetch_by_assembly_location_and_logic_name{

  my($self, $start,$end,$chr,$type, $logic_name) = @_;

  
   if( !defined $type ) {
       $self->throw("Assembly location must be start,end,chr,type");
   }
   if(!defined $logic_name){
     $self->throw("must provide a logic_name to fetch by logic name: $!\n");
   } 

   if( $start !~ /^\d/ || $end !~ /^\d/ ) {
       $self->throw("start/end must be numbers not $start,$end (have you typed the location in the right way around - start,end,chromosome,type");
   }

   my $mapper = $self->db->get_AssemblyMapperAdaptor->fetch_by_type($type);
   $mapper->register_region($start,$end,$chr);

   my @cids = $mapper->list_contig_ids($start,$end,$chr);

   # build the SQL

   my $cid_list = join(',',@cids);

   my $sth = $self->prepare("select p.contig_id,p.contig_start,p.contig_end,p.contig_strand,p.hit_start,p.hit_end,p.hit_name, p.cigar_line,p.analysis_id from dna_align_feature p, analysis a where p.contig_id in($cid_list) and a.analysis_id = p.analysis_id and a.logic_name = '$logic_name'");
   $sth->execute();

   my ($contig_id,$start,$end,$strand,$hstart,$hend,$hname, $cigar,$analysis_id);

   $sth->bind_columns(undef,\$contig_id,\$start,\$end,\$strand,\$hstart,\$hend,\$hname,\$cigar,\$analysis_id);

   my @f;
   
   my %ana;

   while( $sth->fetch ) {
       # we whether this is sensible to use or not
       my @coord_list = $mapper->map_coordinates_to_assembly($start,$end,$strand,$contig_id,"rawcontig");
       
       # coord list > 1 - means does not cleanly map. At the moment, skip
       if( scalar(@coord_list) > 1 ) {
	   next;
       }

       if( !defined $ana{$analysis_id} ) {
	   $ana{$analysis_id} = $self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id);
       }

       # ok, ready to build a sequence feature: do we want this relative or not?

       my $out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();;
       $out->start($coord_list[0]->start);
       $out->end($coord_list[0]->end);
       $out->strand($coord_list[0]->strand);
       $out->seqname($coord_list[0]->seqname);

       $out->hstart($hstart);
       $out->hend($hend);
       $out->hseqname($hname);
       $out->cigar($cigar);

       $out->analysis($ana{$analysis_id});
       
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

   my $sth = $self->prepare("insert into dna_align_feature (contig_id,contig_start,contig_end,contig_strand,hit_start,hit_end,hit_strand,hit_name,cigar_line,analysis_id,score) values (?,?,?,?,?,?,?,?,?,?,?)");

   foreach my $sf ( @sf ) {
       if( !ref $sf || !$sf->isa("Bio::EnsEMBL::FeaturePair") ) {
	   $self->throw("Simple feature must be an Ensembl DnaAlignFeature, not a [$sf]");
       }

       if( !defined $sf->analysis ) {
	   $self->throw("Cannot store sequence features without analysis");
       }
       if( !defined $sf->analysis->dbID ) {
	   # maybe we should throw here. Shouldn't we always have an analysis from the database?
	   $self->throw("I think we should always have an analysis object which has originated from the database. No dbID, not putting in!");
       }
      
       $sth->execute($contig_id,$sf->start,$sf->end,$sf->strand,$sf->hstart,$sf->hend,$sf->hstrand,$sf->hseqname,$sf->cigar_string,$sf->analysis->dbID,$sf->score);
   }


}



sub fetch_featurepair_list_by_contig_id{
  my($self, $contig_id) = @_;

  my @cigar_feats = $self->fetch_by_contig_id($contig_id);
  my @fps;
  foreach my $cigar_feat(@cigar_feats){
    my $f1 = Bio::EnsEMBL::SeqFeature->new();
    my $f2 = Bio::EnsEMBL::SeqFeature->new();
    
    $f1->start($cigar_feat->start);
    $f1->end($cigar_feat->end);
    $f1->score($cigar_feat->score);
    $f1->seqname($cigar_feat->seqname);
    $f1->strand($cigar_feat->strand);
    
    $f2->start($cigar_feat->hstart);
    $f2->end($cigar_feat->hend);
    $f2->strand($cigar_feat->hstrand);
    $f2->seqname($cigar_feat->hseqname);

    my $cigar = $cigar_feat->cigar;
    
    my $dnadna =  Bio::EnsEMBL::DnaDnaAlignFeature->new(-feature1 => $f1,
							-feature2 => $f2,
							-cigar_string    => $cigar);
    my @parsed_fps = $dnadna->_parse_cigar;
    push(@fps, @parsed_fps);
  }

  return @fps;
}

sub fetch_featurepair_list_by_contig_id_and_logic_name{
  my($self, $contig_id, $logic_name) = @_;

  my @cigar_feats = $self->fetch_by_contig_id_and_logic_name($contig_id, $logic_name);
  my @fps;
  foreach my $cigar_feat(@cigar_feats){
    my $f1 = Bio::EnsEMBL::SeqFeature->new();
    my $f2 = Bio::EnsEMBL::SeqFeature->new();
    
    $f1->start($cigar_feat->start);
    $f1->end($cigar_feat->end);
    $f1->score($cigar_feat->score);
    $f1->seqname($cigar_feat->seqname);
    $f1->strand($cigar_feat->strand);
    
    $f2->start($cigar_feat->hstart);
    $f2->end($cigar_feat->hend);
    $f2->strand($cigar_feat->hstrand);
    $f2->seqname($cigar_feat->hseqname);

    my $cigar = $cigar_feat->cigar;
    
    my $dnadna =  Bio::EnsEMBL::DnaDnaAlignFeature->new(-feature1 => $f1,
							-feature2 => $f2,
							-cigar_string    => $cigar);
    my @parsed_fps = $dnadna->_parse_cigar;
    push(@fps, @parsed_fps);
  }

  return @fps;
}


sub fetch_featurepair_list_by_contig_id_and_dbname{
  my($self, $contig_id, $db_name) = @_;

  my @cigar_feats = $self->fetch_by_contig_id_and_dbname($contig_id, $db_name);
  my @fps;
  foreach my $cigar_feat(@cigar_feats){
    my $f1 = Bio::EnsEMBL::SeqFeature->new();
    my $f2 = Bio::EnsEMBL::SeqFeature->new();
    
    $f1->start($cigar_feat->start);
    $f1->end($cigar_feat->end);
    $f1->score($cigar_feat->score);
    $f1->seqname($cigar_feat->seqname);
    $f1->strand($cigar_feat->strand);
    
    $f2->start($cigar_feat->hstart);
    $f2->end($cigar_feat->hend);
    $f2->strand($cigar_feat->hstrand);
    $f2->seqname($cigar_feat->hseqname);

    my $cigar = $cigar_feat->cigar;
    
    my $dnadna =  Bio::EnsEMBL::DnaDnaAlignFeature->new(-feature1 => $f1,
							-feature2 => $f2,
							-cigar_string    => $cigar);
    my @parsed_fps = $dnadna->_parse_cigar;
    push(@fps, @parsed_fps);
  }

  return @fps;
}

sub fetch_featurepair_list_by_dbID{
  my($self, $dbId) = @_;

  my $cigar_feat = $self->fetch_by_dbId($dbId);
  my @fps;
  
    my $f1 = Bio::EnsEMBL::SeqFeature->new();
    my $f2 = Bio::EnsEMBL::SeqFeature->new();
    
    $f1->start($cigar_feat->start);
    $f1->end($cigar_feat->end);
    $f1->score($cigar_feat->score);
    $f1->seqname($cigar_feat->seqname);
    $f1->strand($cigar_feat->strand);
    
    $f2->start($cigar_feat->hstart);
    $f2->end($cigar_feat->hend);
    $f2->strand($cigar_feat->hstrand);
    $f2->seqname($cigar_feat->hseqname);

    my $cigar = $cigar_feat->cigar;
    
    my $dnadna =  Bio::EnsEMBL::DnaDnaAlignFeature->new(-feature1 => $f1,
							-feature2 => $f2,
							-cigar_string    => $cigar);
    my @parsed_fps = $dnadna->_parse_cigar;
    push(@fps, @parsed_fps);
  

  return @fps;
}


sub fetch_featurepair_list_by_assembly_location{
  my($self, $start, $end, $chr, $type) = @_;

  my @cigar_feats = $self->fetch_by_assembly_location($start, $end, $chr, $type);
  my @fps;
  foreach my $cigar_feat(@cigar_feats){
    my $f1 = Bio::EnsEMBL::SeqFeature->new();
    my $f2 = Bio::EnsEMBL::SeqFeature->new();
    
    $f1->start($cigar_feat->start);
    $f1->end($cigar_feat->end);
    $f1->score($cigar_feat->score);
    $f1->seqname($cigar_feat->seqname);
    $f1->strand($cigar_feat->strand);
    
    $f2->start($cigar_feat->hstart);
    $f2->end($cigar_feat->hend);
    $f2->strand($cigar_feat->hstrand);
    $f2->seqname($cigar_feat->hseqname);

    my $cigar = $cigar_feat->cigar;
    
    my $dnadna =  Bio::EnsEMBL::DnaDnaAlignFeature->new(-feature1 => $f1,
							-feature2 => $f2,
							-cigar_string    => $cigar);
    my @parsed_fps = $dnadna->_parse_cigar;
    push(@fps, @parsed_fps);
  }

  return @fps;
}

sub fetch_featurepair_list_by_assembly_location_and_dbname{
  my($self, $start, $end, $chr, $type, $db_name) = @_;

  my @cigar_feats = $self->fetch_by_assembly_location($start, $end, $chr, $type, $db_name);
  my @fps;
  foreach my $cigar_feat(@cigar_feats){
    my $f1 = Bio::EnsEMBL::SeqFeature->new();
    my $f2 = Bio::EnsEMBL::SeqFeature->new();
    
    $f1->start($cigar_feat->start);
    $f1->end($cigar_feat->end);
    $f1->score($cigar_feat->score);
    $f1->seqname($cigar_feat->seqname);
    $f1->strand($cigar_feat->strand);
    
    $f2->start($cigar_feat->hstart);
    $f2->end($cigar_feat->hend);
    $f2->strand($cigar_feat->hstrand);
    $f2->seqname($cigar_feat->hseqname);

    my $cigar = $cigar_feat->cigar;
    
    my $dnadna =  Bio::EnsEMBL::DnaDnaAlignFeature->new(-feature1 => $f1,
							-feature2 => $f2,
							-cigar_string    => $cigar);
    my @parsed_fps = $dnadna->_parse_cigar;
    push(@fps, @parsed_fps);
  }

  return @fps;
}

sub fetch_featurepair_list_by_assembly_location_nad_logic_name{
  my($self, $start, $end, $chr, $type, $logic_name) = @_;

  my @cigar_feats = $self->fetch_by_assembly_location_and_logic_name($start, $end, $chr, $type, $logic_name);
  my @fps;
  foreach my $cigar_feat(@cigar_feats){
    my $f1 = Bio::EnsEMBL::SeqFeature->new();
    my $f2 = Bio::EnsEMBL::SeqFeature->new();
    
    $f1->start($cigar_feat->start);
    $f1->end($cigar_feat->end);
    $f1->score($cigar_feat->score);
    $f1->seqname($cigar_feat->seqname);
    $f1->strand($cigar_feat->strand);
    
    $f2->start($cigar_feat->hstart);
    $f2->end($cigar_feat->hend);
    $f2->strand($cigar_feat->hstrand);
    $f2->seqname($cigar_feat->hseqname);

    my $cigar = $cigar_feat->cigar;
    
    my $dnadna =  Bio::EnsEMBL::DnaDnaAlignFeature->new(-feature1 => $f1,
							-feature2 => $f2,
							-cigar_string    => $cigar);
    my @parsed_fps = $dnadna->_parse_cigar;
    push(@fps, @parsed_fps);
  }

  return @fps;

}

1;


