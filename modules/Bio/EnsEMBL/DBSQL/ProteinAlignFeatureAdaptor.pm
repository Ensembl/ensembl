

#
# BioPerl module for Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor - Adaptor for ProteinAlignFeatures

=head1 SYNOPSIS

    $pfadp = $dbadaptor->get_ProteinAlignFeatureAdaptor();

    my @feature_array = $pfadp->fetch_by_contig_id($contig_numeric_id);

    my @feature_array = $pfadp->fetch_by_assembly_location($start,$end,$chr,'UCSC');
 
    $pfadp->store($contig_numeric_id,@feature_array);


=head1 DESCRIPTION


This is an adaptor for protein features on DNA sequence. Like other
feature getting adaptors it has a number of fetch_ functions and a
store function.


=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::FeatureFactory;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::PepDnaAlignFeature;
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

   my $sth = $self->prepare("select p.contig_id,p.contig_start,p.contig_end,p.contig_strand,p.hit_start,p.hit_end,p.hit_name,p.cigar_line,p.analysis_id, p.score from protein_align_feature p where p.protein_align_feature_id = $id");
   $sth->execute();

   my ($contig_id,$start,$end,$strand,$hstart,$hend,$hname,$cigar,$analysis_id, $score) = $sth->fetchrow_array();

   if( !defined $contig_id ) {
       $self->throw("No simple feature with id $id");
   }

   my $contig = $self->db->get_RawContigAdaptor->fetch_by_dbID($contig_id);
   my $f1 = Bio::EnsEMBL::SeqFeature->new();
   my $f2 = Bio::EnsEMBL::SeqFeature->new();
   
   $f1->start($start);
   $f1->end($end);
   $f1->score($score);
   $f1->seqname($contig->name);
   $f1->strand($strand);
   
   $f2->start($hstart);
   $f2->end($hend);
   $f2->seqname($hname);
   
   
    
   my $dnapep =  Bio::EnsEMBL::DnaPepAlignFeature->new(-feature1 => $f1,
						       -feature2 => $f2,
						       -cigar_string => $cigar);
   $dnapep->analysis($self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id));
   $dnapep->seqname($contig->name);
   $dnapep->attach_seq($contig->seq);
  # my $out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();;
#   $out->start($start);
#   $out->end($end);
#   $out->strand($strand);
#   $out->score($score);
#   $out->hstart($hstart);
#   $out->hend($hend);
#   $out->hseqname($hname);
#   $out->cigar($cigar);
   
#   $out->analysis($self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id));

#   $out->seqname($contig->id);
#   $out->attach_seq($contig->seq);

   return $dnapep;

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

   my $sth = $self->prepare("select p.contig_id,p.contig_start,p.contig_end,p.contig_strand,p.hit_start,p.hit_end,p.hit_name, p.cigar_line,p.analysis_id, p.score from protein_align_feature p where p.contig_id = '$cid'");
   $sth->execute();

   my ($contig_id,$start,$end,$strand,$hstart,$hend,$hname, $cigar,$analysis_id, $score);

   $sth->bind_columns(undef,\$contig_id,\$start,\$end,\$strand,\$hstart,\$hend,\$hname,\$cigar,\$analysis_id, \$score);

   my @f;
   my $contig = $self->db->get_RawContigAdaptor->fetch_by_dbID($cid);
   my %ana;

   while( $sth->fetch ) {
     if( !defined $ana{$analysis_id} ) {
       $ana{$analysis_id} = $self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id);
     }
     my $f1 = Bio::EnsEMBL::SeqFeature->new();
     my $f2 = Bio::EnsEMBL::SeqFeature->new();
     
     $f1->start($start);
     $f1->end($end);
     $f1->score($score);
     $f1->seqname($contig->name);
     $f1->strand($strand);
       
     $f2->start($hstart);
     $f2->end($hend);
     $f2->seqname($hname);
   
       
       
     my $dnapep =  Bio::EnsEMBL::DnaPepAlignFeature->new(-feature1 => $f1,
							 -feature2 => $f2,
							 -cigar_string => $cigar);
     $dnapep->analysis($ana{$analysis_id});
     $dnapep->seqname($contig->name);
  
     push(@f,$dnapep);
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

   my $sth = $self->prepare("select p.contig_id,p.contig_start,p.contig_end,p.contig_strand,p.hit_start,p.hit_end,p.hit_name, p.cigar_line,p.analysis_id, p.score from protein_align_feature p, analysis a where p.contig_id = $cid and a.analysis_id = p.analysis_id and a.logic_name = '$logic_name'");
   $sth->execute();

   my ($contig_id,$start,$end,$strand,$hstart,$hend,$hname, $cigar,$analysis_id, $score);

   $sth->bind_columns(undef,\$contig_id,\$start,\$end,\$strand,\$hstart,\$hend,\$hname,\$cigar,\$analysis_id, \$score);

   my @f;
   my $contig = $self->db->get_RawContigAdaptor->fetch_by_dbID($cid);
   my %ana;

   while( $sth->fetch ) {
     if( !defined $ana{$analysis_id} ) {
       $ana{$analysis_id} = $self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id);
     }
     my $f1 = Bio::EnsEMBL::SeqFeature->new();
     my $f2 = Bio::EnsEMBL::SeqFeature->new();
     
     $f1->start($start);
     $f1->end($end);
     $f1->score($score);
     $f1->seqname($contig->name);
     $f1->strand($strand);
       
     $f2->start($hstart);
     $f2->end($hend);
     $f2->seqname($hname);
   
       
       
     my $dnapep =  Bio::EnsEMBL::DnaPepAlignFeature->new(-feature1 => $f1,
							 -feature2 => $f2,
							 -cigar_string => $cigar);
     $dnapep->analysis($ana{$analysis_id});
     $dnapep->seqname($contig->name);
  
     push(@f,$dnapep);
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
   my $sth = $self->prepare("select p.contig_id,p.contig_start,p.contig_end,p.contig_strand,p.hit_start,p.hit_end,p.hit_name, p.cigar_line,p.analysis_id, p.score from protein_align_feature p where p.contig_id = $cid and p.analysis_id in($analysis_idlist)");
   $sth->execute();

   my ($contig_id,$start,$end,$strand,$hstart,$hend,$hname, $cigar,$analysis_id, $score);

   $sth->bind_columns(undef,\$contig_id,\$start,\$end,\$strand,\$hstart,\$hend,\$hname,\$cigar,\$analysis_id, \$score);

   my @f;
   my $contig = $self->db->get_RawContigAdaptor->fetch_by_dbID($cid);
   my %ana;

   while( $sth->fetch ) {
     if( !defined $ana{$analysis_id} ) {
       $ana{$analysis_id} = $self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id);
     }
     my $f1 = Bio::EnsEMBL::SeqFeature->new();
     my $f2 = Bio::EnsEMBL::SeqFeature->new();
     
     $f1->start($start);
     $f1->end($end);
     $f1->score($score);
     $f1->seqname($contig->name);
     $f1->strand($strand);
       
     $f2->start($hstart);
     $f2->end($hend);
     $f2->seqname($hname);
   
       
       
     my $dnapep =  Bio::EnsEMBL::DnaPepAlignFeature->new(-feature1 => $f1,
							 -feature2 => $f2,
							 -cigar_string => $cigar);
     $dnapep->analysis($ana{$analysis_id});
     $dnapep->seqname($contig->name);
  
     push(@f,$dnapep); 
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

   my @cids = $mapper->list_contig_ids($chr, $start,$end);

   # build the SQL

   my $cid_list = join(',',@cids);
   my $sth = $self->prepare("select p.contig_id,p.contig_start,p.contig_end,p.contig_strand,p.hit_start,p.hit_end,p.hit_name,p.cigar_line, p.analysis_id, p.score from protein_align_feature p where p.contig_id in ($cid_list)");
   $sth->execute();


   my ($contig_id,$contig_start,$contig_end,$strand,$hstart,$hend,$hname,$cigar,$analysis_id, $score);

   $sth->bind_columns(undef,\$contig_id,\$contig_start,\$contig_end,\$strand,\$hstart,\$hend,\$hname,\$cigar,\$analysis_id, \$score);


   my @f;
   my %ana;

   while( $sth->fetch ) {
       # we whether this is sensible to use or not
       my @coord_list = $mapper->map_coordinates_to_assembly($contig_id,$contig_start,$contig_end,$strand,"rawcontig");
       
       # coord list > 1 - means does not cleanly map. At the moment, skip
       if( scalar(@coord_list) > 1 ) {
	   next;
       }

       if( !defined $ana{$analysis_id} ) {
	   $ana{$analysis_id} = $self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id);
       }
       if($coord_list[0]->isa("Bio::EnsEMBL::Mapper::Gap")){
	 #print STDERR "the gap is from ".$coord_list[0]->start." to ".$coord_list[0]->end." on contig ".$contig_id."\n";
	 $self->warn("feature is on a piece of contig not on golden path or in a gap skipping as not needed\n");
	 
	 next;
       }       
       # ok, ready to build a sequence feature: do we want this relative or not?

     
       my $f1 = Bio::EnsEMBL::SeqFeature->new();
       my $f2 = Bio::EnsEMBL::SeqFeature->new();
       
       $f1->start($coord_list[0]->start);
       $f1->end($coord_list[0]->end);
       $f1->score($score);
       $f1->seqname($coord_list[0]->id);
       $f1->strand($coord_list[0]->strand);
       
       $f2->start($hstart);
       $f2->end($hend);
       $f2->seqname($hname);
   
       
       
       my $dnapep =  Bio::EnsEMBL::DnaPepAlignFeature->new(-feature1 => $f1,
							 -feature2 => $f2,
							 -cigar_string => $cigar);
       $dnapep->analysis($ana{$analysis_id});
       $dnapep->seqname($coord_list[0]->id);
       push(@f,$dnapep);
   }

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
   $mapper->register_region($chr, $start,$end);

   my @cids = $mapper->list_contig_ids($chr, $start,$end);

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
   my $sth = $self->prepare("select p.contig_id,p.contig_start,p.contig_end,p.contig_strand,p.hit_start,p.hit_end,p.hit_name, p.cigar_line,p.analysis_id, p.score from protein_align_feature p where p.contig_id in($cid_list) and p.analysis_id in($analysis_idlist)");
   $sth->execute();

   my ($contig_id,$contig_start,$contig_end,$strand,$hstart,$hend,$hname, $cigar,$analysis_id, $score);

   $sth->bind_columns(undef,\$contig_id,\$contig_start,\$contig_end,\$strand,\$hstart,\$hend,\$hname,\$cigar,\$analysis_id, \$score);

   my @f;
    my %ana;

   while( $sth->fetch ) {
       
     # we whether this is sensible to use or not
     my @coord_list = $mapper->map_coordinates_to_assembly($contig_id,$contig_start,$contig_end,$strand,"rawcontig");
     
     # coord list > 1 - means does not cleanly map. At the moment, skip
     if( scalar(@coord_list) > 1 ) {
       next;
     }
     
     if( !defined $ana{$analysis_id} ) {
       $ana{$analysis_id} = $self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id);
     }
     if($coord_list[0]->isa("Bio::EnsEMBL::Mapper::Gap")){
       #print STDERR "the gap is from ".$coord_list[0]->start." to ".$coord_list[0]->end." on contig ".$contig_id."\n";
       $self->warn("feature is on a piece of contig not on golden path or in a gap skipping as not needed\n");
       
       next;
     }       
     # ok, ready to build a sequence feature: do we want this relative or not?
     
     
     my $f1 = Bio::EnsEMBL::SeqFeature->new();
     my $f2 = Bio::EnsEMBL::SeqFeature->new();
     
     $f1->start($coord_list[0]->start);
     $f1->end($coord_list[0]->end);
     $f1->score($score);
     $f1->seqname($coord_list[0]->id);
     $f1->strand($coord_list[0]->strand);
     
     $f2->start($hstart);
     $f2->end($hend);
     $f2->seqname($hname);
     
     
     
     my $dnapep =  Bio::EnsEMBL::DnaPepAlignFeature->new(-feature1 => $f1,
							 -feature2 => $f2,
							 -cigar_string => $cigar);
     $dnapep->analysis($ana{$analysis_id});
     $dnapep->seqname($coord_list[0]->id);
     push(@f,$dnapep);
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
   $mapper->register_region($chr, $start,$end);

   my @cids = $mapper->list_contig_ids($chr, $start,$end);

   # build the SQL

   my $cid_list = join(',',@cids);

   my $sth = $self->prepare("select p.contig_id,p.contig_start,p.contig_end,p.contig_strand,p.hit_start,p.hit_end,p.hit_name, p.cigar_line,p.analysis_id, p.score from protein_align_feature p, analysis a where p.contig_id in($cid_list) and a.analysis_id = p.analysis_id and a.logic_name = '$logic_name'");
   $sth->execute();

   my ($contig_id,$contig_start,$contig_end,$strand,$hstart,$hend,$hname, $cigar,$analysis_id, $score);

   $sth->bind_columns(undef,\$contig_id,\$contig_start,\$contig_end,\$strand,\$hstart,\$hend,\$hname,\$cigar,\$analysis_id, \$score);

   my @f;
   
   my %ana;

   while( $sth->fetch ) {
     # we whether this is sensible to use or not
     my @coord_list = $mapper->map_coordinates_to_assembly($contig_id,$contig_start,$contig_end,$strand,"rawcontig");
     
     # coord list > 1 - means does not cleanly map. At the moment, skip
     if( scalar(@coord_list) > 1 ) {
       next;
     }
     
     if( !defined $ana{$analysis_id} ) {
       $ana{$analysis_id} = $self->db->get_AnalysisAdaptor->fetch_by_dbID($analysis_id);
     }
     if($coord_list[0]->isa("Bio::EnsEMBL::Mapper::Gap")){
       #print STDERR "the gap is from ".$coord_list[0]->start." to ".$coord_list[0]->end." on contig ".$contig_id."\n";
       $self->warn("feature is on a piece of contig not on golden path or in a gap skipping as not needed\n");
       
       next;
     }       
     # ok, ready to build a sequence feature: do we want this relative or not?
     
     
     my $f1 = Bio::EnsEMBL::SeqFeature->new();
     my $f2 = Bio::EnsEMBL::SeqFeature->new();
     
     $f1->start($coord_list[0]->start);
     $f1->end($coord_list[0]->end);
     $f1->score($score);
     $f1->seqname($coord_list[0]->id);
     $f1->strand($coord_list[0]->strand);
     
     $f2->start($hstart);
     $f2->end($hend);
     $f2->seqname($hname);
     
     
     
     my $dnapep =  Bio::EnsEMBL::DnaPepAlignFeature->new(-feature1 => $f1,
							 -feature2 => $f2,
							 -cigar_string => $cigar);
     $dnapep->analysis($ana{$analysis_id});
     $dnapep->seqname($coord_list[0]->id);
     push(@f,$dnapep);
     
     
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

   my $sth = $self->prepare("insert into protein_align_feature (contig_id,contig_start,contig_end,contig_strand,hit_start,hit_end,hit_name,cigar_line,analysis_id,score, evalue, perc_ident) values (?,?,?,?,?,?,?,?,?,?, ?, ?)");

   foreach my $sf ( @sf ) {
       if( !ref $sf || !$sf->isa("Bio::EnsEMBL::FeaturePair") ) {
	   $self->throw("Simple feature must be an Ensembl ProteinAlignFeature, not a [$sf]");
       }

       if( !defined $sf->analysis ) {
	   $self->throw("Cannot store sequence features without analysis");
       }
       if( !defined $sf->analysis->dbID ) {
	   # maybe we should throw here. Shouldn't we always have an analysis from the database?
	   $self->throw("I think we should always have an analysis object which has originated from the database. No dbID, not putting in!");
       }
       #print STDERR "storing ".$sf->gffstring."\n";
       $sth->execute($contig_id,$sf->start,$sf->end,$sf->strand,$sf->hstart,$sf->hend,$sf->hseqname,$sf->cigar_string,$sf->analysis->dbID,$sf->score, $sf->p_value, $sf->percent_id);
   }


}




1;
