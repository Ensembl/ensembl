

#
# BioPerl module for Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Simple Feature Adaptor - database access for simple features 

=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::SimpleFeature;

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

   my $sth = $self->prepare("select s.contig_id,s.contig_start,s.contig_end,s.contig_strand,s.display_label,a.gff_source,a.gff_feature from simple_feature s,analysis a,contig c where s.simple_feature_id = $id and s.analysis_id = a.analysis_id");
   $sth->execute();

   my ($contig_id,$start,$end,$strand,$display,$gff_source,$gff_feature) = $sth->fetchrow_array();

   if( !defined $contig_id ) {
       $self->throw("No simple feature with id $id");
   }

   my $contig = $self->db->get_RawContigAdaptor->fetch_by_dbID($contig_id);
   my $out = Bio::EnsEMBL::SimpleFeature->new();
   $out->start($start);
   $out->end($end);
   $out->strand($strand);
   $out->primary_tag($gff_feature);
   $out->source_tag($gff_source);
   $out->display_text($display);
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

   my $sth = $self->prepare("select s.contig_id,s.contig_start,s.contig_end,s.contig_strand,s.display_label,a.gff_source,a.gff_feature from simple_feature s,analysis a where s.contig_id = $cid and s.analysis_id = a.analysis_id");
   $sth->execute();

   my ($contig_id,$start,$end,$strand,$display,$gff_source,$gff_feature);
   $sth->bind_columns(undef,\$contig_id,\$start,\$end,\$strand,\$display,\$gff_source,\$gff_feature);

   my @f;
   my $contig = $self->db->get_RawContigAdaptor->fetch_by_dbID($cid);
   while( $sth->fetch ) {
       my $out = Bio::EnsEMBL::SimpleFeature->new();
       $out->start($start);
       $out->end($end);
       $out->strand($strand);
       $out->primary_tag($gff_feature);
       $out->source_tag($gff_source);
       $out->display_text($display);
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
   $mapper->register_region($start,$end,$chr);

   my @cids = $mapper->list_contig_ids($start,$end,$chr);

   # build the SQL

   my $cid_list = join(',',@cids);
   my $sth = $self->prepare("select s.contig_id,s.contig_start,s.contig_end,s.contig_strand,s.display_label,a.gff_source,a.gff_feature from simple_feature s,analysis a where s.contig_id in ($cid_list) and s.analysis_id = a.analysis_id");
   $sth->execute();

   my ($contig_id,$start,$end,$strand,$display,$gff_source,$gff_feature);
   $sth->bind_columns(undef,\$contig_id,\$start,\$end,\$strand,\$display,\$gff_source,\$gff_feature);

   my @f;

   while( $sth->fetch ) {
       # we whether this is sensible to use or not
       my @coord_list = $mapper->map_coordinates($start,$end,$strand,$contig_id,"rawcontig");
       
       # coord list > 1 - means does not cleanly map. At the moment, skip
       if( scalar(@coord_list) > 1 ) {
	   next;
       }

       # ok, ready to build a sequence feature: do we want this relative or not?

       my $out = Bio::EnsEMBL::SimpleFeature->new();
       $out->start($coord_list[0]->start);
       $out->end($coord_list[0]->end);
       $out->seqname($coord_list[0]->seqname);
       $out->strand($coord_list[0]->strand);
       $out->primary_tag($gff_feature);
       $out->source_tag($gff_source);
       $out->display_text($display);
       
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

   my $sth = $self->prepare("insert into simple_feature (contig_id,contig_start,contig_end,contig_strand,display_label,analysis_id,score) values (?,?,?,?,?,?,?)");

   foreach my $sf ( @sf ) {
       if( !ref $sf || !$sf->isa("Bio::EnsEMBL::SimpleFeature") ) {
	   $self->throw("Simple feature must be an Ensembl SimpleFeature, not a [$sf]");
       }

       if( !defined $sf->analysis ) {
	   $self->throw("Cannot store sequence features without analysis");
       }
       if( !defined $sf->analysis->dbID ) {
	   # maybe we should throw here. Shouldn't we always have an analysis from the database?
	   $self->throw("I think we should always have an analysis object which has originated from the database. No dbID, not putting in!");
       }

       $sth->execute($contig_id,$sf->start,$sf->end,$sf->strand,$sf->display_text,$sf->analysis->dbID,$sf->score);
   }


}



1;
