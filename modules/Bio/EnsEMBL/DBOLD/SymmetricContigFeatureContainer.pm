
#
# Ensembl module for Bio::EnsEMBL::DBOLD::SymmetricContigFeatureContainer
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBOLD::SymmetricContigFeatureContainer - Binds to SymmetricContigFeature table

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

This module is a container for symmetric contig feature pairs, i.e. pairs of features between contigs that have identical sequence in two versions of a database. The pairs are stored symmetrically, i.e. each feature on each contig is stored in the contig_feature table, and each pair is stored with an id in a separate table.

The method which fetches teh feature pairs breaks the symmetry by asking for all the feature pairs with a certain version of the clone (on which the contig is sitting). The crosmmatching at the moment relies on the sv version of the clones.

=head1 AUTHOR - Ewan Birney

This module is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBOLD::SymmetricContigFeatureContainer;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI


use Bio::EnsEMBL::DBOLD::BaseAdaptor;
use Bio::EnsEMBL::FeatureFactory;

@ISA = qw(Bio::EnsEMBL::DBOLD::BaseAdaptor);

=head2 get_FeaturePair_list_by_rawcontig_id

 Title   : get_FeaturePair_list_by_rawcontig_id
 Usage   : $scfc->get_FeaturePair_list_by_rawcontig_id($rid,15)
 Function: gets all the feature pairs for a specific rawcontig id 
           and clone version
 Example : $scfc->get_FeaturePair_list_by_rawcontig_id('AC000043.12',5)
 Returns : array of Bio::EnsEMBL::FeaturePair
 Args    : id of the rawcontig,sv version of the clone it is on


=cut

sub get_FeaturePair_list_by_rawcontig_id{
   my ($self,$id,$version) = @_;

   if( !defined $id ) {
       $self->throw("Must have a raw contig id");
   }
   if( !defined $version ) {
       $self->throw("Must have a raw contig version");
   }

   my $sth = $self->prepare("select a.seq_start,a.seq_end,a.strand,b.seq_start,b.seq_end,b.strand,b.rawcontigid,p.score  from symmetric_contig_feature a, symmetric_contig_pair_hit p,symmetric_contig_feature b where a.symchid = p.symchid and p.symchid = b.symchid and a.symcfid != b.symcfid and a.rawcontigid = '$id' and a.rawversion=$version");
   #print STDERR "SQL: select a.seq_start,a.seq_end,a.strand,b.seq_start,b.seq_end,b.strand,b.rawcontigid,p.score  from symmetric_contig_feature a, symmetric_contig_pair_hit p,symmetric_contig_feature b where a.symchid = p.symchid and p.symchid = b.symchid and a.symcfid != b.symcfid and a.rawcontigid = '$id'";
   $sth->execute;
   my @out;
   while( my $aref = $sth->fetchrow_arrayref ) {
       my ($start,$end,$strand,$hstart,$hend,$hstrand,$hname,$score) = @{$aref};
       my $out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();
       $out->set_all_fields($start,$end,$strand,$score,$id,'symmetric',$id,
			    $hstart,$hend,$hstrand,$score,$hname,'symmetric',$hname);
       push(@out,$out);
   }

   return @out;
}


=head2 write_FeaturePair_List

 Title   : write_FeaturePair_List
 Usage   : $scfc->write_FeaturePair_List(@fp)
 Function: Writes an array of feature pairs to the db
 Example : $scfc0>write_FeaturePair_List(@fp)
 Returns : nothing
 Args    : array of Bio::EnsEMBL::FeaturePair


=cut

sub write_FeaturePair_List{
   my ($self,@fp) = @_;

   foreach my $fp ( @fp ) {
       my $score = $fp->score;
       my $sth = $self->prepare("INSERT INTO symmetric_contig_pair_hit (symchid,score) VALUES (NULL,$score)");
       $sth->execute;
       $sth = $self->prepare("select LAST_INSERT_ID()");
       $sth->execute;
       my ($hitid) = $sth->fetchrow_array;

       my $seqname = $fp->feature1->seqname;
       $seqname =~ /(\S+)\.(\d+)\.(\S+)/ || $self->throw("Feature pair name does not conform to acc.version.number sequence");
       my $version = $2;
       my $contigid = "$1.$3";
       my $clone=$1;
       $sth = $self->prepare("INSERT INTO symmetric_contig_feature (symcfid,symchid,rawcontigid,rawversion,clone,seq_start,seq_end,strand) VALUES (NULL,$hitid,'".$contigid."',".$version.",'".$clone."',".$fp->feature1->start.",".$fp->feature1->end.",".$fp->feature1->strand.")");
       $sth->execute;

       $seqname = $fp->feature2->seqname;
       $seqname =~ /(\S+)\.(\d+)\.(\S+)/ || $self->throw("Feature pair name does not conform to acc.version.number sequence");
       $version = $2;
       $contigid = "$1.$3";
       $clone=$1;
       $sth = $self->prepare("INSERT INTO symmetric_contig_feature (symcfid,symchid,rawcontigid,rawversion,clone,seq_start,seq_end,strand) VALUES (NULL,$hitid,'".$contigid."',".$version.",'".$clone."',".$fp->feature2->start.",".$fp->feature2->end.",".$fp->feature2->strand.")");
       $sth->execute;
   }

}


1;


