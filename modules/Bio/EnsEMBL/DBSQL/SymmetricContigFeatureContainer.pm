
#
# Ensembl module for Bio::EnsEMBL::DBSQL::SymmetricContigFeatureContainer
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::SymmetricContigFeatureContainer - Binds to SymmetricContigFeature table

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::SymmetricContigFeatureContainer;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI


use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::FeatureFactory;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 get_FeaturePair_list_by_rawcontig_id

 Title   : get_FeaturePair_list_by_rawcontig_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_FeaturePair_list_by_rawcontig_id{
   my ($self,$id) = @_;

   if( !defined $id ) {
       $self->throw("Must have a raw contig id");
   }

   my $sth = $self->prepare("select a.seq_start,a.seq_end,a.strand,b.seq_start,b.seq_end,b.strand,b.rawcontigid,p.score  from symmetric_contig_feature a, symmetric_contig_pair_hit p,symmetric_contig_feature b where a.symchid = p.symchid and p.symchid = b.symchid and a.symcfid != b.symcfid and a.rawcontigid = $id");
   
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


1;


