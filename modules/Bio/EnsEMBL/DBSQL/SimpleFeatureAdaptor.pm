

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

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::SimpleFeature;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);


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

       $sth->execute($contig_id,$sf->start,$sf->end,$sf->strand,$sf->display_label,$sf->analysis->dbID,$sf->score);
   }


}

sub _tablename {
  my $self = shift;
  
  return "simple_feature";
}

sub _columns {
  my $self = shift;

  return qw( simple_feature_id contig_id contig_start contig_end contig_strand
	     display_label analysis_id score );
}


sub _obj_from_hashref {
  my ($self, $hashref) = @_;

  my $aa = $self->db()->get_AnalysisAdaptor();
  my $analysis = $aa->fetch_by_dbID($hashref->{'analysis_id'});
  
  my $rca = $self->db()->get_RawContigAdaptor();
  my $contig = $rca->fetch_by_dbID($hashref->{'contig_id'});

  my $out = Bio::EnsEMBL::SimpleFeature->new();
  $out->start($hashref->{'contig_start'});
  $out->end($hashref->{'contig_end'});
  $out->strand($hashref->{'contig_strand'});
  $out->analysis($analysis);
  $out->display_label($hashref->{'display_label'});
  $out->seqname($contig->name());
  $out->attach_seq($contig); 

  if($hashref->{'score'}) {
    $out->score($hashref->{'score'});
  }

  return $out;
}

1;
