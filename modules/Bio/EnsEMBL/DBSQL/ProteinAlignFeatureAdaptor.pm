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

Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor - 
Adaptor for ProteinAlignFeatures

=head1 SYNOPSIS

    $pfadp = $dbadaptor->get_ProteinAlignFeatureAdaptor();

    my @features = $pfadp->fetch_by_contig_id($contig_numeric_id);

    my @features = $pfadp->fetch_by_assembly_location($start,$end,$chr,'UCSC');
 
    $pfadp->store($contig_numeric_id, @features);


=head1 DESCRIPTION


This is an adaptor for protein features on DNA sequence. Like other
feature getting adaptors it has a number of fetch_ functions and a
store function.


=head1 AUTHOR - Ewan Birney

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root

use Bio::EnsEMBL::DBSQL::BaseAlignFeatureAdaptor;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::SeqFeature;
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAlignFeatureAdaptor);


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
       if( !ref $sf || !$sf->isa("Bio::EnsEMBL::DnaPepAlignFeature") ) {
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
       $sf->dbID($sth->{'mysql_insertid'});
   }


}


# 
# Internal functions not to called be anyone else
#

sub _obj_from_hashref {
  my ($self, $hashref) = @_;

  my $rca = $self->db()->get_RawContigAdaptor();
  my $contig = $rca->fetch_by_dbID($hashref->{'contig_id'});

  my $aa = $self->db()->get_AnalysisAdaptor();

  my $analysis = $aa->fetch_by_dbID($hashref->{'analysis_id'});

  my $f1 = Bio::EnsEMBL::SeqFeature->new();
  my $f2 = Bio::EnsEMBL::SeqFeature->new();
  
  $f1->start($hashref->{'contig_start'});
  $f1->end($hashref->{'contig_end'});
  $f1->strand($hashref->{'contig_strand'});
  
  $f2->start($hashref->{'hit_start'});
  $f2->end($hashref->{'hit_end'});
  $f2->strand(1);
  $f2->seqname($hashref->{'hit_name'});

  $f1->score($hashref->{'score'});
  $f1->p_value($hashref->{'evalue'});
  $f1->percent_id($hashref->{'perc_ident'});

  $f2->score($hashref->{'score'});
  $f2->p_value($hashref->{'evalue'});
  $f2->percent_id($hashref->{'perc_ident'});

  my $cigar = $hashref->{'cigar_line'};

  my $dnapep =  Bio::EnsEMBL::DnaPepAlignFeature->new(-feature1 => $f1,
						      -feature2 => $f2,
						      -cigar_string => $cigar);

  $dnapep->analysis($analysis);
  $dnapep->seqname($contig->name());
  $dnapep->attach_seq($contig);
  #set the 'id' of the feature to the hit name
  $dnapep->id($hashref->{'hit_name'});

  $dnapep->dbID($hashref->{'protein_align_feature_id'});
  
  return $dnapep;
}



sub _tablename {
  my $self = shift;

  return "protein_align_feature";
}

sub _columns {
  my $self = shift;
  
  return qw( protein_align_feature_id contig_id contig_start contig_end
	     analysis_id contig_strand hit_start hit_end hit_name cigar_line
	     evalue perc_ident score );
}

1;
