
#
# Ensembl module for Bio::EnsEMBL::DBSQL::ExternalWrapper
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::ExternalWrapper - Makes a standard Ensembl database a ExternalFeatureFactoryI implementing object

=head1 SYNOPSIS

    # check out DB::ExternalFeatureFactoryI

=head1 DESCRIPTION

This class wraps a standard Ensembl database as if it is an
ExternalFeatureFactory interface, allowing it to serve up Genes and 
(assumming someone gets to write this as well) features.

The idea here is that this database will contain different data (eg, 
data on EMBL CDS from the original entry) which is updated at a different
cycle from the main stuff.

=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::ExternalWrapper;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::RootI
use Bio::EnsEMBL::Gene;
use Bio::Root::RootI;
use Bio::EnsEMBL::DB::ExternalFeatureFactoryI;

@ISA = qw(Bio::EnsEMBL::DB::ExternalFeatureFactoryI Bio::Root::RootI);

sub new {
  my($class,$dbobj) = @_;
  

  my $self = {};
  bless $self,$class;

  if( !defined $dbobj  ) {
      $self->throw("No dbobj or not a dbobj [$dbobj]");
  }
  
  $self->dbobj($dbobj);

  return $self;
}


=head2 get_Ensembl_Genes_clone

 Title   : get_Ensembl_Genes_clone
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Ensembl_Genes_clone{
   my ($self,$cloneid) = @_;
   my $clone;

   my $dbobj = $self->dbobj;

   #print STDERR "Got dbobj $dbobj connected to ",$dbobj->dbname,"\n";

   eval {
       $clone = $self->dbobj->get_Clone($cloneid);
   };

   if( $@ ) {
       # return nothing
       return ();
   }

   my @genes=$clone->get_all_Genes();
   #foreach my $gene ( @genes ) {
   #    print STDERR "got ",$gene->id,"\n";
   #}

   return $clone->get_all_Genes();
}

=head2 get_Ensembl_Genes_contig_list

 Title   : get_Ensembl_Genes_contig_list
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Ensembl_Genes_contig_list{
   my ($self,@contigs) = @_;

   if( scalar @contigs == 0 ) {
       return ();
   }
   #Hash of contigs by gene id
   my %gc;

   #Hash of array of gene objects by contig
   my %cg;
   my $list;
   my @genes;
   my @todocontigs;
   if (my $cgr = $self->cg) {
       #Get them from the cache
       #print STDERR "Getting external genes from the cache\n";
       my %cgh = %$cgr;
       foreach my $c (@contigs) {
	   if ($cgh{$c}) {
	       foreach my $g (@{$cgh{$c}}) {
		   $g->refresh();
	       }
	       push (@genes,@{$cgh{$c}});
	   }
	   else {
	       push (@todocontigs,$c);
	   }
       }
       @genes;
   }
   else {
       foreach my $c (@contigs) {
	   $cg{$c}=[];
       }
       push (@todocontigs,@contigs);
   }
   if( scalar @todocontigs == 0 ) {
       #print STDERR "Returning ".scalar(@genes)." $genes[0] genes...\n";
       return @genes;
   }
   else {
       foreach my $c ( @todocontigs ) {
	   $list .= "'$c',";
       }
       chop $list;
       $list = "($list)";
       
       my $sth = $self->dbobj->prepare("select t.gene,c.id from transcript t,exon_transcript et,exon e,contig c where c.id in $list and c.internal_id = e.contig and e.id = et.exon and t.id = et.transcript");
   
       $sth->execute();
       my @geneids;
       
       my %seen;
       while( my ($id,$contig) = $sth->fetchrow_array ) {
	   if (!$seen{$id}) {
	       $gc{$id}=$contig;
	       $seen{$id}++;
	   }
       }
       push(@geneids,keys(%gc));
       if( scalar(@geneids) == 0 ) {
	   #print STDERR "No ids here...\n";
	   return();
       }
           
      
       
       #print STDERR "Getting external genes normally\n";
       @genes = $self->dbobj->gene_Obj->get_array_supporting('none',@geneids);
       foreach my $gene (@genes) {
	   if (my $contig = $gc{$gene->id}) {
	       push(@{$cg{$contig}},$gene);
	   }
       }
       $self->cg(\%cg);
       
       #print STDERR "Returning ".scalar(@genes)." $genes[0] genes...\n";
       return @genes;
   }
}

=head2 cg

 Title   : cg
 Usage   : $obj->cg($newval)
 Function: Getset for cg value
 Returns : value of cg
 Args    : newvalue (optional)


=cut

sub cg{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'cg'} = $value;
    }
    return $obj->{'cg'};

}

=head2 get_Ensembl_SeqFeatures_contig

 Title   : get_Ensembl_SeqFeatures_contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Ensembl_SeqFeatures_contig{
   my ($self,$contigid) = @_;
   return;
}



=head2 get_Ensembl_SeqFeatures_clone

 Title   : get_Ensembl_SeqFeatures_clone
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut




sub get_Ensembl_SeqFeatures_clone{
   my ($self,$contigid) = @_;
   return;
}










=head2 dbobj

 Title   : dbobj
 Usage   : $obj->dbobj($newval)
 Function: 
 Returns : value of dbobj
 Args    : newvalue (optional)


=cut

sub dbobj{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'dbobj'} = $value;
    }
    return $obj->{'dbobj'};

}







