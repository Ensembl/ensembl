=head1 NAME

GeneCluster

=head1 SYNOPSIS


=head1 DESCRIPTION

This object holds one or more genes which has been clustered according to 
comparison criteria external to this class (for instance, in the 
methods compare and _compare_Genes methods of the class GeneComparison).
Each GeneCluster object holds the IDs of the genes clustered and the beginning and end coordinates
of each one (taken from the start and end coordinates of the first and last exon in the correspondig
each_unique_Exon aray)

=head1 CONTACT

eae@sanger.ac.uk

=cut

# Let the code begin ...

package Bio::EnsEMBL::Utils::GeneCluster;
use Bio::EnsEMBL::Gene;
use strict;

=head1 METHODS

=cut

#########################################################################


=head2 new()

new() initializes the attributes:
_gene_array
_geneID_array
_start
_end

=cut

sub new {
  my ($class,@args)=@_;

  if (ref($class)){
    $class = ref($class);
  }
  my $self = {};
  bless($self,$class);
  $self->{'_gene_array'}= \@args; # array reference that holds the list of genes in the cluster
  $self->{'_geneID_array'}=();    # array that holds the IDs of the genes
  $self->{'_start'}={};           # hash that holds the start position of each gene
  $self->{'_end'}={};             # hash that holds the end position of each gene
  
  foreach my $gene (@args){
    my ( $start , $end)=( _get_start($gene) , _get_end($gene) );
    push ( @ {$self->{'_geneID_array'} }, $gene->id );
    $self->{'_start'}->{$gene->id}=$start;
    $self-> {'_end'} ->{$gene->id}=$end;
  }
  return $self;
}

#########################################################################

=head2 put_Genes()

  function to include one or more genes in the cluster.
  Useful when creating a cluster. It takes as argument an array of genes, it returns nothing.

=cut

sub put_Genes {
  my ($self, @args)= @_;
  my @new_genes = @args;
  push ( @{ $self->{'_gene_array'} }, @new_genes );
  foreach my $new_gene (@new_genes){
    push ( @{ $self->{'_geneID_array'} }, $new_gene->id );
    my ( $start , $end)=( _get_start($new_gene) , _get_end($new_gene) );
    $self->{'_start'}->{$new_gene->id}=$start;
    $self-> {'_end'} ->{$new_gene->id}=$end;
  }
}

#########################################################################

=head2 get_Genes()

  it returns the array of genes in the GeneCluster object

=cut

sub get_Genes {
  my $self = shift @_;
  my @genes = @{ $self->{'_gene_array'} };
  return @genes;
}

#########################################################################

=head2 string()

  it returns a string containing the information from a GeneCluster object,
  i.e. the geneID, the start position and the end position.

=cut

sub string {
  my $self = shift @_;
  my $data='';
  foreach my $gene ( @{ $self->{'_gene_array'} } ){
    my $id = $gene->id;
    while (length($id)<16){
      $id .=' ';
    }
    $data .= $id."  "._get_start($gene)."  "._get_end($gene)."\n";
  }
  return $data;
}

#########################################################################

=head2 _get_start()

 function to get the start position of a gene - it reads the gene object and it returns
 the start position of the first exon

=cut

sub _get_start {
  my $gene = shift @_;
  my @exons = $gene->each_unique_Exon;
  my $st;
  
  if ($exons[0]->strand == 1) {
    @exons = sort {$a->start <=> $b->start} @exons;
    $st = $exons[0]->start;
  } else {
    @exons = sort {$b->start <=> $a->start} @exons; # they're read in opposite direction (from right to left)
    $st = $exons[0]->end;                           # the start is the end coordinate of the right-most exon
  }                                                 # which is here the first of the list @exons
  return $st;
}

#########################################################################

=head2 _get_end()

 function to get the end position of a gene - it reads the gene object and it returns
 the end position of the last exon

=cut

sub _get_end {
  my $gene = shift @_;
  my @exons = $gene->each_unique_Exon;
  my $end;
  
  if ($exons[0]->strand == 1) {
    @exons = sort {$a->start <=> $b->start} @exons;
    $end = $exons[$#exons]->end;
  } else {
    @exons = sort {$b->start <=> $a->start} @exons; # they're read in opposite direction (from right to left)
    $end = $exons[$#exons]->start;                  # the end is the start coordinate of the left-most exon  
  }                                                 # which is here the last of the list @exons
  return $end;
}
#########################################################################


1;
