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

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Gene;
use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);

=head1 METHODS

=cut

#########################################################################


=head2 new()

new() initializes the attributes:

$self->{'_benchmark_types'}
$self->{'_prediction_types'}
$self->{'_benchmark_genes'}
$self->{'_prediction_genes'}

=cut

sub new {
  my ($class,$whatever)=@_;

  if (ref($class)){
    $class = ref($class);
  }
  my $self = {};
  bless($self,$class);

if ($whatever){
    $self->throw( "Can't pass an object to new() method. Use put_Genes() to include Bio::EnsEMBL::Gene in cluster");
  }

  $self->{'_ys_gene'}= {};
  %{$self->{'_statistics'}}=();
  
  return $self;
}

#########################################################################

=head2 put_Genes()

  function to include one or more genes in the cluster.
  Useful when creating a cluster. It takes as argument an array of genes, it returns nothing.

=cut

sub put_Genes {
  my ($self, @new_genes)= @_;
  if ( !defined( $self->{'_benchmark_types'} ) || !defined(  $self->{'_prediction_types'} ) ){
    $self->throw( "Cluster lacks references to gene-types, unable to put the gene");
  }

 GENE:
  foreach my $gene (@new_genes){
    foreach my $type ( @{ $self->{'_benchmark_types'} } ){
      if ($gene->type eq $type){
	push ( @{ $self->{'_benchmark_genes'} }, $gene );
	next GENE; 
      }
    }
    foreach my $type ( @{ $self->{'_prediction_types'} } ){
      if ($gene->type eq $type){
	push ( @{ $self->{'_prediction_genes'} }, $gene );
	next GENE;
      }
    }
  }
}

#########################################################################

=head2 get_Genes()

  it returns the array of genes in the GeneCluster object

=cut

sub get_Genes {
  my $self = shift @_;
  my @genes;
  if ( !defined( $self->{'_benchmark_genes'} ) && !defined( $self->{'_prediction_genes'} ) ){
    $self->warn("The gene array you try to retrieve is empty");
    @genes = ();
  }
  if ( $self->{'_benchmark_genes'} ){
    push( @genes, @{ $self->{'_benchmark_genes'} } );
  }
  if ( $self->{'_prediction_genes'} ){
    push( @genes, @{ $self->{'_prediction_genes'} } );
  }
  return @genes;
}


=head2 get_separated_Genes()

  Handy method to get the genes in the genes in the cluster separated by type.
  It returns two arrayrefs.

=cut


sub get_separated_Genes {
  my ($self) = @_;
  return ( $self->{'_benchmark_genes'}, $self->{'_prediction_genes'} );
}

#########################################################################

=head2 get_Gene_Count()

  it returns the number of genes in the GeneCluster object

=cut

sub get_Gene_Count {
  my $self = shift @_;
  my $count =0;
  if ( $self->{'_benchmark_genes'} ){
    $count += scalar( @{ $self->{'_benchmark_genes'}  } ); 
  }
  if ( $self->{'_prediction_genes'} ){
    $count += scalar( @{ $self->{'_prediction_genes'} } );
  }
  #print STDERR "In GeneCluster.get_Gene_Count(), Count = ".$count."\n";
  return $count;
}

#########################################################################

=head2 gene_Types()
  
  It accepts two array references to set the types. One array holds the gene-types for the 
  benchmark genes and the other on for the predicted genes. 
  It can also be used to get the two type-arrays: ($types1, $types2) = $cluster->gene_Types;
  The conventions throughout are (first entry: benchmark, second entry: prediction)

=cut

sub gene_Types {
  my ($self, $benchmark_types, $prediction_types) = @_;
  if ( $benchmark_types && $prediction_types ) {
    $self->{'_benchmark_types'}  = $benchmark_types;
    $self->{'_prediction_types'} = $prediction_types;
  }
  return ($self->{'_benchmark_types'},$self->{'_prediction_types'});
}

#########################################################################

=head2 get_Genes_by_Type()

  We can get the genes in each cluster of a given type. 
  We pass an arrayref containing the types we want to retrieve.

=cut
  
  sub get_Genes_by_Type() {
    my ($self,$types) = @_;
    unless ($types){
      $self->throw( "must provide a type");
    }
    my @genes = $self->get_Genes;  # this should give them in order, but we check anyway
    my @selected_genes;
    foreach my $type ( @{ $types } ){
      push ( @selected_genes, grep { $_->type eq $type } @genes );
    }
    return @selected_genes;
  }
    


#########################################################################

=head2 get_first_Gene()

  it returns the first gene in the cluster, which usually would be the benchmark gene 

=cut

sub get_first_Gene {
  my $self = shift @_;
  return @{$self->{'_benchmark_genes'}}[0];
}

#########################################################################

=head2 to_String()

  it returns a string containing the information about the genes contained in the
  GeneCluster object

=cut

sub to_String {
  my $self = shift @_;
  my $data='';
  foreach my $gene ( $self->get_Genes ){
    my @exons = $gene->each_unique_Exon;
     
    $data .= sprintf "Id: %-16s"      , $gene->id;
    $data .= sprintf "Contig: %-20s"  , $exons[0]->contig_id;
    $data .= sprintf "Exons: %-3d"    , scalar(@exons);
    $data .= sprintf "Start: %-9d"    , $self->_get_start($gene);
    $data .= sprintf "End: %-9d"      , $self->_get_end  ($gene);
    $data .= sprintf "Strand: %-2d\n" , $exons[0]->strand;
  }
  return $data;
}

#########################################################################

=head2 _get_start()

 function to get the start position of a gene - it reads the gene object and it returns
 the start position of the first exon

=cut

sub _get_start {
  my ($self,$gene) = @_;
  my @exons = $gene->each_unique_Exon;
  my $st;
  
  if ($exons[0]->strand == 1) {
    @exons = sort {$a->start <=> $b->start} @exons;
    $st = $exons[0]->start;
  } else {
    @exons = sort {$b->start <=> $a->start} @exons; # they're read in opposite direction (from right to left)
    $st = $exons[0]->end;                           # the start is the end coordinate of the right-most exon
  }                                                 # which is here the first of the list of sorted @exons
  return $st;
}

#########################################################################

=head2 _get_end()

 function to get the end position of a gene - it reads the gene object and it returns
 the end position of the last exon

=cut

sub _get_end {
  my ($self,$gene) = @_;
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
#Adding new methods to calculate the prediction accuracies of the genes in this cluster
#########################################################################

=head2 nucleotide_level_accuracy()

 function that calculates the difference between the annotated and predicted genes at a nucleotide level
 returns the average sensitivity and specificity of the predictions

=cut

sub nucleotide_level_accuracy {
    my ($self) = @_;	
	my @genes = $self->get_Genes;
    my %statistics;	
	shift @genes; # the first gene in the array should be the yardstick gene

	GENE: foreach my $gene (@genes){
	    my $count =0;	
		my ($sum_sn,$sum_sp) = (0,0);
		my ($gene_sn,$gene_sp) ;
 		TRANS: foreach my $trans ($gene->each_Transcript) {
			
	    	my @results=$self->evaluate_Transcripts($trans);
	   		next TRANS unless @results; 
			$count ++;	#provides a running counter of transcripts which overlap.
			my ($trans_sn, $trans_sp) = @results;
			$sum_sn += $trans_sn;
			$sum_sp += $trans_sp;
			
	       
	    }
		$gene_sn = $sum_sn/$count;	
		$gene_sp = $sum_sp/$count;	
		my @gene_stats = ($gene_sn, $gene_sp);	
		$statistics{$gene->id} = [@gene_stats];
	}
	$self->statistics(%statistics);
}


#########################################################################

=head2 statistics()

 returns a hash containing the statistics of each gene in this cluster

=cut

sub statistics {
  my ($self,%stats) = @_;
  if (%stats){
    %{$self->{'_statistics'}}=%stats;
  }
  return  %{$self->{'_statistics'}};
}


#########################################################################

=head2 evaluate_Transcripts()

 function that compares a transcript with each transcript of the annotated gene in this cluster

=cut

sub evaluate_Transcripts {

    my ($self,$trans) = @_;

    my $yardstick = $self->get_first_Gene();

	my $count=0;
	my ($sum_sn, $sum_sp,$sn,$sp) = (0,0);	# The sums and average sensitivity and specificity of this particular transcript WRT 
				 	 						# all transcripts in the yardstick gene in this cluster.

    TRANS: foreach my $ys_trans($yardstick->each_Transcript) {
      my $trans_tp = 0;
	  my @ys_exons = $ys_trans->translateable_exons;
	
	  foreach my $ys_exon (@ys_exons){
	    
	  my @exons = $trans->translateable_exons;
	    EXON: while (@exons){
	    
	        my $exon = shift @exons;
		
			next EXON unless ($exon->overlaps($ys_exon) && ($exon->strand eq $ys_exon->strand));
       		my $overlap;
			my ($exon_start,$exon_end) =  ($exon->start,$exon->end);
			my ($ys_start,$ys_end) =  ($ys_exon->start,$ys_exon->end);
			my ($start,$end);

			if ($exon_start > $ys_start) { $start = $exon_start;} 
				else {$start = $ys_start;}
			if ($exon_end  < $ys_end) { $end = $exon_end;} 
				else {$end = $ys_end;}

			$overlap = $end - $start;

			$trans_tp += $overlap;
	    }	
		
	  }
    next TRANS unless ($trans_tp ne 0);

	$count ++;	#provides a running counter of transcripts which overlap.
    my $trans_ap = _translateable_exon_length($ys_trans);
    my $trans_pp = _translateable_exon_length($trans);
	
	$sum_sn += sprintf("%.2f",($trans_tp)/($trans_ap)*100);
	$sum_sp += sprintf("%.2f",($trans_tp)/($trans_pp)*100);
    #$sum_sn += $trans_tp/$trans_ap; # sensitivity
	#$sum_sp += $trans_tp/$trans_pp; #specificity 
    }     	
	if ($count eq 0){
	print STDERR "Count eq 0. is there an error?\n";
	return 0 ; 
	}
    $sn = $sum_sn/$count;    
    $sp = $sum_sp/$count;    

	my 	@statistics = (	$sn,$sp);
	return @statistics; 
}

#########################################################################

=head2 _translateable_exon_length()

 internal function that returns the length of the translateable exons 

=cut

sub _translateable_exon_length {
	my ($trans)= @_;
	my @exons = $trans->translateable_exons;
    my $length = 0;
    foreach my $ex (@exons) {
        $length += $ex->length;
    }
    return $length;

}

1;
