=head1 NAME

GeneCluster

=head1 SYNOPSIS


=head1 DESCRIPTION

This object holds one or more genes which has been clustered according to 
comparison criteria external to this class (for instance, in the 
methods compare and _compare_Genes methods of the class GeneComparison).
Each GeneCluster object holds the IDs of the genes clustered and the beginning and end coordinates
of each one (taken from the start and end coordinates of the first and last exon in the correspondig
get_all_Exons array)

=head1 CONTACT

eae@sanger.ac.uk

=cut

# Let the code begin ...

package Bio::EnsEMBL::Utils::GeneCluster;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Root;

@ISA = qw(Bio::EnsEMBL::Root);

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

#########################################################################

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

=head2 pair_Transcripts()

  Title   : pair_Transcripts()
  Usage   : my @array_of_pairs = $gene_cluster->pair_Transcripts
  Function: This method make pairs of transcripts according to maximum reciprocal exon overlap. 
            Transcripts can eventually be repeated because the maximum exon_overlap may coincide 
            in two different pairs. 
  Example : look for instance in the method find_missing_Exons
  Returns : three arrayrefs =  
            1.- a list of Bio::EnsEMBL::Utils::Transcripts, each holding a pair of transcripts, 
            2.- a list with the unpaired transcripts, and 
            3.- a list those transcripts which have been paired up twice or more
  Args    : nothing

=cut
  
sub pair_Transcripts {
  my ($self) = @_;
  
  # get the genes separated by type list 'benchmark'-like or 'prediction'-like
  my ( $genes2, $genes1 ) = $self->get_separated_Genes;
  my (@transcripts1,@transcripts2);
  my (@trans1,@trans2);
  foreach my $gene ( @$genes1 ){
    push( @trans1, $gene->each_Transcript );
  }
  foreach my $gene ( @$genes2 ){
    push( @trans2, $gene->each_Transcript );
  }
  
  # first sort the transcripts by their start position coordinate
  my %start_table1;
  my %start_table2;
  my $i=0;
  foreach my $tran ( @trans1 ) {
    $start_table1{$i} = $tran->start_exon->start;
    $i++;
  }
  my $j=0;
  foreach my $tra ( @trans2  ) {
    $start_table2{$j} = $tra->start_exon->start;
    $j++;
  }
  foreach my $pos ( sort { $start_table1{$a} <=> $start_table1{$b} } keys %start_table1 ){
    push (@transcripts1, $trans1[$pos]);
  }
  foreach my $pos ( sort { $start_table2{$a} <=> $start_table2{$b} } keys %start_table2 ){
    push (@transcripts2, $trans2[$pos]);
  }

  # pair the transcripts, but first, some variable definition

  my %seen1;           # these keep track of those transcript linked and with how much overlap
  my %seen2;           # ditto, for @transcripts2
  my @pairs;           # list of (Bio::EnsEMBL::Utils::TranscriptCluster) transcript-pairs being created 
  my @unpaired;        # list of Bio::EnsEMBL::Transcript which are left unpaired
  my @doubled;         # those which have been paired up twice
  my $overlap_matrix;  # matrix holding the number of exon overaps for each pair of transcripts
  my $link;            # matrix with 1 for the pairs linked and undef otherwise 
  my @overlap_pairs;   # each entry holds an array with the overlap and the two transcripts being compared
  my %repeated;        # to keep track of repeated transcripts

  # first calculate all possible overlaps
  foreach my $tran1 ( @transcripts1 ){
    foreach my $tran2 ( @transcripts2 ){
      $$overlap_matrix{ $tran1 }{ $tran2 } = _compare_Transcripts( $tran1, $tran2 );
      my @list = ( $$overlap_matrix{ $tran1 }{ $tran2 }, $tran1, $tran2 );
      push ( @overlap_pairs, \@list );
      #print STDERR "Overlap( ".$tran1->stable_id.",".$tran2->stable_id." ) = "
      #.$$overlap_matrix{ $tran1 }{ $tran2 }."\n";
    }
  }
  # sort the list of @overlap_pairs on the overlap
  my @sorted_pairs = sort { $$b[0] <=> $$a[0] } @overlap_pairs;
 
  # take the first pair of the list
  my $first = shift @sorted_pairs;
  my ($max_overlap,$tran1,$tran2) =  @$first;
  $seen1{ $tran1 } = $max_overlap;
  $seen2{ $tran2 } = $max_overlap;
  $$link{ $tran1 }{ $tran2 } = 1;
  
  # scan through each overlap
 PAIR:
  foreach my $list ( @sorted_pairs ){
    # each list contains @$list = ( overlap, transcript1, transcript2 )
    
    # first of all, if the overlap is zero, ditch it
    if ( $$list[0] == 0 ){
      next PAIR;
    }

    # if we've seen both transcripts already reject them
    if ( $$link{ $$list[1] }{ $$list[2] } && defined( $seen1{ $$list[1] } ) && defined( $seen2{ $$list[2] } ) ){
      next PAIR;
    }

    # if the same score...
    if ( $$list[0] == $max_overlap ) {
      
      # if we've seen both transcripts already, check they have the highest score
      if ( defined( $seen1{ $$list[1] } ) && defined( $seen2{ $$list[2] } ) ){
	if ( $$list[0] == $seen1{ $$list[1] } && $$list[0] == $seen2{ $$list[2] } ){
	  $$link{ $$list[1] }{ $$list[2] } = 1;
	}
	 next PAIR;
      }
	
      # if the pair is entirely new, we accept it
      if ( !defined( $seen1{ $$list[1] } ) && !defined( $seen2{ $$list[2] } ) ){
	$$link{ $$list[1] }{ $$list[2] } = 1;
	$seen1{ $$list[1] } = $$list[0];
	$seen2{ $$list[2] } = $$list[0];
	next PAIR;
      }

      # we accept repeats only if this is their maximum overlap as well
      if ( !defined( $seen2{$$list[2]} ) && defined( $seen1{$$list[1]} ) && $$list[0] == $seen1{$$list[1]} ){
	$$link{ $$list[1] }{ $$list[2] } = 1;
	$seen2{ $$list[2] } = $$list[0];
	if ( !defined( $repeated{ $$list[1] } ) ){
	  push( @doubled, $$list[1] );
	  $repeated{ $$list[1] } = 1;
	}
	next PAIR;
      }
      if ( !defined( $seen1{$$list[1]} ) && defined( $seen2{$$list[2]} ) && ($$list[0] == $seen2{$$list[2]}) ){ 
	$$link{ $$list[1] }{ $$list[2] } = 1;
	$seen1{ $$list[1] } = $$list[0];
	if ( !defined( $repeated{ $$list[2] } ) ){
	  push( @doubled, $$list[2] );
	  $repeated{ $$list[2] } = 1;
	}
	next PAIR;
      }
    }

    # if the score is lower, only accept if the pair is completely new
    if ( $$list[0] < $max_overlap ){
      if ( !defined( $seen1{ $$list[1] } ) && !defined( $seen2{ $$list[2] } ) ){
	$$link{ $$list[1] }{ $$list[2] } = 1;
	$seen1{ $$list[1] } = $$list[0];
	$seen2{ $$list[2] } = $$list[0];
	$max_overlap = $$list[0];
	next PAIR;
      }
    }
  }
  
  # create a new cluster for each pair linked
  foreach my $tran1 ( @transcripts1 ){
    foreach my $tran2 ( @transcripts2 ){
      if ( $$link{ $tran1 }{ $tran2} && $$link{ $tran1 }{ $tran2 } == 1 ){
	my $pair = Bio::EnsEMBL::Utils::TranscriptCluster->new();
	$pair->put_Transcripts( $tran1, $tran2 );
	push( @pairs, $pair );
      }
    }
  }

  # finally, check for the unseen ones
  foreach my $tran1 ( @transcripts1 ){
    if ( !defined( $seen1{ $tran1 } ) ){
      push( @unpaired, $tran1 );
    }
  }
  foreach my $tran2 ( @transcripts2 ){
    if ( !defined( $seen2{ $tran2 } ) ){
      push( @unpaired, $tran2 );
    }
  }
  print STDERR scalar(@pairs)." transcript pairs created\n";
  #my $count2=1;
  #foreach my $pair ( @pairs ){
  #  print STDERR "pair $count2:\n".$pair->to_String;
  #  $count2++;
  #}
  #$count2=1;
  #print STDERR scalar(@unpaired)." unpaired transcripts\n";
  #foreach my $unpaired ( @unpaired ){
  #  print STDERR "unpaired $count2: ".$unpaired->stable_id."\n";
  #}
  # return the pairs, the unpaired transcripts, and those transcript which have been taken twice or more
  return (\@pairs,\@unpaired,\@doubled);
}

#########################################################################
   

=head2 _compare_Transcripts()

 Title: _compare_Transcripts()
 Usage: this internal function compares the exons of two transcripts according to overlap
        and returns the number of overlaps
=cut

sub _compare_Transcripts {         
  my ($transcript1,$transcript2) = @_;
  my @exons1   = $transcript1->get_all_Exons;
  my @exons2   = $transcript2->get_all_Exons;
  my $overlaps = 0;
  
  foreach my $exon1 (@exons1){
    
    foreach my $exon2 (@exons2){

      if ( ($exon1->overlaps($exon2)) && ($exon1->strand == $exon2->strand) ){
	$overlaps++;
      }
    }
  }
  return $overlaps;  # we keep track of the number of overlaps to be able to choose the best match
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
    my @exons = $gene->get_all_Exons;
     
    $data .= sprintf "Id: %-16s"      , $gene->stable_id;
    $data .= sprintf "Contig: %-20s"  , $exons[0]->contig->id;
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
  my @exons = $gene->get_all_Exons;
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
  my @exons = $gene-get_all_Exons;
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
		$statistics{$gene->stable_id} = [@gene_stats];
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
