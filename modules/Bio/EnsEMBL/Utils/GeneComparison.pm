

=head1 NAME

GeneComparison

=head1 SYNOPSIS


=head1 DESCRIPTION

Perl Class for comparison of two sets of genes.
It can read two references to two arrays of genes, e.g. EnsEMBL built genes and human annotated genes,
and it compares them using methods...

  cluster_Genes (ready)
  get_unmatched_Genes (ready)
  get_fragmented_Genes (ready)
  cluster_Trasncripts_by_Gene  (ready)
    (first cluster all genes and then cluster the transcripts within each gene)   
  cluster_Transcripts  (ready)
    (cluster all transcripts without going through gene-clustering)
  get_Exon_Statistics (ready)
  get_Coding_Exon_Statistics (ready)
  

 compare()
 get_exon_overlaps()
 get_3prime_overlaps()
 get_5prime_overlaps()

The object can be created without passing any genes. The genes to be compared can instead be passed to the
appropriate methods (see below), however, this may lose (in the current version) some info about the gene
types involved in the comparison when there are more than three.

Each GeneComparison object can contain data fields specifying the arrays of genes to be compared.
There are also data fields in the form of two arrays and 
associated to the gene types present in those both arrays.

=head1 CONTACT

eae@sanger.ac.uk

=cut

# Let the code begin ...

package GeneComparison;

use strict;
use GeneCluster;
use TranscriptCluster;
use lib '/nfs/acari/eae/ensembl/modules/';

####################################################################################

=head2 new()

the new() method accepts two array references

=cut

sub new {
  
  my ($class,$gene_array1,$gene_array2) = @_;

  if (ref($class)){
    $class = ref($class);
  }
  my $self = {};
  bless($self,$class);

  if ( $gene_array1 && $gene_array2 ){
    $self->{'_gene_array1'}= $gene_array1;
    $self->{'_gene_array2'}= $gene_array2;
  
    # we keep track of the gene types present in both arrays

    my %types1;
    foreach my $gene ( @{ $gene_array1 } ){
      $types1{$gene->type}=1;
    }
    foreach my $k ( keys(%types1) ){            # put in the array type_array1 
      push( @{ $self->{'_type_array1'} }, $k);  # the gene types present in gene_array1
    }

    my %types2;
    foreach my $gene ( @{ $gene_array2 } ){
      $types2{$gene->type}=1;
    }
    foreach my $k ( keys(%types2) ){            # put in the array _type_array2 
      push( @{ $self->{'_type_array2'} }, $k);  # the gene types present in gene_array2
    }
  }
  
  return $self;
}

######################################################################################
## one could make this more useful by holding the number of genes of each type in each GeneCluster
## object. This could be achieved by, after creating all the clusters in cluster_Genes,
## before returning the array of GeneCluster objects, going through each cluster
## and checking howmany of each type there are and holding this into a couple of variables within
## the GeneCluster object which could then be retrieved later on.
#########################################################################

=head2 gene_Types()

this function sets labels to the arrays of genes passed to new().
It can also be used to retrieve the gene types

=cut

sub gene_Types {
   my ($self,$type1,$type2) = @_;
   if ( $type1 && $type2 ){
     @{ $self->{'_type_array1'} }= $type1;
     @{ $self->{'_type_array2'} }= $type2;
   }
   return ($type1,$type2);
}
  
######################################################################################

=head2 cluster_Genes

  This method takes an array of genes and cluster them
  according to their exon overlap. As a default it takes the genes stored in the GeneComparison object 
  as data fields (or attributes) '_gene_array1' and '_gene_array2'. It can also accept instead as argument
  an array of genes to be clustered, but then information about their gene-type is lost. 
  This method returns an array of GeneCluster objects.

=cut

sub cluster_Genes {

  my ($self, @array) = @_;
  my @genes;
  if (@array){
    @genes = @array;
  }
  else{
    @genes = ( @{ $self->{'_gene_array1'} }, @{ $self->{'_gene_array2'} } );
  }

  #### first sort the genes by the left-most position coordinate ####
  my %start_table;
  my $i=0;
  foreach my $gene (@genes){
    $start_table{$i}=_get_start_of_Gene($gene);
    $i++;
  }
  my @sorted_genes=();
  foreach my $k ( sort { $start_table{$a} <=> $start_table{$b} } keys %start_table ){
    push (@sorted_genes, $genes[$k]);
  }
  @genes = @sorted_genes;

  #### cluster the genes ####
  my @clusters; # this will hold an array of GeneCluster objects
  my $found;
  
  foreach my $gene (@genes){
    $found=0;
  LOOP:
    foreach my $cluster (@clusters){
      foreach my $gene_in_cluster ( $cluster->get_Genes ){  # read genes from GeneCluster object
	
	if ( _compare_Genes($gene,$gene_in_cluster) ){	
	  $cluster->put_Genes($gene);                       # put gene into GeneCluster object
	  $found=1;
	  last LOOP;
	}
      }
    }
    
    if ($found==0){
      my $new_cluster=GeneCluster->new($gene);   # create a GeneCluser object
      push(@clusters,$new_cluster);
    }
  }
  return @clusters;
}
#########################################################################

=head2 _get_start_of_Gene()

  function to get the left-most coordinate of the exons of the gene (start position of a gene).
  For genes in the strand 1 it reads the gene object and it returns the start position of the 
  first exon. For genes in the strand -1 it picks the last exon and reads the start position; this
  is not the proper start of the gene but it is the left-most coordinate.

=cut

sub _get_start_of_Gene {
  my $gene = shift @_;
  my @exons = $gene->each_unique_Exon;
  my $st;
  
  if ($exons[0]->strand == 1) {
    @exons = sort {$a->start <=> $b->start} @exons;
    $st = $exons[0]->start;
  } else {
    @exons = sort {$b->start <=> $a->start} @exons;
    $st = $exons[$#exons]->start;
  }
  return $st;
}

####################################################################################

=head2 cluster_Transcripts_by_Gene()

it first clusters all genes using cluster_Genes and then cluster the transcripts within each gene
cluster. If one has already made a gene-clustering, one can pass the array of clusters
to cluster_Transcripts_by_Gene() to avoid repetition of work

See cluster_Transcripts for a different way of clustering transcripts.

=cut

sub cluster_Transcripts_by_Gene {
  my ($self,@array) = @_;
  my @gene_clusters;      # array of GeneCluster objects
  my @transcript_clusters;
  if (@array){
    @gene_clusters = @array;
  }
  else{
    @gene_clusters = $self->cluster_Genes;
  }
  foreach my $cluster (@gene_clusters){
    
    my @transcripts;   
    foreach my $gene ( $cluster->get_Genes ){
      push ( @transcripts, $gene->each_Transcript );
    }
    push ( @transcript_clusters, $self->cluster_Transcripts(@transcripts) );
    # we pass the array of transcripts to be clustered to the method cluster_Transcripts
  }
  return @transcript_clusters;
}

####################################################################################

=head2 cluster_Transcripts()

  This method cluster all the transcripts in the gene arrays passed to the GeneComparison constructor new().
  It also accepts an array of transcripts to be clustered. The clustering is done according to exon
  overlaps, which is implemented in the function _compare_Transcripts.
  This method returns an array of TranscriptCluster objects.

=cut
  
sub cluster_Transcripts {
  my ($self,@array) = @_;
  my @transcripts;
  if (@array){                 # transcripts to be clustered can be either provided as argument
    @transcripts = @array;
  }
  else{                        # or else taken from the gene_array data fields
    my @genes = ( @{ $self->{'_gene_array1'} }, @{ $self->{'_gene_array2'} } );
    foreach my $gene (@genes){
      my @more_transcripts = $gene->each_Transcript;
      push ( @transcripts, @more_transcripts );
    }
  }
  # we do the clustering with the array @transcripts
  
  # first sort the transcripts by their start position coordinate
  my %start_table;
  my $i=0;
  foreach my $transcript (@transcripts){
    $start_table{$i} = $transcript->start_exon->start;
    $i++;
  }
  my @sorted_transcripts=();
  foreach my $pos ( sort { $start_table{$a} <=> $start_table{$b} } keys %start_table ){
    push (@sorted_transcripts, $transcripts[$pos]);
  }
  @transcripts = @sorted_transcripts;

  #### cluster the sorted transcripts
  my @transcript_clusters;
  my $found;
  
  foreach my $transcript (@transcripts){
    $found=0;
  LOOP:
    foreach my $cluster (@transcript_clusters){
      foreach my $transcript_in_cluster ( $cluster->get_Transcripts ){  
	                                  # read transcripts from TranscriptCluster object
	
	if ( _compare_Transcripts($transcript, $transcript_in_cluster) ){	
	  $cluster->put_Transcripts ($transcript);  # put gene into Transcript Cluster object
	  $found=1;
	  last LOOP;
	}
      }
    }
    
    if ($found==0){
      my $new_cluster=TranscriptCluster->new($transcript);   # create a new TranscriptCluser object
      push(@transcript_clusters,$new_cluster);
    }
  }
  return @transcript_clusters;
}

####################################################################################

=head2 get_Exon_Statistics()

  This method produces an histogram with the number of exon pairs per percentage overlap.
  The percentage overlap of two exons, say e1 and e2, is calculated in _compare_Transcripts method as
  100*intersection(e1,e2)/max_length(e1,e2). It takes whole exons, i.e. without chopping out the 
  non-coding part. This method returns a hash containing the number of occurences as values and
  the integer percentage overlap as keys

=cut
  
sub get_Exon_Statistics{

  my ($self,$array1,$array2) = @_;
  my (@transcripts1,@transcripts2);
  my %stats;
  if ($array1 && $array2){         # we can pass two arrays references (to transcripts) as argument
    @transcripts1 = @{$array1};
    @transcripts2 = @{$array2};
    
  }
  else{                            # or else read the transcripts from the gene_array data fields
    my @genes1 = @{ $self->{'_gene_array1'} };
    foreach my $gene (@genes1){
      my @more_transcripts = $gene->each_Transcript;
      push ( @transcripts1, @more_transcripts );
    }
    my @genes2 = @{ $self->{'_gene_array2'} };
    foreach my $gene (@genes2){
      my @more_transcripts = $gene->each_Transcript;
      push ( @transcripts2, @more_transcripts );
    }
  }
  foreach my $t1 (@transcripts1){
    foreach my $t2 (@transcripts2){
      if ( _compare_Transcripts($t1,$t2) ){
	my %partial_stats = _exon_Statistics($t1,$t2); # produces the stats for the current two transcripts
	                                               # The hash holds number of occurences as values and
	                                               # integer percentage overlap as keys
	foreach my $k ( keys %partial_stats ) {
	  $stats{$k} += $partial_stats{$k};            # Here it adds up to the overall value
	}
      }
    }
  }
  return %stats;
}

#########################################################################

=head2 _exon_Statistics()

This internal function reads two transcripts and calculates the percentage of overlap 
between their exons = (INTERSECT($exon1,$exon2))/MAX($exon1,$exon2)

=cut

sub _exon_Statistics {
  my ($transcript1,$transcript2) = @_;
  my @exons1 = $transcript1->each_Exon; # transcripts get their exons in order
  my @exons2 = $transcript2->each_Exon;

  my %stats;

  foreach my $exon1 (@exons1){
  
    foreach my $exon2 (@exons2){

      if ( ($exon1->overlaps($exon2)) && ($exon1->strand == $exon2->strand) ){
	
	# calculate the percentage of exon overlap = (INTERSECT($exon1,$exon2))/MAX($exon1,$exon2)
	my $max;
	if ( ($exon1->length) < ($exon2->length) ){
	  $max = $exon2->length;
	}
	else{
	  $max = $exon1->length;
	}
	# compute the overlap extent
	
	my ($s,$e);  # start and end coord of the overlap
	my ($s1,$e1) = ($exon1->start,$exon1->end);
	my ($s2,$e2)= ($exon2->start,$exon2->end);
	if ($s1<=$s2 && $s2<$e1){
	  $s=$s2;
	}
	if ($s1>=$s2 && $e2>$s1){
	  $s=$s1;
	}
	if ($e1<=$e2 && $e1>$s2){
	  $e=$e1;
	}
	if ($e2<=$e1 && $e2>$s1){
	  $e=$e2;
	}
	my $common = ($e - $s + 1);
	my $percent = int( (100*$common)/$max );
	$stats{$percent}++;	
      }
    }
  }
  return %stats;
}    

####################################################################################

=head2 get_Coding_Exon_Statistics()

  The same aim as the get_Exon_Statistics method but chopping the non-coding part from the exons.
  This method returns a hash containing the number of occurences as values and
  the integer percentage overlap as keys

=cut
  
sub get_Coding_Exon_Statistics{

  my ($self,$array1,$array2) = @_;
  my (@transcripts1,@transcripts2);
  my %stats;
  if ($array1 && $array2){       # we can pass two arrays references of transcripts as argument
    @transcripts1 = @{$array1};
    @transcripts2 = @{$array2};
    
  }
  else{                          # or else read the transcripts from the gene_array data fields
    my @genes1 = @{ $self->{'_gene_array1'} };
    foreach my $gene (@genes1){
      my @more_transcripts = $gene->each_Transcript;
      push ( @transcripts1, @more_transcripts );
    }
    my @genes2 = @{ $self->{'_gene_array2'} };
    foreach my $gene (@genes2){
      my @more_transcripts = $gene->each_Transcript;
      push ( @transcripts2, @more_transcripts );
    }
  }
  # in order to get info about the non-coding regions of the transcripts we have
  # to get translation objects for them, however, some of them do not have it defined.
  # In those cases we just ignore the comparison.
  
  foreach my $t1 (@transcripts1){
    foreach my $t2 (@transcripts2){
      if ( _compare_Transcripts($t1,$t2) ){
	if ($t1->translation && $t2->translation){
	  my %partial_stats = _coding_Exon_Statistics($t1,$t2);
	  foreach my $k ( keys %partial_stats ) {
	    $stats{$k} += $partial_stats{$k};
	  }
	}
	else{
	  print "Transcript without translation:\n";
	  if (!$t1->translation){
	    print $t1->id."\n";
	  }
	  if (!$t2->translation){
	    print $t2->id."\n";
	  }
	  print "\n";
	}
      }
    }
  }
  return %stats;
}

####################################################################################  

=head2 _coding_Exon_Statistics()

This internal function reads two transcripts and calculates the percentage 
of overlap between the coding exons = (INTERSECT($exon1,$exon2))/MAX($exon1,$exon2)

=cut

sub _coding_Exon_Statistics {         

  # the coding region may start in any exon, not necessarily the first one
  # and may end in any one as well

  my ($transcript1,$transcript2) = @_;

  my @exons1 = $transcript1->each_Exon;
  my @exons2 = $transcript2->each_Exon;

  my $translation1 = $transcript1->translation;
  my $translation2 = $transcript2->translation;

  # IDs of the exons where the coding region starts and ends
  my ($s_id1,$e_id1) = ($translation1->start_exon_id,$translation1->end_exon_id);
  my ($s_id2,$e_id2) = ($translation2->start_exon_id,$translation2->end_exon_id);

  # identify those exons in each transcript
  my ($s_exon1,$e_exon1);  # these will be the exons where the coding region (starts,ends) in the 1st transcript
  foreach my $exon1 (@exons1){
    if ($exon1->id eq $s_id1){
      $s_exon1 = $exon1;
    }
    if ($exon1->id eq $e_id1){
      $e_exon1 = $exon1;
    }
  }
  my ($s_exon2,$e_exon2);  # these will be the exons where the coding region (starts,ends) in the 2nd transcript
  foreach my $exon2 (@exons2){
    if ($exon2->id eq $s_id2){
      $s_exon2 = $exon2;
    }
    if ($exon2->id eq $e_id2){
      $e_exon2 = $exon2;
    }
  }
  # take these exons and those in between these exons (since only these exons contain coding region)
  my @coding_exons1;
  foreach my $exon1 (@exons1){
    if ( $exon1->start >= $s_exon1->start && $exon1->start <= $e_exon1->start ){
      push ( @coding_exons1, $exon1 );
    }
  }
  my @coding_exons2;
  foreach my $exon2 (@exons2){
    if ( $exon2->start >= $s_exon2->start && $exon2->start <= $e_exon2->start ){
      push ( @coding_exons2, $exon2 );
    }
  }
  @exons1 = @coding_exons1;
  @exons2 = @coding_exons2;
    
  # start and end of coding regions 
  # relative to the origin of the first and last exons where the coding region starts
  my ($s_code1,$e_code1) = ($translation1->start,$translation1->end);
  my ($s_code2,$e_code2) = ($translation2->start,$translation2->end);

  # here I hold my stats results
  my %stats;
  
  foreach my $exon1 (@exons1){

    foreach my $exon2 (@exons2){

      if ( ($exon1->overlaps($exon2)) && ($exon1->strand == $exon2->strand) ){

	# start and end coord of the overlap
	my ($s,$e)=(0,0);  

	# exon positions
	my ($s1,$e1) = ($exon1->start,$exon1->end);
	my ($s2,$e2)= ($exon2->start,$exon2->end);
	
	#exon lengths
	my ($l1,$l2) = ($exon1->length, $exon2->length);

	# One has to chop off the non-coding part out of the first and last exons in 
	# the coding region : ($s_id1,$e_id1) and ($s_id2,$e_id2)
	# Recall that ($s_code1,$e_code1) and ($s_code2,$e_code2) are the start/end of the coding regionS
	
	my ($fs1,$fe1,$fs2,$fe2)=(0,0,0,0); # put a flag on the first and last exons to print them out

	if ($exon1->id eq $s_id1){
	  $s1 = $s1 + $s_code1 - 1;
	  $l1 = $e1 - $s1 + 1;
	  $fs1=1;
	}
	if ($exon1->id eq $e_id1){
	  $e1 = $s1 + $e_code1 - 1;
	  $l1 = $e1 - $s1 + 1;
	  $fe1=1;
	}
	if ($exon2->id eq $s_id2){
	  $s2 = $s2 + $s_code2 - 1;
	  $l2 = $e2 - $s2 + 1;
	  $fs2=1;
	}
	if ($exon2->id eq $e_id2){
	  $e2 = $s2 + $e_code2 - 1;
	  $l2 = $e2 - $s2 + 1;
	  $fe2=1;
	}

	# calculate the percentage of exon overlap = (INTERSECT($exon1,$exon2))/MAX($exon1,$exon2)
	my $max;
	if ( $l1 < $l2 ){
	  $max = $l2;
	}
	else{
	  $max = $l1;
	}
	
	# compute the overlap extent
	if ($s1<=$s2 && $s2<$e1){
	  $s=$s2;
	}
	if ($s1>=$s2 && $e2>$s1){
	  $s=$s1;
	}
	if ($e1<=$e2 && $e1>$s2){
	  $e=$e1;
	}
	if ($e2<=$e1 && $e2>$s1){
	  $e=$e2;
	}
	my $common = ($e - $s + 1);
	if ($common != 0){ # we check since after cutting the non-coding piece we might lose the overlap
	  my $percent = int( (100*$common)/$max );
	  $stats{$percent}++;	
	  if ($fs1 || $fs2 ) { 
	    print "(".$s1.",".$e1.")";
	    if ($fs1){
	      print "-> start coding exon";
	    }
	    print "\t".$exon1->id."\n";
	    
	    print "(".$s2.",".$e2.")";
	    if ($fs2) { 
	      print "-> start coding exon";
	    }
	    print "\t".$exon2->id."\n";
	  }
	  if ( $fe1 || $fe2 ) { 
	    print "(".$s1.",".$e1.")";
	    if ($fe1){
	      print "-> end coding exon";
	    }
	    print   "\t".$exon1->id."\n";
	    
	    print "(".$s2.",".$e2.")";
	    if ($fe2) { print "-> end coding exon";
		      }  
	    print "\t".$exon2->id."\n";
	  }
	  if ($fs1 || $fs2 || $fe1 || $fe2){ print "(".$s.",".$e.") Overlap --> ".$percent."\n\n";}
	}
      }
    }
  }
  return %stats;
}    


####################################################################################  
  
=head2 get_unmatched_Genes()

  This function returns those genes that have not been identified with any other gene in the
  two arrays (of type1 and type2) passed to GeneComparison->new(). It derives the unmatched genes
  directly from the gene_arrays passed to the GeneComaprison object. It returns two references to 
  two arrays of GeneCluster objects. The first one corresponds to the genes of type1 which have 
  not been identified in type2, and the second one holds the genes of type2 that have not been 
  identified with any of the genes of type1. The actual Gene objects within it can be retrieved 
  using the method get_Genes on the GeneCluster objects

  In order to recover the type one could as well use the get_Genes method on the array elements
  and then the type method on the Bio::EnsEMBL::Gene objects
  
=cut

sub get_unmatched_Genes {
  my $self = shift @_;
  my %found1;
  my %found2;

  foreach my $gene1 ( @{ $self->{'_gene_array1'} } ){
    foreach my $gene2 ( @{ $self->{'_gene_array2'} } ){
      if ( _compare_Genes($gene1,$gene2)){
	$found1{$gene1->id} = 1;
	$found2{$gene2->id} = 1;
      }
    }
  }
  my @unmatched1;
  foreach my $gene1 ( @{ $self->{'_gene_array1'} } ){
    unless ( $found1{$gene1->id} ){
      my $new_cluster = GeneCluster->new( $gene1 );
      push (@unmatched1, $new_cluster);
    }
  }
  my @unmatched2;
  foreach my $gene2 ( @{ $self->{'_gene_array2'} } ){
    unless ( $found2{$gene2->id} ){
      my $new_cluster = GeneCluster->new( $gene2 );
      push (@unmatched2, $new_cluster);
    }
  }
  return (\@unmatched1,\@unmatched2);
}

#########################################################################

=head2 get_fragmented_Genes()

it returns an array of GeneCluster objects, where these objects contain
fragmented genes: genes of a given type which are overlapped with more than one
gene of the other type. In order to avoid repetition of work, if a gene-clustering
has been already performed, one can pass the array of clusters as an argument.
This, however, is only allowed if _type_array1 and _type_array2 are defined, since the method
makes use of them.

=cut

sub get_fragmented_Genes {
  my ($self,@array) = @_;
  my @clusters;
  my @fragmented;

  if (@array && $self->{'_type_array1'} && $self->{'_type_array2'} ){
    @clusters = @array;
  }
  else{
    @clusters = $self->cluster_Genes;
  }
  
  foreach my $cluster (@clusters){
    
    my @genes = $cluster->get_Genes;
    my (@type1,@type2);
    
    foreach my $gene ( @genes ){
     
      my $type = $gene->type;
      push( @type1, grep /$type/, @{ $self->{'_type_array1'} } );
      push( @type2, grep /$type/, @{ $self->{'_type_array2'} } );
      # @type1 and @type2 hold all the occurrences of gene-types 1 and 2, respectively
      if ( ( @type1 && scalar(@type1)>1 ) || ( @type2 && scalar(@type2) >1 ) ) {
	push (@fragmented, $cluster);
      }
    }
  }
  return @fragmented;
}

#########################################################################

=head2 get_3prime_overlaps()

=cut

#########################################################################

=head2 get_5prime_overlaps()

=cut

#########################################################################


=head2 _compare_Genes()

 Title: _compare_Genes
 Usage: this internal function compares the exons of two genes on overlap and returns the number of overlaps

=cut

sub _compare_Genes {         
  my ($gene1,$gene2) = @_;
  my @exons1 = $gene1->each_unique_Exon;
  my @exons2 = $gene2->each_unique_Exon;
  
  my $overlap_count=0;
  foreach my $exon1 (@exons1){
  
    foreach my $exon2 (@exons2){

      if ( ($exon1->overlaps($exon2)) && ($exon1->strand == $exon2->strand) ){
	$overlap_count++;
      }
    }
  }
  return $overlap_count;
}      
#########################################################################

=head2 _compare_Transcripts()

 Title: _compare_Transcripts()
 Usage: this internal function compares the exons of two transcripts according to overlap
        and returns the number of overlaps

=cut

sub _compare_Transcripts {         
  my ($transcript1,$transcript2) = @_;
  my @exons1 = $transcript1->each_Exon;
  my @exons2 = $transcript2->each_Exon;
  my $overlap_count=0;
  
  foreach my $exon1 (@exons1){
    
    foreach my $exon2 (@exons2){

      if ( ($exon1->overlaps($exon2)) && ($exon1->strand == $exon2->strand) ){
	$overlap_count++;
      }
    }
  }
  return $overlap_count;
}    

#########################################################################

1;

