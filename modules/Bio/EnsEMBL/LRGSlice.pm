=head1 LICENSE

  Copyright (c) 1999-2010 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::LRGSlice - Arbitary Slice of a genome

=head1 SYNOPSIS

  $sa = $db->get_SliceAdaptor;

  $slice =
    $sa->fetch_by_region( 'LRG', 'LRG3');

  # get some attributes of the slice
  my $seqname = $slice->seq_region_name();
  my $start   = $slice->start();
  my $end     = $slice->end();

  # get the sequence from the slice
  my $seq = $slice->seq();

  # get some features from the slice
  foreach $gene ( @{ $slice->get_all_Genes } ) {
    # do something with a gene
  }

  foreach my $feature ( @{ $slice->get_all_DnaAlignFeatures } ) {
    # do something with dna-dna alignments
  }

=head1 DESCRIPTION

A LRG Slice object represents a region of a genome.  It can be used to retrieve
sequence or features from an area of interest.

=head1 METHODS

=cut

package Bio::EnsEMBL::LRGSlice;
use vars qw(@ISA);
use strict;

use Bio::PrimarySeqI;


#Test for LRG tag in cvs

my $reg = "Bio::EnsEMBL::Registry";

use Bio::EnsEMBL::Slice;

sub new{
  my $class = shift;

  my $self = bless {}, $class ;

  my $slice = Bio::EnsEMBL::Slice->new(@_);
#  my $self = $class->SUPER::new( @_);

  my $max=-99999999999;
  my $min=9999999999;
  my $chrom;
  my $strand;

  foreach my $segment (@{$slice->project('chromosome')}) {
    my $from_start = $segment->from_start();
    my $from_end    = $segment->from_end();
    my $to_name    = $segment->to_Slice->seq_region_name();
    $chrom = $to_name;

    my $to_start    = $segment->to_Slice->start();
    my $to_end    = $segment->to_Slice->end();
    if($to_start > $max){
      $max = $to_start;
    }
    if($to_start < $min){
      $min = $to_start;
    }
    if($to_end > $max){
      $max = $to_end;
    }
    if($to_end <  $min){
      $min = $to_end;
    }
    my $ori        = $segment->to_Slice->strand();
    $strand = $ori;   
    
#   print "$from_start-$from_end  => $to_name $to_start-$to_end ($ori) \n";
  }
  if(!defined($chrom)){
    die "Could not project to chromosome for ".$slice->name."??\n";
  }
  my $sa = $slice->adaptor;
#  print "creating chrom slice from $min to $max\n";
  my $chrom_slice = $sa->fetch_by_region("chromosome",$chrom, $min, $max, $strand);

#  print "chrom slcie start = ".$chrom_slice->start."  end = ".$chrom_slice->end."\n";

#  print $chrom_slice."\n";
  $self->{'_orig_slice'} = $slice;
  $self->{'_chrom_slice'} = $chrom_slice;


#  print "CHROM : ".$chrom_slice->seq_region_name."\t".$chrom_slice->start."\t".$chrom_slice->end."\n";
#  print "LRG   : ".$slice->seq_region_name."\t".$slice->start."\t".$slice->end."\n";


  my $asma = "Bio::EnsEMBL::Registry"->get_adaptor($sa->db->species,"core","assemblymapper");
  my $csa = "Bio::EnsEMBL::Registry"->get_adaptor($sa->db->species,"core","coordsystem");
  
  
  
#  my $cs1 = $csa->fetch_by_name("Chromosome","GRCh37");
#  my $cs1 = $csa->fetch_by_name("Chromosome","NCBI36");
  my $cs1 = $chrom_slice->coord_system;
  my $cs2 = $slice->coord_system;
  
  
  my $asm = $asma->fetch_by_CoordSystems($cs1,$cs2);


#  print "mapper to be used for lrg is ".ref($asm)."\n";
 $self->{'_asm'} = $asm;

 return $self;
}

use vars '$AUTOLOAD';


sub AUTOLOAD {
  my $self = shift;

  my $method = $AUTOLOAD;
  $method =~ s/.*:://;


  if($method =~ /^get_all_Attribute/){
    print STDERR "get_all_Attribbutes called\n";
    return  $self->{'_orig_slice'}->$method(@_);    
  }
  elsif($method =~ /^get_all_/ ){
    my $features = $self->{'_chrom_slice'}->$method(@_);
    my @new_features;
    foreach my $ft (@{$features}){
      if($ft->start > $ft->end){
	my $temp = $ft->start;
        $ft->start($ft->end);
        $ft->end($temp);
      }	
      if(($ft->start+$ft->slice->start) > $self->{'_chrom_slice'}->end or ($ft->end+$ft->slice->start) < $self->{'_chrom_slice'}->start){
	print STDERR "start before orig start???\n";
	next;
      }
      print STDERR "FT: ".$ft->dbID."\t(".$ft->start.") ".($ft->start+$ft->slice->start)."\t(".$ft->end.") ".($ft->end+$ft->slice->start)."   ".$ft->slice->seq_region_name."\n";
      my $new_ft = $ft->transfer($self->{'_orig_slice'});  
      if(defined($new_ft)){
#	print "NEW FT: ".$new_ft."\t".($new_ft->start+$new_ft->slice->start)."\t".($new_ft->end+$new_ft->slice->start)."\n";
	push @new_features, $new_ft;
      }
      else{
      # DO i want to give a message here or just ignore them???
	print STDERR "problem transfering $ft start =".($ft->start+$ft->slice->start)." end= ".($ft->end+$ft->slice->end)."\n";
      }
   }
    
    return \@new_features;
  }
#  print "CAlling $method on lrg slice\n";
  return  $self->{'_orig_slice'}->$method(@_);
}

sub DESTROY{
}


1;
