=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Utils::Slice - Utility functions for slices

=head1 SYNOPSIS

  use Bio::EnsEMBL::Utils::Slice qw(split_Slices);

  # ...

  # get all chromosomes in the database
  my $slices = $slice_adaptor->fetch_all('chromosome');

  # split the chromosomes into equal chunks of size less than 1MB
  # with an overlap of 1kb
  $slices = split_Slices( $slices, 1e6, 1e3 );

=head1 METHODS

=cut


package Bio::EnsEMBL::Utils::Slice;

use strict;
use warnings;

use Exporter;

use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(&split_Slices);

use Bio::EnsEMBL::Utils::Exception qw(throw);
use POSIX;

=head2 split_Slices

  Arg [1]    : ref to list of slices
  Arg [2]    : int maxlength of sub slices
  Arg [3]    : int overlap length (optional)
  Example    : my $sub_slices = split_Slices($slices,$maxlen,$overlap)
  Description: splits a slice into smaller slices 
  Returntype : ref to list of slices
  Exceptions : maxlen <1 or overlap < 0

=cut

sub split_Slices{
  my ($slice_big,$max_length,$overlap)=@_;

  if(!defined($max_length) or $max_length < 1){
    throw("maxlength needs to be set and > 0");
  }

  if(!defined($overlap)){
    $overlap = 0;
  }
  elsif($overlap < 0){
    throw("negative overlaps not allowed");
  }

  my @out=();

  foreach my $slice (@$slice_big){

    my $start = $slice->start;
    my $end;
    my $multiple;
    my $number;
    my $length = $slice->length;

    if($max_length && ($length > $overlap)) {
      #No seq region may be longer than max_length but we want to make
      #them all similar size so that the last one isn't much shorter.
      #Divide the seq_region into the largest equal pieces that are shorter
      #than max_length

      #calculate number of slices to create
      $number = ($length-$overlap) / ($max_length-$overlap);
      $number = ceil($number); #round up to int

      #calculate length of created slices
      $multiple = $length / $number;
      $multiple   = floor($multiple); #round down to int
    } else {
      #just one slice of the whole seq_region
      $number = 1;
      $multiple = $length;
    }

    my $i;
    for(my $i=0; $i < $number; $i++) {
      $end = $start + $multiple + $overlap;

      #any remainder gets added to the last slice of the seq_region
      $end = $slice->end if($i == $number-1);
      push @out, Bio::EnsEMBL::Slice->new
        (-START             => $start,
         -END               => $end,
         -STRAND            => 1,
         -SEQ_REGION_NAME   => $slice->seq_region_name,
         -SEQ_REGION_LENGTH => $slice->seq_region_length,
         -COORD_SYSTEM      => $slice->coord_system,
         -ADAPTOR           => $slice->adaptor);
      $start += $multiple + 1;
    }
  }

  return \@out;
}



1;
