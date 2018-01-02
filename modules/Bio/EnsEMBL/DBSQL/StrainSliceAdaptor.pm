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

Bio::EnsEMBL::DBSQL::StrainSliceAdaptor - adaptor/factory for MappedSlices
representing alternative assemblies

=head1 SYNOPSIS

  my $slice =
    $slice_adaptor->fetch_by_region( 'chromosome', 14, 900000, 950000 );

  my $msc = Bio::EnsEMBL::MappedSliceContainer->new(-SLICE => $slice);

  # create a new strain slice adaptor and attach it to the MSC
  my $ssa = Bio::EnsEMBL::DBSQL::StrainSliceAdaptor->new($sa->db);
  $msc->set_StrainSliceAdaptor($ssa);
  
  # now attach strain
  $msc->attach_StrainSlice('Watson');

=head1 DESCRIPTION

NOTE: this code is under development and not fully functional nor tested
yet.  Use only for development.

This adaptor is a factory for creating MappedSlices representing
strains and attaching them to a MappedSliceContainer. A mapper will be created
to map between the reference slice and the common container slice coordinate
system.

=head1 METHODS

  new
  fetch_by_name

=head1 REALTED MODULES

  Bio::EnsEMBL::MappedSlice
  Bio::EnsEMBL::MappedSliceContainer
  Bio::EnsEMBL::AlignStrainSlice
  Bio::EnsEMBL::StrainSlice

=cut

package Bio::EnsEMBL::DBSQL::StrainSliceAdaptor;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use Bio::EnsEMBL::MappedSlice;
use Bio::EnsEMBL::Mapper;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

our @ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 new

  Example     : my $strain_slice_adaptor =
                  Bio::EnsEMBL::DBSQL::StrainSliceAdaptor->new;
  Description : Constructor.
  Return type : Bio::EnsEMBL::DBSQL::StrainSliceAdaptor
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub new {
  my $caller = shift;
  deprecate("new is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::StrainSliceAdaptor::new instead");
  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);
  
  return $self;
}


=head2 fetch_by_name

  Arg[1]      : Bio::EnsEMBL::MappedSliceContainer $container - the container
                  to attach MappedSlices to
  Arg[2]      : String $name - the name of the strain to fetch
  Example     : my ($mapped_slice) = @{ $msc->fetch_by_name('Watson') };
  Description : Creates a MappedSlice representing a version of the container's
                reference slice with variant alleles from the named strain
  Return type : listref of Bio::EnsEMBL::MappedSlice
  Exceptions  : thrown on wrong or missing arguments
  Caller      : general, Bio::EnsEMBL::MappedSliceContainer
  Status      : At Risk
              : under development

=cut

sub fetch_by_name {
  my $self = shift;
  deprecate("fetch_by_name is deprecated and will be removed in e88. Please use Bio::EnsEMBL::Variation::DBSQL::StrainSliceAdaptor::fetch_by_name instead");
  my $container = shift;
  my $name = shift;

  # argueent check
  unless ($container and ref($container) and
          $container->isa('Bio::EnsEMBL::MappedSliceContainer')) {
    throw("Need a MappedSliceContainer.");
  }

  unless ($name) {
    throw("Need a strain name.");
  }

  my $slice = $container->ref_slice;

  # get a connection to the variation DB
  my $variation_db = $self->db->get_db_adaptor('variation');

  unless($variation_db) {
    warning("Variation database must be attached to core database to retrieve variation information");
    return '';
  }
  
  # now get an allele feature adaptor
  my $af_adaptor = $variation_db->get_AlleleFeatureAdaptor;
  
  # check we got it
  unless(defined $af_adaptor) {
    warning("Not possible to retrieve AlleleFeatureAdaptor from variation database");
    return '';
  }
  
  # now get a sample adaptor
  my $sample_adaptor = $variation_db->get_SampleAdaptor;
  
  # check we got it
  unless(defined $sample_adaptor) {
    warning("Not possible to retrieve SampleAdaptor from variation database");
    return '';
  }
  
  # fetch sample object for this strain name
  my $sample = shift @{$sample_adaptor->fetch_all_by_name($name)};
  
  # check we got a result
  unless(defined $sample) {
    warn("Strain ".$name." not found in the database");
    return '';
  }
  
  
  ## MAP STRAIN SLICE TO REF SLICE
  ################################
  
  # create a mapper
  my $mapper = Bio::EnsEMBL::Mapper->new('mapped_slice', 'ref_slice');
  
  # create a mapped_slice object  
  my $mapped_slice = Bio::EnsEMBL::MappedSlice->new(
    -ADAPTOR   => $self,
    -CONTAINER => $container,
    -NAME      => $slice->name."\#strain_$name",
  );
  
  # get the strain slice
  my $strain_slice = $slice->get_by_strain($sample->name);
  
  # get all allele features for this slice and sample
  #my @afs = sort {$a->start() <=> $b->start()} @{$af_adaptor->fetch_all_by_Slice($slice, $sample)};
  
  # get allele features with coverage info
  my $afs = $strain_slice->get_all_AlleleFeatures_Slice(1);
  
  # check we got some data
  #warning("No strain genotype data available for slice ".$slice->name." and strain ".$sample->name) if ! defined $afs[0];
  
  
  my $start_slice = $slice->start;
  my $start_strain = 1;
  my $sr_name = $slice->seq_region_name;
  #my $sr_name = 'ref_slice';
  my ($end_slice, $end_strain, $allele_length);
  
  my $indel_flag = 0;
  my $total_length_diff = 0;
  
  # check for AFs
  if(defined($afs) && scalar @$afs) {
  
    # go through each AF
    foreach my $af(@$afs) {
      
      # find out if it changes the length
      if($af->length_diff != 0) {
        
        $indel_flag = 1;
        $total_length_diff += $af->length_diff;
        
        # get the allele length
        $allele_length = $af->length + $af->length_diff();
        
        $end_slice = $slice->start + $af->start() - 2;
  
        if ($end_slice >= $start_slice){
            $end_strain = $end_slice - $start_slice + $start_strain;
            
            #add the sequence that maps
            $mapper->add_map_coordinates('mapped_slice', $start_strain, $end_strain, 1, $sr_name, $start_slice, $end_slice);
            
            #add the indel
            $mapper->add_indel_coordinates('mapped_slice', $end_strain+1, $end_strain+$allele_length, 1, $sr_name,$end_slice+1,$end_slice + $af->length);
            
            $start_strain = $end_strain + $allele_length + 1;
        }
        
        else{
          
            #add the indel
            $mapper->add_indel_coordinates('mapped_slice', $end_strain+1,$end_strain + $allele_length, 1, $sr_name,$end_slice+1,$end_slice + $af->length);
            
            $start_strain += $allele_length;
        }
        
        $start_slice = $end_slice + $af->length+ 1;
      }
    }
  }
  
  # add the remaining coordinates (or the whole length if no indels found)
  $mapper->add_map_coordinates('mapped_slice', $start_strain, $start_strain + ($slice->end - $start_slice), 1, $sr_name, $start_slice, $slice->end);
  
  # add the slice/mapper pair
  $mapped_slice->add_Slice_Mapper_pair($strain_slice, $mapper);  
  
  
  
  ## MAP REF_SLICE TO CONTAINER SLICE
  ###################################
  
  if($total_length_diff > 0) {
  
    # create a new mapper
    my $new_mapper = Bio::EnsEMBL::Mapper->new('ref_slice', 'container');
    
    # get existing pairs
    my @existing_pairs = $container->mapper->list_pairs('container', 1, $container->container_slice->length, 'container');
    my @new_pairs = $mapper->list_pairs('mapped_slice', 1, $strain_slice->length(), 'mapped_slice');
    
    # we need a list of indels (specifically inserts)
    my @indels;
    
    # go through existing first
    foreach my $pair(@existing_pairs) {
      
      if($pair->from->end - $pair->from->start != $pair->to->end - $pair->to->start) {
        my $indel;
        $indel->{'length_diff'} = ($pair->to->end - $pair->to->start) - ($pair->from->end - $pair->from->start);
        
        # we're only interested in inserts here, not deletions
        next unless $indel->{'length_diff'} > 0;
        
        $indel->{'ref_start'} = $pair->from->start;
        $indel->{'ref_end'} = $pair->from->end;
        $indel->{'length'} = $pair->to->end - $pair->to->start;
        
        push @indels, $indel;
      }
    }
    
    # now new ones
    foreach my $pair(@new_pairs) {
      
      if($pair->from->end - $pair->from->start != $pair->to->end - $pair->to->start) {
        my $indel;
        $indel->{'length_diff'} = ($pair->from->end - $pair->from->start) - ($pair->to->end - $pair->to->start);
        
        # we're only interested in inserts here, not deletions
        next unless $indel->{'length_diff'} > 0;
        
        $indel->{'ref_start'} = $pair->to->start;
        $indel->{'ref_end'} = $pair->to->end;
        $indel->{'length'} = $pair->from->end - $pair->from->start;
        
        push @indels, $indel;
      }
    }
    
    # sort them
    @indels = sort {
      $a->{'ref_start'} <=> $b->{'ref_start'} ||  # by position
      $b->{'length_diff'} <=> $a->{'length_diff'} # then by length diff so we only keep the longest
    } @indels;
    
    # clean them
    my @new_indels = ();
    my $p = $indels[0];
    push @new_indels, $indels[0] if scalar @indels;
    
    for my $i(1..$#indels) {
      my $c = $indels[$i];
      
      if($c->{'ref_start'} != $p->{'ref_start'} && $c->{'ref_end'} != $p->{'ref_end'}) {
        push @new_indels, $c;
        $p = $c;
      }
    }
    
    $start_slice = $slice->start;
    $start_strain = 1;
    $sr_name = $slice->seq_region_name;
    #$sr_name = 'ref_slice';
    
    foreach my $indel(@new_indels) {
      
      $end_slice = $indel->{'ref_end'};
      $end_strain = $start_strain + ($end_slice - $start_slice);
      
      $allele_length = $indel->{'length'} + $indel->{'length_diff'};
      
      $new_mapper->add_map_coordinates($sr_name, $start_slice, $end_slice, 1, 'container', $start_strain, $end_strain);
      
      $new_mapper->add_indel_coordinates($sr_name,$end_slice+1,$end_slice + $indel->{'length'}, 1, 'container', $end_strain+1, $end_strain+$allele_length);
      
      $start_strain = $end_strain + $allele_length + 1;
      $start_slice = $end_slice + $indel->{'length'} + 1;
    }
    
    $new_mapper->add_map_coordinates($sr_name, $start_slice, $slice->end, 1, 'container', $start_strain, $start_strain + ($slice->end - $start_slice));
    
    # replace the mapper with the new mapper
    $container->mapper($new_mapper);
    
    # change the container slice's length according to length diff
    $container->container_slice($container->container_slice->expand(undef, $total_length_diff, 1));
  }
  
  return [$mapped_slice];
}


1;

