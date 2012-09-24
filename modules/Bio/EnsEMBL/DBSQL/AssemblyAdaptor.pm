=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::AssemblyAdaptor - Retrieves meta information
related to the assembly, density features/counts per chromosome or if none
provided, all top level seq regions


=head1 SYNOPSIS


=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::AssemblyAdaptor;

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::MetaContainer;
use Bio::EnsEMBL::Attribute;

use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);



=head2 new

  Arg [1]    : Bio::EnsEMBL::DBAdaptor $dbadaptor the adaptor for
               the database this assembly info adaptor is using.
  Example    : my $aia = new Bio::EnsEMBL::AssemblyAdaptor($dbadaptor);
  Description: Creates a new AssemblyAdaptor object
  Returntype : Bio::EnsEMBL::DBSQL::AssemblyAdaptor
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::DBAdaptor
  Status     : Stable

=cut

sub new {
  my($class, $dbadaptor) = @_;

  my $self = $class->SUPER::new($dbadaptor);

  return $self;
}

=head2 fetch_info
    
  Description: Returns a hash containing information about the assembly
               stored in the meta table, such as assembly name, date etc., 
               a reference to array of top level seq_region names and a
               reference to array of all coordinate system versions found
  Returntype : reference to a hash with assembly info key and value pairs
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut


sub fetch_info {
  my $self = shift;

  #fetch assembly information stored in the meta table

  my $meta_container = $self->db()->get_adaptor('MetaContainer');

  my @meta_keys = qw(assembly.name       assembly.date                    genebuild.start_date
                     genebuild.method    genebuild.initial_release_date   genebuild.last_geneset_update);
  my %assembly_info;

  foreach my $meta_key (@meta_keys) {
      my @values = @{ $meta_container->list_value_by_key($meta_key) };     
      if (@values) {
	  $assembly_info{$meta_key} = $values[0];
      }
  }
 
  my $schema_build = $self->db()->_get_schema_build();
  if ($schema_build) {
      $assembly_info{'schema_build'} = $schema_build;
  }
  
 #fetch available coordinate systems

  my $csa = $self->db()->get_adaptor('CoordSystem');
  my %versions;
  foreach my $cs (@{$csa->fetch_all()}) {
    $versions{$cs->version()} = 1;
  }
  my @coord_system_versions = keys %versions;

  $assembly_info{'coord_system_versions'} = \@coord_system_versions;

  #fetch top level seq_region names

  my $sa = $self->db()->get_adaptor('Slice');

  my $slices = $sa->fetch_all('toplevel');
  
  my %unique = map { $_->seq_region_name(), 0 } @{$slices};
  my $names = [sort { $a cmp $b } keys %unique];
  $assembly_info{'top_level_seq_region_names'} = $names;

  return \%assembly_info;
}


=head2 fetch_stats

  Arg [1]    : string $seq_region_name (optional)
               The name of the toplevel seq_region for which statistics should be fetched

  Description: Returns a reference to a hash containing density features/ density related 
               seq_region attributes for a toplevel seq_region provided or if none
               provided - all top level seq regions
  Returntype : hashref
  Exceptions : throw if the toplevel slice with seq_region_name provided does not exist
  Caller     : general
  Status     : Stable

=cut


sub fetch_stats {
  my $self = shift;

  my $seq_region_name = shift;

  my @slices;

  my %assembly_stats;

  my $sa = $self->db()->get_adaptor('Slice');
  
  if ($seq_region_name) {
      my $slice = $sa->fetch_by_region('toplevel',$seq_region_name);
      if (!$slice) {
	  throw("Top level slice $seq_region_name not found");
      }
      push(@slices, $slice);
      $assembly_stats{'seq_region_name'} = $seq_region_name;
  } else {
      @slices = @{$sa->fetch_all('toplevel')};
  }

  my @density_types  = qw(genedensity knowngenedensity snpdensity percentgc);

  my @attrib_types = qw(GeneNo% SNPCount);

  my $aa = $self->db()->get_adaptor('Attribute');

  my $dfa = $self->db()->get_adaptor('DensityFeature');

  #used to calculate the average density value for density types represented as ratios
 
  my %density_ft_count = ();

  foreach my $slice (@slices) {

     $assembly_stats{'Length (bps)'} += $slice->length();

     foreach my $density_type (@density_types) {
	      
	  my $density_features = $dfa->fetch_all_by_Slice($slice,$density_type);
	  
	  foreach my $density_feature (@$density_features) {

	      if ($density_feature->density_type()->value_type() eq 'ratio') {
		  $density_ft_count{$density_feature->density_type()->analysis()->display_label()} += 1;
	      }

	       $assembly_stats{$density_feature->density_type()->analysis()->display_label()} += $density_feature->density_value(); 
	  }
     }

     foreach my $attrib_type (@attrib_types) {

	  my $attribs = $aa->fetch_all_by_Slice($slice,$attrib_type);
	      
	  foreach my $attrib (@$attribs) {
		 $assembly_stats{$attrib->description()} += $attrib->value(); 
	  }
     }
  }

  foreach my $density_analysis (keys %density_ft_count) {

      if ($density_ft_count{$density_analysis} > 1) {
	  $assembly_stats{$density_analysis} /= $density_ft_count{$density_analysis};
	  $assembly_stats{$density_analysis} = sprintf "%.2f", $assembly_stats{$density_analysis}; 
	  $assembly_stats{$density_analysis} .= '%';
      }
  }

  return \%assembly_stats;
}



1;

