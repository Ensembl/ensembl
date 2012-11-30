package Bio::EnsEMBL::Pipeline::Production::DensityGenerator;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Pipeline::Base/;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;

## Default run method
## For a given species, generates the required density features in the core database
sub run {
  my ($self) = @_;
  my $species = $self->param('species');
  my $dba        = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');
  my $logic_name = $self->param('logic_name');
  my $helper     = $dba->dbc()->sql_helper();
  my $dfa        = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'DensityFeature');
  my $analysis   = $self->get_analysis($logic_name);
  my $max_run    = $self->param('max_run');

  $self->delete_old_features($dba, $logic_name);
  $self->check_analysis($dba);

  my $density_type = $self->get_density_type($analysis);
  Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'DensityType')->store($density_type);
  my $slices = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'slice')->fetch_all('toplevel');
  my $option = $self->get_option();
  my $total  = $self->get_total($option);
  my $count  = 0;
  my @features;
  my $iteration = 0;

  my @sorted_slices =
	sort({$a->coord_system()->rank() <=> $b->coord_system()->rank() || $b->seq_region_length() <=> $a->seq_region_length()} @$slices);
  while (my $slice = shift @sorted_slices) {
	$iteration++;
	$count += $self->get_density($slice, $option);
	if ($slice->has_karyotype) {
	  my @blocks = $self->generate_blocks($slice);
	  for my $block (@blocks) {
		my $feature = $self->get_density($block, $option);
		my $df = Bio::EnsEMBL::DensityFeature->new(-seq_region    => $slice,
												   -start         => $block->start - $slice->start + 1,
												   -end           => $block->end - $slice->start + 1,
												   -density_type  => $density_type,
												   -density_value => $feature);
		if ($feature > 0) {
		  push(@features, $df);
		}
	  }
	} else {
	  my $feature = $self->get_density($slice, $option);
	  my $df = Bio::EnsEMBL::DensityFeature->new(-seq_region    => $slice,
												 -start         => $slice->start,
												 -end           => $slice->end,
												 -density_type  => $density_type,
												 -density_value => $feature);
	  if ($feature > 0) {
		push(@features, $df);
	  }
	}
	if ($count >= $total || $iteration == $max_run) {
	  last;
	}
  } ## end while (my $slice = shift ...)
  if (scalar(@features) > 0) {
	$dfa->store(@features);
  }
  return;
} ## end sub run

## Creating density type for fixed number of blocks
sub get_density_type {
  my ($self, $analysis) = @_;
  my $density_type = Bio::EnsEMBL::DensityType->new(-analysis        => $analysis,
													-region_features => $self->param('bin_count'),
													-value_type      => $self->param('value_type'));
  return $density_type;
}

## Creating density type for a fixed block size
sub get_density_type_block {
  my ($self, $analysis) = @_;
  my $density_type = Bio::EnsEMBL::DensityType->new(-analysis   => $analysis,
													-block_size => $self->param('block_size'),
													-value_type => $self->param('value_type'));
  return $density_type;
}

## Generates a fixed number of blocks for a given slice
sub generate_blocks {
  my ($self, $slice) = @_;
  my $block_size = $slice->length()/$self->param('bin_count');
  my ($current_end, $current_start, @blocks);

  for (my $i = 0; $i < $self->param('bin_count'); $i++) {
	$current_start = int($block_size*$i + 1);
	$current_end   = int($block_size*($i + 1));
	if ($current_end > $slice->end()) {
	  $current_end = $slice->end();
	}
	push(@blocks, $slice->sub_Slice($current_start, $current_end));
  }
  return @blocks;
}

## Generates blocks of fixed size for a given slice
sub generate_blocks_fixed {
  my ($self, $slice) = @_;
  my $nb_block = $slice->length()/$self->param('block_size');
  my ($current_end, $current_start, @blocks);
  for (my $i = 0; $i < $nb_block; $i++) {
	$current_start = int($self->param('block_size')*$i + 1);
	$current_end   = int($self->param('block_size')*($i + 1));
	if ($current_end > $slice->end()) {
	  $current_end = $slice->end();
	}
	push(@blocks, $slice->sub_Slice($current_start, $current_end));
  }
  return @blocks;
}

## Deletes all entries associated to a given analysis
sub delete_old_features {
  my ($self, $dba, $logic_name) = @_;
  my $helper = $dba->dbc()->sql_helper();
  my $sql    = q{
     DELETE df, dt
     FROM   density_feature df, density_type dt, analysis a, seq_region s, coord_system cs
     WHERE  df.seq_region_id = s.seq_region_id 
     AND    s.coord_system_id = cs.coord_system_id
     AND    cs.species_id = ? AND a.analysis_id = dt.analysis_id
     AND    dt.density_type_id = df.density_type_id
     AND    a.logic_name = ? };
  $helper->execute_update(-SQL    => $sql,
						  -PARAMS => [$dba->species_id(), $logic_name]);
  
  #Cleanup anything which didn't have a feature linked. Happened during production
  my $left_join_sql = <<'SQL';
DELETE    dt 
FROM      density_type dt 
LEFT JOIN density_feature df using (density_type_id) 
JOIN      analysis a using (analysis_id)
WHERE     df.density_type_id IS NULL
AND       a.logic_name =?
SQL
  $helper->execute_update(-SQL => $left_join_sql, -PARAMS => [$logic_name]);
}

## Checks if the analysis already exists in the database
## If yes, update the last update date
## If not, add a new analysis entry
sub check_analysis {
  my ($self, $dba) = @_;
  my $logic_name = $self->param('logic_name');
  my $aa         = Bio::EnsEMBL::Registry->get_adaptor($self->param('species'), 'core', 'analysis');
  my $analysis   = $aa->fetch_by_logic_name($logic_name);
  if (!defined($analysis)) {
	$analysis = $self->get_analysis($logic_name);
	$aa->store($analysis);
  } else {
	my $support = 'Bio::EnsEMBL::Utils::ConversionSupport';
	$analysis->created($support->date());
	$aa->update($analysis);
  }
}

## Creates a new analysis object using the associated information from the production database
sub get_analysis {
  my ($self, $logic_name) = @_;
  my $prod_dba = $self->get_production_DBAdaptor();
  my $helper   = $prod_dba->dbc()->sql_helper();
  my $sql      = q{
     SELECT distinct display_label, description
     FROM analysis_description 
     WHERE is_current = 1 and logic_name = ? };
  my $anals = $helper->execute(-SQL => $sql, -PARAMS => [$logic_name])->[0];
  throw "Could not find analysis $logic_name" unless (defined $anals);
  my ($display_label, $description) = @{$anals};
  my $analysis = new Bio::EnsEMBL::Analysis(-logic_name    => $logic_name,
											-display_label => $display_label,
											-description   => $description,
											-displayable   => 1);
  return $analysis;
}

sub get_biotype_group {
  my ($self, $group) = @_;
  my $helper   = $self->get_production_DBAdaptor()->dbc()->sql_helper();
  my $sql      = q{
     SELECT name
     FROM biotype
     WHERE object_type = 'gene'
     AND is_current = 1
     AND biotype_group = ?
     AND db_type like '%core%' };
  return $helper->execute_simple(-SQL => $sql, -PARAMS => [$group]) || [];
}

# Empty method if no specific option is needed
sub get_option {
  my $self = @_;
}

sub get_density {
  my ($self, $block, $biotypes) = @_;
  my $species = $self->param('species');
  my $ga = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'gene');
  return $ga->count_all_by_Slice($block, $biotypes);
}

sub get_total {
  my ($self, $option) = @_;
  my $species = $self->param('species');
  my $ga = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'gene');
  return $ga->count_all_by_biotype($option);
}

1;
