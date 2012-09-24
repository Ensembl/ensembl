package Bio::EnsEMBL::Pipeline::FASTA::Base;

use strict;
use warnings;
use base qw/Bio::EnsEMBL::Pipeline::Base/;

use File::Spec;

sub fasta_path {
  my ( $self, @extras ) = @_;
  return $self->get_dir('fasta', $self->param('species'), @extras);
}

sub old_path {
  my ($self, $species) = @_;
  my $base = $self->param('ftp_dir');
  my $prod = $self->production_name($species);
  my $release = $self->param('previous_release');
  my $dir = File::Spec->catdir($base, "release-$release", 'fasta', $prod, 'dna');
}

# Filter a FASTA dump for reference regions only and parse names from the
# headers in the file e.g. NNN:NNN:NNN:1:80:1
sub filter_fasta_for_nonref {
  my ($self, $source_file, $target_fh) = @_;

  my $allowed_regions = $self->allowed_regions();
  
  $self->info('Filtering from %s', $source_file);
  
  open my $src_fh, '-|', "gzip -c -d $source_file" or $self->throw("Cannot decompress $source_file: $!");
  my $transfer = 0;
  while(my $line = <$src_fh>) {
    #HEADER
    if(index($line, '>') == 0) {
      #regex is looking for NNN:NNN:NNN:1:80:1 i.e. the name
      my ($name) = $line =~ />.+\s(.+:.+:.+:\d+:\d+:\d+)/; 
      $transfer = ($allowed_regions->{$name}) ? 1 : 0;
      if($transfer) {
        $self->info('%s was an allowed Slice', $name);
      }
      else {
        $self->info('%s will be skipped; not a reference Slice', $name);
      }
    }
    print $target_fh $line if $transfer;
  }
  close($src_fh);
  
  return;
}

sub allowed_regions {
  my ($self) = @_;
  my $filter_human = 1;
  my @slices = grep { $_->is_reference() } @{$self->get_Slices('core', $filter_human)};
  my %hash = map { $_->name() => 1 } @slices;
  return \%hash;
}

sub has_non_refs {
  my ($self) = @_;
  my $sql = <<'SQL';
select count(*)
from attrib_type at
join seq_region_attrib sra using (attrib_type_id)
join seq_region sr using (seq_region_id)
join coord_system cs using (coord_system_id)
where cs.species_id =?
and at.code =?
SQL
  my $dba = $self->get_DBAdaptor();
  return $dba->dbc()->sql_helper()->execute_single_result(
    -SQL => $sql, -PARAMS => [$dba->species_id(), 'non_ref']);
}


1;
