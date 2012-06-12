package Bio::EnsEMBL::Pipeline::FASTA::Indexer;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Pipeline::FASTA::Base/;

use File::Copy qw/copy/;
use File::Spec;
use Bio::EnsEMBL::Utils::Exception qw/throw/;

sub decompress {
  my ($self) = @_;
  my $source = $self->param('file');
  my $target_dir = $self->target_dir();
  my ($vol, $dir, $file) = File::Spec->splitpath($source);
  my $target = File::Spec->catdir($target_dir, $file);
  my $gunzipped_target = $target;
  $gunzipped_target =~ s/.gz$//;
  $self->info('Copying from %s to %s', $source, $target);
  copy($source, $target) or throw "Cannot copy $source to $target: $!";
  $self->info('Decompressing %s to %s', $source, $gunzipped_target);
  system("gunzip -f $target") and throw sprintf('Could not gunzip. Exited with code %d', ($? >>8));
  return $gunzipped_target;
}

sub repeat_mask_date {
  my ($self) = @_;
  my $res = $self->get_DBAdaptor()->dbc()->sql_helper()->execute_simple(
    -SQL => <<'SQL',
select max(date_format( created, "%Y%m%d"))
from analysis a join meta m on (a.logic_name = lower(m.meta_value))
where meta_key =?
SQL
    -PARAMS => ['repeat.analysis']
  );
  return $res->[0] if @$res;
  return q{};
}

sub run {
  my ($self) = @_;
  my $decompressed = $self->decompress();
  $self->index_file($decompressed);
  $self->cleanup_DBAdaptor();
  return;
}

sub index_file {
  die "Implement";
}

sub target_dir {
  die "Implement";
}


1;
