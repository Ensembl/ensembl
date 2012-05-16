package Bio::EnsEMBL::Pipeline::Flatfile::Base;

use strict;
use warnings;
use base qw/Bio::EnsEMBL::Pipeline::Base/;

sub data_path {
  my ($self) = @_;
  return $self->get_dir($self->param('type'), $self->param('species'));
}

1;
