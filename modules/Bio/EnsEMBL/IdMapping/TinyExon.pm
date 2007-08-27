package Bio::EnsEMBL::IdMapping::TinyExon;

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS


=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut


# internal data structure (array indices):
#
#  0  dbID
#  1  stable_id
#  2  start
#  3  end
#  4  strand
#  5  seq_region_name
#  6  coord_system_name
#  7  coord_system_version
#  8  seq
#  9  phase
# 10  need_project
# 11  common_start
# 12  common_end
# 13  common_strand
# 14  common_sr_name


use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::TinyFeature;
our @ISA = qw(Bio::EnsEMBL::IdMapping::TinyFeature);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


sub start {
  my $self = shift;
  $self->[5] = shift if (@_);
  return $self->[5];
}


sub end {
  my $self = shift;
  $self->[6] = shift if (@_);
  return $self->[6];
}


sub strand {
  my $self = shift;
  $self->[7] = shift if (@_);
  return $self->[7];
}


sub seq_region_name {
  my $self = shift;
  $self->[8] = shift if (@_);
  return $self->[8];
}


sub coord_system_name {
  my $self = shift;
  $self->[9] = shift if (@_);
  return $self->[9];
}


sub coord_system_version {
  my $self = shift;
  $self->[10] = shift if (@_);
  return $self->[10];
}


sub seq {
  my $self = shift;
  $self->[11] = shift if (@_);
  return $self->[11];
}


sub phase {
  my $self = shift;
  $self->[12] = shift if (@_);
  return $self->[12];
}


sub need_project {
  my $self = shift;
  $self->[13] = shift if (@_);
  return $self->[13];
}


sub common_start {
  my $self = shift;

  # when used as a setter, always set a value
  $self->[14] = shift if (@_);

  # when used as a getter
  if (scalar(@$self) > 14) {
    # return value for common coord_system if available (but avoid
    # autovivification gotcha!)
    return $self->[14];
  } elsif ($self->need_project) {
    # return undef if common value expected but not there (e.g. no projection
    # found
    return undef;
  } else {
    # return native value
    return $self->start;
  }
}


sub common_end {
  my $self = shift;

  # when used as a setter, always set a value
  $self->[15] = shift if (@_);

  # when used as a getter
  if (scalar(@$self) > 14) {
    # return value for common coord_system if available (but avoid
    # autovivification gotcha!)
    return $self->[15];
  } elsif ($self->need_project) {
    # return undef if common value expected but not there (e.g. no projection
    # found
    return undef;
  } else {
    # return native value
    return $self->end;
  }
}


sub common_strand {
  my $self = shift;

  # when used as a setter, always set a value
  $self->[16] = shift if (@_);

  # when used as a getter
  if (scalar(@$self) > 14) {
    # return value for common coord_system if available (but avoid
    # autovivification gotcha!)
    return $self->[16];
  } elsif ($self->need_project) {
    # return undef if common value expected but not there (e.g. no projection
    # found
    return undef;
  } else {
    # return native value
    return $self->strand;
  }
}


sub common_sr_name {
  my $self = shift;

  # when used as a setter, always set a value
  $self->[17] = shift if (@_);

  # when used as a getter
  if (scalar(@$self) > 14) {
    # return value for common coord_system if available (but avoid
    # autovivification gotcha!)
    return $self->[17];
  } elsif ($self->need_project) {
    # return undef if common value expected but not there (e.g. no projection
    # found
    return undef;
  } else {
    # return native value
    return $self->seq_region_name;
  }
}


sub length {
  my $self = shift;
  return ($self->end - $self->start + 1);
}


sub is_known {
  return 1;
}


1;

