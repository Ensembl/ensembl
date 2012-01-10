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

Bio::EnsEMBL::IdMapping::TinyExon - lightweight exon object

=head1 SYNOPSIS

  # fetch an exon from the db and create a lightweight exon object
  # from it
  my $exon = $exon_adaptor->fetch_by_stable_id('ENSE000345437');
  my $lightweight_exon = Bio::EnsEMBL::IdMapping::TinyExon->new_fast( [
      $exon->dbID,
      $exon->stable_id,
      $exon->version,
      $exon->created_date,
      $exon->modified_date,
      $exon->start,
      $exon->end,
      $exon->strand,
      $exon->slice->seq_region_name,
      $exon->slice->coord_system_name,
      $exon->slice->coord_system->version,
      $exon->slice->subseq( $exon->start, $exon->end, $exon->strand ),
      $exon->phase,
      $need_project,
  ] );

=head1 DESCRIPTION

This is a lightweight exon object for the stable Id mapping. See the
documentation in TinyFeature for general considerations about its
design.

=head1 METHODS

  start
  end
  strand
  seq_region_name
  coord_system_name
  coord_system_version
  seq
  phase
  need_project
  common_start
  common_end
  common_strand
  common_sr_name
  length
  is_known

=cut


package Bio::EnsEMBL::IdMapping::TinyExon;

# internal data structure (array indices):
#
#  0-4 see TinyFeature
#  5  start
#  6  end
#  7  strand
#  8  seq_region_name
#  9  coord_system_name
# 10  coord_system_version
# 11  seq
# 12  phase
# 13  need_project
# 14  common_start
# 15  common_end
# 16  common_strand
# 17  common_sr_name


use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::TinyFeature;
our @ISA = qw(Bio::EnsEMBL::IdMapping::TinyFeature);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);


=head2 start

  Arg[1]      : (optional) Int - the exon's start coordinate
  Description : Getter/setter for the exon's start coordinate.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub start {
  my $self = shift;
  $self->[5] = shift if (@_);
  return $self->[5];
}


=head2 end

  Arg[1]      : (optional) Int - the exon's end coordinate
  Description : Getter/setter for the exon's end coordinate.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub end {
  my $self = shift;
  $self->[6] = shift if (@_);
  return $self->[6];
}


=head2 strand

  Arg[1]      : (optional) Int - the exon's strand
  Description : Getter/setter for the exon's strand.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub strand {
  my $self = shift;
  $self->[7] = shift if (@_);
  return $self->[7];
}


=head2 seq_region_name

  Arg[1]      : (optional) String - seq_region name
  Description : Getter/setter for the seq_region name of the slice the exon is
                on.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub seq_region_name {
  my $self = shift;
  $self->[8] = shift if (@_);
  return $self->[8];
}


=head2 coord_system_name

  Arg[1]      : (optional) String - coord_system name
  Description : Getter/setter for the coord_system name of the slice the exon is
                on.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub coord_system_name {
  my $self = shift;
  $self->[9] = shift if (@_);
  return $self->[9];
}


=head2 coord_system_version

  Arg[1]      : (optional) String - coord_system version
  Description : Getter/setter for the coord_system version of the slice the
                exon is on.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub coord_system_version {
  my $self = shift;
  $self->[10] = shift if (@_);
  return $self->[10];
}


=head2 seq

  Arg[1]      : (optional) String - the exon's sequence
  Description : Getter/setter for the exon's sequence.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub seq {
  my $self = shift;
  $self->[11] = shift if (@_);
  return $self->[11];
}


=head2 phase

  Arg[1]      : (optional) Int - the exon's phase
  Description : Getter/setter for the exon's phase.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub phase {
  my $self = shift;
  $self->[12] = shift if (@_);
  return $self->[12];
}


=head2 need_project

  Arg[1]      : (optional) Boolean - attribute to set
  Description : Getter/setter for the attribute determining whether an exon
                needs to be projected onto a common coord_system. You don't need
                to do so if the native coord_system is common to the source and
                target assemblies, or if no common coord_system is found (the
                Cache object has methods to determine this).
  Return type : Boolean
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub need_project {
  my $self = shift;
  $self->[13] = shift if (@_);
  return $self->[13];
}


=head2 common_start

  Arg[1]      : (optional) Int - the exon's start in common coord_system
                coordinates
  Description : Getter/setter for the exon's start in common coord_system
                coordinates. Will return $self->start if no projection to a
                common coord_system is required.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

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


=head2 common_end

  Arg[1]      : (optional) Int - the exon's end in common coord_system
                coordinates
  Description : Getter/setter for the exon's end in common coord_system
                coordinates. Will return $self->end if no projection to a
                common coord_system is required.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

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


=head2 common_strand

  Arg[1]      : (optional) Int - the exon's strand in common coord_system
                coordinates
  Description : Getter/setter for the exon's strand in common coord_system
                coordinates. Will return $self->strand if no projection to a
                common coord_system is required.
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

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


=head2 common_sr_name

  Arg[1]      : (optional) String - seq_region name of the exon's slice on the
                common coord_system
  Description : Getter/setter for the seq_region name of the exon's slice on the
                common coord_system coordinates. Will return
                $self->seq_region_name if no projection to a common coord_system
                is required.
  Return type : String
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

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


=head2 length

  Description : Returns the exon length (distance between start and end).
  Return type : Int
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub length {
  my $self = shift;
  return ($self->end - $self->start + 1);
}


=head2 is_known

  Description : Determine whether an exon is known. In the context of stable Id
                mapping, this is true for all exons.
  Return type : Boolean
  Exceptions  : none
  Caller      : general
  Status      : At Risk
              : under development

=cut

sub is_known {
  return 1;
}


1;

