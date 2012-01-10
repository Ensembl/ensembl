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

Bio::EnsEMBL::SNP

=head1 SYNOPSIS

    $snp = new Bio::EnsEMBL::SNP(
      -start   => 10,
      -end     => 10,
      -strand  => 1,
      -source  => 'The SNP Consortium',
      -score   => 99,                     # new meaning
      -status  => 'suspected',            # new
      -alleles => 't|c'                   # new
    );

   # add it to an annotated sequence

  $annseq->add_SeqFeature($feat);

=head1 DESCRIPTION

This class was written because the EnsEMBL::ExternalData::Variation
object is way too slow.  There was simply too much chaining to bioperl
methods many, many layers deep.  This object behaves like a Variation
but has a much faster constructor, and faster accessors for the relevant
methods needed by the web.

=head1 METHODS

=cut


package Bio::EnsEMBL::SNP;
use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::ExternalData::Variation;
use Scalar::Util qw(weaken isweak);

@ISA = qw( Bio::EnsEMBL::ExternalData::Variation );


sub new_fast {
  my $class = shift;
  my $hashref = shift;
  my $self = bless $hashref, $class;
  weaken($self->{adaptor})  if ( ! isweak($self->{adaptor}) );
  return $self;
}

sub dbID {
  my $self = shift;
  
  if(@_) {
    $self->{'dbID'} = shift;
  }

  return $self->{'dbID'};
}

sub position {
  my ($self, $arg) = @_;

  if(defined $arg) {
    $self->{_gsf_start} = $arg;
    $self->{_gsf_end}   = $arg;
  }

  return $self->{_gsf_start};
}

sub start {
  my ($self, $arg) = @_;

  if(defined $arg) {
    $self->{_gsf_start} = $arg;
  }
  
  return $self->{_gsf_start};
}

sub end {
  my ($self, $arg) = @_;

  if(defined $arg) {
    $self->{_gsf_end} = $arg;
  }

  return $self->{_gsf_end};
}


sub source {
   my ($self, $arg) = @_;

   if(defined $arg) {
     $self->{_source} = $arg;
   }

   return $self->{_source};
 }

sub score {
  my ($self, $arg) = @_;

  if(defined $arg) {
    $self->{_gsf_score} = $arg;
  }

  return $self->{_gsf_score};
}


sub source_tag {
  my ($self, $arg) = @_;

  if(defined $arg) {
    $self->{_source_tag} = $arg;
  }

  return $self->{_source_tag};
}

sub source_version {
  my ($self, $arg) = @_;

  if(defined $arg) {
    $self->{_source_version} = $arg;
  }

  return $self->{_source_version};
}


=head2 display_name

  Arg [1]    : none
  Example    : print $snp->display_name();
  Description: This method returns a string that is considered to be
               the 'display' identifier.  For snps this is the
               returns the same thing as the id method.
  Returntype : string
  Exceptions : none
  Caller     : web drawing code

=cut

sub display_name {
  my $self = shift;
  return $self->id();
}



1;
