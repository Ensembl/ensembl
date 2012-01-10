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

Bio::EnsEMBL::SeqRegionSynonym -
Object representing an alternatice name.

=head1 SYNOPSIS

=head1 DESCRIPTION

This object holds information about alternative name to
Ensembl seq regions.

=head1 METHODS

=cut

package Bio::EnsEMBL::SeqRegionSynonym;

use strict;
use warnings;
no warnings qw(uninitialized);

use Bio::EnsEMBL::Storable;
use Bio::Annotation::DBLink;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(deprecate);

our @ISA = qw(Bio::EnsEMBL::Storable);

=head2 new

  Args [...] : list of named parameters 
  Example    : my $srs = new Bio::EnsEMBL::SeqRegionSynonym(
                    -adaptor => $adaptor,
                    -synonym => $alt_name,
                    -external_db_id => 1234
                    -seq_region_id  => 12);
  Description: Creates a new SeqRegionSynonym object
  Returntype : Bio::EnsEMBL::SeqRegionSynonym
  Exceptions : none
  Caller     : Bio::EnsEMBL::SeqRegionSynonymAdaptor
  Status     : At Risk
=cut

sub new {
  my ($class, @args) = @_;
  
  my $self = bless {},$class;

  my ( $adaptor, $synonym, $ex_db, $seq_region_id, $dbid) =
    rearrange ( ['ADAPTOR','SYNONYM','EXTERNAL_DB_ID','SEQ_REGION_ID','DBID'], @args );

  $self->adaptor($adaptor);

  if( defined $ex_db ) { $self->external_db_id( $ex_db ) }
  if( defined $seq_region_id ) { $self->seq_region_id( $seq_region_id ) }
  if (defined $dbid) { $self->{'dbID'} = $dbid}

  if( defined $synonym ) { 
    $self->name( $synonym ) ;
  } else {
    warn "No alternative name given\n";
    return undef;
  }

  return $self;
}

sub name{
 my $self = shift;
  $self->{'name'} = shift if(@_);
  return $self->{'name'};
}

sub external_db_id{
 my $self = shift;
  $self->{'ex_db'} = shift if(@_);
  return $self->{'ex_db'};
}

sub seq_region_id{
  my $self = shift;
  $self->{'seq_region_id'} = shift if(@_);
  return $self->{'seq_region_id'};
}

1;
