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

Bio::EnsEMBL::OntologyXref

=head1 DESCRIPTION

This class extends the DBEntry in order to associate Evidence Tags
to the relationship between EnsEMBL objects and ontology accessions
(primarily GO accessions).

The relationship to GO that is stored in the database is actually
derived through the relationship of EnsEMBL peptides to SwissProt
peptides, i.e. the relationship is derived like this:

  ENSP -> SWISSPROT -> GO

And the evidence tag describes the relationship between the SwissProt
Peptide and the GO entry.

In reality, however, we store this in the database like this:

  ENSP -> SWISSPROT
  ENSP -> GO

and the evidence tag hangs off of the relationship between the ENSP and
the GO identifier.  Some ENSPs are associated with multiple closely
related Swissprot entries which may both be associated with the same GO
identifier but with different evidence tags.  For this reason a single
'OntologyXref' can have multiple evidence tags.

=head1 SYNOPSIS

  my $ontology_xref = Bio::EnsEMBL::OntologyXref->new();
  $ontology_xref->add_linkage_type('IEA');

  foreach my $evtag ( @{ $ontology_xref->get_all_linkage_types() } ) {
    print "$evtag\n";
  }

=head1 METHODS

=cut

package Bio::EnsEMBL::OntologyXref;

use strict;

use base qw( Bio::EnsEMBL::DBEntry );

=head2 add_linkage_type

  Arg [1]    : string $value
               allowed values:
               'IC', 'IDA', 'IEA', 'IEP', 'IGI', 'IMP', 'IPI',
               'ISS', NAS', 'ND', 'TAS', 'NR', 'RCA'
  Arg [2]    : (optional) Bio::EnsEMBL::DBEntry $source
  Example    : $ontology_xref->add_linkage_type('IGI');
  Description: Associates a linkage type and source DBEntry with
               this ontology_xref
  Returntype : integer; number of linkages
  Exceptions : thrown if $linkage_type argument not supplied or
               the optional DBEntry is not a DBEntry object.
  Caller     : DBEntryAdaptor
  Status     : Experimantal

=cut

sub add_linkage_type {
  my ( $self, $lt, $source_dbentry ) = @_;

  if ( !defined($lt) ) {
    $self->throw("linkage type argument required");
  }

  if ( defined($source_dbentry)
       && !$source_dbentry->isa('Bio::EnsEMBL::DBEntry') )
  {
    $self->throw("source_dbentry must be a Bio::EnsEMBL::DBEntry");
  }

  $self->{'linkage_types'} ||= [];

  push @{ $self->{'linkage_types'} },
    [ $lt, ( $source_dbentry || () ) ];
}


=head2 get_all_linkage_info

  Arg [1]    : none
  Example    :

    foreach ( @{ $ontology_xref->get_all_linkage_info() } ) {
      print "evidence: $_->[0] via $_->[1]->display_id";
    }

  Description: Retrieves a list of evidence-tag/source-DBEntry pairs
               associated with this ontology_xref
  Returntype : listref of listrefs
  Exceptions : none
  Caller     : geneview? general.
  Status     : Experimental

=cut

sub get_all_linkage_info {
  my ($self) = @_;

  return $self->{'linkage_types'} || [];
}


=head2 get_all_linkage_types

  Arg [1]    : none
  Example    :

    print( join( ' ', @{ $ontology_xref->get_all_linkage_types() } ),
           "\n" );

  Description: Retrieves a unique list of evidence tags associated with
               this ontology_xref
  Returntype : none
  Exceptions : none
  Caller     : geneview? general
  Status     : Stable

=cut

sub get_all_linkage_types {
  my ($self) = @_;

  my %seen;
  return [ grep  { !$seen{$_}++ }
             map { $_->[0] } @{ $self->{'linkage_types'} } ];

  #return [ map{ $_->[0]} @{ $self->{'linkage_types'} || [] } ];
}


=head2 flush_linkage_types

  Arg [1]    : none
  Example    : $ontology_xref->flush_linkage_types();
  Description: Removes any associated evidence tags
  Returntype : none
  Exceptions : none
  Caller     : general 
  Status     : Stable

=cut

sub flush_linkage_types {
  my ($self) = @_;

  $self->{'linkage_types'} = [];
}

1;
