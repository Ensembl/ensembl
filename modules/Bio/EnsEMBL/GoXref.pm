
#
# EnsEMBL module for GoXref.pl
#
# Cared for by Arne Stabenau <stabenau@ebi.ac.uk>
#
# Copyright EnsEMBL 2000-2003
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::GoXref

=head1 DESCRIPTION

This class extends the DBEntry in order to associate Evidence Tags to the 
relationship between EnsEMBL objects and GO identifiers.  The relationship
to GO that is stored in the database is actually derived through the 
relationship of EnsEMBL peptides to SwissProt peptides.

I.e. The relationship is derived like this:

ENSP  -> SWISSPROT -> GO
                   
An the evidence tag describes the relationship between the SWISSPROT Peptide 
and the GO entry.

In reality, however, we store this in the database like this:

ENSP -> SWISSPROT
ENSP -> GO

and the evidence tag hangs off of the relationship between the ENSP and
the GO identifier. Some ENSPs are associated with multiple closely related 
Swissprot entries which may both be associated with the same GO identifier but
with different evidence tags.  For this reason a single 'GoXref' can have 
multiple evidence tags.


=head1 SYNOPSIS

my $goxref = Bio::EnsEMBL::GoXref->new;
$goxref->add_linkage_type('IEA')

foreach my $evtag (@{$goxref->get_all_linkage_types()}) {
  print "$evtag\n";
}

=head1 CONTACT

Post questions to the ensembl development list: <ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut

package Bio::EnsEMBL::GoXref;
use vars qw(@ISA);
use strict;


@ISA = qw( Bio::EnsEMBL::DBEntry );


=head2 add_linkage_type

  Arg [1]    : string $value
    allowed are "IC", "IDA", "IEA", "IEP", "IGI", "IMP", "IPI", "ISS",
                "NAS", "NS", "TAS", "NR"
  Arg [2]    : (optional) Bio::EnsEMBL::DBEntry $source
  Example    : $go_xref->add_linkage_type('IGI');
  Description: Associates a linkage type and source DBEntry with this go_xref
  Returntype : integer; number of linkages
  Exceptions : thrown if $linkage_type argument not supplied or
               the optional DBEntry is not a DBEntry object.
  Caller     : DBEntryAdaptor
  Status     : Experimantal

=cut

sub add_linkage_type {
  my $self = shift;
  my $lt = shift;
  my $source_dbentry = shift || undef();
  
  $self->throw("linkage type argument required") if(!$lt);
  $self->throw("source_xref must be a Bio::EnsEMBL::DBEntry")
      if($source_dbentry 
	 and ! UNIVERSAL::isa($source_dbentry,'Bio::EnsEMBL::DBEntry'));

  $self->{'linkage_types'} ||= [];
  push @{$self->{'linkage_types'}}, [$lt, ( $source_dbentry || () )];
}


=head2 get_all_linkage_info

  Arg [1]    : none
  Example    : foreach (@{$gox->get_all_linkage_info})
                { print "evidence: $_->[0] via $_->[1]->display_id" }
  Description: Retrieves a list of evidence-tag/source-DBEntry pairs
               associated with this go_xref
  Returntype : listref of listrefs
  Exceptions : none
  Caller     : geneview? general.
  Status     : Experimental

=cut

sub get_all_linkage_info {
  my $self = shift;

  return $self->{'linkage_types'} || [];
}


=head2 get_all_linkage_types

  Arg [1]    : none
  Example    : foreach my $et (@{$goxr->get_all_linkage_types}){ print "$et ";}
  Description: Retrieves a unique list of evidence tags associated with 
               this go_xref
  Returntype : none
  Exceptions : none
  Caller     : geneview? general
  Status     : Stable

=cut

sub get_all_linkage_types {
  my $self = shift;
  my %seen;
  return[ grep{!$seen{$_}++} map{$_->[0]} @{$self->{'linkage_types'}} ];
  #return [ map{ $_->[0]} @{ $self->{'linkage_types'} || [] } ];
}


=head2 flush_linkage_types

  Arg [1]    : none
  Example    : $goxr->flush_linkage_types
  Description: Removes any associated evidence tags
  Returntype : none
  Exceptions : none
  Caller     : general 
  Status     : Stable

=cut

sub flush_linkage_types {
  my $self = shift;
  
  $self->{'linkage_types'} = [];
}

1;
