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

Bio::EnsEMBL::AltAlleleGroup

=head1 SYNOPSIS

  use Bio::EnsEMBL::AltAlleleGroup;
  use Bio::EnsEMBL::DBSQL::AltAlleleGroupAdaptor;
  
  my $aag_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor("Human","core","AltAlleleGroup");
  
  $aag_adaptor->fetch_all_Groups_by_type('MANUAL');

=head1 DESCRIPTION

    Alt allele groups keep track of which alleles are tied to a particular Gene
    They allow related genes to be located.
=cut

package Bio::EnsEMBL::AltAlleleGroup;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(warning);

use base qw/Bio::EnsEMBL::Storable/;

=head2 new

  Arg [-MEMBERS]: A list reference of [gene_id,is_ref,type]
             :gene_id is a dbID for Gene (consistent only within one release)
             :is_ref is a boolean flag denoting the reference allele
             :type is a string for descriptive purposes
  Example    : $aag = Bio::EnsEMBL::AltAlleleGroup->new(...);
  Description: Creates a new alt-allele group object
  Returntype : Bio::EnsEMBL::AltAlleleGroup
  Exceptions : none
  Caller     : general

=cut

sub new {
    my $caller = shift;

    my $class = ref($caller) || $caller;
    my $self = $class->SUPER::new(@_);
    my ( $list ) = rearrange( [ 'MEMBERS'], @_ );
    
    $self->{'MEMBERS'} = $list;
    
    return $self;
}

sub ref_Gene_id {
    my $self = shift;
    my $list = $self->{'MEMBERS'};
    
    foreach my $allele (@$list) {
        my ($gene_id,$is_ref,$type) = @$allele;
        if ($is_ref) {
            return $gene_id;
        }
    }
    warning("No reference allele set for this group");
    return;
}

sub get_all_Gene_ids {
    my $self = shift;    
    my $list = $self->{'MEMBERS'};
    
    my @gene_ids;
    
    foreach my $allele (@$list) {
        my ($gene_id,$is_ref,$type) = @$allele;
        push @gene_ids,$gene_id;
    }
    return \@gene_ids;
}

1;